# Performance Assessment — INTOMO v1.0

> Analysis date: 2026-05-08  
> Reference configuration: `exp3_20260101` (24-epoch single-day run, `conf.m` as committed)  
> Method: full static code review of all first-party `.m` sources

---

## 1. Executive Summary

INTOMO is a compute-intensive pipeline whose runtime is dominated by **3-D iterative ray tracing**. On a representative 24-epoch, ~100-station run with integrated RO, the theoretical lower bound on wall-clock time is several hours using a single CPU core. The primary levers for improvement — in order of estimated impact — are:

| Rank | Intervention | Est. Speedup | Effort |
|------|-------------|-------------|--------|
| 1 | Parallelize ray tracing with `parfor` | 4–16× | Low |
| 2 | Remove per-ray console I/O inside loops | 1.3–3× | Low |
| 3 | Cache epoch-invariant covariance matrices | 1.1–1.5× | Low |
| 4 | Switch A matrix to sparse storage | 1.2–2× | Medium |
| 5 | Replace `pinv` with `\` in linear solvers | 1.1–1.5× | Low |
| 6 | Vectorize per-ray-point inner loops | 1.1–1.4× | Medium |
| 7 | Use faster MAT-file format for amtrix cache | 1.05–1.2× | Low |
| 8 | Merge triple `excessphase` calls into one pass | 1.05–1.15× | Medium |

Parallelization alone, with MATLAB's Parallel Computing Toolbox and 8 workers, could reduce a multi-hour run to under an hour without any algorithmic changes.

---

## 2. Pipeline Architecture and Timing Context

The pipeline has two structurally different cost centres.

### 2.1 Configuration & Setup (RUNINTOMO.m, runs once)

These stages run once before the epoch loop and their cost is negligible relative to the per-epoch work:
- Station read and filter (`readBLh`, `boundingTOMOLAB`, `deleteStat`)
- ERA5 grid construction (`ERA5GridTOMO`) — cached in `model.mat`
- SP3 orbit read and interpolation (`readSP3dat`, `interSP3`)
- ZTD/RO data read (`readtxtOBS`, `findRO`)
- `station.mat` construction (`construct_station_LAB`)

### 2.2 Per-Epoch Processing (intomolab.m, the dominant cost)

The epoch loop in `intomolab.m` calls the following in sequence for each of **E** epochs:

```
initialOBS   →   covarianceAprRT   →   matrices_epochRT
             →   weightObs         →   getGain / getGainR
             →   (save KAL .mat)
```

`initialOBS` itself calls `spaceRT` (for RO) and a nested loop over stations/satellites that calls `groundRT` → `voxel_dist_3D_combined` for each ray. This is where nearly all CPU time is spent.

### 2.3 Scale Numbers for Reference Configuration

| Quantity | Value | Source |
|----------|-------|--------|
| Epochs | 24 | 1-day run, 3600 s interval |
| Tomographic voxels (inner) | 6 × 14 × 11 = **924** | `model.lat_TOMO`, `model.lon_TOMO`, `model.levels_TOMO` |
| RT grid size | 129 × 193 × ~95 ≈ **2.4 M cells** | `model.res = 0.25°`, `levels_TOMO_RT` |
| GNSS rays per epoch (estimated) | ~800 (100 stations × 8 visible sats) | typical dense European network |
| RT outer iterations per GNSS ray | up to **5** | `n_iterend = 5` (ground mode) |
| RT outer iterations per RO sample | up to **30** | `n_iterend = 30` (default) |
| RO samples per occultation | ~50–250 | depends on `switches.ROres = 2` |

For a single epoch: ~800 `voxel_dist_3D_combined` calls for GNSS plus potentially thousands for RO.

---

## 3. Bottleneck Analysis

### B-1. Iterative Ray Tracing — `voxel_dist_3D_combined` ★★★★★

**Files:** `RayTracer/Tomo_Fcn/voxel_dist_3D_combined.m`, called from `groundRT.m` and `spaceRT.m`

This is the single most expensive function in the entire codebase. It implements a shooting-method ray tracer with an outer correction loop and a dense inner propagation loop.

**Structure:**
```
while diff_dist > dt                    % outer correction: up to n_iterend iterations
    for each propagation step:          % inner: ~200–20,000 steps per ray
        cspice_recgeo(...)              %   geodetic coordinate conversion (×2–3)
        refr_interp_3D(...)             %   3-D trilinear interpolation on RT grid
        dp_alt / dpacalc / dhrefcalc   %   gradient computation (3 functions)
        gradrec(...)                   %   spherical-to-Cartesian transform
    end
    savevar(...)                        % full ray snapshot saved per iteration
end
excessphase(...)  ×3                   % post-loop: called for dLs, dLw, dLh
```

**Cost drivers:**
- `cspice_recgeo` is a compiled MEX function but called 2–3 times per step. For a GNSS ray with `ds = dsL/50000`, `dsL ≈ 20,200 km`, and 5 outer iterations: roughly **5 × (steps_per_trace) × 3 cspice calls**. At moderate elevations this is tens of thousands of cspice calls per single ray.
- `refr_interp_3D` is called once per step when inside the atmosphere and performs a nearest-layer lookup plus interpolation in a 2.4 M-cell grid.
- `savevar` copies the entire growing ray coordinate arrays on every outer iteration — this is O(n²) in the number of steps.

**No parallelism.** The loop over rays is purely sequential.

---

### B-2. Sequential Ray Loop in `initialOBS` ★★★★★

**File:** `Tomography/ENGINE/initialOBS.m`, lines 84–102

```matlab
for nr = 1:stac                        % stations — sequential
    sat_number = length(station.h(nr).satellite);
    while sat < sat_number             % satellites — sequential
        sat = sat + 1;
        ...
        groundRT(...)                  % calls voxel_dist_3D_combined
    end
end
```

All station–satellite pairs are processed one by one. Each pair is **fully independent**: no shared mutable state that prevents parallelisation. With `parfor` (Parallel Computing Toolbox) over the outer `nr` loop (or a flattened index loop), all available CPU cores could be exploited with zero algorithmic change.

---

### B-3. Sequential RO Tracing Loop in `spaceRT` ★★★★☆

**File:** `Tomography/ENGINE/spaceRT.m`, lines 47–260

```matlab
for roit = 1:size(station.ro,2)        % occultations — sequential
    for i = idex:step:limit            % along-track samples — sequential
        voxel_dist_3D_combined(...)    % calls RT per sample
        interateModelsA(...)
    end
end
```

Each occultation (outer `roit`) is independent. The inner `i` loop processes samples of a single occultation sequentially; samples within one occultation *are* independent and could also be parallelised. With `switches.ROres = 2`, a 50 Hz occultation file sampled at every 2nd point still yields ~200–500 ray-trace calls per occultation.

---

### B-4. Dense Kalman Gain Computation — `getGain` / `getGainR` ★★★☆☆

**Files:** `Tomography/ENGINE/getGain.m`, `getGainR.m`

The Kalman gain is:

```
K = P⁻ · Aᵀ · (A·P⁻·Aᵀ + R)⁻¹
```

Computed via full SVD:
```matlab
Anew = (A * Pminus * A' + R_SWD.^2);   % (m×m) dense matrix, m = observations
[U,S,V] = svd(Anew, 'econ');           % O(m² × min(m,n)) cost
```

With ~800 GNSS + ~200 RO observations per epoch (`m ≈ 1000`) and `n ≈ 924` state variables:
- `A * Pminus * A'` forms a 1000×1000 dense matrix from a 1000×924 and 924×924 matrix multiplication: **O(m² × n)** ≈ 8.5 × 10⁸ FLOPs
- `svd` on 1000×1000: **O(m³)** ≈ 10⁹ FLOPs

This is called every epoch (24 times). The covariance matrix `Pminus` (924×924) is also dense, stored and updated as a full matrix.

**Note:** `Pminus` in `intomolab.m` (line 428) is updated as:
```matlab
Pminus = Phi * Pplus * Phi' + Q;
```
`Pplus` starts dense from `Pzero = diag(Pzero.^2)` (line 337) and stays dense. Using diagonal or banded approximations of the state covariance is a standard technique in large-scale Kalman filtering that is not yet exploited.

---

### B-5. `weightObs` — `pinv` on Observation Normal Matrix ★★★☆☆

**File:** `Tomography/ENGINE/weightObs.m`, lines 231 and 236/264

```matlab
C = pinv(var_obs);             % O(n³) pseudo-inverse, dense full matrix
...
x = pinv(A'*C*A)*A'*C*SWD_obs;  % further pinv of normal matrix
```

`pinv` computes a full SVD internally. For a dense `var_obs` of size (m × m) where m ≈ number of observations, this is **O(m³)** per epoch. The `\` (backslash) operator would be at minimum 3–5× faster for well-conditioned systems. The outer `while` loop (lines 261–275) iterates up to 3 times with a `pinv` each iteration.

---

### B-6. Per-Ray Console Output Inside Loops ★★★☆☆

The following print statements execute **on every ray step or ray**:

| File | Line | Statement | Called how often |
|------|------|-----------|-----------------|
| `groundRT.m` | 38, 42 | `tic` / `toc` pair | Once per GNSS ray — prints elapsed time every ray |
| `initialOBS.m` | 89 | `disp(['Slant number: ...'])` | Once per GNSS ray |
| `spaceRT.m` | 115 | `fprintf('\n spaceRT: ...')` | Once per RO sample |
| `interateModelsA.m` | 46–49 | `fprintf('interateModelsA: ...')` | Once per ray |

MATLAB's `disp`/`fprintf` calls are synchronous and involve OS-level I/O. Inside a tight loop of 800 rays per epoch × 24 epochs = 19,200 calls, this is a significant source of wall-clock overhead, especially when running in batch mode (`matlab -batch`). Measured overhead is typically 5–30 ms per print call depending on terminal buffer speed.

---

### B-7. Epoch-Invariant Covariance Recomputed Every Epoch ★★☆☆☆

**Files:** `Tomography/ENGINE/covarianceNwRT.m`, `covarianceAprRT.m`

Both functions compute:
```matlab
for i = 1:size(Nw_apr, 1)             % loop over ALL epochs every call
    NWPone = reshape(Nw_apr(i,:), num_lat*num_lon, num_lev-1);
    NWPext(i,:) = mean(NWPone);
end
modQ = cov(NWPext);                    % covariance over all epochs
```

This loops over all **E** epoch rows of `Nw_apr` every time it is called. With `DETER` apriori, `Nw_apr` is fixed before the epoch loop, so `modQ` (and consequently `Q = diag(modQ.^2)`) does not change between epochs. The computation is currently repeated 24 times unnecessarily.

---

### B-8. amtrix Files Use HDF5 Format ★★☆☆☆

**File:** `Tomography/ENGINE/initialOBS.m`, line 192

```matlab
save(save_filename, 'A', 'Avec', ..., '-v7.3');
```

The `-v7.3` flag selects HDF5/MAT7.3 format, which has significantly higher read/write overhead than the legacy `-v7` (zlib-compressed) or `-v6` (uncompressed) formats for matrices below ~2 GB. For small- to medium-sized `A` matrices (924 × ~1000), the HDF5 overhead on every new epoch can add several seconds per file.

---

### B-9. Array Growth by Concatenation in Inner Loops ★★☆☆☆

**Files:** `groundRT.m` (lines 63–65), `spaceRT.m` (lines 185–186, 195–196)

```matlab
iposEM1 = [iposEM1 iposver];    % appends to growing array in each loop pass
iposEM2 = [iposEM2 distver];
```

MATLAB does not optimise in-place concatenation; each `[A B]` allocates a new array and copies both inputs. Inside loops over `rayRT.voxEM.iposEm` (which can have 10–100 rows for a single ray), this causes **O(rows²)** total memory operations. Pre-allocating to maximum size and using an index counter eliminates this.

---

### B-10. `excessphase` Called Three Times per Ray ★★☆☆☆

**File:** `Tomography/ENGINE/interateModelsA.m`, lines 61–65 (and equivalently in `voxel_dist_3D_combined.m`)

```matlab
switches.refopt = 1; [rayIND, dLs] = excessphase(rayIND, switches);
switches.refopt = 2; [rayIND, dLw] = excessphase(rayIND, switches);
switches.refopt = 3; [rayIND, dLh] = excessphase(rayIND, switches);
```

`excessphase` iterates over the ray-point arrays once per call. The three calls differ only in which refractivity component (`refr`, `refw`, or `refh`) is used for the trapezoidal integration. A single combined pass computing all three integrals simultaneously would reduce this overhead by ~3×.

---

### B-11. Full Dense State Covariance `Pplus` Stored and Updated ★★☆☆☆

**File:** `Tomography/ENGINE/intomolab.m`, lines 347–349, 428, 454

```matlab
Pplus  = Pminus - K * A_k * Pminus;   % (924×924) dense matrix operation
Pminus = Phi * Pplus * Phi' + Q;       % (924×924) matrix triple-product
```

With 924 state variables, `Pplus` and `Pminus` are 924×924 dense matrices (~6.5 MB each in double). The triple product `Phi * Pplus * Phi'` where `Phi = diag(...)` is a full O(n³) = 924³ ≈ 7.9 × 10⁸ FLOPs. Since `Phi` is diagonal (set on line 411), the update simplifies to an elementwise row/column rescaling: `Pplus_new[i,j] = phi[i] * phi[j] * Pplus[i,j] + Q[i,j]`, which is O(n²). This optimisation is not currently applied.

---

### B-12. `deleteStat` pairwise-distance matrix ★☆☆☆☆

**File:** `Tomography/PPROCESS/deleteStat.m`, line 18

```matlab
distances = squareform(pdist(coords));   % O(S²) with S = #stations
```

`squareform(pdist(...))` builds a full S×S distance matrix and stores it. For S = 500 stations this is a 500×500 = 250,000-element matrix built once at startup. Cost is negligible for typical station counts (< 1000) and the function is called only once.

---

## 4. Recommendations

### R-1. Parallelize ray tracing with `parfor` (★★★★★ impact, Low effort)

**`initialOBS.m`:** Replace the sequential station loop with a flat `parfor` over all (station, satellite) pairs:

```matlab
% Build a flat list of all (station, satellite) pairs first
ray_jobs = [];
for nr = 1:stac
    for s = 1:length(station.h(nr).satellite)
        ray_jobs(end+1,:) = [nr, s];
    end
end

parfor k = 1:size(ray_jobs, 1)
    nr  = ray_jobs(k, 1);
    sat = ray_jobs(k, 2);
    [A_col, SWD_k, ...] = groundRT(model, station, epoch, nr, sat, switches);
    % accumulate into output arrays (requires careful indexing)
end
```

Because each GNSS ray writes to a unique column of `A`, there is no write conflict. MATLAB's `parfor` requires explicit handling of slice-parallel outputs; using a cell array indexed by `k` and assembling after the loop is the standard pattern.

**`spaceRT.m`:** Similarly, parallelize over `roit` (occultation index). The inner `i` (sample) loop per occultation can remain sequential unless individual occultations are too short to justify the overhead.

With MATLAB Parallel Computing Toolbox and 8 workers, the expected speedup from parallelising the ray loop alone is **4–8×** on the ray-tracing stage (which dominates runtime).

---

### R-2. Remove per-ray console output from inner loops (★★★☆☆ impact, Low effort)

Remove or gate the following with a verbosity flag:

| Location | Change |
|----------|--------|
| `groundRT.m` lines 38, 42 | Remove `tic`/`toc` per ray; move outside loop or delete |
| `initialOBS.m` line 89 | Suppress or throttle: print only every N-th ray |
| `spaceRT.m` line 115 | Suppress or gate with `switches.verbose` |
| `interateModelsA.m` lines 46–49 | Remove `fprintf` per ray |

Suggested pattern for any per-ray logging:
```matlab
if mod(b, 50) == 0          % print every 50th ray
    fprintf('Ray %d / %d\n', b, total_rays);
end
```

---

### R-3. Cache epoch-invariant Q matrix (★★☆☆☆ impact, Low effort)

In `intomolab.m`, move the `covarianceNwRT` call outside the epoch loop and reuse the result:

```matlab
% Compute Q once before the epoch loop (Nw_apr is fixed in DETER mode)
Q_cached = covarianceNwRT(Nw_apr, A(1).A, num_lat, num_lon, ...
                          size(BLh_pudel_proj,3), BLh_pudel_proj_num);

for epoch = 2 : size(solution_epochs, 2)
    ...
    Q = Q_cached;        % reuse instead of recomputing
    ...
end
```

Guard this with a check: if `switches.aprModel == 'ERA5'` the apriori changes per epoch and the cache should be invalidated per epoch.

---

### R-4. Use sparse A matrix (★★★☆☆ impact, Medium effort)

The observation matrix `A` (voxels × observations) is extremely sparse: each ray traverses at most 20–50 of the 924 voxels, giving a fill density of ~3–5 %. Declaring `A` as sparse from the start and passing `sparse(A)` into the linear solvers would:
- Reduce memory usage by ~20×
- Speed up `A * Pminus * A'` by exploiting sparsity (MATLAB's sparse BLAS routines)
- Speed up `A' * C * A` in `weightObs` similarly

In `initialOBS.m`, change the pre-allocation:
```matlab
A = sparse((model.num_lat_TOMO-1)*(model.num_lon_TOMO-1)*(model.num_levels_TOMO-1), assum_sat*stac);
```

The existing `svds` branch in `getGain`/`getGainR` already handles sparse input:
```matlab
if issparse(Anew) == 1
    [U,S,V] = svds(Anew, ile_w);
```
so the downstream SVD path is already prepared for sparse input.

---

### R-5. Replace `pinv` with `\` in `weightObs` (★★★☆☆ impact, Low effort)

```matlab
% Current (weightObs.m line 231, 236, 264)
C = pinv(var_obs);
x = pinv(A'*C*A) * A'*C*SWD_obs;

% Faster equivalents
C = var_obs \ eye(size(var_obs));      % if var_obs is full rank
x = (A'*C*A) \ (A'*C*SWD_obs);        % normal equations
```

For well-conditioned systems `\` is 3–5× faster than `pinv` because it avoids the full SVD. If the matrix is possibly rank-deficient, use `lsqminnorm` (MATLAB R2017b+) which is faster than `pinv` while still handling rank deficiency.

---

### R-6. Pre-allocate arrays in per-ray accumulation loops (★★☆☆☆ impact, Low effort)

In `groundRT.m` and `spaceRT.m`, replace:
```matlab
iposEM1 = [];
iposEM2 = [];
for it = 1:size(rayRT.voxEM.iposEm, 1)
    iposEM1 = [iposEM1 iposver];   % grows every iteration
    iposEM2 = [iposEM2 distver];
end
```

With:
```matlab
max_em = size(rayRT.voxEM.iposEm, 1) * size(rayRT.voxEM.iposEm, 2);
iposEM1 = zeros(1, max_em);
iposEM2 = zeros(1, max_em);
ptr = 0;
for it = 1:size(rayRT.voxEM.iposEm, 1)
    n = numel(iposver);
    iposEM1(ptr+1:ptr+n) = iposver;
    iposEM2(ptr+1:ptr+n) = distver;
    ptr = ptr + n;
end
iposEM1 = iposEM1(1:ptr);
iposEM2 = iposEM2(1:ptr);
```

---

### R-7. Use faster MAT-file format for amtrix cache (★★☆☆☆ impact, Low effort)

In `initialOBS.m` line 192, change:
```matlab
save(save_filename, 'A', 'Avec', 'SWD', 'SIWV', 'R_SWD', ..., '-v7.3');
```
to:
```matlab
save(save_filename, 'A', 'Avec', 'SWD', 'SIWV', 'R_SWD', ...);  % defaults to -v7
```

The `-v7` default format uses zlib compression but is pure MAT (not HDF5), with significantly lower metadata overhead for modest-sized matrices. Only switch back to `-v7.3` if any variable exceeds 2 GB.

---

### R-8. Optimise diagonal Phi update in Kalman prediction (★★☆☆☆ impact, Low effort)

In `intomolab.m` (line 428), `Phi = diag(phi_vec)` is a diagonal matrix, so `Phi * Pplus * Phi'` is equivalent to an outer product rescaling:

```matlab
% Current: O(n³) triple matrix product
Pminus = Phi * Pplus * Phi' + Q;

% Optimised: O(n²) elementwise scaling
phi_vec = diag(Phi);               % extract diagonal
Pminus = (phi_vec * phi_vec') .* Pplus + Q;
```

This turns a matrix triple-product into two vectorised operations, reducing FLOPs from O(n³) to O(n²).

---

### R-9. Combine triple `excessphase` calls into one pass (★★☆☆☆ impact, Medium effort)

Refactor `excessphase.m` to accept a flag vector and return all three integrals in a single loop over ray-point arrays:

```matlab
% Instead of three calls:
[rayIND, dLs] = excessphase(rayIND, switches, 1);   % total
[rayIND, dLw] = excessphase(rayIND, switches, 2);   % wet
[rayIND, dLh] = excessphase(rayIND, switches, 3);   % dry

% Single combined call:
[rayIND, dLs, dLw, dLh] = excessphase_combined(rayIND, switches);
```

---

### R-10. Profile before committing to further optimisation

Before investing effort beyond R-1–R-7, run MATLAB's profiler to measure actual timings:

```matlab
profile on
RUNINTOMO
profile viewer
```

Focus profiling on a single epoch (modify `intomolab` to break after `epoch == 1`) to identify which function accounts for the largest fraction of time in your actual data configuration.

---

## 5. Memory Considerations

| Object | Approximate size | Notes |
|--------|-----------------|-------|
| `Pplus` / `Pminus` (dense 924×924 double) | ~6.5 MB each | Created/updated every epoch |
| `A` matrix per epoch (924 × ~1000, dense) | ~7 MB | Sparse would be ~0.3 MB |
| `model.refrRT{epoch}` grid (129×193×95) | ~190 MB | Held in memory for all epochs |
| `model.mat` (full model structure) | >> 1 GB | Saved with `-v7.3`; loading at startup may take 10–60 s |

If the full `model` structure cannot be held in RAM (e.g., on machines with < 8 GB), demand-loading individual epoch fields would reduce peak memory at the cost of extra I/O.

---

## 6. Summary Priority Table

| ID | Component | File(s) | Type | Impact | Effort | Suggested Action |
|----|-----------|---------|------|--------|--------|-----------------|
| R-1 | Ray loop parallelism | `initialOBS.m`, `spaceRT.m` | Parallelism | ★★★★★ | Low | `parfor` over ray index |
| R-2 | Per-ray console I/O | `groundRT.m`, `spaceRT.m`, `initialOBS.m`, `interateModelsA.m` | I/O | ★★★☆☆ | Low | Remove / throttle prints |
| R-3 | Covariance Q cache | `intomolab.m` | Algorithmic | ★★☆☆☆ | Low | Move outside epoch loop |
| R-4 | Sparse A matrix | `initialOBS.m`, `getGain.m`, `weightObs.m` | Data structure | ★★★☆☆ | Medium | `sparse(...)` allocation |
| R-5 | `pinv` → `\` | `weightObs.m` | Algorithmic | ★★★☆☆ | Low | Replace `pinv` calls |
| R-6 | Array pre-allocation | `groundRT.m`, `spaceRT.m` | Memory | ★★☆☆☆ | Low | Pre-allocate + index ptr |
| R-7 | MAT-file format | `initialOBS.m` | I/O | ★★☆☆☆ | Low | Remove `-v7.3` flag |
| R-8 | Diagonal Phi update | `intomolab.m` | Algorithmic | ★★☆☆☆ | Low | Elementwise rescaling |
| R-9 | Triple `excessphase` | `interateModelsA.m`, `voxel_dist_3D_combined.m` | Algorithmic | ★★☆☆☆ | Medium | Merge into one pass |
| R-10 | Profile first | — | Process | — | Low | `profile on` + viewer |
