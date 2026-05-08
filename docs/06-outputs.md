# Output Files

INTOMO produces several MATLAB `.mat` files during a run. They fall into three categories: workspace state files (saved incrementally), per-epoch observation matrices, and the final tomography result.

> **Important:** The output sub-directories `OUT/` and `KAL/` under `PATH_EXTERNALSAVE/<PROJECT_NAME>/` are **not created automatically** by INTOMO. They must exist before running or the corresponding `save` calls will fail silently (they are wrapped in `try/catch` blocks). See [09-issues.md](09-issues.md).

---

## 1. Model state — `model.mat`

**Path:** `Tomography/DATA/<PROJECT_NAME>/WORK/model.mat`  
**Format:** MATLAB binary, `-v7.3` (HDF5)  
**Written by:** `RUNINTOMO.m` (twice)

### First write (after Stage 2)

Contains the `model` structure populated by `ERA5GridTOMO.m` for all epochs. This represents the ray tracing and tomography grids together with the apriori refractivity fields. If this file already exists on startup, Stage 2 is skipped and this file is loaded instead.

Key variables: `model` (see [07-data-structures.md](07-data-structures.md) for field details).

### Second write (after Stage 7)

Overwrites the first write and adds:

| Variable | Description |
|----------|-------------|
| `model` | Updated `model` structure (now includes `BLh`, `NAME`, `X/Y/Z`, station metadata) |
| `station` | Full `station` structure array (see [07-data-structures.md](07-data-structures.md)) |
| `values` | Apriori refractivity index/value arrays for inner and outer model |
| `pathTOMO` | Resolved path to WORK directory |
| `PATH_SAVE` | Output save path |
| `obs_set` | Epoch matrices for TOMO, SP3, and METEO processing |
| `switches` | Full `switches` structure |
| `apriori` | Apriori refractivity values from DETER/ERA5 |
| `PROJECT_NAME` | Project identifier string |

---

## 2. Station observations — `station.mat`

**Path:** `Tomography/DATA/<PROJECT_NAME>/WORK/station.mat`  
**Format:** MATLAB binary, default version  
**Written by:** `RUNINTOMO.m` Stage 6 (`construct_station_LAB`)

Contains the `station` structure array indexed by epoch. Each element holds:

- Ground GNSS delays (ZTD, SWD, STD) and gradients for each receiver/satellite pair
- Satellite coordinates (ECEF)
- RO excess phase data and LEO/GPS coordinates (integrated mode)

If this file already exists when `RUNINTOMO.m` runs, it is loaded rather than recomputed, allowing the ray tracing step to restart from Stage 8 without repeating data reading.

For the complete field-by-field reference see [07-data-structures.md](07-data-structures.md).

---

## 3. Per-epoch observation matrices — `amtrix_<n>.mat`

**Path:** `Tomography/DATA/<PROJECT_NAME>/WORK/<switches.amtrix>/amtrix_<epoch>.mat`  
**Format:** MATLAB binary, `-v7.3`  
**Written by:** `initialOBS.m` (Stage 8a)

One file is produced per epoch. The sub-folder name is `const` for voxel parameterisation and `node` for bilinear-h parameterisation (set via `switches.amtrix`).

| Variable | Size / Type | Description |
|----------|-------------|-------------|
| `A` | `[n_obs × n_vox]` sparse/dense | Observation matrix with refractivity path derivatives [m/TECU or dimensionless] |
| `Avec` | `[n_vox × n_obs × 3]` | Ray path direction vectors in each voxel |
| `SWD` | `[n_obs × 1]` | Observed slant delays or excess phase [mm] |
| `SIWV` | `[n_obs × 1]` | Slant Integrated Water Vapour [kg/m²] |
| `R_SWD` | `[n_obs × 1]` | SWD observation errors [mm] |
| `R_SIWV` | `[n_obs × 1]` | SIWV observation errors |
| `not_cro` | vector | Indices of voxels not crossed by any ray path in this epoch |
| `elev` | `[n_obs × 1]` | Satellite elevation angles [deg]; RO observations are assigned `elev = 0` |
| `SAT` | `[n_obs × 1]` | GPS satellite PRN numbers; RO observations are assigned `SAT = 100` |
| `station_name` | cell array | GNSS station name strings; RO rows contain `'RO'` |
| `RelDist` | `[n_obs × 1]` | Total signal path length [km] |
| `dSWD` | `[n_obs × 1]` | Difference between observed and apriori-simulated slant delay [mm] |
| `coord` | cell array | Voxel boundary intersection coordinates for RO paths |
| `SWD_nodes_integ` | vector | Integral value of SWD from apriori data: `A * Nw_apr(epoch,:)'` |

If a file already exists for a given epoch it is loaded without recomputation, enabling partial re-runs.

---

## 4. Final tomography output — `OUTPUT_TOMO.mat`

**Path:** `PATH_EXTERNALSAVE/<PROJECT_NAME>/OUT/OUTPUT_TOMO.mat`  
**Format:** MATLAB binary  
**Written by:** `RUNINTOMO.m` after `intomolab` returns

| Variable | Size | Description |
|----------|------|-------------|
| `output.xP` | `[n_epochs × n_vox]` | Estimated refractivity values per voxel per epoch [N-units] |
| `output.mxP` | `[n_epochs × n_vox]` | Standard deviation of refractivity estimates (sqrt of Kalman `Pplus` diagonal) |
| `output.n_obs` | `[1 × n_epochs]` | Number of observations used in each epoch |
| `output.n_c` | `[1 × n_epochs]` | Number of apriori constraint rows used in each epoch |

For the LSQ method, `output` additionally contains `dxP` (difference between estimated and apriori values).

---

## 5. Per-epoch Kalman diagnostics — `kalman_mat_<n>.mat`

**Path:** `PATH_EXTERNALSAVE/<PROJECT_NAME>/KAL/kalman_mat_<epoch>.mat`  
**Format:** MATLAB binary  
**Written by:** `intomolab.m` per epoch (inside `try/catch`)

Contains the `thin` structure for Kalman filter diagnostics:

| Field | Description |
|-------|-------------|
| `thin.Q` | Diagonal of process noise covariance matrix Q |
| `thin.P2` | Observation weight matrix P2 (diagonal or full) |
| `thin.R_SWD_k` | Stacked observation error variances used in the gain calculation |
| `thin.Pplus` | Diagonal of corrected state covariance (posterior) |
| `thin.Pminus` | Diagonal of predicted state covariance (prior) |
| `thin.id_del` | Indices of observations deleted by the outlier filter |
| `thin.theta` | Log of eigenvalues from SVD (first epoch only) |
| `thin.theta_ind` | Sorted `theta` values |

---

## Console output

During processing INTOMO prints progress information to the MATLAB command window:

- ERA5 file name and epoch number as each grid is computed
- `'Tomography and ray tracing domain data(model.mat file) calculated'` or `'...found'`
- `'Slant number: NNN/MMMM   <datetime>   Station: N   Satellite: N'` for each ray traced
- `'Epoch: N / M'` for each inversion epoch
- `'Tomography processing completed'`

Warnings are emitted when data are missing or optional steps fail (e.g., `'findRO: No matching radiooccultation found'`, `'intomolab: Unable to save kalman filtering errors. Missing path'`).

---

## Output dependency diagram

```
RUNINTOMO.m
│
├── Stage 2 ──────────────────── WORK/model.mat  (first write)
├── Stage 6 ──────────────────── WORK/station.mat
├── Stage 7 ──────────────────── WORK/model.mat  (second write, full)
│
└── intomolab.m
    ├── Per epoch ─────────────── WORK/<amtrix>/amtrix_<n>.mat
    ├── Per epoch ─────────────── KAL/kalman_mat_<n>.mat
    │
    └── Final ─────────────────── OUT/OUTPUT_TOMO.mat
```
