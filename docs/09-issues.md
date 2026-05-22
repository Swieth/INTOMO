# Known Issues (Current)

This document lists only unresolved issues in INTOMO v1.0.

---

## 1. Portability

### P-04 — `downloadORB.m` uses `uncompress` (Unix only)

**File:** [`Tomography/READ/downloadORB.m`](../Tomography/READ/downloadORB.m), line 40

```matlab
system(['uncompress ' pathORB gps_names{j} '.Z']);
```

`uncompress` is a Unix command. On Windows this system call fails. The function also depends on legacy FTP host assumptions that may no longer serve authoritative IGS combined SP3 products, so users may still need manual SP3 placement.

---

## 2. Logical Bugs

### L-02 — First-epoch Kalman path always calls `getGain` regardless of robust flag

**File:** [`Tomography/ENGINE/intomolab.m`](../Tomography/ENGINE/intomolab.m), lines 339–343

```matlab
if switches.run_KALMAN==0 || switches.run_KALMAN==1
    ...
    [K, con(epoch,1), Y(epoch).Y, ~, thin] = getGain(Pminus, A_k, R_SWD_k);
end
```

The condition is always true when `switches.run_KALMAN` is either `0` or `1`, so first epoch always uses `getGain`, even when robust filtering is selected.

### L-06 — `TomoLSQ.m` references `switches.cutModel` which is never set

**File:** [`Tomography/ENGINE/TomoLSQ.m`](../Tomography/ENGINE/TomoLSQ.m), line 28

```matlab
if switches.cutModel
```

`switches.cutModel` is not defined in `conf.m` and not initialized elsewhere; LSQ mode can fail with an undefined field error.

---

## 3. Robustness

### R-04 — `findRO.m` may return undefined `coordRO`

**File:** [`Tomography/READ/findRO.m`](../Tomography/READ/findRO.m), lines 111–114

```matlab
if ~isempty(names)
    ...
    if ~exist('coordRO','var')
        coordRO = [];
        warning('findRO: No matching radiooccultation found')
    end
end
% No fallback when names is empty
```

If `pathRO` is empty or contains no files, outer `if ~isempty(names)` is skipped and `coordRO` may remain unassigned.

### R-05 — Silent error suppression with bare `try/catch` or empty `catch`

Several locations wrap critical operations in low-information `try/catch` blocks, reducing diagnosability:

| Location | Risk |
|----------|------|
| `RUNINTOMO.m` — station coordinate read | Re-throws a generic error message that can hide root cause |
| `intomolab.m` — consistency check block | Incomplete/ineffective error handling path |
| `intomolab.m` — epoch processing | Per-epoch exceptions collapsed into generic message, stack lost |
| `intomolab.m` — robust Kalman gain inversion | Warning-only path without recovery |
| `intomolab.m` / `RUNINTOMO.m` — matrix save paths | Warning-only on failures in critical outputs |

---

## 4. New issues identified (2026-05-22)

### L-07 — `covarianceAprRT.m` produces NaN in `modQ`, triggering destructive "Cutting observations"

**File:** [`Tomography/ENGINE/covarianceAprRT.m`](../Tomography/ENGINE/covarianceAprRT.m)

Two paths produce NaN values in the output `modQ` vector:

1. `interp1(href, modQref, h, 'spline')` where `href = [150, ..., 11750]` m. The TOMO grid includes `h = 0 m` and `h = 14500 m`, both outside this table; `'spline'` extrapolation produces large or NaN values for out-of-range inputs.

2. When `switches.totalN` is true, `Pmean1(lvl) = abs(mean(values.Nw_apr_num - values.aprioriEra))` propagates NaN if either array contains NaN at that level.

The NaNs flow into `P_apriori(epoch).P_apriori` and then into the stacked `R_SWD_k` in `matrices_epochRT.m`. The existing handler there deletes NaN rows and truncates `R_SWD_k` to `size(A_k,1)` without tracking which rows are dropped, potentially silently removing real GNSS observations.

**Diagnostics added:** `matrices_epochRT.m` now reports per-source NaN counts and rewords the warning to make the destructive truncation explicit.

**Suggested fix (not yet applied):** Replace `'spline'` with `'pchip'` and clamp `h` to `[min(href), max(href)]` before the call; or impute NaN-free values using the empirical polynomial fallback already present in the function.

---

### L-08 — `coordVector.m` dead interpolation block and unsafe NaN-fill indexing (fixed)

**File:** [`Tomography/PPROCESS/coordVector.m`](../Tomography/PPROCESS/coordVector.m)

Two defects were present and have been corrected:

1. **Dead interpolation loop (lines 13–18):** A loop filled intermediate NaN-coordinate rows by linear interpolation between valid boundary points, but the very next line (`coordV = diff(coordV(id_coord,:))`) selected only the original valid rows, discarding all interpolated output. The loop had zero effect on results.

2. **Unsafe NaN-fill (`coordV(min(id)-1,:)`):** When the sentinel NaN row was appended at the end, `min(id)` was the index of that row. If `min(id) == 1` (first row is NaN), the expression `min(id)-1 == 0` caused an indexing error, silently swallowed by the `try/catch` in `groundRT.m` and `spaceRT.m`, resulting in a silently zeroed `Avec` column.

Both issues are now fixed. The dead loop is removed; the NaN-fill is guarded with an explicit error if no valid direction row precedes the first NaN.

**Note:** The `i_pos` argument is accepted but not used; it is retained for backward compatibility with call sites in `groundRT.m` and `spaceRT.m`.

---

### R-06 — Post-run `save` crash: `free(): chunks in smallbin corrupted` (deferred)

**File:** [`Tomography/RUNINTOMO.m`](../Tomography/RUNINTOMO.m), line 291

```matlab
save([PATH_SAVE, PROJECT_NAME, '/OUT/OUTPUT_TOMO.mat'], 'output')
```

MATLAB terminates with a glibc heap-corruption error immediately after `Tomography processing completed` is printed. The `save` call is the first operation after `intomolab` returns. Heap corruption of this kind is almost always caused by a MEX file writing past an allocated buffer earlier in the run (prime suspects: cspice MEX functions called during ray tracing, or HDF5/`-v7.3` serialisation of a large nested struct).

**To diagnose:**
- Run with `save(..., '-v6')` (MAT4/5 format, no HDF5) to determine whether the crash is format-specific.
- Selectively save sub-fields of `output` to isolate which struct triggers the corruption.
- Run under `valgrind --tool=memcheck` with a small synthetic dataset if MATLAB's `-batch` mode supports it.
