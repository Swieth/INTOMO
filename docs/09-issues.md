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
