# Configuration

All configuration is done inside [`Tomography/conf.m`](../Tomography/conf.m). When `conf.m` is executed (automatically at the start of `RUNINTOMO.m` via the `conf;` call), it saves all workspace variables to `Tomography/CONF/<PROJECT_NAME>.mat`, which the driver then loads back as needed throughout the run.

---

## 1. Time mode

```matlab
switches.time_mode = 'POSTPROCESSING';
```

Only `'POSTPROCESSING'` is implemented. Real-time processing is not available in v1.0.

---

## 2. Paths

| Variable | Description |
|----------|-------------|
| `PATH_INSTALL` | Absolute path to the repository root (trailing separator required) |
| `PATH_EXTERNALSAVE` | Absolute path to the external output directory |
| `PROJECT_NAME` | Project identifier string — used as the sub-folder name for both input data and results |

Derived paths set automatically by `conf.m`:

| Variable | Value |
|----------|-------|
| `PATH_SAVE` | Same as `PATH_EXTERNALSAVE` |
| `pathCONF` | `PATH_INSTALL + 'Tomography/CONF'` |
| `pathTOMO` | `PATH_INSTALL + 'Tomography/DATA/' + PROJECT_NAME + '/WORK/'` |
| `pathORB` | `PATH_INSTALL + 'Tomography/DATA/' + PROJECT_NAME + '/ORB/'` |
| `pathMETEO` | `PATH_INSTALL + 'Tomography/DATA/' + PROJECT_NAME + '/METEO'` |
| `pathATM` | `PATH_INSTALL + 'Tomography/DATA/' + PROJECT_NAME + '/ATM/'` |
| `pathRO` | `PATH_INSTALL + 'Tomography/DATA/' + PROJECT_NAME + '/RO'` (only when `switches.integrated = "yes"`) |
| `pathEXPORT` | `PATH_SAVE + PROJECT_NAME + '/OUT/'` |

---

## 3. Time slots and data intervals

All times use the format `[yyyy mm dd hh mm ss]`.

| Variable | Description | Example |
|----------|-------------|---------|
| `observation_start_TOMO` | First tomography epoch | `[2026 01 01 00 00 00]` |
| `observation_end_TOMO` | Last tomography epoch | `[2026 01 01 23 59 59]` |
| `estimation_interval_TOMO` | Spacing between tomography epochs [s] | `3600` |
| `observation_interval_SP3` | SP3 data interval [s] | `3600` |
| `observation_interval_ZTD` | ZTD observation interval [s] | `3600` |
| `observation_interval_METEO` | ERA5 file interval [s] | `3600` |
| `interpolation_interval_METEO` | Refractivity interpolation interval [s] | `3600` |
| `observation_interval_NWP` | NWP model update interval [s] | `3600` |
| `observation_interval_APRIORI` | A priori model update interval [s] | `3600` |

> The current `readtxtOBS.m` reader assumes an hourly cadence (`estimation_interval_TOMO = 3600`) and does not support other intervals. See [09-issues.md](09-issues.md).

---

## 4. Ray tracing domain

The ray tracing (RT) model covers a larger area than the tomography model and is used as the medium through which the 3-D ray paths travel.

| Variable | Description | Example |
|----------|-------------|---------|
| `model.lat_TOMO_RT` | Latitude grid [deg] | `[33:0.25:65]` |
| `model.lon_TOMO_RT` | Longitude grid [deg] | `[-12:0.25:36]` |
| `model.levels_TOMO_RT` | Altitude layers [km] | multi-segment array from 0.01 to 86 km |
| `model.res` | Horizontal resolution [deg] | `0.25` (matches ERA5) |
| `model.ds` | Ray step size [km] | `2` |
| `model.radii` | Earth semi-axes [km] | `[6378.137, 6356.752314245]` |
| `model.GRIDboundaries` | `[lon1 lon2 res lat1 lat2 res]` | `[-12 36 0.25 33 65 0.25]` |
| `unduFile` | Undulation `.mat` file name | `'undu.mat'` |

---

## 5. Tomography domain

The tomography (TOMO) model is the inner grid where refractivity values are estimated.

| Variable | Description | Example |
|----------|-------------|---------|
| `model.lat_TOMO` | Latitude faces [deg] | `[43:2.0:55]` |
| `model.lon_TOMO` | Longitude faces [deg] | `[-2:2.0:26]` |
| `model.levels_TOMO` | Altitude layer boundaries [m] | `[0 500 1000 ... 14500]` |
| `model.west_limit_TOMO` | Western bounding box [deg] | `-1.5` |
| `model.east_limit_TOMO` | Eastern bounding box [deg] | `25.5` |
| `model.south_limit_TOMO` | Southern bounding box [deg] | `43.5` |
| `model.north_limit_TOMO` | Northern bounding box [deg] | `54.5` |
| `model.cut_off_angle` | Minimum satellite elevation [deg] | `10` |

GNSS stations outside the bounding box are removed by `boundingTOMOLAB.m`. Stations closer than `switches.stat_range` km to each other are thinned by `deleteStat.m`.

---

## 6. The `switches` structure

`switches` controls all processing options. It is saved into the project `.mat` file and threaded through every function.

### Observation and model options

| Field | Values | Tested | Description |
|-------|--------|--------|-------------|
| `switches.integrated` | `"yes"` / `"no"` | yes | Include RO data (integrated tomography) |
| `switches.totalN` | `true` / `false` | yes | Process total refractivity (true) or wet-only (false) |
| `switches.aprModel` | `'DETER'` / `'ERA5'` | yes | A priori refractivity source |
| `switches.observations` | `{'SWD',''}` | yes | Type of observations; `'SWD'` is the only supported value |
| `switches.solution` | `{'REAL'}` / `{'SYNTHETIC'}` | yes | Use measured or synthetic delays |
| `switches.parametrization` | `{'constant'}` / `{'bilinear-h'}` | yes / partial | Voxel constant or bilinear node parameterisation |
| `switches.regular` | `{'yes'}` | yes | Regular node parameterisation only in v1.0 |
| `switches.method` | `{'KALMAN'}` / `{'LSQ'}` | yes | Kalman filter or Least Squares inversion |

### Kalman filter options

| Field | Values | Tested | Description |
|-------|--------|--------|-------------|
| `switches.filter` | `{'KALMAN'}` / `{'ROBUST'}` | partial | Standard or robust (IGGIII downweighting) Kalman filter |
| `switches.phi` | `{'identity'}` | yes | State transition matrix; `'identity'` uses a diagonal transition |
| `switches.stacking` | `{'NO'}` / `{n}` | no | Accumulate `n` epochs before solving (multi-epoch stacking) |

### Apriori and constraint options

| Field | Values | Tested | Description |
|-------|--------|--------|-------------|
| `switches.apriori` | `{'INNER','OUTER','TOP','BOTTOM'}` | partial | Which model regions receive apriori values |
| `switches.constraints` | `{'',''}` or `{'HORIZONTAL','VERTICAL'}` | no | Smoothness constraints (not operational in v1.0) |
| `switches.decorelation` | `{'NO'}` | no | Remove linearly dependent rows in A matrix |

### Input/output options

| Field | Values | Tested | Description |
|-------|--------|--------|-------------|
| `switches.coord` | `{'FORMATTED'}` | yes | Coordinate input format; only `'FORMATTED'` is supported |
| `switches.stat_range` | number [km] | yes | Minimum inter-station distance; `0` disables thinning |
| `switches.saveAtmParam` | `true` / `false` | yes | Save temperature, pressure, humidity to `model.mat` output |
| `switches.amtrix` | `'const'` / `'node'` | yes | Sub-folder name for observation matrix files (set automatically) |
| `switches.ROres` | integer | yes | RO position subsampling factor (1 = 50 Hz, 10 = 5 Hz) |
| `switches.refron` | `true` / `false` | yes | Compute ray path from refractivity values (true) or from T/p/e (false) |
| `switches.time_mode` | `'POSTPROCESSING'` | yes | Processing time mode |

### How `switches.amtrix` is set

```matlab
if strcmp(switches.parametrization,'bilinear-h')
    switches.amtrix = 'node';
else
    switches.amtrix = 'const';
end
```

Observation matrix files are saved to `WORK/<switches.amtrix>/amtrix_<epoch>.mat`.

---

## 7. Saving the configuration

`conf.m` ends with:

```matlab
save([PATH_INSTALL 'Tomography/CONF/' PROJECT_NAME '.mat']);
clearvars -except PATH_INSTALL PROJECT_NAME time_mode;
```

This persists the full workspace (all path variables, `model`, `switches`, time settings) to `Tomography/CONF/<PROJECT_NAME>.mat`. `RUNINTOMO.m` loads individual variables from this file at each stage with `load(confFILE, 'var1', 'var2', ...)`.

After running, all settings (including `obs_set`, `station`, `values`, and `apriori`) are also appended to `WORK/model.mat` for archival purposes.
