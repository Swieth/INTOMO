# Data Structures

This document describes the key MATLAB structures used throughout INTOMO. The primary source is [`INTOMO_variables.md`](../INTOMO_variables.md); where code inspection revealed additional or corrected information it has been incorporated.

---

## `model`

Populated by `ERA5GridTOMO.m` (Stages 2 and 6). Holds everything needed to describe both the ray tracing domain and the tomography inversion grid.

### Grid geometry

| Field | Description |
|-------|-------------|
| `model.lat_TOMO` | Latitude faces of TOMO grid [deg] |
| `model.lon_TOMO` | Longitude faces of TOMO grid [deg] |
| `model.levels_TOMO` | Altitude layer boundaries of TOMO grid [m] |
| `model.lat_TOMO_RT` | Latitude grid of RT domain [deg] |
| `model.lon_TOMO_RT` | Longitude grid of RT domain [deg] |
| `model.levels_TOMO_RT` | Altitude layers of RT domain [km] |
| `model.num_lat_TOMO` | Number of TOMO latitude faces |
| `model.num_lon_TOMO` | Number of TOMO longitude faces |
| `model.num_levels_TOMO` | Number of TOMO altitude faces |
| `model.num_lat_TOMO_RT` | Number of RT latitude nodes |
| `model.num_lon_TOMO_RT` | Number of RT longitude nodes |
| `model.num_levels_TOMO_RT` | Number of RT altitude layers |
| `model.res` | Horizontal resolution [deg] |
| `model.radii` | Earth semi-axes [km]: `[a, b]` |
| `model.rWGS` | Mean Earth radius at each RT node [km] |
| `model.GRIDboundaries` | `[lon1 lon2 res lat1 lat2 res]` |
| `model.BLh_pudel` | Square coordinates of TOMO model nodes |
| `model.BLh_pudel_num` | Square coordinates of TOMO model voxel centres |
| `model.BLh_pudel_rad` | Geographical coordinates of TOMO model nodes [rad] |
| `model.BLh_pudel_num_rad` | Geographical coordinates of TOMO model voxel centres [rad] |
| `model.BLh_outer` | Square coordinates of outer TOMO model nodes |
| `model.BLh_outer_num` | Square coordinates of outer TOMO model voxel centres |
| `model.BLh_outer_rad` | Geographical coordinates of outer TOMO nodes [rad] |
| `model.BLh_outer_num_rad` | Geographical coordinates of outer TOMO voxel centres [rad] |
| `model.mid_levels_TOMO` | Altitudes of TOMO voxel centres [m] |
| `model.mid_lat_TOMO` | Latitudes of TOMO voxel centres [deg] |
| `model.mid_lon_TOMO` | Longitudes of TOMO voxel centres [deg] |
| `model.LAT` | Latitude matrix for TOMO nodes |
| `model.LON` | Longitude matrix for TOMO nodes |
| `model.num_inner` | Voxel indices inside the inner model |
| `model.num_outer` | Voxel indices inside the outer model |

### Atmospheric fields (cell arrays indexed by epoch)

| Field | Description |
|-------|-------------|
| `model.temp{epoch}` | Temperature at RT grid nodes [K] |
| `model.pres{epoch}` | Pressure at RT grid nodes [hPa] |
| `model.wvpr{epoch}` | Specific humidity at RT grid nodes [g/kg] |
| `model.refr{epoch}` | Wet refractivity at RT grid nodes [N-units] |
| `model.refrRT{epoch}` | Total refractivity at RT grid nodes [N-units] |
| `model.temp_num_apr{epoch}` | Temperature at RT voxel centres |
| `model.pres_num_apr{epoch}` | Pressure at RT voxel centres |
| `model.wvpr_num_apr{epoch}` | Specific humidity at RT voxel centres |
| `model.refr_num_apr{epoch}` | Wet refractivity at RT voxel centres |
| `model.refrRT_num{epoch}` | Total refractivity at RT voxel centres |
| `model.refr_apr{epoch}` | Wet refractivity at TOMO nodes/voxel centres |
| `model.refr_aprF{epoch}` | Total refractivity at TOMO nodes/voxel centres |
| `model.wvpr_apr{epoch}` | Specific humidity at TOMO nodes/voxel centres |
| `model.temp_apr{epoch}` | Temperature at TOMO nodes/voxel centres |
| `model.pres_apr{epoch}` | Pressure at TOMO nodes/voxel centres |

### Station metadata (added in Stage 6)

| Field | Description |
|-------|-------------|
| `model.BLh` | `[n_sta × 4]` — `[id lat lon h]` after filtering |
| `model.BLH` | `[n_sta × 4]` — `[id lat lon H]` (with orthometric height) |
| `model.BLh_ori` | Full `[id lat lon H h]` from `stations.txt` |
| `model.NAME` | Cell array of station name strings |
| `model.X`, `model.Y`, `model.Z` | ECEF coordinates of retained stations [km] |
| `model.cut_off_angle` | Minimum satellite elevation angle [deg] |
| `model.ds` | Ray step size [km] |

### Bilinear-h mode extras

| Field | Description |
|-------|-------------|
| `model.Xbil` | Latitude mesh for bilinear interpolation |
| `model.Ybil` | Longitude mesh for bilinear interpolation |
| `model.bilmeshLatLon` | Reshaped lat/lon coordinate matrix |

---

## `station(t)`

Created by `construct_station_LAB.m`. A MATLAB structure array of length `n_epochs`. Each element `station(t)` has:

### Ground GNSS sub-structure `station(t).h(nr)`

`nr` runs from 1 to the number of retained GNSS stations.

| Field | Description |
|-------|-------------|
| `.nazwa` | 4-character station name string |
| `.parameters` | `[id lat lon h H]` |
| `.ZTD` | ZTD, gradients, and conversion factor |

#### Satellite sub-structure `station(t).h(nr).satellite(sat)`

`sat` runs from 1 to 32 (GPS satellites). Fields are empty when the satellite is not visible.

| Field | Description |
|-------|-------------|
| `.PRN` | GPS PRN number |
| `.elevation` | Satellite elevation angle [deg] |
| `.SWD` | Slant Wet Delay [mm] |
| `.M_SWD` | SWD measurement error [mm] |
| `.STD` | Slant Total Delay [mm] |
| `.M_STD` | STD measurement error [mm] |
| `.SIVW` | Slant Integrated Water Vapour [kg/m²] |
| `.M_SIVW` | SIWV error |
| `.vmf1w` | Vienna Mapping Function — wet component |
| `.vmf1h` | Vienna Mapping Function — hydrostatic component |
| `.coord` | Satellite ECEF position `[X Y Z]` [km] |

### RO sub-structure `station(t).ro`

| Field | Description |
|-------|-------------|
| `.coordT` | Transmitter (GPS) ECEF coordinates [m] |
| `.coordR` | Receiver (LEO) ECEF coordinates [m] |
| `.date` | Occultation start date `[yyyy mm dd hh mm ss]` |
| `.exL2` | L2 excess phase [m] |
| `.exLC` | Ionosphere-free excess phase LC [m] |

---

## `values`

Holds the apriori refractivity data mapped to voxel index arrays. Field names have suffix `_num` for `'constant'` parameterisation and no suffix for `'bilinear-h'`.

| Field (`constant`) | Field (`bilinear-h`) | Description |
|--------------------|---------------------|-------------|
| `Nw_apr` | `Nw_apr` | `[n_epochs × n_vox]` apriori wet refractivity |
| `Nw_out_num` | `Nw_out` | Apriori wet refractivity for outer model voxels |
| `num_Nw_num(t).number` | `num_Nw(t).number` | Voxel indices for inner model, per epoch |
| `Nw_obs_num(t,:)` | `Nw_obs(t,:)` | Apriori values at inner model voxels, per epoch |
| `num_Nw_out_num(t).number` | `num_Nw_out(t).number` | Voxel indices for outer model, per epoch |
| `Nw_obs_out_num(t,:)` | `Nw_obs_out(t,:)` | Apriori values at outer model voxels, per epoch |
| `WV_apr` | `WV_apr` | Apriori water vapour [g/m³] |
| `rem_REC` | `rem_REC` | Receiver IDs to exclude from processing |
| `rem_SAT` | `rem_SAT` | Satellite PRNs to exclude from processing |

---

## `switches`

See [03-configuration.md](03-configuration.md) for the full list of `switches` fields and valid values. During `intomolab.m` the raw cell/string fields are translated into boolean runtime flags:

| Flag | Source | Meaning |
|------|--------|---------|
| `switches.run_SWD` | `observations{1,1} == 'SWD'` | Process slant wet delays |
| `switches.run_INNER` | `apriori{1,1} == 'INNER'` | Use inner model apriori |
| `switches.run_OUTER` | `apriori{1,2} == 'OUTER'` | Use outer model apriori |
| `switches.run_HORIZONTAL` | `constraints{1,1} == 'HORIZONTAL'` | Apply horizontal smoothness |
| `switches.run_VERTICAL` | `constraints{1,2} == 'VERTICAL'` | Apply vertical smoothness |
| `switches.run_STACKING` | `stacking{1,1}` | Number of epochs to stack |
| `switches.run_KALMAN` | `filter{1,1}`: `0`=KALMAN, `1`=ROBUST | Kalman mode |
| `switches.run_DECOR` | `decorelation{1,1} == 'YES'` | Decorrelation mode |

---

## `obs_set`

Created by `time_listing.m`. All matrices have the format `[JD DOY YEAR GPS_WEEK GPS_SECONDS DOW MM DD HH MM flag]`.

| Field | Description |
|-------|-------------|
| `obs_set.observation_set` | Tomography processing epochs |
| `obs_set.observation_set_SP3` | SP3 file epochs (with file-boundary flags in column 11) |
| `obs_set.interpolation_set` | METEO interpolation epochs |

Column mapping (1-indexed):

| Col | Content |
|-----|---------|
| 1 | Julian Date |
| 2 | Day of Year |
| 3 | Year |
| 4 | GPS Week |
| 5 | GPS seconds since epoch |
| 6 | Day of Week (0=Sunday) |
| 7 | Month |
| 8 | Day of Month |
| 9 | Hour (UTC) |
| 10 | Minute |
| 11 | File boundary flag (1 at daily boundaries) |

---

## `ray`

Computed inside `voxel_dist_3D_combined.m` and passed through `groundRT`/`spaceRT`. **Not saved** to any output file by default.

| Field | Description |
|-------|-------------|
| `ray.de_lat_ray_fin` | Bent ray latitude coordinates for each ray point [rad] |
| `ray.de_lon_ray_fin` | Bent ray longitude coordinates for each ray point [rad] |
| `ray.de_alt_ray_fin` | Bent ray altitude coordinates for each ray point [km] |
| `ray.de_lat_ray_fins` | Straight ray latitude coordinates [rad] |
| `ray.de_lon_ray_fins` | Straight ray longitude coordinates [rad] |
| `ray.de_alt_ray_fins` | Straight ray altitude coordinates [km] |
| `ray.refr` | Refraction index at each ray point |
| `ray.refh` | Hydrostatic refraction index at each ray point |
| `ray.refw` | Wet refraction index at each ray point |
| `ray.X_ray_b` | Bent ray Cartesian coordinates [km] |
| `ray.X_ray` | Straight ray Cartesian coordinates [km] |
| `ray.dL` | Final phase delay [m] |
| `ray.step` | Straight ray segment length [km] |
| `ray.stepb` | Bent ray segment length [km] |
| `ray.nstepb` | Bent ray segment length × refractivity |
| `ray.i_pos` | Indices of voxels containing ray points |
| `ray.vox` | Meteorological values at ray points |
| `ray.voxEM` | Meteorological values in voxels traversed but without ray points |
| `ray.voxIN` | Meteorological values in voxels with ray points |
| `ray.grad_n` | Cartesian gradient at each ray point (below model top) |
| `ray.grad_ng` | Spherical gradient at each ray point |
| `ray.t` | Tangent vector at each ray point |
| `ray.dpa` | Vertical gradient at each ray point |
| `ray.h` | Second gradient derivative |
| `ray.g` | Vector update for ray position |
| `ray.de_b` | Correction of initial tangent vector |
| `ray.diff_dist` | Distance from last ray point to receiver |
| `ray.toc` | Wall-clock time for one ray iteration |

---

## `thin`

Produced per epoch by `intomolab.m` and saved to `KAL/kalman_mat_<epoch>.mat`.

| Field | Description |
|-------|-------------|
| `thin.Q` | Diagonal of process noise covariance Q |
| `thin.P2` | Observation weights (diagonal or full matrix) |
| `thin.R_SWD_k` | Stacked observation error variances |
| `thin.Pplus` | Posterior state covariance diagonal |
| `thin.Pminus` | Prior state covariance diagonal |
| `thin.id_del` | Row indices removed by `filterOBSRT` |
| `thin.theta` | Log-eigenvalues from SVD (first epoch only) |
| `thin.theta_ind` | Sorted `theta` values |

---

## `output`

Returned by `intomolab.m` and saved to `OUT/OUTPUT_TOMO.mat`.

| Field | Size | Description |
|-------|------|-------------|
| `output.xP` | `[n_epochs × n_vox]` | Estimated refractivity [N-units] |
| `output.mxP` | `[n_epochs × n_vox]` | Estimation uncertainty (std dev) [N-units] |
| `output.n_obs` | `[1 × n_epochs]` | Number of observations per epoch |
| `output.n_c` | `[1 × n_epochs]` | Number of apriori constraint rows per epoch |

For the LSQ method `output` also contains `output.dxP` (difference from apriori).
