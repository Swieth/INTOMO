# INTOMO — Important variables

## Variable reference

| Variable | Description |
|----------|-------------|
| `A` | Observational matrix with signal derivatives |
| `A_apriori` | The A apriori matrix with 0 or 1 values |
| `A_constr` | Constraints given to A matrix |
| `A_k` | Stacked A matrices |
| `A_RO` | Observational matrix with signal derivatives for RO |
| `Avec` | The ray path vectors of each signal in tomography voxel |
| `azi` | Azimuth of the satellite |
| `BLh_pudel` | Coordinates of voxel edges/centers |
| `coord` | Coord of intersections between RO and voxel model |
| `dSWD` | Difference between observed and simulated slant delay |
| `elev` | Elevation angles of ray traced signals |
| `elevat_k` | Stacked elevation angle vector |
| `epoch` | Processing epoch number |
| `exPh` | The values of excess phase |
| `K` | Kalman gain matrix |
| `levels` | Altitude layers of tomography model |
| `model` | Parameters of tomography and ray tracing models (structure) |
| `num_lat` | Number of voxels in latitudinal direction |
| `num_lev` | Number of voxels in vertical direction |
| `num_lon` | Number of voxels in longitudinal direction |
| `Nw_apr` | Apriori values of refractivities |
| `P_apriori` | The apriori weights for the tomography grid |
| `P_constr` | Constraints given to P |
| `P2` | Weights of observations |
| `path_save` | Path to save the output matrices |
| `Pminus` | Kalman filtering predicted covariance matrix of estimation uncertainties |
| `Q` | The first guess for covariance matrix |
| `R_SIWV` | The errors of Slant Integrated Water Vapour |
| `R_SWD` | The errors of the slant delays |
| `R_SWD_k` | Stacked R_SWD vector |
| `ray` | Parameters of ray traced signal (structure) |
| `RelDist` | Total distance traversed by GNSS signal |
| `rem_REC` | Id of receivers to delete from apriori data |
| `rem_SAT` | Id of satellites to delete from apriori data |
| `sat` | Id of GPS satellite |
| `SAT_k` | Stacked id of GPS vector |
| `satcoord` | Coordinates of the GNSS satellite |
| `SIWV` | The values of Slant Integrated Water Vapour |
| `stat_name_k` | Stacked names of GNSS vector |
| `station` | Structural matrix containing all informations regarding GNSS observations (structure) |
| `station_name` | Names of GNSS stations |
| `SWD` | The values of Slant Delay (SD) or excess phase (EP) |
| `SWD_apriori` | The a priori values of Slant Delay (SD) |
| `SWD_constr` | Constraints given to SD |
| `SWD_k` | Stacked SWD vector |
| `SWD_obs` | The values of slant delays |
| `switches` | Structural matrix of settings of the INTOMO processing (structure) |
| `thin` | Kalman filtering and selected tomgoraphy parameters (structure) |
| `values` | Structural matrix containing inner/outer model data |
| `xPminus` | Kalman filtering predicted covariance matrix of estimation uncertainties |

## Important structures

### Ray

- `ray.de_lat_ray_fin` — bent ray latitude coordinates for each ray point (RP) [rad]
- `ray.de_lon_ray_fin` — bent ray longitude coordinates for each RP [rad]
- `ray.de_alt_ray_fin` — bent ray altitude coordinates for each RP [km]
- `ray.de_lat_ray_fins` — straight ray latitude coordinates for each RP [rad]
- `ray.de_lon_ray_fins` — straight ray longitude coordinates for each RP [rad]
- `ray.de_alt_ray_fins` — straight ray altitude coordinates for each RP [km]
- `ray.refr` — refraction index for each RP
- `ray.refh` — hydrostatic refraction index for each RP
- `ray.refw` — wet refraction index for each RP
- `ray.X_ray_b` — bent ray cartesian coordinates for each RP
- `ray.X_ray` — straight ray cartesian coordinates for each RP
- `ray.diff_dist` — last ray point distance to receiver
- `ray.grad_n` — cart. grad. in each RP below upper boundary of the model
- `ray.grad_ng` — spherical grad. in each RP below upper boundary of the model
- `ray.t` — tangent vector in each RP below upper boundary of the model
- `ray.dpa` — vert. grad. in each RP below upper boundary of the model
- `ray.toc` — time of one iteration calculation
- `ray.h` — 2nd grad. deriv. in each RP below upper boundary of the model
- `ray.g` — vector update for ray position
- `ray.de_b` — correction of initial tangent vector
- `ray.dL` — final phase delay [m]
- `ray.nstepb` — bent ray segment length multiplied by refractivity
- `ray.step` — straight ray segment length
- `ray.stepb` — bent ray segment length
- `ray.i_pos` — indexes of voxels containing ray points
- `ray.vox` — selected meteorological values at ray points positions
- `ray.voxEM` — structure for voxels parameters without raypoints but traversed by raypath (e.g. pressure, temperature, water vapor)
- `ray.voxIN` — structure for voxels parameters with raypoints (e.g. pressure, temperature, water vapor)

**Warning:** The structure is not saved as an output variable in default settings!

### values

- `Nw_apr` — apriori wet refractivity
- `Nw_out` — apriori wet refractivity for outer model
- `num_Nw_out` — id number of wet refractivity for outer model voxel
- `Nw_obs_out` — additional wet refractivity for outer model voxel
- `num_Nw` — id number of wet refractivity for inner model voxel
- `Nw_obs` — additional wet refractivity for inner model voxel
- `WV_apr` — apriori water vapour
- `WV_out` — apriori water vapour for outer model
- `WV_obs_out` — additional water vapour for outer model voxel
- `num_WV` — id number of water vapour for inner model voxel
- `WV_obs` — additional water vapour for inner model voxel

### Switches

- `observations{'SWD'}` — type of observations
- `priori{'INNER','OUTER'}` — apriori model
- `constraints{'HORIZONTAL','VERTICAL'}` — apriori model constraints (%not tested)
- `stacking{'NO'/number}` — %not tested
- `filter{'KALMAN'/'ROBUST'}` — Kalman filtering type
- `solution{'REAL'/'SYNTHETIC'}` — solution type
- `decorelation{'NO'/'YES'}` — decorelation type
- `aprModel{'DETER'/'ERA5'}` — a priori model source
- `amtrix{folder_name}` — name for work folder to contain observation matrices
- `parameterization{'constant'}` — type of parameterization
- `phi{'identity'}` — include identity matrix
- `method{'KALMAN'/'LSQ'}` — Kalman or LSQ processing
- `totalN{true/false}` — processing in wet refractivity or (total) refractivity
- `saveAtmParam{true/false}` — save atmospheric apriori param (e.g. pressure) to tomography output
- `coord{'Formatted'}` — input format for coordinates of GNSS stations
- `ROres{(resolution of RO processing [Hz])}` — e.g. 1 for 50Hz, 10 for 5Hz
- `refron{true/false}` — ray tracing from apriori refractivities (if available)/from atmospheric parameters
- `integrated{'yes'/'no'}` — integrated tomography processing

### station(t)

`t` — number of epochs

- `h(nr)` — `nr` — number of stations
  - `ZTD` — initial ZTD, gradients and conversion factor if present
  - `parameters` — id and coordinates of GNSS receiver [id lat lon h H]
  - `name` — station 4 letters abbrev
  - `satellite(sat)` — `sat` — satellite number currently 1:32
    - `PRN` — PRN number
    - `elevation` — satellite elevation [deg]
    - `SWD` — slant wet delay
    - `M_SWD` — SWD error
    - `STD` — slant total delay
    - `M_STD` — STD error
    - `SIVW` — slant integrated water vapour
    - `M_SIVW` — SIWV error
    - `vmf1w` — Vienna Mapping Function wet
    - `vmf1h` — Vienna Mapping Function hydrostatic
    - `coord` — satellite position (ECEF)
- `ro` — radio occultation (RO) data
    - `coordT` — coordinates of the transmitter (ECEF)
    - `coordR` — coordinates of the receiver (ECEF)
    - `date` — date of RO start
    - `exL2` — RO excess phase on L2
    - `exLC` — RO excess phase on LC

### output

- `xP` — estimated values of refractivities
- `mxP` — errors of refractivities estimation

### obs_set

- `obervations_set` — tomography processing epochs
- `obervations_set_SP3` — SP3 processing epochs
- `interpolation_set` — ray tracing processing epochs

All above in the format: `[time_since:DOY:YEAR:GPS_WEEK:seconds_since:DOW:MM:DD:HH:MM]`

### model

- `temp` — temperature in raytracing (RT) model nodes
- `pres` — pressure in RT model nodes
- `wvpr` — specific humidity in RT model nodes
- `refr` — wet refractivity in RT model nodes
- `refrRT` — total refractivity in RT model nodes
- `temp_num_apr` — temperature in RT voxel centers
- `pres_num_apr` — pressure refractivity in RT voxel centers
- `wvpr_num_apr` — specific humidity refractivity in RT voxel centers
- `refr_num_apr` — wet refractivity in RT voxel centers
- `refrRT_num` — total refractivity in RT voxel centers
- `refr_apr` — wet refr. in tomography (TOMO) model nodes\voxel centers
- `refr_aprF` — total refractivity in TOMO model nodes\voxel centers
- `wvpr_apr` — specific humidity in TOMOmodel nodes\voxel centers
- `temp_apr` — temperature in TOMO model nodes\voxel centers
- `pres_apr` — pressure in TOMO model nodes\voxel centers
- `BLh_pudel` — square coordinates of TOMO model nodes
- `BLh_pudel_num` — square coordinates of TOMO model voxel centers
- `BLh_pudel_rad` — geographical coordinates of TOMO model nodes
- `BLh_pudel_num_rad` — geographical coordinates of TOMO model voxel centers
- `BLh_outer` — square coordinates TOMO model nodes of outer model
- `BLh_outer_num` — square coord. of TOMO model voxel centers of outer model
- `BLh_outer_rad` — geographical coordinates of TOMO nodes of outer model
- `BLh_outer_num_rad` — geo.. coord. of TOMO model voxel centers of outer model
- `levels_mid_TOMO` — altitudes of TOMO model voxel centers
- `num_inner` — indices of voxels in TOMO inner model
- `num_outer` — indices of voxels in TOMO outer model
- `rWGS` — mean Earth radius in RT model nodes
- `LAT` — latitude matrix for TOMO model nodes
- `LON` — longitude matrix for TOMO model nodes
- `mid_levels_TOMO` — altitudes of TOMO model voxel centers
- `mid_lat_TOMO` — latitudes of TOMO model voxel centers
- `mid_lon_TOMO` — longitudes of tomography model voxel centers

### Thin

- `Q` — the first guess for covariance matrix
- `P2` — weights of observations
- `R_SWD_k` — stacked errors of the slant delays and tomography model
- `id_del` — observations deleted from the processing (see Subsec. 1…)
- `theta` — log of eigenvalues from singular value decomposition
- `theta_ind` — sorted theta values

## Paths

### Path to be set manually

- `PATH_INSTALL` — installation folder
- `PATH_EXTERNALSAVE` — tomography output save folder

### Paths set automatically

- `pathRO` — path to folder with RO.nc files
- `pathTOMO` — path for WORK directory
- `pathORB` — path to SP3 data
- `pathCONF` — path to configuration folder
- `pathMETEO` — path to meteorological data
- `pathATM` — path to ground-based data
- `pathEXPORT` — tomography output save folder

## Others

**Amtrix** variable saved as tomography output contains data for each epoch of the processing with both ray tracing data and slant delays/excess phase values. The variables included are: `A`, `Avec`, `coord`, `dSWD`, `elev`, `not_cro` (indexes of voxels not crossed by any ray path), `R_SIWV`, `R_SWD`, `RelDist`, `SAT`, `SIWV`, `station_name`, `SWD` and `SWD_nodes_integ` (integral value of SWD from apriori data)
