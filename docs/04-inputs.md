# Input Files

INTOMO requires five categories of input data. Each category is stored in a designated sub-folder under `Tomography/DATA/<PROJECT_NAME>/`. One additional static file lives under `Tomography/CONF/`.

---

## 1. GNSS station coordinates

**Location:** `Tomography/CONF/stations.txt`  
**Reader:** [`Tomography/READ/readBLh.m`](../Tomography/READ/readBLh.m)  
**Required when:** Always (both `REAL` and `SYNTHETIC` modes)

### Format

Tab-delimited text file, no header row. Columns:

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `Name` | string | 4-character station abbreviation |
| 2 | `num` | integer | Sequential station number |
| 3 | `lat` | float | Geodetic latitude [decimal degrees, North positive] |
| 4 | `lon` | float | Geodetic longitude [decimal degrees, East positive] |
| 5 | `h` | float | Ellipsoidal height [m] |
| 6 | `H` | float | Orthometric height (above geoid) [m] |

### Example rows

```
AAER    1    49.7824740000000    4.6426460000000    197.9280000000000    151.5270000000000
BOR1    8    52.2769600000000    17.0734630000000   124.3520000000000    89.0050000000000
```

### Notes

- The reader uses `formatSpec = '%s%f%f%f%f%f'` with `\t` as the delimiter.
- After reading, stations outside the bounding box defined in `conf.m` are removed by `boundingTOMOLAB.m`.
- Stations closer than `switches.stat_range` km to any other retained station are further thinned by `deleteStat.m`.
- A template `stations.txt` is included in the repository.

---

## 2. Ground GNSS ZTD observations

**Location:** `Tomography/DATA/<PROJECT_NAME>/ATM/<YYYYMMDD>.txt`  
**Reader:** [`Tomography/READ/readtxtOBS.m`](../Tomography/READ/readtxtOBS.m)  
**Required when:** `switches.solution = {'REAL'}`

### Format

Tab-delimited text file with a one-line header. The file name encodes the date as `YYYYMMDD.txt`. One row per station per epoch.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `ID` | string | 4-character station abbreviation |
| 2 | `DOMES` | string | DOMES number (read but not used internally) |
| 3 | `Date` | datetime string | `YYYY-MM-DD HH:MM:SS.s` |
| 4 | `B` | float | Latitude [decimal degrees] |
| 5 | `L` | float | Longitude [decimal degrees] |
| 6 | `Helips` | float | Ellipsoidal height [m] |
| 7 | `Hnorm` | float | Orthometric height [m] |
| 8 | `ZTD` | float | Zenith Total Delay [mm] |
| 9 | `mZTD` | float | ZTD uncertainty [mm] |
| 10 | `GradN` | float | North gradient of ZTD [mm] |
| 11 | `GradE` | float | East gradient of ZTD [mm] |
| 12 | `mGradN` | float | North gradient uncertainty [mm] |
| 13 | `mGradE` | float | East gradient uncertainty [mm] |

### Example rows

```
ID      DOMES           Date                    B           L           Helips      Hnorm       ZTD         mZTD    GradN   GradE   mGradN  mGradE
BART    12234M001       2020-07-08 00:00:00.0   54.25069    20.81544    93.233      65.259      2397.2      1.6     0.05    -0.65   0.23    0.22
```

### Notes

- The reader strips comma-as-thousands-separator from numeric fields (e.g., `2,397.2` → `2397.2`).
- The reader currently supports exactly **one file per day** and assumes **24 hourly epochs**. Processing periods spanning multiple days or using sub-hourly intervals are not supported. See [09-issues.md](09-issues.md).
- Only stations whose names appear in `model.NAME` (i.e., stations that survived the bounding-box and density filters) are retained.
- The Zenith Hydrostatic Delay (ZHD) is computed separately from pressure via `pBLh2ZHDRT.m` and subtracted to obtain the Slant Wet Delay.
- A template observation file is included in `Tomography/DATA/project_name/ATM/`.

---

## 3. SP3 satellite orbits

**Location:** `Tomography/DATA/<PROJECT_NAME>/ORB/`  
**Reader:** [`Tomography/READ/readSP3dat.m`](../Tomography/READ/readSP3dat.m), [`Tomography/READ/downloadORB.m`](../Tomography/READ/downloadORB.m)  
**Required when:** Always

### Format

Standard IGS SP3-c/d format (ECEF coordinates, GPS constellation). File naming convention: `igs<GPS_WEEK><DOW>.sp3` (e.g., `igs23993.sp3`).

INTOMO tries to use the following priority order when looking for a given GPS week/day:

1. `igs<week><dow>.sp3` — final IGS combined
2. `igr<week><dow>.sp3` — rapid IGS
3. `igu<week><dow>_18.sp3` — ultra-rapid (18h prediction)
4. `igu<week><dow>_12.sp3`
5. `igu<week><dow>_06.sp3`
6. `igu<week><dow>_00.sp3`
7. `igu<prev_week><prev_dow>_18.sp3`
8. `igu<prev_week><prev_dow>_12.sp3`

### Automatic download

`downloadORB.m` attempts to download missing SP3 files via anonymous FTP from `igs.bkg.bund.de`. The compressed `.Z` files are fetched and decompressed with `uncompress`. See [09-issues.md](09-issues.md) for portability caveats.

If SP3 files are placed manually they do not need to be downloaded. Files for the processing period (including the day before if the first epoch is before 06:00 UTC) must be present.

### Processing

`readSP3dat.m` reads all `.sp3` files in `pathORB` that cover the processing period. `interSP3.m` then interpolates satellite ECEF coordinates to match the tomography epoch times.

---

## 4. Radio Occultation (RO) data

**Location:** `Tomography/DATA/<PROJECT_NAME>/RO/`  
**Reader:** [`Tomography/READ/findRO.m`](../Tomography/READ/findRO.m)  
**Required when:** `switches.integrated = "yes"`

### Format

NetCDF Level-1b files in UCAR atmPhs format. Two naming patterns are recognised:

- `atmPhs_<sat>.<year>.<doy>.<HH>.<MM>.<GNSS_prn>_0001.0001_nc` — COSMIC-2/UCAR format (no `.nc` extension)
- `<prefix>_atmPhs_*_nc` or `spire_gnss-ro_L1B_atmPhs_*.nc` — SPIRE format

A file is recognised as `atmPhs` if its name contains the string `atmPhs`; otherwise it is treated as a `conPhs` file.

### NetCDF variables read

| Variable | `atmPhs` field | `conPhs` field | Description |
|----------|----------------|----------------|-------------|
| GPS satellite position X | `xGps` | `xGnssLR` | ECEF X [m] |
| GPS satellite position Y | `yGps` | `yGnssLR` | ECEF Y [m] |
| GPS satellite position Z | `zGps` | `zGnssLR` | ECEF Z [m] |
| LEO satellite position X | `xLeo` | `xLeoLR` | ECEF X [m] |
| LEO satellite position Y | `yLeo` | `yLeoLR` | ECEF Y [m] |
| LEO satellite position Z | `zLeo` | `zLeoLR` | ECEF Z [m] |
| L1 excess phase | `exL1` | `exL1` | [m] |
| L2 excess phase | `exL2` | `exL2` | [m] |

### NetCDF global attributes read

| Attribute | Description |
|-----------|-------------|
| `year`, `month`, `day`, `hour`, `minute`, `second` | Occultation start time |
| `occfreq2` | L2 frequency [Hz], optional — defaults to 1227.60 MHz if absent |

### Derived quantity

The ionosphere-free excess phase is computed internally:

```
f1 = 1575.42 MHz
f2 = occfreq2 (or 1227.60 MHz default)
exLC = exL1 + (f2² / (f1² - f2²)) × (exL1 - exL2)
```

If the sum of `exLC` is negative the array is flipped and its absolute value is used.

### Epoch matching

`findRO.m` matches RO files to tomography epochs by comparing the string `YYYYMMDDH` derived from the file's global attributes against the epoch hour string from `obs_set`.

---

## 5. Numerical Weather Model (ERA5)

**Location:** `Tomography/DATA/<PROJECT_NAME>/METEO/`  
**Reader:** [`Tomography/ENGINE/ERA5GridTOMO.m`](../Tomography/ENGINE/ERA5GridTOMO.m)  
**Required when:** `switches.aprModel = 'ERA5'` (not needed for `'DETER'` mode)

### Format

NetCDF files covering the processing domain. File naming convention: `ERA5_<yyyy>-<m>-<d>.nc` (e.g., `ERA5_2026-1-1.nc`).

The file must contain pressure-level data on a regular lat/lon grid matching `model.GRIDboundaries`. The variables expected inside the file are standard ERA5 pressure-level variables: temperature, specific humidity, geopotential (or pressure levels). The exact variable names are resolved inside `ERA5GridTOMO.m`.

ERA5 data can be downloaded from the Copernicus Climate Data Store:  
<https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels>

---

## 6. Geoid undulation

**Location:** `Tomography/CONF/undu.mat`  
**Used by:** [`Tomography/PPROCESS/Undulation.m`](../Tomography/PPROCESS/Undulation.m), called from `ERA5GridTOMO.m`  
**Required when:** Always (if not present, `ERA5GridTOMO.m` generates it from `EGM2008` via `geoidheight2.m`)

The file stores a pre-computed undulation grid for the RT domain. If it does not exist on the first call, `ERA5GridTOMO.m` creates and saves it automatically. On subsequent epochs the cached file is reloaded.

An alternative undulation file `unduera5.mat` is also present in `Tomography/CONF/`; it is not referenced directly in the driver but may be used as the `unduFile` variable if set accordingly in `conf.m`.

---

## Summary table

| Input | Folder | Format | Required |
|-------|--------|--------|---------|
| Station coordinates | `CONF/stations.txt` | Tab-delimited TXT | Always |
| ZTD observations | `ATM/<YYYYMMDD>.txt` | Tab-delimited TXT | `solution=REAL` |
| SP3 orbits | `ORB/*.sp3` | IGS SP3 | Always |
| RO excess phase | `RO/*_nc` or `RO/*.nc` | UCAR/SPIRE NetCDF | `integrated=yes` |
| ERA5 NWM | `METEO/ERA5_*.nc` | NetCDF | `aprModel=ERA5` |
| Undulation grid | `CONF/undu.mat` | MATLAB `.mat` | Auto-generated if absent |
