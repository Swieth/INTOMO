# INTOMO v1.0 - Integrated GNSS Tomography Library

## General Information

**INTOMO** (Integrated Tomography Model) is a MATLAB-based GNSS tomography tool based on Radon transform theory for the troposphere monitoring using both ground-based and space-based GNSS observations. 

This software was developed as part of the PhD dissertation by **Adam Cegla** at **Wroclaw University of Environmental and Life Sciences**, Poland. It is part of the NAWA OPUS project titled:

> *"3D Integrated Sensing of the Troposphere Using Ground and Space-Based GNSS Observations"*  
> (Project No. 2020/37/B/ST10/03703)

> **Note:** This is an experimental release (v1.0). Some features may not work as expected. For issues, contact: [adam.cegla@upwr.edu.pl](mailto:adam.cegla@upwr.edu.pl)

### Requirements

- MATLAB 2018a or newer
- NAIF CSPICE (MATLAB/MEX interface)
- MATLAB Aerospace Toolbox
- MATLAB Mapping Toolbox
- MATLAB Statistics and Machine Learning Toolbox

---

## References

Please cite the following if using this library:

- Cegla, Adam, et al. *"Application of integrated GNSS tomography in observation study over the area of southern Poland."* **Advances in Space Research**, 74(8), 2024: 3654-3667.

For methodology:

- Cegla, Adam, et al. *"INTOMO operator for GNSS multi-source tomography based on 3D ray tracing technique."* **Journal of Geodesy**, 98(11), 2024: 101.

---

## INTOMO Processing Options

You can configure processing options using the `conf.m` file:

| Option | Values | Description |
|--------|--------|-------------|
| `switches.integrated` | `'yes'` / `'no'` | Use integrated or ground-based only tomography |
| `switches.parametrization` | `'constant'` / `'bilinear-h'` | Type of parameterization |
| `switches.solution` | `'real'` / `'synthetic'` | Type of input data |
| `switches.method` | `'kalman'` / `'lsq'` | Filtering method |
| `switches.totalN` | `'true'` / `'false'` | Total or wet-only refractivity mode |
| `aprModel` | `'DETER'` / `'ERA5'` | A priori model source |

Additional options are included but untested. See `INTOMO_variables.pdf` and the `switches` structure for details.

---

## Getting Started

Before running `RUNINTOMO.m`, ensure `conf.m` is set up:

1. Set `PATH_INSTALL` and `PATH_EXTERNALSAVE`.
2. Define `PROJECT_NAME` and date/time settings.
3. Configure tomography and ray tracing domain settings.
4. Set processing options (see above).
5. Create the following folder structure:
Tomography/
      └── DATA/
      └── project_name/
      ├── ATM # Ground GNSS data
      ├── METEO # Numerical Weather Model data (e.g., ERA5)
      ├── ORB # SP3 data
      ├── RO # Radio Occultation data (Level 1b, COSMIC)
      └── WORK # Workspace for intermediate files

6. Prepare GNSS station coordinates and data in `.txt` format (templates included).
7. Place all data files in the correct folders.
8. Run processing.

---

## Required Input Data

- **RO data**: NetCDF Level1b files in `RO/`
- **GNSS receiver positions**: `.txt` in `CONF/`
- **GNSS observations**: `.txt` in `ATM/`
- **NWM data**: NetCDF files in `METEO/`
- **SP3 files**: `.sp3` format in `ORB/`

> ⚠️ INTOMO does not currently generate `.txt` GNSS files automatically; templates are included.

---

## Processing Stages

1. **Data Preparation**
   - Read GNSS stations and a priori data
   - Load SP3 and RO observations
   - Compute refractivity and slant delays
2. **RO Ray Tracing**
3. **GNSS Ray Tracing**
4. **Epoch-wise Iteration**
5. **A priori V-Cov Matrices**
6. **Stack Matrices**
7. **Observation Filtering & Weighting**
8. **Kalman Filtering**
9. **Iteration of Processing and Filtering**

---

## Output Files

Saved to: `Tomography/DATA/project_name/WORK`

| File | Description |
|------|-------------|
| `model.mat` | Model parameters |
| `station.mat` | Observational structures |
| `amtrix.mat` | Ray tracing observation matrices |
| `OUTPUT_TOMO.mat` | Final tomography output in `OUT/` folder |

---

## MATLAB Functions by Stage

Only key functions are listed here. See the source code for details.

### Data Preparation
- `construct_station_LAB.m`, `ERA5GridTOMO.m`, `gridcalcRT.m`, ...
- `readBLh.m`, `findRO.m`, `interSP3.m`, ...

### RO Ray Tracing
- `spaceRT.m`, `A_row_bilin_regRT.m`, `coordVector.m`, ...

### GNSS Ray Tracing
- `groundRT.m`, `bilinearRT.m`, `interateModelsA.m`, ...

### Other Stages
- Kalman: `getGain.m`, `intomolab.m`
- LSQ: `TomoLSQ.m`
- Filtering: `filterOBSRT.m`, `weightObs.m`
- Covariance: `covarianceAprRT.m`, `covarianceNwRT.m`

### Ray Tracing Utilities
- `select_id.m`, `vox_distance.m`, `refcalc.m`, `excessphase.m`, ...

For the full list, refer to the comments in each `.m` file.

---

## Contact

**Adam Cegla**  
Email: [adam.cegla@upwr.edu.pl](mailto:adam.cegla@upwr.edu.pl)  
Affiliation: Wroclaw University of Environmental and Life Sciences

---

## Acknowledgment

This research was supported by the NAWA OPUS project titled:

> *"3D Integrated Sensing of the Troposphere Using Ground and Space-Based GNSS Observations"* (2020/37/B/ST10/03703)

