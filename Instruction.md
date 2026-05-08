# INTOMO — Instruction

## 1. General Info

This is the Integrated Tomography Model (INTOMO) v.1.0 library. INTOMO is an integrated GNSS tomography tool based on radio occultation and ground-based GNSS observations. It was developed as part of the PhD dissertation by Adam Cegla, a PhD candidate at Wroclaw University of Environmental and Life Sciences, Poland, and funded by the NAWA OPUS project entitled “3D Integrated Sensing of the Troposphere Using Ground and Space-Based GNSS Observations” (2020/37/B/ST10/03703).

Please note that this is the first experimental version of INTOMO v1.0, and some processing options may not work properly. In case of any issues, don't hesitate to contact: adam.cegla@upwr.edu.pl.

The library is written for MATLAB 2018a or newer. It requires the NAIF CSPICE toolkit (MATLAB/MEX), the MATLAB Aerospace Toolbox and the MATLAB Mapping Toolbox.

## 2. References

If you use this library please cite the work below:

Cegla, Adam, et al. "Application of integrated GNSS tomography in observation study over the area of southern Poland." Advances in Space Research 74.8 (2024): 3654-3667.

If you would like to know more about the processing methodology please refer to:

Cegla, Adam, et al. "INTOMO operator for GNSS multi-source tomography based on 3D ray tracing technique." Journal of Geodesy 98.11 (2024): 101.

or

Cegla, Adam, et al. "Application of integrated GNSS tomography in observation study over the area of southern Poland." Advances in Space Research 74.8 (2024): 3654-3667.

It is a constantly expanding project, and new functionalities are added as you read that text.

## 3. INTOMO processing options

Current functionalities of the INTOMO library include:

- Processing GNSS tomography in ground-based only and integrated tomography: `switches.integrated` (`'yes'`/`'no'`)
- Processing GNSS tomography with voxel-based or node-based parameterization: `switches.parametrization` (`'constant'`/`'bilinear-h'`)
- Processing real or synthetic data: `switches.solution` (`'real'`/`'synthetic'`)
- Processing with the Kalman and Least Squares methods: `switches.method` (`'kalman'`/`'lsq'`)
- Processing in total or wet-only refractivity mode: `switches.totalN` (`'true'`/`'false'`)
- A priori model source: `aprModel` (`'DETER'`/`'ERA5'`)

These settings can be configured in the `conf.m` file. Other settings available in the `conf.m` file are planned and already included in the setting file, however were not tested yet. For details please refer to `INTOMO_variables.md` and structure `switches`.

## 4. Getting started

Before running the `RUNINTOMO.m` script, please set the appropriate variables in the `conf.m` script.

- Set the path to the INTOMO function as `PATH_INSTALL` and the path to save the results as `PATH_EXTERNALSAVE`.
- Define `PROJECT_NAME` as the name of the folders containing the initial tomography data and `Tomography/DATA/project_name`
- Set the dates for processing in the format `[yyyy mm dd hh mm ss]`.
- Define the tomography domain settings. Since INTOMO uses two different grid domains for ray tracing and GNSS processing, it is advised to design the ray tracing grid to cover a larger area than the one used for tomography processing. The horizontal and vertical resolution depends on the Numerical Weather Model used (currently, only ERA5 has been tested).
- Select the settings from the "INTOMO processing options".
- Create following folders:
  - `Tomography/DATA/project_name/ATM` — contains GNSS ground-based data
  - `Tomography/DATA/project_name/METEO` — contains NWM data
  - `Tomography/DATA/project_name/ORB` — contains SP3 data
  - `Tomography/DATA/project_name/RO` — contains RO data (Level 1b, COSMIC)
  - `Tomography/DATA/project_name/WORK` — work directory to save model data and observation matrices A
- Prepare GNSS station coordinates in a proper format (more details in point 5).
- Prepare the data and copy them to the right folders (more details in point 5).
- Launch processing.

Make sure all these settings are configured properly in the `conf.m` file before running the script.

## 5. Data reading

The following a priori data is needed for the processing:

- RO data in the form of Phs NetCDF files (Level1b processing from sources such as UCAR), moved to the `Tomography/DATA/project_name/RO` folder.
- Positions of ground GNSS receivers in the form of `.txt` files, moved to the `Tomography/CONF` folder.
- Note that the current version of INTOMO does not provide an algorithm to create a `.txt` file in the proper format; a template `.txt` file is included in the library.
- Ground GNSS receiver observations in the form of `.txt` files, moved to the `Tomography/DATA/project_name/ATM` folder.
- As with the receiver positions, INTOMO does not currently provide an algorithm to create these `.txt` files, but a template is included in the library.
- NetCDF file with NWM data, moved to the `Tomography/DATA/project_name/METEO` folder. This version has only been tested with ERA5 data.
- SP3 data in the form of `.sp3` files, moved to the `tomography/DATA/project_name/ORB` folder.

## 6. Processing stages

1. Preparation of data for tomography processing
   1. reading GNSS stations coordinates,
   2. calculating a priori refractivity fields in ray tracing and tomography domain,
   3. reading `.sp3` files and searching for available ground-based observations
   4. reading zenith delays and excess phases with LEO satellites coordinates,
   5. calculating gradients and slant delays,
   6. preparing data structures,
2. Ray Tracing of RO,
3. Ray Tracing of GNSS signals,
4. Points 2 and 3 iterated by every epoch,
5. Preparing apriori variance-covariance matrices,
6. Building stacked matrices _k (see INTOMO_variables file),
7. Calculating the weights and filtering of the observations,
8. Kalman filtering
9. Points 4 and 8 iterated by every epoch,

## 7. Processing outputs

The algorithm saves the following output matrices:

- `model.mat` in the `Tomography/DATA/project_name/WORK` folder, containing the model's parameters.
- `station.mat` in the `Tomography/DATA/project_name/WORK` folder, containing a structure with observations from ground and space-based sources.
- `amtrix.mat` in the `Tomography/DATA/project_name/WORK/obsmatrices` folder, containing observational matrices from ray tracing.
- `OUTPUT_TOMO.mat` in the `save_folder/project_name/OUT` folder, containing refractivity estimates.

For more details on the variables please refer to `INTOMO_variables.md`.

## 8. Matlab functions used

The functions were described by assigning them into processing stages from point 6. Detailed description can be found in each of the `.m` files. Only more important functions were described.

### ad.1.

- `construct_station_LAB.m` — convert the specific ZTD data into the station structural matrix
- `ERA5GridTOMO.m` — generate TOMO and RT models
- `gridcalcRT.m` — calculate values of refractivity in models nodes/centers
- `boundingTOMOLAB.m` — cut off GNSS stations outside tomography model
- `deleteStat.m` — filter out station based in case of too dense network
- `pBLh2ZHD.m` — calculate ZHD in GNSS stations locations
- `pudel2.m` — create 2D matrix of coordinates of TOMO models nodes/voxel centers
- `Undulation.m` — create undulation grid for RT
- `time_listing.m` — create obs_set variable with processing epochs
- `interSP3.m` — find satellite coordinates matching the GNSS stations
- `downloadORB.m` — download SP3 data
- `findRO.m` — read RO data (excess phase and LEO-GNSS orbits)
- `RUNINTOMO.m` — define input variables for tomography processing
- `readBLh.m` — read coordinates of GNSS stations

### ad.2.

- `A_row_bilin_regRT.m` — calculate the deriveratives of spline/bilinear functions
- `groundRT.m` — calculate ray tracing and derivatives in A matrix for ground-based obs
- `interateModelsA.m` — find the distances traversed by ray path inside TOMO model
- `bilinearRT.m` — Bilinear/spline interpolation function
- `coordVector.m` — calculates local vectors of ray paths segments in voxels
- `perler_T.m` — calculate derivatives for spline interpolation
- `screenZTD.m` — filter out ZTD data
- `ECEF2LLA.m` — recalculate coordinates from ECEF to LLA

### ad.3.

- `A_row_bilin_regRT.m` — calculate the deriveratives of spline/bilinear functions
- `interateModelsA.m` — find the distances traversed by ray path inside TOMO model
- `bilinearRT.m` — Bilinear/spline interpolation function according
- `spaceRT.m` — calculate ray tracing and derivatives in A matrix for space-based obs
- `coordVector.m` — calculates local vectors of ray paths segments in voxels
- `perler_T.m` — calculate derivatives for spline interpolation
- `ECEF2LLA.m` — recalculate coordinates from ECEF to LLA

### ad.4.

- `initialOBS.m` — prepare and stack A matrices from space and ground - based obs

### ad.5.

- `covarianceAprRT.m` — get first estimate of dynamic disturbance of refractivity
- `covarianceNwRT.m` — get variance-covariance of dynamic disturbance of refractivity
- `apConstRT.m` — setup a priori information in form of a matrices

### ad.6.

- `matrices_epochRT.m` — stack matrices for Kalman filtering

### ad.7

- `filterOBSRT.m` — filter our observations exceeding defined threshold
- `weightObs.m` — calculate weights of the observations
- `stackedNwcalc.m` — find zenith delay values for RO ray path-voxel boundary inters.

### ad.8.

- `getGain.m` — get estimation of Kalman gain matrix
- `getGainR.m` — get robust estimation of Kalman gain matrix
- `svd.m` / `svdecon.m` — matrix Singular Value Decomposition

### ad.9.

- `intomolab.m` — calculate GNSS tomo estimates of refractivities with Kalman filtering
- `TomoLSQ.m` — calculate GNSS tomo estimates of refractivities with least-square method

In case of ray tracing algorithm, following functions are available:

- `select_id.m` — select id of the next voxel in ray path direction
- `savevar.m` — save variables to structure
- `dhrefcalc.m` — Function for voxel horizontal gradient computation
- `dp_alt.m` — Function to obtain mean voxel difference of refraction value and height
- `dpacalc.m` — Voxel vertical gradient computation and it's interpolation to the ray point
- `excessphase.m` — Calculate excess phase based on ray path points coordinates.
- `find_num.m` — Find ID of the voxel where ray point is located
- `gradrec.m` — transform gradients between spherical and cartesian coordinates
- `IDW_atom.m` — IDW interpolation on ray point coordinates
- `recgrad.m` — transform grad. between cart. and sphere. coordinates
- `refcalc.m` — calculates refractivity
- `ref_interp_3D.m` — interpolate refractivity for selected coordinates
- `statprofile.m` — calculates zenith delays for selected profile
- `vox_distance.m` — calculates distance traversed by ray path in voxel
- `voxel_dist_3D_combined` — reconstructs the 3D ray path through the voxel model
- `vdist` — distance between points on the WGS-84 ellipsoid Earth
