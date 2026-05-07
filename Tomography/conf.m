%% This is a configuration file for INTOMO.1.0 The major purpose of this file is to set up:
%0.Time mode of the solution (postprocessing)
%1.Data sources locations and associated names for source files, as well as destination folders for exported data,
%2.Time slots and time interval of various data sources,
%3.Ray tracing and tomography model domain settings

%All information will be stored as a 'mat' file in working directory of INTOMO software
%!!!Please, make sure you have all the data in required folders and file names are consistent.
%% The list of minimum required data is as follows:
%I.   Coordinates file in a format consistent with stations.txt file format
%in INTOMO/Tomography/CONF
%II.  Orbits file in a format consistent with IGS ECEF format (SP3) folder: 
%III. RO Level1b UCAR format files in .net (https://cdaac-www.cosmic.ucar.edu/cdaac/doc/formats.html)
%IV. Apriori/Initial conditions for tomography models
%    a) DETER - deterministic models such as gpt/gpt2/UNB3mm are used to estabilish crude initial conditions in the voxels
%    b) ERA5 - reanalysis data format in the form of .net file
%    (https://cds.climate.copernicus.eu/datasets/reanalysis-era5-pressure-levels?tab=overview)
% GNSS ground-based observations according to the format in
% INTOMO/DATA/project_name/ATM
% Switches to be set:
%        .observations{'SWD'} - type of observations
%        .apriori{'INNER','OUTER'} - apriori model 
%        .constraints{'HORIZONTAL','VERTICAL'}  - apriori model cotnraints %not tested
%        .stacking{'NO'/number} %not tested
%        .filter{'KALMAN'/'ROBUST'} - Kalman filtering type
%        .solution{'REAL'/'SYNTHETIC'} - solution type
%        .decorelation{'NO'/'YES'} - decorelation type %not tested
%        .aprModel{'DETER'/'ERA5'} - apriori model source
%        .amtrix{folder_name} - name for work folder to contain observation matrices 
%        .parameterization{'constant'} - type of parameterization
%        .phi{'identity'} - include identity matrix
%        .method{'KALMAN'/'LSQ'} - Kalman or LSQ processing
%        .totalN{true/false} - processing in wet refractivity or (total) refracitivty
%        .saveAtmParam{true/false} - save atmospheric apriori param (e.g. pressure) to tomography output
%        .coord{'Formatted'} - input format for coordiates of GNSS stations
%        .ROres{(resolution of RO processing [Hz])} - e.g. 1 for 50Hz, 10 for 5Hz 
%        .refron{true/false} - ray tracing from apriori refractivities (if available)/from atmospheric parameters
%        .integrated{'yes'/'no'} - integrated tomography processing
%Path to be set manually:
% PATH_INSTALL......  installation folder
% PATH_EXTERNALSAVE.  tomography output save folder
%Paths set automatically:
% pathRO............  path to folder with RO.nc files
% pathTOMO..........  path for WORK directory
% pathORB...........  path to SP3 data
% pathCONF..........  path to configuration folder
% pathMETEO.........  path to meteorological data
% pathATM...........  path to ground-based data
% pathEXPORT........  tomography output save folder

%% Other variables
%           radii......   earth radius
%           GRIDboundaries  boumndaries of the RT model
%           res........   horizontal resolution of RT model
%           lat/lon/leveld_TOMO_RT   latitude/longitude/altitude of RT model nodes
%% ad 0. TIME MODE
%The selection of time mode is as follows:
switches.time_mode = 'POSTPROCESSING';

%% Data sources locations and associated names for source files, as well as destination folders for exported/visualised data,
PATH_INSTALL = 'C:\Users\Tanja\Documents\INTOMO\INTOMO-main\'; 
PATH_EXTERNALSAVE = 'C:\Users\Tanja\Documents\INTOMO\INTOMO-output\'; 
addpath([PATH_INSTALL 'Tomography\CONF']);
addpath([PATH_INSTALL 'Tomography\ENGINE']);
addpath(genpath([PATH_INSTALL 'Tomography\EXTERNAL']));
addpath([PATH_INSTALL 'Tomography\PPROCESS']);
addpath([PATH_INSTALL 'Tomography\READ']);
addpath([PATH_INSTALL 'Tomography\WRITE']);
addpath([PATH_INSTALL 'Tomography\VISUAL']);
addpath([PATH_INSTALL 'Tomography\DATA']);
addpath([PATH_INSTALL 'RayTracer\Frgn_Fcn']);
addpath([PATH_INSTALL 'RayTracer\Other_Fcn']);
addpath([PATH_INSTALL 'RayTracer\Tomo_Fcn']);
addpath([PATH_INSTALL 'RayTracer\UndulationFiles']);

% Select directory for saving the output 
PATH_SAVE = PATH_EXTERNALSAVE;
pathCONF = [PATH_INSTALL 'Tomography/CONF'];

%MISSION DATA
PROJECT_NAME = 'exp3_20260101';

%EPOCH RAY TRACING DATA FOLDER
pathTOMO = [PATH_INSTALL 'Tomography/DATA/'  PROJECT_NAME '/WORK/'];
addpath(pathTOMO);

%% Time slot for tomography processing and time interval of data sources

if strcmp(switches.time_mode,'POSTPROCESSING') == 1
    %TOMO2 TIME settings
    observation_start_TOMO = [2026 01 01 00 00 00]; %[year month day hour minute second] 
    observation_end_TOMO = [2026 01 01 23 59 59]; %[year month day hour minute second] 
    estimation_interval_TOMO = 3600; %[sec]
    %ORBITS Data
    observation_start_SP3 =observation_start_TOMO; 
    observation_end_SP3 = observation_end_TOMO;
    observation_interval_SP3 = 3600;
    %ZTD Observations
    observation_start_ZTD = observation_start_TOMO; 
    observation_end_ZTD = observation_end_TOMO;
    observation_interval_ZTD = 3600; 
end

%METEO Observations
observation_start_METEO = observation_start_TOMO;
observation_end_METEO = observation_end_TOMO;
observation_interval_METEO = 3600;
interpolation_interval_METEO = 3600;
observation_interval_NWP = 3600;

%% Processing parameters
% Set if integrated tomography [yes/no]
switches.integrated = "yes";
if switches.integrated == "yes"
    pathRO = [PATH_INSTALL 'Tomography/DATA/'  PROJECT_NAME '/RO' ];
end

%% INTOMO domains settings 
%Ray Tracing Model
model.radii = [6378.137, 6356.752314245]; %Earth radius
lam1 =  -12; %1st longitude of ERA5 data apriori file [deg]
lam2 = 36; %2nd longitude of ERA5 data apriori file  [deg] 
model.res = 0.25; %horizontal resolution of ERA5 grid [deg]
lat1 = 33; %1st latitude of ERA5 data apriori file  [deg]
lat2 = 65; %2nd latitude of ERA5 data apriori file [deg] 
unduFile= 'undu.mat'; 
model.levels_TOMO_RT = unique([(0.01:0.1:2.6), (2.6:0.2:6), (6:0.5:16), (16:1:36), (36:5:86)]);
model.ds = 2; %stepsize for RO ray tracing [km]
switches.refron = false; %simulate ray path from only refraction values ['true'], simulate ray path from temp, pres, wvpres ['true']
switches.ROres = 2; %frequency of LEO/GPS RO positions (1 is 50Hz, 10 is 5Hz etc)
model.lat_TOMO_RT = [lat1:model.res:lat2];
model.lon_TOMO_RT = [lam1:model.res:lam2];
model.num_lat_TOMO_RT = size(model.lat_TOMO_RT,2);
model.num_lon_TOMO_RT = size(model.lon_TOMO_RT,2);
model.num_levels_TOMO_RT = size(model.levels_TOMO_RT,2);
% Tomography Model
model.lat_TOMO = [43:2.0:55];  %specifying the location of model faces [deg],
model.lon_TOMO = [-2:2.0:26]; % the location of model faces [deg],
model.levels_TOMO = [0 500 1000 1500 2000 2500 3000 4500 6000 7500 9000 14500]; %specifying the location of model faces [m],
model.GRIDboundaries = [lam1 lam2 model.res lat1 lat2 model.res];
model.num_lat_TOMO = size(model.lat_TOMO,2);
model.num_lon_TOMO = size(model.lon_TOMO,2);
model.num_levels_TOMO = size(model.levels_TOMO,2);
clear lam1 lam2 lat1 lat2 res

% bounding box - all stations outside bounding box defined as [west_limit_TOMO, east_limit_TOMO, south_limit_TOMO, north_limit_TOMO] are going to be removed from processing.
model.west_limit_TOMO = -1.5;
model.east_limit_TOMO = 25.5;
model.south_limit_TOMO = 43.5;
model.north_limit_TOMO = 54.5;

% cut off angle for slant measurements
model.cut_off_angle = 10;

% GNSS station coordinates format
switches.coord = {'FORMATTED'};

% Randomly filter out GNSS station in selected range from each other ('0' to switch off)
switches.stat_range = 0; %[km]

% Save pressure, temperature and humidity grid information to model output
switches.saveAtmParam = true;

%% Orbits data path
pathORB = [PATH_INSTALL 'Tomography/DATA/'  PROJECT_NAME '/ORB/'];
addpath(pathORB);

%% Troposphere model path
pathMETEO = [PATH_INSTALL 'Tomography/DATA/'  PROJECT_NAME '/METEO'];
addpath(pathMETEO);  

%% Apriori/Initial conditions for tomography models {DETER/ERA5}
% >DETER - deterministic gpt/gpt2/UNB3mm models
% >ERA5 - ERA5 model
switches.aprModel='DETER';
observation_interval_APRIORI = 3600;

%% Ground-based data observations path
pathATM = [ PATH_INSTALL 'Tomography/DATA/'  PROJECT_NAME '/ATM/'];

%% Switches that refer to use (or not) of specific property of INTOMO model
% observations{'SWD'} chose between refractivity and or water
switches.observations = {'SWD',''};
% processing with total (true) or wet only (false) refractivity
switches.totalN = false;
% apriori{'INNER','OUTER','TOP','BOTTOM'} chose if you will use/provide apriori (initial conditions) for INNER (within bounding box) or/and OUTER (outside the bounding box) model
%it is highly recommend to provide for OUTER model apriori data
switches.apriori = {'INNER','OUTER','TOP','BOTTOM'}; %% not tested in other settings
% constraints{'HORIZONTAL','VERTICAL'} chose if HORIZONTAL/VERTICAL constraints are going to be used (as of 15/07/2015  VERTICAL constraints not operational)
%constraints might help to stabilise equation system in case of low horizontal troposphere variation, however it may damp horizontal variability in case of severe weather
switches.constraints = {'',''}; %% not tested
% stacking{'NO'} chose if the solution is going to be based on single observation epoch in the interval specified in estimation_interval_TOMO
%(in agreement with ZTD interval) or is going to be estimated every n*estimation_interval_TOMO using data from all epochs, you need to specify 'n - number', this
% will be helpful in constraint solution in case of low number of observations might smooth out active weather.
switches.stacking = {'NO'}; 
% method{'KALMAN'/'LSQ'} chose the method for the solution - Kalman filter ('KALMAN') or Least Squares ('LSQ')
switches.method = {'KALMAN'}; 
% filter{'KALMAN'/'ROBUST'} chose the level of Kalman filter robustness, standard KF - 'KALMAN', robust KF - 'ROBUST' for details please refere to the paper:
% Rohm W., Zhang K., Bosy J. Limited constraint, robust Kalman filtering for GNSS troposphere tomography Atmospheric Measurement Techniques, Vol. 7 No. 5, 2014, pp. 1475-1486
% DOI: 10.5194/amt-7-1475-2014
switches.filter = {'ROBUST'}; %% not tested in other settings
% solution{'REAL'/'SYNTHETIC'} chose whether you use REAL data (SWD/SIWV) hence the observation uncertainty is set realistic or the data provided are SYNTHETIC (SWD/SIWV)
% and the uncertainties are set to minimum (0.001m)
switches.solution = {'SYNTHETIC'}; 
% decorelation{'NO'/'YES'} chose whether you need to use the A matrix decorrelation procedure (will remove linearly dependent rows) for details please consult:
% Rohm W., Zhang K., Bosy J. Limited constraint, robust Kalman filtering for GNSS troposphere tomography Atmospheric Measurement Techniques, Vol. 7 No. 5, 2014, pp. 1475-1486
% DOI: 10.5194/amt-7-1475-2014
switches.decorelation = {'NO'}; %% not available yet
% phi{'identity'} chose whether you want to use the identity matrix for the propagation of wet refractivity ('identity')
% or define changes based on WRF data ('WRF', this option works only if WRF data are used as a priori information)
switches.phi = {'identity'}; 
% parametrization{'constant'/'bilinear-h'} chose whether you want to use the constant parametrization using constant value of Nw/WV in each voxel ('constant'), trilinear parametrization
% based on bilinear parametrization with height changes according to PWV
% changes ('bilinear-h') [Perler et al. 2011, Ding et al.2018]
switches.parametrization = {'constant'};
% grid parameterization
switches.regular = {'yes'}; %only regular node parametetriaztion is available in this version of INTOMO 

%Folder to save observation matrix A
if strcmp(switches.parametrization,'bilinear-h')
    switches.amtrix ='node';
else
    switches.amtrix ='const';
end

%% Export format
pathEXPORT = [PATH_SAVE PROJECT_NAME '/OUT/'];
%% Save conf file settings
save([PATH_INSTALL 'Tomography/CONF/' PROJECT_NAME '.mat']);
clearvars -except PATH_INSTALL PROJECT_NAME time_mode;
