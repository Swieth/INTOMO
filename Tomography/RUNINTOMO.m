%RUNINTOMO - this file is a script that distribute the processes in INTOMO v. 1.0 software
%it takes into account the conf -file MISSION specific. Please have in mind
%that this version of INTOMO is in experimental phase and some of the
%processing options may not work properly.
%THIS FILE SHOULD NOT BE modified unless you really know what are you doing in case of
%questions and remarks please consult FAQ, if this won't help contact me on:
%adam.cegla@upwr.edu.pl

clear
clc
tic
conf;
confFILE = [PATH_INSTALL 'Tomography/CONF/' PROJECT_NAME '.mat'];
load(confFILE,'switches','pathTOMO'); %loading the configuration for further processing

%% Find coordinates of GNSS stations 
% For now there is no universal function to read the coordiantes and
% observation data for GNSS station. The users need to write it by
% themselves
if strcmp(switches.coord,'FORMATTED')
    load(confFILE,'model');
    try
        [NAME,BLh_ori] = readBLh([PATH_INSTALL 'Tomography/CONF/stations.txt']); 
    catch
        error('RUNINTOMO: Wrong format of input data for GNSS stations coordinates')
    end
    %here other functions to read GNSS stations coordiantes may be added
end

%% Recalculate to XYZ
for i = 1:size(BLh_ori,1)
    Xsta(:,i) = cspice_georec(BLh_ori(i,3)*pi()/180, BLh_ori(i,2)*pi()/180, BLh_ori(i,5)/1000, model.radii(1),(model.radii(1)-model.radii(2))/model.radii(1));
end
% Cut stations outisde tomography model
[X,Y,Z,lat,lon,h,H,NAME] = boundingTOMOLAB(model.east_limit_TOMO,model.west_limit_TOMO,model.north_limit_TOMO,model.south_limit_TOMO,Xsta,BLh_ori,NAME);
% Filter out stations to selected spatial density
[X,Y,Z,lat,lon,h,H,NAME] = deleteStat(X,Y,Z,lat,lon,h,H,NAME,switches.stat_range);
BLh = [1:size(lat,1);lat'*180/pi();lon'*180/pi(); h'];
BLh = BLh';
BLh_ori = [1:size(lat,1);lat'*180/pi();lon'*180/pi(); H'; h'];
NAME_stats = NAME; 

%% Calculate refractivities for ray tracing and tomography model
% load initial variables
load(confFILE,'observation_start_TOMO','observation_end_TOMO','estimation_interval_TOMO', 'observation_interval_SP3','observation_interval_ZTD','observation_interval_METEO','observation_interval_NWP','interpolation_interval_METEO','observation_interval_APRIORI','model','pathMETEO','pathCONF','pathTOMO','unduFile');
model.BLh = BLh;
model.BLH = BLh_ori(1:4,:)';
model.BLh_ori = BLh_ori;
model.NAME = NAME;
model.X = X;
model.Y = Y;
model.Z = Z;
% find time of the epochs
[obs_set] = time_listing(switches,observation_start_TOMO,observation_end_TOMO, estimation_interval_TOMO,observation_interval_SP3,observation_interval_ZTD,observation_interval_NWP,interpolation_interval_METEO);
% Calculate variables
if ~exist([pathTOMO,'model.mat'],'file')
    for epoch = 1:size(obs_set.observation_set,1)
        date0 = [num2str(obs_set.observation_set(epoch,3)),'-',num2str(obs_set.observation_set(epoch,7)),'-',num2str(obs_set.observation_set(epoch,8)),' ',num2str(obs_set.observation_set(epoch,9)),':00:00'];
        t0 = datenum(date0, 'yyyy-mm-dd HH:MM:SS');
        ERAname = ['ERA5_',num2str(obs_set.observation_set(epoch,3)),'-',num2str(obs_set.observation_set(epoch,7)),'-',num2str(obs_set.observation_set(epoch,8)),'.nc'];
        model =  ERA5GridTOMO(ERAname, model,model.GRIDboundaries,date0,pathMETEO,pathCONF,observation_start_TOMO,epoch,unduFile,switches);
        disp([ERAname,'_',num2str(epoch)])
    end   
    save([pathTOMO,'model.mat'],'model');
    disp('Tomography and ray tracing domain data(model.mat file) calculated')
else
    load([pathTOMO,'model.mat'],'model');
    disp('Tomography and ray tracing domain data(model.mat file) found')
end
save(confFILE,'obs_set','-append');

%% Download and read sp3 data
load(confFILE, 'pathORB','pathTOMO','pathORB','switches');
downloadORB(pathORB,obs_set.observation_set_SP3);
%SP3 processing
SP3data = readSP3dat(obs_set,pathORB);
[SP3Xn,SP3Yn,SP3Zn,PRNn,obs_set.observation_set_SP3] = interSP3(SP3data,obs_set.observation_set_SP3);
clear SP3data   

%% Read ZTD values from observation files
if strcmp(switches.solution,'REAL')
    load(confFILE,'pathATM');
    [ZTDA, MZTDA, DGNA, MDGNA, DGEA, MDGEA, NAMES, MISS_STAT] = readtxtOBS(pathATM, model.NAME, obs_set.observation_set);
    [ZTDA, MZTDA, DGNA, MDGNA, DGEA, MDGEA, model] = screenZTD(ZTDA, MZTDA, DGNA, MDGNA, DGEA, MDGEA, MISS_STAT, model);
end

%% Calculate apriori Nw values of GPT2 model
if strcmp(switches.parametrization,'constant')
    for epoch = 1:size(obs_set.observation_set,1)
         %Temperature
         [T,p,undu] = distr_T_gpt2RT(model.BLh_pudel_num_rad,jd2mjd(obs_set.observation_set(epoch,1)));
         %Water Vapour Pressure
         [E,~]=distr_e_unb3RT(model.BLh_pudel_num_rad,undu,obs_set.observation_set(epoch,1));
         %%%%REFRACTIVITY
         if switches.totalN
            [NhT, NwT ]= refcalc(p', T, E,'a');
            NwT = NhT+NwT;
            NwT = real(NwT)';
         else
            [NwT,~] = eT2Nw(E,T);
            NwT = NwT';
         end
         %d)forming model for all observation hours
         NwT = NwT';
         synthetic.Nw_num(epoch,:) = NwT;
         synthetic.Nw_out_num = [];
         apriori.Nw_DETER_num(epoch,:) = NwT;
         apriori.Nw_DETER_out_num = [];
         %%%%%%%%%WATER VAPOUR%%%%%%%%%%%%
         R = 461.525;
         WVT = 100*E./(R*T)*1000;
         WVT = WVT';
         synthetic.WV_num(epoch,:) = WVT;
         synthetic.WV_out_num = [];
         clear NwT Nw3DT WVT WV3DT R T E T3D E3D
    end
elseif  strcmp(switches.parametrization,'bilinear-h')  
    for epoch = 1:size(obs_set.observation_set,1)
         %Temperature
         [T,p,undu] = distr_T_gpt2RT(model.BLh_pudel_rad,jd2mjd(obs_set.observation_set(epoch,1)));
         %Water Vapour Pressure
         [E,~]=distr_e_unb3RT(model.BLh_pudel_rad,undu,obs_set.observation_set(epoch,1));
         %%%%REFRACTIVITY
         if switches.totalN
            [NhT, NwT ]= refcalc(p', T, E,'a');
            NwT = NhT+NwT;
            NwT = real(NwT)';
         else
            [NwT,~] = eT2Nw(E,T);
            NwT = NwT';
         end
         %d)forming model for all observation hours
         NwT = NwT';
         synthetic.Nw(epoch,:) = NwT;
         synthetic.Nw_out = [];
         apriori.Nw_DETER(epoch,:) = NwT;
         apriori.Nw_DETER_out = [];
         %%%%%%%%%WATER VAPOUR%%%%%%%%%%%%
         R = 461.525;
         WVT = 100*E./(R*T)*1000;
         WVT = WVT';
         synthetic.WV(epoch,:) = WVT;
         synthetic.WV_out = [];
         clear NwT Nw3DT WVT WV3DT R T E T3D E3D
    end
end
 
%% Collect RO data
if strcmp(switches.integrated,'yes')
    load(confFILE,'pathRO');
    [cordRO] = findRO(pathRO,obs_set);
else
    cordRO = [];
end

%% Assing apriori model values to tomography variables
for epoch = 1:size(obs_set.observation_set,1)  
    if strcmp(switches.parametrization,'constant')
        if switches.totalN
            NwT1_num = model.refrRT_num{epoch};
        else
            NwT1_num = model.refr_num_apr{epoch};
        end
        E1_num = model.wvpr_num_apr{epoch};
        T1_num = model.temp_num_apr{epoch};
        NwT_num = reshape(NwT1_num,1,size(NwT1_num,1)*size(NwT1_num,2)*size(NwT1_num,3));
        E_num = reshape(E1_num,1,size(NwT1_num,1)*size(NwT1_num,2)*size(NwT1_num,3));
        T_num = reshape(T1_num,1,size(NwT1_num,1)*size(NwT1_num,2)*size(NwT1_num,3));
        NwT_num = NwT_num';
        apriori.Nw_ERA5_num(epoch,:) = NwT_num;
        R = 461.525;
        WVT_num = 100*E_num./(R*T_num)*1000;
        WVT_num = WVT_num';
        apriori.WV_ERA5_num(epoch,:) = WVT_num;
    elseif strcmp(switches.parametrization,'bilinear-h')      
        if switches.totalN
            NwT1 = model.refr_aprF{epoch};
        else
            NwT1 = model.refr_apr{epoch};
        end
        E1 = model.wvpr_apr{epoch};
        T1 = model.temp_apr{epoch};
        NwT = reshape(NwT1,1,size(NwT1,1)*size(NwT1,2)*size(NwT1,3));
        E = reshape(E1,1,size(NwT1,1)*size(NwT1,2)*size(NwT1,3));
        T = reshape(T1,1,size(NwT1,1)*size(NwT1,2)*size(NwT1,3));
        NwT = NwT';
        apriori.Nw_ERA5(epoch,:) = NwT;
        R = 461.525;
        WVT = 100*E./(R*T)*1000;
        WVT = WVT';
        apriori.WV_ERA5(epoch,:) = WVT;
    end   
end

%% Calculate the dry part of Zenith Total Delay
if strcmp(switches.solution,'REAL')
    ZHD = pBLh2ZHDRT(model,NAMES);
    clear p
end

%% Create station structure with GNSS and RO coordinates and delays for use in RT 
station_filename = [pathTOMO 'station.mat'];
if exist(station_filename,'file')==0
   if strcmp(switches.solution,'SYNTHETIC')
       [station] = construct_station_LAB(model.BLh,model.BLH,0,0,0,model.NAME,PRNn,SP3Xn,SP3Yn,SP3Zn,obs_set.observation_set_SP3,model.cut_off_angle,cordRO,switches,0,0,0,0);
   elseif strcmp(switches.solution,'REAL')
       [station] = construct_station_LAB(model.BLh,model.BLH,ZTDA,MZTDA,ZHD,model.NAME,PRNn,SP3Xn,SP3Yn,SP3Zn,obs_set.observation_set_SP3,model.cut_off_angle,cordRO,switches,DGNA,MDGNA,DGEA,MDGEA);
   end
   save(station_filename,'station');
else
   load(station_filename,'station');
end

%% Set up apriori values
if strcmp(switches.aprModel,'DETER') 
    if strcmp(switches.observations{1,1},'SWD') 
        if strcmp(switches.parametrization,'constant')
            [num_Nw_out_num, Nw_obs_out_num] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_outer,size(obs_set.observation_set,1),1))',apriori.Nw_DETER_num);
            [num_Nw_num, Nw_obs_num] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_inner',size(obs_set.observation_set,1),1))',apriori.Nw_DETER_num);
        elseif strcmp(switches.parametrization,'nodes') || strcmp(switches.parametrization,'bilinear-h')  
            [num_Nw_out, Nw_obs_out] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_outer,size(obs_set.observation_set,1),1))',apriori.Nw_DETER);
            [num_Nw, Nw_obs] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_inner',size(obs_set.observation_set,1),1))',apriori.Nw_DETER);
        end
   end
elseif strcmp(switches.aprModel,'ERA5')
   if strcmp(switches.observations{1,1},'SWD') 
        if strcmp(switches.parametrization,'constant')
            [num_Nw_out_num, Nw_obs_out_num] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_outer,size(obs_set.observation_set,1),1))',apriori.Nw_ERA5_num);
            [num_Nw_num, Nw_obs_num] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_inner',size(obs_set.observation_set,1),1))',apriori.Nw_ERA5_num);        
        elseif strcmp(switches.parametrization,'nodes') || strcmp(switches.parametrization,'bilinear-h')  
            [num_Nw_out, Nw_obs_out] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_outer,size(obs_set.observation_set,1),1))',apriori.Nw_ERA5);
            [num_Nw, Nw_obs] = aprioriCONSTRNWP([],[],[],[],(1:size(obs_set.observation_set,1))',(repmat(model.num_inner',size(obs_set.observation_set,1),1))',apriori.Nw_ERA5);
        end
    end
end
  
%% Assign apriori values to tomography variables
if strcmp(switches.aprModel,'DETER') 
    if strcmp(switches.observations{1,1},'SWD') 
        if strcmp(switches.parametrization,'constant')
            values.Nw_apr = apriori.Nw_DETER_num;
            values.num_Nw_num = num_Nw_num;
            values.Nw_obs_num = Nw_obs_num;
            values.num_Nw_out_num = num_Nw_out_num;
            values.Nw_obs_out_num = Nw_obs_out_num;
            values.Nw_out_num = apriori.Nw_DETER_num; 
        elseif strcmp(switches.parametrization,'bilinear-h')          
            values.Nw_apr = apriori.Nw_DETER;
            values.num_Nw = num_Nw;
            values.Nw_obs = Nw_obs;
            values.num_Nw_out = num_Nw_out;
            values.Nw_obs_out = Nw_obs_out;
            values.Nw_out = apriori.Nw_DETER;
        end
    end
elseif strcmp(switches.aprModel,'ERA5')
    if strcmp(switches.observations{1,1},'SWD') 
        if strcmp(switches.parametrization,'constant')
            values.Nw_apr = apriori.Nw_ERA5_num;
            values.num_Nw_num = num_Nw_num;
            values.Nw_obs_num = Nw_obs_num;
            values.num_Nw_out_num = num_Nw_out_num;
            values.Nw_obs_out_num = Nw_obs_out_num;
            values.Nw_out_num = apriori.Nw_ERA5_num;  
        elseif strcmp(switches.parametrization,'bilinear-h')          
            values.Nw_apr = apriori.Nw_ERA5;
            values.num_Nw = num_Nw;
            values.Nw_obs = Nw_obs;
            values.num_Nw_out = num_Nw_out;
            values.Nw_obs_out = Nw_obs_out;
            values.Nw_out = apriori.Nw_ERA5;
        end
    end
end

%% Add satellites or receiver to exclude from processing
values.rem_REC = {''};
values.rem_SAT = [];
model.add_bias = 0;
model.add_var = 0;
switches.nearest_n = 8;
if strcmp(switches.parametrization,'bilinear-h')  
    model.Xbil= repmat(model.lat_TOMO,length(model.lon_TOMO),1);
    model.Ybil = repmat(model.lon_TOMO',1,length(model.lat_TOMO));
    model.bilmeshLatLon = reshape(model.BLh_pudel_rad,3,size(model.BLh_pudel_rad,2)*size(model.BLh_pudel_rad,3)); 
end
toc
load(confFILE,'PATH_SAVE');
save([pathTOMO,'\model.mat'],'model','station','values','pathTOMO','PATH_SAVE','obs_set','switches','apriori','PROJECT_NAME')

%% Tomography processing start
switches.project_name = PROJECT_NAME;
output = intomolab(station,model,values,pathTOMO,PATH_SAVE,obs_set.observation_set,switches);
disp('Tomography processing completed')

%% Save the results
save([PATH_SAVE,'\',PROJECT_NAME,'\OUT\OUTPUT_TOMO.mat'],'output')



