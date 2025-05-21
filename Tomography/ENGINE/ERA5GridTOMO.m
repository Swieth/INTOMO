  function model = ERA5GridTOMO(ERAname, model,boundRT,date0,pathERA,pathCONF,obs_start,epoch,unduFile,switches)
% Function to generate .mat refraction grid file from .net cdf ERA5. Some
% of the parameters like pressure, temperature and humidity grids are not saved in output 
% 'model' structure to save the space, however user may unlock it if
% needed with switches.saveAtmParam set to true.
%---------------------------------------------------------------------
%   INPUT:
%       ERAname...   ERA5 net cdf file name
%       model.....   structure with tomography and raytracing model variables
%       pathERA...   path to net cdf ERA models folder
%       pathCONF..   path to configuration folder
%       boundRT...   boundaries of ray tracing model [deg]
%       obs_start.   matrix with epochs of the processing
%       epoch.....   processing epoch number
%       unduFile..   name of the file with undulation grid
%       switches..   structural matrix of settings of the INTOMO processing
% OUTPUT:
%       model
%           temp.......... temperature in RT grid nodes
%           pres.......... pressure in RT grid nodes
%           wvpr.......... specific humidity in RT grid nodes
%           refr.......... wet refractivity in RT grid nodes
%           refrRT........ total refractivity in RT grid nodes
%           temp_num_apr.. temperature in RT voxel centers 
%           pres_num_apr.. pressure refractivity in RT voxel centers 
%           wvpr_num_apr.. specific humidity refractivity in RT voxel centers 
%           refr_num_apr.. wet refractivity in RT voxel centers 
%           refrRT_num.... total refractivity in RT voxel centers
%           refr_apr...... wet refractivity in tomography model nodes\voxel centers
%           refr_aprF..... total refractivity in tomography model nodes\voxel centers 
%           wvpr_apr...... specific humidity in tomography model nodes\voxel centers
%           temp_apr...... temperature in tomography model nodes\voxel centers 
%           pres_apr...... pressure in tomography model nodes\voxel centers
%           BLh_pudel..... square coordinates of tomography model nodes
%           BLh_pudel_num. square coordinates of tomography model voxel centers
%           BLh_pudel_rad. geopgraphical coordinates of tomography model nodes  
%           BLh_pudel_num_rad. geopgraphcial coordinates of tomography model voxel centers
%           BLh_outer ..... square coordinates tomography model nodes of outer model
%           BLh_outer_num. square coordinates of tomography model voxel  centers of outer model
%           BLh_outer_rad . geopgraphical coordinates of tomography model nodes of outer model
%           BLh_outer_num_rad geopgraphcial coordinates of tomography model voxel centers of outer model
%           levels_mid_TOMO   altitudes of tomography model voxel centers
%           num_inner..... indices of voxels in tomography inner model
%           num_outer..... indices of voxels in tomography outer model 
%           rWGS.......... mean Earth radius in RT model nodes
%           LAT........... latitude matrix for tomography model nodes
%           LON........... longitude matrix for tomography model nodes
%           mid_levels_TOMO altitudes of tomography model voxel centers
%           mid_lat_TOMO    latitudes of tomography model voxel centers
%           mid_lon_TOMO    longitudes of tomography model voxel centers

        %% Part for Blh pudel
        if epoch == 1
            BLHstruc = pudel2(model,model.lat_TOMO,model.lon_TOMO,model.levels_TOMO,switches);
            BLHstrucRT = pudel2(model,model.lat_TOMO_RT,model.lon_TOMO_RT,model.levels_TOMO_RT,switches);
        else
            BLHstruc = [];
        end
        %% Generate unulation grid if not exists
        if ~exist([pathCONF,'/',unduFile])
           if license('test','map_toolbox')==1 
                unduera = Undulation(boundRT(4),boundRT(5),boundRT(1),boundRT(2),0.25,pathCONF,'unduera5');
           else
               start_time = obs_start(1,4)/24 + obs_start(1,5)/60/24 + obs_start(1,6)/60/60/24 ;
               jd_start = cal2jd(obs_start(1,1),obs_start(1,2),obs_start(1,3)+start_time);
               [~,~,unduera] = distr_T_gpt2RT(BLHstrucRT.BLh_pudel_rad,jd2mjd(jd_start));
               unduera = reshape(unduera,size(model.lat_TOMO_RT,2),size(model.lon_TOMO_RT,2),size(model.levels_TOMO_RT,2));
               unduera = flip(squeeze(unduera(:,:,1))',2);
           end
           save([pathCONF,'/',unduFile], 'unduera');
        else
           load([pathCONF,'/',unduFile], 'unduera'); % x  Kelut's area
        end
        
        %% Generate refracticion values for ray tracing use (ray tracing model)
        % layers spacing
        lam_N = [boundRT(1):0.25:boundRT(2)];
        if isempty(lam_N)
            lam_N = [boundRT(1):0.25:360-0.25 0:0.25:boundRT(2)];
        end
        phi_N = [boundRT(4):0.25:boundRT(5)];
        g0 = 9.80665;
      
        hpre = model.levels_TOMO_RT;
        % Starting date
        t0 = datenum(date0, 'yyyy-mm-dd HH:MM:SS');
        % Create timeseries
            qEpoch=datevec(t0);
        % find ERA files foe datetime series
        [eNwp, ~ ] = fera5Epoch(qEpoch);
        % find unique files and epoch
        [ueNwp, ~, ~] = unique(eNwp, 'rows');
        epoch = ueNwp(1,5); %1 za ep
        % LOAD NWP FILES AND GENERATE GRIDS WITH METEOROLOGICAL PARAMETERS
        ERA_ID = [pathERA,'/',ERAname];
        eradim = ncinfo(ERA_ID, 'z');
        size_era = eradim.Size;
        % create lon lat grid
        lat_var = ncread(ERA_ID,'latitude');
        lon_var = ncread(ERA_ID,'longitude');
        [LAT, LON] = meshgrid(lat_var, lon_var);
        LAT = LAT.*pi/180;
        LON = LON.*pi/180;  
        % specific humidity
        Spechum = ncread(ERA_ID,'q', [1 1 1 epoch], [size_era(1) size_era(2) size_era(3) 1], [1 1 1 1]);
        % temperature
        tempGrid = ncread(ERA_ID,'t', [1 1 1 epoch], [size_era(1) size_era(2) size_era(3) 1], [1 1 1 1]);
        % geopotential
        Geop = ncread(ERA_ID,'z', [1 1 1 epoch], [size_era(1) size_era(2) size_era(3) 1], [1 1 1 1]);
        % Flip geopotential heights structure according to .nc ECMWF data format
        try
            iso = ncread(ERA_ID,'level');
            Spechum = flip(Spechum,3);
            tempGrid = flip(tempGrid,3);
            iso = flip(iso);
        catch
            iso = ncread(ERA_ID,'pressure_level');
            Geop  = flip(Geop,3);
        end
        % flip all data since ERA provides layers from top to bottom
        Geoph = flip(Geop,3)./(1000*g0);
        % Convert to hPa and resize to 3-dimensional array - dim: (26, 361, 720)
        pGrid = repmat(iso,1,size_era(1),size_era(2));
        pGrid = double(permute(pGrid,[2 3 1]));
        if model.GRIDboundaries(3) ~= 0.25
            n = model.GRIDboundaries(3)/0.25;
        else
            n = 1;
        end
        LAT = LAT(1:n:end,1:n:end);
        LON = LON(1:n:end,1:n:end);
        lam_N = lam_N(1:n:end);
        phi_N = phi_N(1:n:end);
        [pGrid3D, eGrid3D, tempGrid3D, RTstruc.rWGS ] = gridcalcRT(LAT,LON, unduera(1:n:end,1:n:end), pGrid(1:n:end,1:n:end,:), tempGrid(1:n:end,1:n:end,:), Spechum(1:n:end,1:n:end,:), Geoph(1:n:end,1:n:end,:), hpre,phi_N,lam_N,[],model);
        RTstruc.pGrid3D = pGrid3D(:,:,1:size(hpre,2));
        RTstruc.eGrid3D = eGrid3D(:,:,1:size(hpre,2));
        RTstruc.tempGrid3D = tempGrid3D(:,:,1:size(hpre,2));
        [LAT, LON] = meshgrid(flip(phi_N), lam_N);
        RTstruc.LAT = LAT.*pi/180;
        RTstruc.LON = LON.*pi/180;  
        %% Convert meteorological parameters to refraction grid
        [Nh, Nw ] = refcalc(RTstruc.pGrid3D, RTstruc.tempGrid3D, RTstruc.eGrid3D,'b');
        RTstruc.N3D = Nw;  
        RTstruc.N3D_RT = Nw + Nh;
        
        %% Generate refracticion values for tomography use (tomography domain)
        clearvars pGrid3D_num eGrid3D_num tempGrid3D_num  pGrid3D eGrid3D tempGrid3D N_num Nw_num 
        [LAT, LON] = meshgrid(lat_var, lon_var);
        LAT = LAT(1:n:end,1:n:end);
        LON = LON(1:n:end,1:n:end);
        lam_N = model.lon_TOMO;
        phi_N = model.lat_TOMO;
        [pGrid3D, eGrid3D, tempGrid3D, BLHRTstruc.rWGS ] = gridcalcRT(LAT*pi/180,LON*pi/180, unduera(1:n:end,1:n:end), pGrid(1:n:end,1:n:end,:), tempGrid(1:n:end,1:n:end,:), Spechum(1:n:end,1:n:end,:), Geoph(1:n:end,1:n:end,:), model.levels_TOMO./1000,phi_N,lam_N,RTstruc,model);
        BLHRTstruc.pGrid3D = pGrid3D(:,:,1:size(model.levels_TOMO,2));
        BLHRTstruc.eGrid3D = eGrid3D(:,:,1:size(model.levels_TOMO,2));
        BLHRTstruc.tempGrid3D = tempGrid3D(:,:,1:size(model.levels_TOMO,2));
        [LAT, LON] = meshgrid(flip(phi_N), lam_N);
        BLHRTstruc.LAT = LAT.*pi/180;
        BLHRTstruc.LON = LON.*pi/180;  
        %% Convert meteorological parameters to refraction grid
        [Nh, Nw] = refcalc(BLHRTstruc.pGrid3D, BLHRTstruc.tempGrid3D, BLHRTstruc.eGrid3D,'c');
        BLHRTstruc.N3D = Nw;  
        BLHRTstruc.N3D_RT = Nw + Nh;    
        %% Calculate refractivity values in the voxel centers (only for tomography domain)
        latu = LAT(1,1:end-1)'+diff(LAT(1,:))'./2;
        lonu = LON(1:end-1,1)+diff(LON(:,1))./2;
        altu = model.levels_TOMO(1:end-1)'+diff(model.levels_TOMO)'./2;
        if epoch == 1
            BLHRTstrucnum = pudel2(model,flip(latu)',lonu',altu',switches);
        end
        for z = 1:size(altu,1) 
            for y = 1:size(latu,1)
                for x = 1:size(lonu,1)
                    N = IDW_atom(LAT*pi()/180, LON*pi()/180, BLHRTstruc.tempGrid3D, BLHRTstruc.eGrid3D, BLHRTstruc.pGrid3D,BLHRTstruc.N3D_RT, BLHRTstruc.N3D ,BLHRTstruc.rWGS,[(90-latu(y))*pi()/180 lonu(x)*pi()/180 altu(z)],model.levels_TOMO,switches);
                    BLHRTstrucnum.N3D_num(x,y,z) = N.Nw;
                    BLHRTstrucnum.tempGrid3D_num(x,y,z) = N.T;
                    BLHRTstrucnum.eGrid3D_num(x,y,z) = N.h;
                    BLHRTstrucnum.pGrid3D_num(x,y,z) = N.p;
                    BLHRTstrucnum.N3D_RT_num(x,y,z) = N.Nt;
                end
            end
        end
        %% Save variables
        if switches.saveAtmParam
            model.temp{epoch,:,:,:} = RTstruc.tempGrid3D;
            model.pres{epoch,:,:,:} = RTstruc.pGrid3D;
            model.wvpr{epoch,:,:,:} = RTstruc.eGrid3D;
            model.refr{epoch,:,:,:} = RTstruc.N3D;
        end
        model.refrRT{epoch,:,:,:} = RTstruc.N3D_RT;
        if strcmp(switches.parametrization,'constant')
            if switches.saveAtmParam
                model.temp_num_apr{epoch,:,:,:} = flip(BLHRTstrucnum.tempGrid3D_num,2);
                model.pres_num_apr{epoch,:,:,:} = flip(BLHRTstrucnum.pGrid3D_num,2);
                model.wvpr_num_apr{epoch,:,:,:} = flip(BLHRTstrucnum.eGrid3D_num,2);
                model.refrRT_num{epoch,:,:,:} = flip(BLHRTstrucnum.N3D_RT_num,2);
            end
            model.refr_num_apr{epoch,:,:,:} = flip(BLHRTstrucnum.N3D_num,2);
        else
            model.refr_apr{epoch,:,:,:} = flip(Nw,2);
            model.refr_aprF{epoch,:,:,:} = flip(Nw,2) + flip(Nh,2); 
            if switches.saveAtmParam
                model.wvpr_apr{epoch,:,:,:} = flip(BLHRTstruc.eGrid3D,2);
                model.temp_apr{epoch,:,:,:} = flip(BLHRTstruc.tempGrid3D,2);
                model.pres_apr{epoch,:,:,:} = flip(BLHRTstruc.pGrid3D,2);
            end
        end
        if epoch == 1
            model.BLh_pudel = BLHstruc.BLh_pudel;
            model.BLh_pudel_num = BLHstruc.BLh_pudel_num;
            model.BLh_pudel_rad = BLHstruc.BLh_pudel_rad;
            model.BLh_pudel_num_rad = BLHstruc.BLh_pudel_num_rad;
            model.BLh_outer = BLHstruc.BLh_outer;
            model.BLh_outer_num = BLHstruc.BLh_outer_num;
            model.BLh_outer_rad = BLHstruc.BLh_outer_rad;
            model.BLh_outer_num_rad = BLHstruc.BLh_outer_num_rad;
            model.levels_mid_TOMO  = [0;squeeze(model.BLh_pudel_num_rad(3,1,:))];
            model.num_inner = BLHRTstrucnum.num_inner;
            model.num_outer = BLHRTstrucnum.num_outer;
            model.rWGS = RTstruc.rWGS;
            model.LAT = RTstruc.LAT;
            model.LON = RTstruc.LON;
            model.mid_levels_TOMO = altu;
            model.mid_lat_TOMO = latu;
            model.mid_lon_TOMO = lonu;
        end
  end
 