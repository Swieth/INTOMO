function [A, SWD_GNSS,SIWV,REC,elev,SAT,station_name,R_SIWV,R_SWD,SWD_nodes_integ,RelDist,dSWD,Avec] = groundRT(A,Avec,model,station,epoch,X_st,satcoord,elevat,ray_num,nr,sat,nodes_columns,levels,num_lev,Tmatrix,Nw_apr,switches)                    
%Function used to build basic matrices from GNSS data based on 3D ray
%tracing
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       A..........   observational matrix with signal deriveratives
%       Avec.......   the ray path vectors of each signal in tomography voxel
%       model......   parameters of tomography and ray tracing models
%       station....   structural matrix containing all informations regarding GNSS observations
%       epoch......   number of processing epoch
%       X_st.......   coordinates of the GNSS receiver
%       satcoord...   coordinates of the GNSS satellite
%       elevat.....   elevation angle to satellite
%       ray_num....   number of row in A matrix
%       sat........   id of GPS satellite
%       nodes_columns segragated id's of voxels for use in bilin interp
%       levels.....   altitude layers of tomography model
%       num_lev....   number of voxels in vertical direction
%       Tmatrix....   second derivates of refractivity for use in bilin interp
%       Nw_apr.....   apriori values of refractivities
%       swtiches...   structural matrix of settings of the INTOMO processing
%%%OUTPUT:
%       A..........   observational matrix with signal deriveratives
%       SWD_GNSS...   the values of Slant Delay (SD)
%       SIWV.......   the values of Slant Integrated Water Vapour
%       REC........   id of GNSS stations to remove
%       elev.......   elevation angles of ray traced signals
%       SAT........   id of the satellites
%       station_name  names of GNSS stations
%       R_SWD......   the errors of the slant delays
%       R_SIWV.....   the errors of Slant Integrated Water Vapour
%       SWD_nodes_integ.  SWD integral based on a priori data
%       RelDist....   total distance travesred by GNSS signal
%       dSWD.......   difference between ray traced and observed SD or EP 
%       Avec.......   the ray path vectors of each signal in tomography voxel

    %% Calculate GNSS signal delay with 3D ray tracing
    tic;
    switches.refron = false;
    switches.ground = true;
    rayRT = voxel_dist_3D_combined(model.refrRT{epoch},model.refr{epoch},model.levels_TOMO_RT'*1000,X_st',satcoord,0,0.005,model.rWGS,elevat*pi()/180,model.LAT,model.LON,model.temp{epoch},model.pres{epoch},model.wvpr{epoch},switches);
    toc
    %% Calculate distances in tomography model
    rayRT = interateModelsA(model,rayRT,switches);               
    SWD_nodes_integ = [];
    % Check if inside the tomography model
    if ~isnan(satcoord(1)) && isfield(rayRT,'voxIN') 
        if  switches.totalN
            SWD = rayRT.dL; 
        else
            SWD = rayRT.dLw; 
        end
        RelDist = rayRT.voxIN.distanceRel;    
        %% Sum distances in each voxel and add to A matrix
        if isfield(rayRT.voxEM, 'distanceCorr')
            iposEM1 = [];
            iposEM2 = [];
            for it = 1:size(rayRT.voxEM.iposEm,1)
                iposver = rayRT.voxEM.iposEm(it,:);
                distver = rayRT.voxEM.distanceCorr(it,:);
                iposver = nonzeros(iposver)';
                distver = distver(1:size(iposver,2));
                iposEM1 = [iposEM1 iposver];
                iposEM2 = [iposEM2 distver];  
                clear iposver distver
            end
            iposEM = [iposEM1;iposEM2];
            iposIN = [rayRT.i_pos; rayRT.vox(:,8)'];
            ind = [iposIN iposEM]'; 
        else
            iposIN = [rayRT.i_pos; rayRT.vox(:,8)];
            ind = iposIN'; 
        end
        ind = sortrows(ind);
        [row_ind,~,row_ori]=unique(ind(:,1));
        ind_sum=[row_ind cell2mat(accumarray(row_ori,1:numel(row_ori), [],@(x) {sum(ind(x,2:end),1)}))];
        ind_sum = ind_sum'; 
        % Calculate vectors of the signal in each voxels
        if strcmp(switches.parametrization,'constant')
            % Index Avec by the IN-voxel ids only (ascending, matching coordVector
            % row order).  EM-voxel rows of Avec remain zero, which is correct
            % because voxIN.coord carries no entry/exit chord for those voxels.
            in_vox = unique(rayRT.i_pos);
            Avec(in_vox, ray_num, :) = coordVector(rayRT.voxIN.coord, rayRT.i_pos);
            A(ind_sum (1,:),ray_num) = ind_sum(2,:);
        elseif strcmp(switches.parametrization,'bilinear-h')
            azi = station.h(nr).satellite(sat).azi;
            if isfield(rayRT.rayIND,'X_ray_b')
                coordIN = rayRT.rayIND.X_ray_b';
                coorIntersect = rayRT.voxIN.coord;
                coorIntersect = [squeeze(coorIntersect(:,:,1));squeeze(coorIntersect(:,:,2))];
                coorIntersect(isnan(coorIntersect(:,1)),:) = [];
                [~,ia,~] = unique(rayRT.i_pos');
            else
                exPh = 0;
                dexPh = [];
                RelDist = 0;
                R_SWD = 0;
                dexPh = 0;
                coord = [];
                return
            end
            coordINpar = coordIN(sort(ia),:);
            if any(any(isnan(coordINpar))) 
                par = [];
            else
                par = [coordINpar;coorIntersect(end,:)];
            end
            for i = 2:size(par,1)
                par(i,4) = vecnorm(par(1,1:3)-par(i,1:3));
            end
            [A,~,~] = bilinearRT(model,A,elevat,azi,nodes_columns,levels,num_lev,Tmatrix,Nw_apr,epoch,ray_num,par,switches);
        end
        %% Save Slant Delay value
        if strcmp(switches.solution,'REAL')
            if switches.totalN
                rayIND = rayRT.rayIND;
                SWD_GNSS = station.h(nr).satellite(sat).STD;
                R_SWD = station.h(nr).satellite(sat).M_STD;
                dSWD = SWD_GNSS - SWD;
                % Check whether the signal is in the model
                boundaryCheck = rayIND.dLs - rayRT.dL;
                if abs(boundaryCheck)> 0.1
                    SWD_GNSS = NaN;
                    dSWD = NaN;
                    R_SWD =  NaN; 
                    warning('groundRT: Simulated ray path is out of the tomography model. Deleting from A matrix.')
                end
            else
                SWD_GNSS = station.h(nr).satellite(sat).SWD;
                R_SWD = station.h(nr).satellite(sat).M_SWD;
                dSWD = SWD_GNSS - SWD;
            end
        else
            SWD_GNSS = SWD;
            dSWD = NaN;
            R_SWD =  NaN; 
        end
        SIWV = NaN;
        REC = nr;
        elev = elevat;
        SAT = station.h(nr).satellite(sat).PRN;
        station_name = {station.h(nr).nazwa};
        R_SIWV =  NaN;
        clear rayRT iposEM  iposIN  ind ind_sum
    else
        SWD_GNSS = 0;
        SIWV = NaN;
        REC = nr;
        elev = elevat;
        SAT = station.h(nr).satellite(sat).PRN;
        station_name = {station.h(nr).nazwa};
        R_SIWV =  NaN;
        R_SWD =  NaN; 
        RelDist = NaN;
        dSWD = NaN;
    end
                