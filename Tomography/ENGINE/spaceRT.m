function [A_RO, exPh,RelDist,R_SWD,dexPh,coord,Avec] = spaceRT(model,A_RO,Avec,station,epoch,nodes_columns,levels,num_lev,Tmatrix,Nw_apr,switches,path_save)   
%Function used to build basic matrices from RO data based on 3D ray
%tracing
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       model......   parameters of tomography and ray tracing models
%       A_RO.......   observational matrix with signal deriveratives for RO
%       Avec.......   the ray path vectors of each signal in tomography voxel
%       station....   structural matrix containing all informations regarding GNSS observations
%       epoch......   .processing epoch number
%       nodes_columns segragated id's of voxels for use in bilin interp
%       levels.....   altitude levels of tomography model
%       num_lev..     number of altitude levels of tomography model
%       Tmatrix....   second derivates of refractivity for use in bilin interp
%       Nw_apr.....   apriori values of refractivities
%       swtiches...   structural matrix of settings of the INTOMO processing
%%%OUTPUT:
%       A_RO.......   observational matrix with signal deriveratives for RO
%       exPh.......   the values of excess phase
%       RelDist....   total distance travesred by GNSS signal
%       R_SWD......   the errors of the slant delays
%       dexPh......   difference between ray traced and observed excess
%       coord......   coord. of interesetcions between RO and voxel model
%       Avec.......    the ray path vectors of each signal in tomography voxel

        %% Check whether RO is available during selected epoch
        ray_num = 0;
        switches.refron = false;
        switches.ground = false;
        radii = model.radii; 
        if ~isfield(station,'ro')
            A_RO = [];
            exPh = [];
            RelDist = [];
            R_SWD = [];
            dexPh = [];
            coord = [];
            return
        end
        %% Prepare RO data for ray tracing
        if ~isempty(station.ro) 
            exPh= [];
            RelDist = [];
            R_SWD = [];
            dexPh = [];
            coord = [];
            for roit = 1:size(station.ro,2)
                % Get coordinates
                GPScoord = station.ro(roit).cordT;
                LEOcoord = station.ro(roit).cordR;
                if ~isempty(GPScoord)
                    dateRO = station.ro(roit).date;
                    if strcmp(switches.solution,'REAL')
                        exLC = station.ro(roit).exLC;
                    end
                    exL2 = station.ro(roit).exL2;
                    % Get positions interval based on RO netCDF structure
                    % (ascending/descending). In order to work properly,
                    % .nc file must contain the same number of excess phase
                    % values and LEO/GPS coordinates.
                    posit = [1:switches.ROres:size(GPScoord,1)];
                    if strcmp(switches.solution,'REAL')
                        % REAL mode: original logic - exL2 and coordinates on same grid,
                        if length(station.ro(roit).exL2)==size(GPScoord,2)
                            if sum(station.ro(roit).exL2) > 0 && sum(diff(station.ro(roit).exL2)) > 0
                                idex = find(exL2>30);
                                idex = round(min(idex)./switches.ROres,0);
                                limit = length(posit);
                                step = 1;
                            elseif sum(station.ro(roit).exL2) < 0
                                idex = find(exL2>(min(exL2)+30));
                                idex = round(max(idex)./switches.ROres,0);
                                limit = 1;
                                step = -1;
                            end
                        else
                            warning('spaceRT: Improper structure of .nc radiooccultation file. Ray tracing all available LEO postions')
                            if sum(station.ro(roit).exL2) > 0 && sum(diff(station.ro(roit).exL2)) > 0
                                posit = [1:switches.ROres:size(GPScoord,1)];
                                idex = 1;
                                limit = length(posit);
                                step = 1;
                            elseif sum(station.ro(roit).exL2) < 0
                                posit = [1:switches.ROres:size(GPScoord,1)];
                                idex = size(GPScoord,1);
                                limit = 1;
                                step = -1;
                            end
                        end
                    else
                        % SYNTHETIC mode: exL2 is never used in computation.
                        % Determine ascending/descending via dot product of
                        % LEO->GPS vector with LEO velocity at mid-point.
                        mid_idx = round(length(posit)/2);
                        leo_mid = LEOcoord(posit(mid_idx),:);
                        gps_mid = GPScoord(posit(mid_idx),:);
                        v_leo   = LEOcoord(posit(mid_idx)+1,:) - LEOcoord(posit(mid_idx)-1,:);
                        leo2gps = gps_mid - leo_mid;
                        dot_C   = dot(leo2gps, v_leo);
                        if dot_C < 0
                            % setting/descending
                            idex = 1;
                            limit = length(posit);
                            step = 1;
                        else
                            % rising/ascending
                            idex = length(posit);
                            limit = 1;
                            step = -1;
                        end
                    end
                    %% Prepare data for ray tracing
                    for i = idex:step:limit 
                        %disp(posit(i))
                        fprintf('\n spaceRT: sample index=%d / %d\n', posit(i), posit(end));
                        % Recalculate coordinates to ECI
                        X_Tinit = [GPScoord(posit(i),1),GPScoord(posit(i),2),GPScoord(posit(i),3)]; 
                        X_Rinit = [LEOcoord(posit(i),1),LEOcoord(posit(i),2),LEOcoord(posit(i),3)];
                        llaT = eci2lla(X_Tinit*1000,double(dateRO)); 
                        llaR = eci2lla(X_Rinit*1000,double(dateRO));
                        % Recalculate to rectangular wgs
                        X_T = cspice_georec(llaT(2)*pi()/180, llaT(1)*pi()/180, llaT(3)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                        X_R = cspice_georec(llaR(2)*pi()/180, llaR(1)*pi()/180, llaR(3)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                        %% Calculate excess phase
                        rayRT = voxel_dist_3D_combined(model.refrRT{epoch},0,model.levels_TOMO_RT'*1000,X_T',X_R',model.ds,5,model.rWGS,0,model.LAT,model.LON,model.temp{epoch},model.pres{epoch},model.wvpr{epoch},switches); 
                        
                        fprintf('RO ray: alt_pass=%d, diff_dist=%.4f km, n_iter=%d\n', ...
                            isfield(rayRT,'alt_pass'), ...
                            rayRT.diff_dist(end), ...
                            length(rayRT.diff_dist));   
                        
                        
                        rayRT2(i).rayRT = rayRT;
                        %% Calculate excess phase in tomography model
                        rayRT = interateModelsA(model,rayRT,switches); 
                        %% Calculate signal distances in tomography voxels
                        if sum(rayRT.refr>1) && isfield(rayRT,'voxIN') && ~isfield(rayRT,'alt_pass')
                            rayIND = rayRT.rayIND;
                            ray_num = ray_num + 1;
                            disp(ray_num)
                            % Find coordinates of intersections
                            try
                                coord(ray_num,1).OU = rayRT.voxEM.coord;
                                coord(ray_num,1).IN = rayRT.voxIN.coord;
                            end
                            %% Calculate psuedo excess phase 
                            if strcmp(switches.solution,'REAL')
                                if switches.totalN
                                    ratio1 = rayIND.dLs/rayRT.dL;
                                    exPh(ray_num,1) =  ratio1*exLC(posit(i)); 
                                    hei(ray_num,1) = min(rayRT.de_alt_ray_fin);
                                    dexPh(ray_num,1) = exPh(ray_num,1) - rayIND.dLs;
                                else
                                    ratio1 = rayIND.dLs/rayRT.dL;
                                    ratio2 = rayIND.dLw/rayIND.dLs;
                                    exPh(count,1) = exLC(posit(i)) - rayRT.dLh;
                                    dLhsynth = rayRT.dLh * exLC(posit(i))/rayRT.dL;
                                    exPh1(ray_num,1) = (exLC(posit(i)) - dLhsynth)*ratio1;
                                    exPh2(ray_num,1) =   exLC(posit(i)) *ratio1*ratio2;
                                    exPh(ray_num,1) = (exPh1(count,1) + exPh2(count,1))./2;
                                    dexPh1(ray_num,1) = exPh1(ray_num,1) - rayIND.dLw;
                                    dexPh2(ray_num,1) = exPh2(ray_num,1) - rayIND.dLw;
                                    dexPh(ray_num,1) = (dexPh1(count,1) + dexPh2(count,1))./2;
                                    hei(ray_num,1) = min(rayRT.de_alt_ray_fin);
                                end
                            else
                                if switches.totalN
                                    exPh(ray_num,1) = rayIND.dLs;
                                else
                                    exPh(ray_num,1) = rayIND.dLw;
                                end
                                hei(ray_num,1) = min(rayRT.de_alt_ray_fin);
                                dexPh = [];
                            end
                            %% Add distances traversed by the signal to A matrix
                            RelDist(ray_num,1) = rayRT.voxIN.distanceRel;
                            R_SWD(ray_num,1) = NaN; 
                            iposEMRO1 = [];
                            iposEMRO2 = [];
                            for it = 1:size(rayRT.voxEM.iposEm,1)
                                iposver = rayRT.voxEM.iposEm(it,:);
                                distver = rayRT.voxEM.distanceCorr(it,:);
                                iposver = nonzeros(iposver)';
                                distver = distver(1:size(iposver,2));
                                iposEMRO1 = [iposEMRO1 iposver];
                                iposEMRO2 = [iposEMRO2 distver];  
                                clear iposver distver
                            end
                             iposROEM = [iposEMRO1;iposEMRO2];
                             iposROIN = [rayRT.i_pos;rayRT.voxIN.distanceVox'];
                             ind = [iposROIN iposROEM]'; 
                             ind = sortrows(ind);
                             [row_ind,~,row_ori]=unique(ind(:,1));
                             ind_sumRO =[row_ind cell2mat(accumarray(row_ori,1:numel(row_ori), [],@(x) {sum(ind(x,2:end),1)}))];
                             ind_sumRO = ind_sumRO';  
                             clear row_ind row out
                             % Calculate vector matrices
                             if strcmp(switches.parametrization,'constant')
                                try
                                    Avec(ind_sumRO(1,:),ray_num,:) = coordVector(rayRT.voxIN.coord,rayRT.i_pos);
                                catch
                                    warning('spaceRT: A matrix vector calculation failed') 
                                end
                                A_RO(ind_sumRO(1,:),ray_num) = ind_sumRO(2,:);
                             elseif strcmp(switches.parametrization,'bilinear-h')
                                if isfield(rayRT.rayIND,'X_ray_b')
                                    coordIN = rayRT.rayIND.X_ray_b';
                                    id_par = rayRT.i_pos';
                                    try
                                        coorIntersect = rayRT.voxIN.coord;
                                        coorIntersect = [squeeze(coorIntersect(:,:,1));squeeze(coorIntersect(:,:,2))];
                                        coorIntersect(isnan(coorIntersect(:,1)),:) = [];
                                        parend = coorIntersect(end,:);
                                        parst= coorIntersect(1,:);
                                    catch
                                        parend = [];
                                        parst = [];
                                    end
                                else % if RO outside the model
                                    exPh = 0;
                                    dexPh = [];
                                    RelDist = 0;
                                    R_SWD = 0;
                                    dexPh = 0;
                                    coord = [];
                                    return
                                end
                                %coordINpar = coordIN(sort(ia),:);
                                coordINpar = coordIN(1:1:end,:);
                                nancoord = isnan(coordINpar);
                                if any(any(nancoord)) 
                                    id = isnan(coordINpar(:,1));
                                    coordINpar(id,:) = [];
                                    ind(id,:) = [];
                                    id_par = ind(1,:);
                               else
                                    par = [parst;coordINpar;parend];
                               end
                               for i = 1:size(par,1)
                                    par(i,4) = vecnorm(par(1,1:3)-par(i,1:3));
                               end
                               % Calculate integrals of bilinear interpolation
                               [A_RO,~,~] = bilinearRT(model,A_RO,0,0,nodes_columns,levels,num_lev,Tmatrix,Nw_apr,epoch,ray_num,par,switches);
                               clearvars he la lo altst idst id1 par par1
                            end 
                        end
                        if isfield(rayRT,'alt_pass') 
                            fid = fopen([path_save 'alt_pass_log.txt'],'a');
                            fprintf(fid, 'epoch=%d | sample_index=%d | alt_range=[%.2f %.2f] m | n_iter=%d\n', ...
                                epoch, ...
                                posit(i), ...
                                min(rayRT.de_alt_ray_fin)*1000, ...
                                max(rayRT.de_alt_ray_fin)*1000, ...
                                length(rayRT.diff_dist));
                            fclose(fid);
                            break
                        end  

                    end
                end
            end
        else % if no satellite
           A_RO = [];
           exPh = [];
           RelDist = [];
           R_SWD = [];
           dexPh = [];
           coord = [];
        end
end