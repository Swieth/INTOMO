function [voxIN,voxEM]  = vox_distance(lat_v, lon_v, alt_v,raylon ,radii, id_in, ray,bound)
% The function vox_distance calculates distance traversed by ray path in each
% voxel. As a results function delivers distances for each 
% ray point as well as coordinates of ray path and voxel model
% intersection.
% vox_distance is based of two 3D lines intersection, first between two
% consecutive ray points, second between ray points projected on voxel
% boundary.

%%%INPUT
% lat_v         ... latidudes voxel model [degrees]
% lon_v         ... longitudes voxel models [degrees]
% alt_p         ... altitude voxel models  [m]
% raylon        ... ray points longitude converted to clockwise positive
% radii         ,.. Earth radius
% id_in         ... id of ray points below atmlimit (deflaut = 86km)
% ray           ... structure to save data

%%%OUTPUT
% voxIN         ... structure for voxels with raypoints
%   distanceCrd ... distance traversed by signal between two ray points
%   distanceRel ... simulated travel raypath length (for comparison)
%   distanceVox ... distance traversed by signal in ray point voxel
%   coord       ... coordinates of projected point
% voxEM         ... structure for voxels without rapoints
%   dist2tar    ... distance to target point
%   idalt       ... index of altitude layer with projected point
%   idlon       ... index of longitude layer with projected point
%   idlat       ... index of latitudelayer with projected point
%   dist2tarCorr... corrected distance to target point
%   distanceCorr... corrected distance traversed by signal in voxel
%   iposEm      ... index of voxel traversed by signal
%   posEm       ... indexes of voxel boundaries
%   coord       ... coordinates of projected point
%   distance    ... distance traversed by signal in voxel
%   idcum       ... index of projected point in reference to ray points



%% In-voxel raypints
coord = zeros(size(id_in,1),3,2);
    for i = 1:size(id_in,1) 
        % Find ray point coordinates
        space = true;
        try
            lat_p1 = ray.de_lat_ray_fin(id_in(i)-1)*180/pi();
            lon_p1 = raylon(id_in(i)-1);
            alt_p1 = ray.de_alt_ray_fin(id_in(i)-1)*1000;
            X1 = ray.X_ray_b(:,id_in(i)-1);
            [pos_X1,~] = find_num(lat_p1,lon_p1,alt_p1,lat_v,lon_v,alt_v,1);
            if length(pos_X1) == 5
                pos_X1(6) = pos_X1(5);
            end
        catch
            space = false;
        end
        lat_p2 = ray.de_lat_ray_fin(id_in(i))*180/pi();
        lon_p2 = raylon(id_in(i));
        alt_p2 = ray.de_alt_ray_fin(id_in(i))*1000;
        [pos_X2,~] = find_num(lat_p2,lon_p2,alt_p2,lat_v,lon_v,alt_v,1);
        X2 = ray.X_ray_b(:,id_in(i));
        if length(pos_X2) == 5
            pos_X2(6) = pos_X2(5);
        end
        lat_p3 = ray.de_lat_ray_fin(id_in(i)+1)*180/pi();
        lon_p3 = raylon(id_in(i)+1);
        alt_p3 = ray.de_alt_ray_fin(id_in(i)+1)*1000;
        % Get index of closest to ray point latitude, longitude, alittude
        [pos_X3,~] = find_num(lat_p3,lon_p3,alt_p3,lat_v,lon_v,alt_v,1);  
        X3 = ray.X_ray_b(:,id_in(i)+1);
        % Fill altitude index if next ray point outside the model
        if length(pos_X3) == 5
            pos_X3(6) = pos_X3(5);
        end     
        % Calculate distance to voxel boundaries
        try
            t_distance(i,1) = vecnorm(X1 - X2);
        end
        t_distance(i,2) = vecnorm(X2 - X3);          
        %% Distance to previus ray point (in ray path direction)
        % Calculate distance if ray points in the same voxel
        if space
            if pos_X1(3) == pos_X2(3) && pos_X1(4) == pos_X2(4) && pos_X1(1)  == pos_X2(1) && pos_X1(2) == pos_X2(2) && pos_X1(5) == pos_X2(5) && pos_X1(6) == pos_X2(6)
                distance(i,1) = vecnorm(X1 - X2);
                coord(i,:,1) = NaN;
            %If not in the same voxel
            else 
                %Select layer index, project in altitude direction
                if alt_p1 <= alt_p2
                    idalt = min(pos_X2(5),pos_X2(6));
                else
                    idalt = max(pos_X2(5),pos_X2(6));
                end
                [X_T1] = cspice_georec(lon_p2*pi()/180, lat_p2*pi()/180, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_T2] = cspice_georec(lon_p1*pi()/180, lat_p1*pi()/180, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_alt,~] = lineIntersect3D([X_T1';X2'],[X_T2';X1']);
                % Calculate distance
                distAlt = vecnorm(X_alt' - X2); 
                 %Select layer index, project in longitude direction
                if lon_p1 <= lon_p2
                    idlon = min(pos_X2(3),pos_X2(4));
                else
                    idlon = max(pos_X2(3),pos_X2(4));
                end
                [X_T1] = cspice_georec(lon_v(idlon)*pi()/180, lat_p2*pi()/180, alt_p2/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_T2] = cspice_georec(lon_v(idlon)*pi()/180, lat_p1*pi()/180, alt_p1/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_lon,~] = lineIntersect3D([X_T1';X2'],[X_T2';X1']);  
                 % Calculate distance
                distLon = vecnorm(X_lon' - X2); 
                %Select layer index, project in latitude direction
                if lat_p1 <= lat_p2
                    idlat = min(pos_X2(1),pos_X2(2));
                else
                    idlat = max(pos_X2(1),pos_X2(2));
                end
                [X_T1] = cspice_georec(lon_p2*pi()/180, lat_v(idlat)*pi()/180, alt_p2/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_T2] = cspice_georec(lon_p1*pi()/180, lat_v(idlat)*pi()/180, alt_p1/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_lat,~] = lineIntersect3D([X_T1';X2'],[X_T2';X1']); 
                % Calculate distance
                distLat = vecnorm(X_lat' - X2); 
                % Select coordinates of intersection points and distance by
                % shortest path to projected point
                coordp = [X_alt;X_lon;X_lat];
                [distance(i,1),ind] = min([distAlt;distLon;distLat]);
                try
                    coord(i,:,1) = coordp(ind,:);
                catch
                    disp(i)
                end
              
            end
        end
        %%  Distance to next ray point (in ray path direction)
        % Calculate distance if ray points in the same voxel
        if pos_X3(3) == pos_X2(3) && pos_X3(4) == pos_X2(4) && pos_X3(1)  == pos_X2(1) && pos_X3(2) == pos_X2(2) && pos_X3(5) == pos_X2(5) && pos_X3(6) == pos_X2(6)
            distance(i,2) = 0; %Calculated in next iteration
            coord(i,:,2) = NaN;
        %If not in the same voxel
        else
            %Select layer index, project in altitude direction
            if alt_p2 >= alt_p3
                idalt = min(pos_X2(5),pos_X2(6));
            else
                idalt = max(pos_X2(5),pos_X2(6));
            end
            [X_T2] = cspice_georec(lon_p2*pi()/180, lat_p2*pi()/180, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
            [X_T1] = cspice_georec(lon_p3*pi()/180, lat_p3*pi()/180, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
            [X_alt,~] = lineIntersect3D([X_T1';X3'],[X_T2';X2']);
            % Calculate distance
            distAlt = vecnorm(X_alt' - X2); 
            %Select layer index, project longitude direction
            if lon_p2 >= lon_p3
                idlon = min(pos_X2(3),pos_X2(4));
            else
                idlon = max(pos_X2(3),pos_X2(4));
            end
            [X_T2] = cspice_georec(lon_v(idlon)*pi()/180, lat_p2*pi()/180, alt_p2/1000, radii(1),(radii(1)-radii(2))/radii(1));
            [X_T1] = cspice_georec(lon_v(idlon)*pi()/180, lat_p3*pi()/180, alt_p3/1000, radii(1),(radii(1)-radii(2))/radii(1));
            [X_lon,~] = lineIntersect3D([X_T1';X3'],[X_T2';X2']);
            % Calculate distance
            distLon = vecnorm(X_lon' - X2);
            %Select layer index, project in latitude direction
            if lat_p2 >= lat_p3
                idlat = min(pos_X2(1),pos_X2(2));
            else
                idlat = max(pos_X2(1),pos_X2(2));
            end
            [X_T2] = cspice_georec(lon_p2*pi()/180, lat_v(idlat)*pi()/180, alt_p2/1000, radii(1),(radii(1)-radii(2))/radii(1));
            [X_T1] = cspice_georec(lon_p3*pi()/180, lat_v(idlat)*pi()/180, alt_p3/1000, radii(1),(radii(1)-radii(2))/radii(1));
            [X_lat,~] = lineIntersect3D([X_T1';X2'],[X_T2';X3']);
            %Calculate distance
            distLat = vecnorm(X_lat' - X2);
            % Select coordinates of intersection points and distance by
            % shortest path to projected point
            [distance(i,2),ind] = min([distAlt;distLon;distLat]);
            coordp = [X_alt;X_lon;X_lat];
            coord(i,:,2) = coordp(ind,:);  
            clearvars -except lat_v lon_v alt_v bound radii id_in ray raylon distance coord t_distance voxIN
        end
    end
    % Save variables
    id = find(distance(1:bound,1) > t_distance(1:bound,1));
    distance(id,1) = t_distance(id,1);
    id = find(distance(1:bound,2) > t_distance(1:bound,2));
    distance(id,2) = t_distance(id,1);
    distLon = [distance(:,1);0];
    distLat = [0;distance(:,2)];
    voxIN.distanceCrd = sum([distLon(2:end),distLat(2:end)],2);
    voxIN.distanceVox = sum(distance,2);
    voxIN.distanceRel = sum(t_distance(1:bound,1))+voxIN.distanceCrd(end);
    voxIN.coord = coord;
    %% Interpolation for voxels without raypoint
    %Check distance between two voxel boundary-ray intersection points and save cooridnates of
    %voxel boundary-ray intersection points
    for i = 1:size(coord,1)-1
        dist(i) = vecnorm(squeeze(coord(i,:,2))-squeeze(coord(i+1,:,1)));
        coordEm(i,:,1) = coord(i,:,2);
        coordEm(i,:,2) = coord(i+1,:,1);
    end
    %Check id of voxel boundary-ray intersection points
    idcum = intersect(find(~isnan(coordEm(:,1,1))),find(dist >0.01)');
    id = find(isnan(coordEm(:,1,1)));
    %Delete NaNs and zeros (for two or more voxel boundary-ray intersection points in one voxel)
    coordEm(id,:,:) = [];
    dist(id) = [];
    id = find(dist <0.01);
    coordEm(id,:,:) = [];
    if size(coordEm,1)  > 1
        %Iterate each pair of voxel boundary-ray intersection points
        for i = 1:size(coordEm,1)-1
            %Counter for iteration
            b = 1;
            %Coordinates of a pair of voxel boundary-ray intersection points
            X_HM = squeeze(coordEm(i,:,1));
            X2 = squeeze(coordEm(i,:,2));
            [lon_ray1,lat_ray1,alt_ray1] = cspice_recgeo(X_HM',radii(1),(radii(1)-radii(2))/radii(1)); 
            [lon_ray2,lat_ray2,alt_ray2] = cspice_recgeo(X2',radii(1),(radii(1)-radii(2))/radii(1)); 
            % Find closest voxel boundaries for target point
            [pos_X2,~] = find_num(lat_ray2*180/pi(),lon_ray2*180/pi(),alt_ray2*1000,lat_v,lon_v,alt_v,1);
            % Check direction of lat/lon/alt increament
            if lon_ray1 < lon_ray2
                switches.direction(1) = true;
            else
                switches.direction(1) = false;
            end
            if lat_ray1 < lat_ray2
                switches.direction(2) = true;
            else
                switches.direction(2) = false;
            end
            if alt_ray1 < alt_ray2
                switches.direction(3) = true;
            else
                switches.direction(3) = false;
            end 
            % Set/Reset change of lat/lon/alt index
            switches.change = [false;false;false];
             % Stop iterating while hit the target
            while vecnorm(X_HM - X2) > 0.01 
                % Stop iprocessing after 100 iterations (ray missed due to projection innacuracies)
                if b == 100 
                    % Cut vectors to minimum distance to target
                    [voxEM,distanceEm,i_pos,pos]  = vox_distance_cut(voxEM,i,distanceEm,i_pos,pos);
                    break
                end
                X1 = X_HM;
                % Find closest voxel boundaries for iterated point
                [lon_ray1,lat_ray1,alt_ray1] = cspice_recgeo(X1',radii(1),(radii(1)-radii(2))/radii(1)); 
                [pos_X1,~] = find_num(lat_ray1*180/pi(),lon_ray1*180/pi(),alt_ray1*1000,lat_v,lon_v,alt_v,1);
                % Calculate id's of voxel boundaries
                if b == 1
                    [idlon,idlat,idalt,switches] = select_id(switches,pos_X1,0,0,0,b);
                else
                    [idlon,idlat,idalt,switches] = select_id(switches,pos_X1,idlon,idlat,idalt,b);               
                end
                % Ray points in the same voxel (only in case innacuracies of projection)
                if pos_X1(3) == pos_X2(3) && pos_X1(4) == pos_X2(4) && pos_X1(1)  == pos_X2(1) && pos_X1(2) == pos_X2(2) && pos_X1(5) == pos_X2(5) && pos_X1(6) == pos_X2(6)
                    voxEM.dist2tar(i,b) = vecnorm(X_HM - X2);
                    distanceEm(i,b) = vecnorm(X_HM - X2);
                    coord(i,:,2) = NaN;
                    [pos(i,:,b),i_pos(i,b)] = find_num(lat_ray1*180/pi(),lon_ray1*180/pi(),alt_ray1*1000,lat_v,lon_v,alt_v,0);
                    break
                end
                %Project point to selected longitude
                [X_M1] = cspice_georec(lon_v(idlon)*pi()/180, lat_ray1, alt_ray1, radii(1),(radii(1)-radii(2))/radii(1));
                [X_M2] = cspice_georec(lon_v(idlon)*pi()/180, lat_ray2, alt_ray2, radii(1),(radii(1)-radii(2))/radii(1));
                [X_lon,~] = lineIntersect3D([X_M2';X2],[X_M1';X1]); 
                distLon = vecnorm(X_lon - X1); 
                %Project point to selected latitude
                [X_M1] = cspice_georec(lon_ray1, lat_v(idlat)*pi()/180, alt_ray1, radii(1),(radii(1)-radii(2))/radii(1));
                [X_M2] = cspice_georec(lon_ray2, lat_v(idlat)*pi()/180, alt_ray2, radii(1),(radii(1)-radii(2))/radii(1));
                [X_lat,~] = lineIntersect3D([X_M2';X2],[X_M1';X1]);
                distLat = vecnorm(X_lat - X1);
                %Project point to selected altitude
                [X_M1] = cspice_georec(lon_ray1, lat_ray1, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_M2] = cspice_georec(lon_ray2, lat_ray2, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_alt,~] = lineIntersect3D([X_M2';X2],[X_M1';X1]);
                distAlt = vecnorm(X_alt - X1); 
                %Select shortest distance
                coordp = [X_lon;X_lat;X_alt];
                [distanceEm(i,b),ind] = min([distLon;distLat;distAlt]);
                coordEM(i,:,b) = coordp(ind,:);
                X_HM = coordEM(i,:,b);
                switches.change(ind) = true;  
                %Calculate voxel id
                halfdist = distanceEm(i,b)/2;
                eX = (X_HM-X1)./(norm(X_HM-X1));
                Xm = X_HM+eX*halfdist;
                [lon_rayM,lat_rayM,alt_rayM] = cspice_recgeo(Xm',radii(1),(radii(1)-radii(2))/radii(1));
                % Stop iterating if missed the target (only if signal goes out of model)
                try
                    [pos(i,:,b),i_pos(i,b)] = find_num(lat_rayM*180/pi(),lon_rayM*180/pi(),alt_rayM*1000,lat_v,lon_v,alt_v,0);
                catch
                    % Cut vectors to minimum distance to target
                    [voxEM,distanceEm,i_pos,pos]  = vox_distance_cut(voxEM,i,distanceEm,i_pos,pos);
                    break
                end
                % Save results
                voxEM.dist2tar(i,b) = vecnorm(X_HM - X2);
                voxEM.id_alt(i,b) = idalt;
                try
                    voxEM.id_lon(i,b) = idlon;
                    voxEM.id_lat(i,b) = idlat;
                catch
                    voxEM.id_lon(i,b) = NaN;
                    voxEM.id_lat(i,b) = NaN;       
                end
                 % Stop iterating if missed the target
                if b>2
                    if voxEM.dist2tar(i,b-1) - voxEM.dist2tar(i,b) < -1 && abs(voxEM.dist2tar(i,b-1)) < 2
                        % Cut vectors to minimum distance to target
                        [voxEM,distanceEm,i_pos,pos]  = vox_distance_cut(voxEM,i,distanceEm,i_pos,pos);
                        break
                    end
                end
                % Next iteration
                b = b + 1;
                clear Xm eX lon_rayM lat_rayM alt_rayM halfdist coordp distLon distLat distAlt X_alt X_lon X_alt distances
            end
            clear dlat1 dlat2 dlon1 dlon2 dalt1 dalt2 n_alt1_1 n_alt1_2 n_alt2_1 n_alt2_2 n_lon1_1 n_lon1_2 n_lon2_1 n_lon2_2 n_lat1_1 n_lat1_2 n_lat2_1 n_lat2_2
            %Check and correct distances
            try
                dX2 = diff(voxEM.dist2tar(i,1:end));
                iddX2 = find(dX2>0);
                dx2val = zeros(length(dX2)+1,1);
                dx2val(iddX2+1) = dX2(iddX2); 
                voxEM.dist2tarCorr{i,:} =  voxEM.dist2tar(i,1:size(dx2val,1)) - dx2val';
                voxEM.distanceCorr{i,:}  = distanceEm(i,1:size(dx2val,1)) - dx2val';
            catch
                disp('Correction of inner voxels ray point failed')
                voxEM.missed(i) = 999;
            end
        end
        % Turn cell array to matrix
         sizeCell = max(cellfun('size',voxEM.dist2tarCorr,2));
         dist2tarCorr = zeros(size(i_pos,1),sizeCell);
         distanceEmCorr = zeros(size(i_pos,1),sizeCell);
         for i = 1:size(voxEM.dist2tarCorr,1)
            cell1 = voxEM.dist2tarCorr{i};
            cell2 = voxEM.distanceCorr{i};
            dist2tarCorr(i,1:size(cell1,2)) = dist2tarCorr(i,1:size(cell1)) + cell1;
            distanceEmCorr(i,1:size(cell1,2)) = distanceEmCorr(i,1:size(cell1)) + cell2;
         end
         % Save results
         voxEM.dist2tarCorr = dist2tarCorr;
         voxEM.distanceCorr = distanceEmCorr;
         voxEM.iposEm = i_pos;
         voxEM.posEm = pos;
         voxEM.distance = distanceEm;
         voxEM.idcum = idcum;
         % There can be no additional points on boundaries
         try
             voxEM.coord = coordEM;
         end
    elseif size(coordEm,1) == 1 % Valid only for rare densified, tomography solution   
         voxEM.distanceCorr = vecnorm(squeeze(coordEm(1,:,2))-squeeze(coordEm(1,:,1)));
         %Calculate voxel id
         halfdist = voxEM.distanceCorr/2;
         eX = (coordEm(1,:,2)-coordEm(1,:,1))./(norm(coordEm(1,:,2)-coordEm(1,:,1)));
         Xm = coordEm(1,:,2)+eX*halfdist;
         [lon_rayM,lat_rayM,alt_rayM] = cspice_recgeo(Xm',radii(1),(radii(1)-radii(2))/radii(1));
         [~,i_pos] = find_num(lat_rayM*180/pi(),lon_rayM*180/pi(),alt_rayM*1000,lat_v,lon_v,alt_v,0);
         voxEM.iposEm = i_pos;           
    else
        voxEM.coord = [];
        voxEM.distanceCorr = [];
        voxEM.iposEm = [];
    end

end 
