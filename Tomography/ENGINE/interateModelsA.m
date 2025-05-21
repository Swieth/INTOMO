function ray = interateModelsA(model,ray,switches)
% Function used to get entries of ray traced signal inside the
% voxel model and calulate slant delay
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       model......   parameters of tomography and ray tracing models
%       ray........   parameters of ray traced signal
%       swtiches...    structural matrix of settings of the INTOMO processing
%%%OUTPUT:        
%       ray.d_voxel.........   ray path in each voxel [m]
%       ray.n_voxel.........   affected voxel [n m]
%       ray.de_lat_ray_fin..   bent ray latitude coordinates [rad]
%       ray.de_lon_ray_fin..   bent ray longitude coordinates [rad]
%       ray.de_alt_ray_fin..   bent ray alttitude coordinates [km]
%       ray.de_lat_ray_fins.   straigt ray latitude coordinates [rad]
%       ray.de_lon_ray_fins.   straight ray longitude coordinates [rad]
%       ray.de_alt_ray_fins.   straight ray alttitude coordinates [km]
%       ray.refr............   refraction index at ray point
%       ray.refh............   hydrostatic refraction index at ray point
%       ray.refw............   wet refraction index at ray point
%       ray.X_ray_b.........   bent ray cooridnates in carthesian
%       ray.X_ray...........   staright ray cooridnates in carthesian
%       ray.diff_dist.......   solution distance to receiver
%       ray.grad_n..........   carthesian gradient in ray point
%       ray.grad_ng.........   spherical gradient in ray point
%       ray.t...............   tangent vector
%       ray.dpa.............   voxel vertical gradient
%       ray.toc.............   time of one iteration calculation
%       ray.h...............   second gradient derivative
%       ray.g...............   vector update for ray position
%       ray.de_b............   correction of inital tangent vector
%       ray.dL..............   final phase delay [m]
%       ray.nstepb..........   bent ray segment length multiplied by refr
%       ray.step............   straight ray segment length
%       ray.stepb...........   bent ray segment length
%       ray.i_pos...........   indexes of voxels containing ray points
%       ray.vox.............   selected meteo values of ray points
%       ray.voxEM...........   see func. vox_distance* 
%       ray.voxIN...........   see func. vox_distance* 

        try
            %% Find ray path fragments intersecting tomography model
            raylon = ray.de_lon_ray_fin.*180/pi;
            id_in  = find(ray.de_lat_ray_fin.*180/pi <= model.lat_TOMO(end) & ray.de_lat_ray_fin.*180/pi >= model.lat_TOMO(1) & raylon <= model.lon_TOMO(end) & raylon >= model.lon_TOMO(1) & ray.de_alt_ray_fin.*1000 <= model.levels_TOMO(end) & ray.de_alt_ray_fin.*1000 >= model.levels_TOMO(1));
            if length(id_in) > 1 && ~switches.refron
                switches.stop = false; %
                rayIND.refr = ray.refr(id_in);
                rayIND.refw = ray.refw(id_in);
                rayIND.refh = ray.refh(id_in);
                rayIND.X_ray_b = ray.X_ray_b(:,id_in); 
                rayIND.X_ray = ray.X_ray(:,id_in); 
                rayIND.de_alt_ray_fin = ray.de_alt_ray_fin(id_in);
                switches.refopt= 1;
                %% Caculate slant delay from parts of signal in tomography model
                [rayIND,dLs] = excessphase(rayIND,switches);
                switches.refopt= 2;
                [rayIND,dLw] = excessphase(rayIND,switches);
                switches.refopt= 3;
                [rayIND,dLh] = excessphase(rayIND,switches);
                rayIND.dLs = dLs; 
                rayIND.dLw = dLw; 
                rayIND.dLh = dLh; 
                ray.rayIND = rayIND;
                %% Find id of voxel intersected by the signal
                for i = 1:length(id_in)
                    [vox(i,:),i_pos(i)] = find_num(ray.de_lat_ray_fin(id_in(i)).*180/pi, raylon(id_in(i)),ray.de_alt_ray_fin(id_in(i))*1000,model.lat_TOMO,model.lon_TOMO,model.levels_TOMO,5); 
                end
                %% Calculate the distance of signal in voxel model
                try
                    [voxIN,voxEM] = vox_distance(model.lat_TOMO, model.lon_TOMO, model.levels_TOMO, raylon,model.radii, id_in, ray,size(id_in,1));
                catch
                    warning('interateModelsA: Failed to calculate signal distances in the tomo model')    
                end
                if size(voxIN.distanceVox,1) > size(ray.refr,1)
                    id_in = id_in(1:size(ray.refr,1));
                    vox = vox(1:size(ray.refr,1),:);
                    i_pos = i_pos(1,1:size(ray.refr,1));
                    warning('interateModelsA: Unequal values of signal distances. Observation might be corrupted')    
                end  
                %% Save variables
                try
                    vox(:,8) = voxIN.distanceVox(1:size(id_in,1));
                    %Additional variables
                    if ~switches.refron
                        vox(:,9) = ray.refr(1:size(id_in,1));
                        vox(:,10) = ray.refh(1:size(id_in,1));
                        vox(:,11) = ray.refw(1:size(id_in,1));
                        vox(:,12) = ray.Temp(1:size(id_in,1));
                        vox(:,13) = ray.pres(1:size(id_in,1));
                        vox(:,14) = ray.hum(1:size(id_in,1));
                    end
                    ray.i_pos = i_pos;
                    ray.vox   = vox;
                    voxIN.id_in = id_in;
                    ray.voxEM = voxEM;
                    ray.voxIN = voxIN;
                catch
                    warning('interateModelsA: Failed to save ray tracing variables... repeating')
                    try
                        vox(:,8) = voxIN.distanceVox(1:size(id_in,1));
                        ray.i_pos = i_pos;
                        voxIN.id_in = id_in;
                        ray.vox   = vox;
                        ray.voxEM = voxEM;
                        ray.voxIN = voxIN; 
                    catch
                        warning('interateModelsA: Failed to save ray tracing variables 2nd time')
                    end
                end
            else
                % Return dummy value '0' for voxel [1 1]
                ray.d_voxel = 0;
                ray.n_voxel = [1 1];
            end
        end
end



