function ray =  voxel_dist_3D_combined(ref_a,ref_w,alt_N,X_T,X_R,ds,dt,rwgs,elev,LAT,LON,Temp,pLevel,WVpres,set)
% Function1 voxel_dist_3D reconstructs the 3D ray path through the voxel model 
% given in lat,lon,height. It is prepared for space and ground -
% based observation. Radio occultation has to start above model maximum
% altitude.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
%---------------------------------------------------------------------------------------------------
%%%Processing modes:
%%Simulation method
% set.refron(true) - quicker simualtion based on refractivity values interpolated from
% grid model / set.refron(false) - simualtion based on temp/pres/wvpres values interpolated from
% grid model and recalculated to refractivity. 
%%Tomography
% set.integrated - simulation for tomography use, if true, calculates no
% voxIN, voxEM (see below)
%
%%%INPUT:
%       ref_a.....   refractivity [ppm], matrix containing all levels
%       ref_w.....   wet refractivity [ppm], matrix containing all levels
%       alt_N.....   ellipsoidal height of refractivity [m]
%       X_T.......   cartesian transmitter coordinates at time of signal transmission [km]
%       X_R.......   cartesian receiver coordinates at time of signal reception [km]
%       ds........   step size, length of one ray segment [km]
%       dt........   threshold which defines a 'close' solution [km]
%       LAT.......   ellipsoidal latitude of refractivity [degrees]
%       LON.......   ellipsoidal longitude of refractivity [degrees]
%       Temp......   temperature [K], matrix containing all levels
%       pLevel....   preasure [hPa], matrix containing all levels
%       WVpres....   water vapour preassure [kg kg-1], matrix containing all levels
%       set.......   settings
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
%       ray.voxEM...........   strucutre for voxels without raypoints but traversed by raypath 
%       ray.voxIN...........   strucutre for voxels with raypoints
if ~isnan(X_R)
    %% Voxel model boundaries
    [lat_v,lon_v,alt_v] = voxel_bound(LAT,LON,alt_N);

    %% Initial parameters for ray-tracing
    % Initial ray position in (x,y,z)
    X_ray(:,1)   = X_T;
    X_ray_b(:,1) = X_T; 
    % Initial settings for ray-tracing
    diff_dist = 1000;        % [km] initial value for closest ray point
    n_iter    = 1;           % Counter for iterations
    n_ray     = 2;           % Counter for ray positions
    n_iterend = 30;          % Deflaut number of iterations
    hlim      = max(alt_N)/1000; % Max altitude to limit refraction calculation
    % Work variables
    inp_data  = [];          % Structure to save work variables
    dh3       = 0;           % Distance between nodal points in colatitude direction
    dh4       = 0;           % Distance between nodal points in longitude direction
    g         = [0;0;0];     % Vector update for ray position
    refr      = 0;           % Vector to store refractivity values
    grad_n    = [0 0 0]';    % Matrix to save cartesian gradient values
    grad_ng   = [0 0 0]';    % Matrix to save spherical gradient values
    ts        = [0 0 0]';    % Matrix to save raypath unit_vectors
    ray       = [];          % Structure to save results
    t0s       = [0 0 0]';    % Matrix to save initital unit_vectors after next iterations
    refrh     = 0;           % hydrostatic refractivity
    refrw     = 0;           % wet refractivity
    Tempprof  = 0;           % temperature along the ray path
    presprof  = 0;           % presure along the ray path
    humprof   = 0;           % humidity along the ray path
    id        = 0;           % ray point id
    ray.error = [];          % vector to collect errors
    k         = 0;
    %% Get radii of reference ellipsoid and spherical coordinates
    radii = [6378.137, 6356.752314245]; 
    % Initial ray position (lon,lat,alt) at time of transmission
    [lon_ray(1),lat_ray(1),alt_ray(1)]       = cspice_recgeo(X_ray(:,1)  ,radii(1),(radii(1)-radii(2))/radii(1));
    %[~,lat_ray2(1),lon_ray(1)] = cspice_recsph(X_ray(:,1));
    [lon_ray_b(1),lat_ray_b(1),alt_ray_b(1)] = cspice_recgeo(X_ray_b(:,1),radii(1),(radii(1)-radii(2))/radii(1));    
    % Targeting ray position (lon,lat,alt), i.e. position of the receiver
    [~,~,alt_end] = cspice_recgeo(X_R',radii(1),(radii(1)-radii(2))/radii(1));
    %% Additional settings
    if alt_ray(1) < hlim
        set.ground = true;
    else
        set.ground = false;
    end
    set.stop = false;

    %% Calculate initial unit vector
    if  ~set.ground    
        eX   = (X_R-X_T)'./(norm(X_R-X_T)); 
        eX_b = (X_R-X_T)'./(norm(X_R-X_T));
        t0   = eX_b;
        t    = eX_b;
    elseif elev*180/pi() < 90 && set.ground 
        eX   = (X_R-X_T)'./(norm(X_R-X_T)); 
        eX_b = (X_R-X_T)'./(norm(X_R-X_T));
        % New step size based on tranmitter receiver position
        dsL = vecnorm(X_R-X_T); 
        ds =  dsL/200000; 
        % Change number of maximum iterations to speed up
        n_iterend = 5;
        % Refraction value at ground station position
        [~, ~, NhRay, NwRay,  ~, ~,  ~, ~, ~] = statprofile(LAT,LON, pLevel, Temp, WVpres, alt_N, [lat_ray_b(1) lon_ray_b(1) alt_ray_b(1)*1000], rwgs);
        Nst = NhRay(1) + NwRay(1);
        % Elevation angle corretcion
        eaCor = 1e-7*Nst*cot(elev); %6
        eaObs = elev - eaCor; 
        % Next ray point position to determine azimuth based on initial unit
        % vector
        X_T2 = X_T + ds*eX_b';
        [~, azimt1,~, ~] = fxyz2ea(X_T, X_T2, 'e');
        wgs84 = wgs84Ellipsoid2;
        % New ray point position to determine new unit vector
        [lat,lon,h] = aer2geodetic2(azimt1*180/pi(),eaObs*180/pi(),(dsL/100000)*1000,lat_ray_b(1)*180/pi(),lon_ray_b(1)*180/pi(),alt_ray_b(1)*1000,wgs84);
        [X_T2bis] = cspice_georec( lon*pi()/180, lat*pi()/180, h/1000, radii(1),(radii(1)-radii(2))/radii(1));
        % New unit_vector
        eX   = (X_T2bis'-X_T)'./(norm(X_T2bis'-X_T)); 
        eX_b = (X_T2bis'-X_T)'./(norm(X_T2bis'-X_T));
        t0   = eX_b;
        t    = eX_b;
    end

    %% Ray-tracing loop
    while diff_dist > dt
       % Compute excess phase and break loop if number of iterations equals k      
       if n_iter > n_iterend     
            ray.alt_pass = 1;
            ray.error = [ray.error 100];
            % Compute excessphase for refraction and meteo-parameters solution
             try
                if set.refron
                    set.refopt= 1;
                    [ray,dLs] = excessphase(ray,set);
                    ray.dL = dLs;
                else
                    set.refopt= 1;
                    [ray,dLs] = excessphase(ray,set);
                    set.refopt= 2;
                    [ray,dLw] = excessphase(ray,set);
                    set.refopt= 3;
                    [ray,dLh] = excessphase(ray,set);
                    ray.dL = dLs;
                    ray.dLw = dLw;
                    ray.dLh = dLh;
                    %fprintf('%9.4f %9.6f\n',ray.altend,ray.dLw)
                end
            catch
                warning('voxel_dist_3D_combined: Failed to calcuate refraction')
                ray.error = [ray.error 400];
            end
            %Break loop
            set.stop = true;
            break
       end
       % Set new step size for ground based solution
       if set.ground && alt_ray(n_ray-1) > hlim
          ds = dsL/4000;
       %elseif set.ground && alt_ray(n_ray-1) < hlim*0.05 
         % ds = dsL/200000; %50000
       elseif set.ground && alt_ray(n_ray-1) < hlim 
          ds = dsL/50000; %50000
       end  
        % Straight-line ray-tracing
        X_ray(:,n_ray) = X_ray(:,n_ray-1) + ds*eX;
        % Convert straight ray coordinates from x,y,z to lat,lon,alt
        [lon_ray(n_ray),lat_ray(n_ray),alt_ray(n_ray)] = cspice_recgeo(X_ray(:,n_ray),radii(1),(radii(1)-radii(2))/radii(1)); 

        %% Step-wise processing outside the atmosphere for approaching rays
        if alt_ray(n_ray) > hlim && (alt_ray(n_ray)-alt_ray(n_ray-1)) <= 0 % [km]
            % New ray position
            X_ray_b(:,n_ray) = X_ray_b(:,n_ray-1) + ds*t0; 
            % Convert bended ray coordinates from x,y,z to lat,lon,alt
            [lon_ray_b(n_ray),lat_ray_b(n_ray),alt_ray_b(n_ray)] = cspice_recgeo(X_ray_b(:,n_ray),radii(1),(radii(1)-radii(2))/radii(1));       
        end 
        dpa(n_ray-1) = 0;
        %% Step-wise processing inside the atmosphere
        if alt_ray(n_ray) <= hlim 
           % Bended ray-tracing
           % Compute refractive index for current ray position
           if set.refron
                [N,inp_data] = refr_interp_3D(ref_a,ref_w,LAT*180/pi,LON*180/pi,alt_N,(pi()/2- lat_ray_b(n_ray-1))*180/pi,lon_ray_b(n_ray-1)*180/pi,alt_ray_b(n_ray-1)*1000,hlim*1000,Temp,WVpres,pLevel,rwgs,set,inp_data); % [ppm]
                N0 = N;
                refrw(n_ray-1) = 1 + N*10^-6;
           else
                [N,inp_data] = refr_interp_3D(ref_a,ref_w,LAT*180/pi,LON*180/pi,alt_N,(pi()/2- lat_ray_b(n_ray-1))*180/pi,lon_ray_b(n_ray-1)*180/pi,alt_ray_b(n_ray-1)*1000,hlim*1000,Temp,WVpres,pLevel,rwgs,set,inp_data); % [ppm]
                if isstruct(N)
                    Nh  = N.Nh;
                    Nw  = N.Nw;
                    try
                        Tempprof(n_ray-1) = N.T;
                        presprof(n_ray-1) = N.p;
                        humprof(n_ray-1) = N.h;
                    end
                else
                    Nw = 0;
                    Nh = 0;
                    Tempprof(n_ray-1) = 0;
                    presprof(n_ray-1) = 0;
                    humprof(n_ray-1) = 0;
                end
                refrw(n_ray-1) = 1 + Nw*10^-6; 
                refrh(n_ray-1) = 1 + Nh*10^-6;
                N0 = Nh + Nw;
           end
              N0stack(n_ray-1) = N0;
              refr(n_ray-1) = 1 + N0stack(n_ray-1)*10^-6 ;
            % Loop over coordinate components
            if N0 == 0
                grad_n(:,n_ray-1) = [0 0 0]';
            else
                % Compute gradient components in x,y,z directions
                 %Processing in spherical coordinates
                [lon,lat,alt] = cspice_recgeo(X_ray_b(:,n_ray-1),radii(1),(radii(1)-radii(2))/radii(1));           
                %Find closest GRID layers below and above ray point 
                [m1, ~] = find(alt_N > alt*1000);
                [m2, ~] = find(alt_N < alt*1000);                       
                %Vertical [alt] gradient interpolation       
                [dref,dref2,dh,dh2,m1,m2] = dp_alt(m1,m2,alt_N,inp_data);
                [dpa,dp] = dpacalc(dh,dh2,dref,dref2,n_ray,dpa,alt,alt_N,m2);
                grad_n(3,n_ray-1) = dp*10^-3; 
                %Horizontal [lon,lat]  gradient interpolation2
                if n_ray>2
                   [drefm1,drefm2,dh3,dh4] = dhrefcalc(pi()/2 - lat,lon,m1,m2,dh3,dh4,lon_ray_b(n_ray-2),pi()/2- lat_ray_b(n_ray-2),set,inp_data);
                else
                   [drefm1,drefm2,dh3,dh4] = dhrefcalc(pi()/2 - lat,lon,m1,m2,dh3,dh4,lon_ray_b(n_ray-1),pi()/2- lat_ray_b(n_ray-1),set,inp_data);
                end
                grad_n(1,n_ray-1) = (drefm1/dh3)*10^-3;      
                grad_n(2,n_ray-1) = (drefm2/dh4)*10^-3; 
                grad_ng(:,n_ray-1) = grad_n(:,n_ray-1);
                %Recalculation to cartesian
                grad_n(:,n_ray-1) = gradrec([lon,pi()/2 - lat,alt],grad_n(:,n_ray-1));              
            end        
            % Compute 2nd derivative of time series expansion of x(s)
            h = (grad_n(:,n_ray-1) - dot(grad_n(:,n_ray-1),t).*t)./refr(n_ray-1);
            if norm(h) > 10^-20
                % Calculate ray position update
                g = t + h*ds/2;           
                % Unit vector of ray position update
                g = g/norm(g);            
                % Update ray position vector
                X_ray_b(:,n_ray) = X_ray_b(:,n_ray-1) + g*ds;            
                % Update the ray tangent vector
                t = t + h*ds;           
                % Unit tangent vector
                t = t/norm(t);
            else
                % New ray position
                X_ray_b(:,n_ray) = X_ray_b(:,n_ray-1) + ds*t;
            end
            ts(:,n_ray) = t;
            ray.h(:,n_ray-1) = h; 
            ray.g(:,n_ray-1) = g; 
            % Convert bended ray coordinates from x,y,z to lat,lon,alt
            [lon_ray_b(n_ray),lat_ray_b(n_ray),alt_ray_b(n_ray)] = cspice_recgeo(X_ray_b(:,n_ray),radii(1),(radii(1)-radii(2))/radii(1));       
        end
        %% Stop processing once the signal is below the receiver
        if alt_ray(n_ray) > hlim && (alt_ray(n_ray)-alt_ray(n_ray-1)) > 0 %
            % Convert bended ray coordinates from x,y,z to lat,lon,alt
             X_ray_b(:,n_ray) = X_ray_b(:,n_ray-1) + ds*eX_b; 
            % Convert bended ray coordinates from x,y,z to lat,lon,alt
             [lon_ray_b(n_ray),lat_ray_b(n_ray),alt_ray_b(n_ray)] = cspice_recgeo(X_ray_b(:,n_ray),radii(1),(radii(1)-radii(2))/radii(1));
            %% Stop processing once the signal is below the receiver 
            if alt_ray_b(end) > alt_end  + dt   
                % Find the closest point to the receiver
                [~,id] = min(vecnorm(X_ray(:,end-1)-X_ray_b,1));
                diff_dist  = vecnorm(X_ray(:,end-1)-X_ray_b(:,id),1);               
                %Store data if n_iter equals 30
                %disp(diff_dist)
                %% Compute correction to initial vector as long the the reconstructed signal is outside the given thresholds.
                %Save variables 
                n_iter = n_iter + 1;
                if ~set.refron
                    ray = savevar(lat_ray_b,lon_ray_b,alt_ray_b,lat_ray,lon_ray,alt_ray,refr,X_ray,X_ray_b,diff_dist,grad_n,grad_ng,ts,ds,ray,n_iter,refrh,refrw,Tempprof,presprof,humprof,id);
                else
                    ray = savevar(lat_ray_b,lon_ray_b,alt_ray_b,lat_ray,lon_ray,alt_ray,refr,X_ray,X_ray_b,diff_dist,grad_n,grad_ng,ts,ds,ray,n_iter,0,0,0,0,0,0);
                end
                if diff_dist > dt % [km]
                    % Corrections of initial vector
                    dX_ray_b = X_ray_b(:,id) - X_ray_b(:,1);
                    % Vector between transmitter and receiver
                    dX_ray   = X_ray(:,end-1)   - X_ray(:,1);
                    de_b(:,n_iter+1) = (dX_ray_b./norm(dX_ray_b)-dX_ray./norm(dX_ray));
                    % New initial unit vector
                    if set.ground
                        t  = t0 - de_b(:,n_iter+1)*100; 
                    else
                       if isfield(ray,'diff_dist') 
                          if n_iter > 3 && ray.diff_dist(n_iter-2) < ray.diff_dist(n_iter-1)
                            k = k+1;
                            t  = t0 - de_b(:,n_iter+1)./(1*k); 
                          elseif k == 0
                            t  = t0 - de_b(:,n_iter+1)./1; 
                          else
                            t  = t0 - de_b(:,n_iter+1)./(1*k); 
                          end
                        else
                         t  = t0 - de_b(:,n_iter+1)./1; 
                        end
                    end
                    t  = t./norm(t);
                    t0 = t;
                    t0s(:,n_iter-1) = t0;             
                    % Control: Check elevation angles after each iteration
                    %[eath(n_iter), azimt1(n_iter), ranget1, ENU1] = fxyz2ea(X_ray_b(:,1)', X_ray_b(:,2)', 'e');
                    % Store variables for next iteration except for final iteration
                    n_ray = 1;
                    %Save variables
                    if ~set.refron
                        ray = savevar(lat_ray_b,lon_ray_b,alt_ray_b,lat_ray,lon_ray,alt_ray,refr,X_ray,X_ray_b,diff_dist,grad_n,grad_ng,ts,ds,ray,n_iter,refrh,refrw,Tempprof,presprof,humprof,id);
                    else
                        ray = savevar(lat_ray_b,lon_ray_b,alt_ray_b,lat_ray,lon_ray,alt_ray,refr,X_ray,X_ray_b,diff_dist,grad_n,grad_ng,ts,ds,ray,n_iter,0,refrw,0,0,0,id);
                    end
                    % Clear variables for next iteration except for final iteration
                    X_ray(:,2:end)   = [];
                    X_ray_b(:,2:end) = [];
                    lon_ray(2:end)   = [];
                    lat_ray(2:end)   = [];
                    alt_ray(2:end)   = [];
                    lon_ray_b(2:end) = [];
                    lat_ray_b(2:end) = [];
                    alt_ray_b(2:end) = [];
                    refr             = [];
                    dpa              = [];
                    dp               = [];
                    ts               = [];
                    refrw        = [];
                    if ~set.refron
                        refrh        = [];
                    end
                % Stop processing. The signal path has been reconstructed
                % within the given thresholds.
                ray.toc(n_iter) = toc;
                else
                    ray.toc(n_iter) = toc;
                end           
            end
        end     
        % Increase counter
        n_ray = n_ray + 1;   
    end
    % Check if bent ray path crossing the Earth surface
    if min(alt_ray_b) < 0
        ray.alt_pass = 1;
    end
    if sum(ray.refr>1)
        %% Atmospheric delay dL [m]
        try
            if set.refron && ~set.stop
                set.refopt= 1;
                [ray,dLs] = excessphase(ray,set);
                ray.dL = dLs;
            elseif ~set.refron && ~set.stop
                set.refopt= 1;
                [ray,dLs] = excessphase(ray,set);
                set.refopt= 2;
                [ray,dLw] = excessphase(ray,set);
                set.refopt= 3;
                [ray,dLh] = excessphase(ray,set);
                ray.dL = dLs;
                ray.dLw = dLw ;
                ray.dLh = dLh;
                %fprintf('%9.4f %9.6f\n',ray.altend,ray.dLw)
            end
        catch
            warning('voxel_dist_3D_combined: Failed to calcuate refraction')
            ray.error = [ray.error 400];
        end 
        ray.t0 = t0s;

        if ~strcmp(set.integrated,'yes')
            %% Path lengths in each traversed voxel
            % Get entries inside the voxel model
            try
                if min(ray.de_lon_ray_fin.*180/pi)<0
                    id1 = find(ray.de_lon_ray_fin.*180/pi<0);
                    id2 = find(ray.de_lon_ray_fin.*180/pi>0);
                    raylon = ray.de_lon_ray_fin(id2).*180/pi;
                    if id1(end) < id2(1)
                        raylon = [360+ray.de_lon_ray_fin(id1).*180/pi;raylon];
                    else
                       id1(end) > id2(1);
                       raylon = [raylon;360+ray.de_lon_ray_fin(id1).*180/pi];
                    end
                else
                    raylon = ray.de_lon_ray_fin.*180/pi;
                end
                %lat = unique(LAT(1,:)*180/pi());
                %lon = unique(LON(:,1)*180/pi());
                %id_in  = find(ray.de_lat_ray_fin.*180/pi <= lat_v(end) & ray.de_lat_ray_fin.*180/pi >= lat_v(1) & raylon <= lon_v(end) & raylon >= lon_v(1) & ray.de_alt_ray_fin.*1000 <= alt_v(end) & ray.de_alt_ray_fin.*1000 >= alt_v(1));
                id_in  = find(ray.de_lat_ray_fin <= LAT(1,1) & ray.de_lat_ray_fin >= LAT(1,end) & raylon <= LON(end,1)*180/pi() & raylon >= LON(1,1)*180/pi() & ray.de_alt_ray_fin.*1000 <= alt_N(end) & ray.de_alt_ray_fin.*1000 >= alt_N(1));

                if ~isempty(id_in)
                    % Loop over all entries inside the voxel model
                    for i = 1:length(id_in)
                        % Get coordinates of ray point
                        %[vox(i,:),i_pos(i)] = find_num(ray.de_lat_ray_fin(id_in(i)).*180/pi, raylon(id_in(i)),ray.de_alt_ray_fin(id_in(i))*1000,lat_v,lon_v,alt_v,0); 
                        [vox(i,:),i_pos(i)] = find_num(ray.de_lat_ray_fin(id_in(i)).*180/pi, raylon(id_in(i)),ray.de_alt_ray_fin(id_in(i))*1000,flip(LAT(1,:)*180/pi()),LON(:,1)'*180/pi(),alt_N,0); 
                    end
                    try
                        bound = size(ray.refr,1);
                        if bound > size(id_in,1)
                            bound = size(id_in,1);
                            id_in = id_in(1:bound);
                            vox = vox(1:bound,:);
                            i_pos = i_pos(1:bound);

                        end
                        %[voxIN,voxEM] = vox_distance(lat_v, lon_v, alt_v, raylon,radii, id_in, ray,bound);
                        [voxIN,voxEM] = vox_distance(flip(LAT(1,:)*180/pi()),LON(:,1)'*180/pi(), alt_N, raylon,radii, id_in, ray,bound);

                    catch
                        disp('Failed to calculate distance')
                        ray.error = [ray.error 200];
                    end
                    if size(voxIN.distanceVox,1) > size(ray.refr,1)
                        id_in = id_in(1:size(ray.refr,1));
                        vox = vox(1:size(ray.refr,1),:);
                        i_pos = i_pos(1,1:size(ray.refr,1));
                    end    
                    try
                        vox(:,8) = voxIN.distanceVox(1:size(id_in,1));
                        %Additional variables
                        vox(:,9) = ray.refr(1:size(id_in,1));
                        vox(:,10) = ray.refh(1:size(id_in,1));
                        vox(:,11) = ray.refw(1:size(id_in,1));
                        vox(:,12) = ray.Temp(1:size(id_in,1));
                        vox(:,13) = ray.pres(1:size(id_in,1));
                        vox(:,14) = ray.hum(1:size(id_in,1));
                        % Save variables
                        voxIN.id_in = id_in;
                        ray.i_pos = i_pos;
                        ray.vox   = vox;
                        ray.voxEM = voxEM;
                        ray.voxIN = voxIN;
                    catch
                        warning('voxel_dist_3D_combined: Failed to save variables... repeating')
                        ray.error = [ray.error 500];
                    end
                    if isfield(ray.error)
                        try
                            vox(:,8) = voxIN.distanceVox(1:size(id_in,1));
                            % Save variables
                            voxIN.id_in = id_in;
                            ray.i_pos = i_pos;
                            ray.vox   = vox;
                            ray.voxEM = voxEM;
                            ray.voxIN = voxIN; 
                        catch
                            warning('voxel_dist_3D_combined: Failed to save variables 2nd time')
                        end
                    end
                else
                    % Return dummy value '0' for voxel [1 1]
                    ray.d_voxel = 0;
                    ray.n_voxel = [1 1];
                end
            %catch
                %ray.error = [ray.error 300];      
            end
        end
        
        
    else
        ray.dL = NaN;
    end
    if set.stop
        return
    end
else
    ray.dL = NaN;
end
