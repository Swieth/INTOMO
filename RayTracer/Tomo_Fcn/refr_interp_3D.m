function  [refr_ray,inp_data] = refr_interp_3D(ref_a,ref_w,lat_N,lon_N,alt_N,lat_ray,lon_ray,alt_ray,hlim, Temp,hum,pres,rwgs,set,inp_data)

% Function refr_interp_3D reads the hydrostatic and wet refractivity values 
% as provided on a lat/lon/height grid and computes refractivity for any point
% of interest (POI) inside or outside the grid

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

%% Input data:
% ref_h........... hydrostatic refractivity [ppm], matrix containing all levels
% ref_w........... wet refractivity [ppm], matrix containing all levels
% lat_N........... ellipsoidal latitude of refractivity [degrees]
% lon_N........... ellipsoidal longitude of refractivity [degrees]
% alt_N........... ellipsoidal height of refractivity [m]
% lat_ray......... latitude of POI [degrees]
% lon_ray......... latitude of POI [degrees]
% alt_ray......... altitude of POI [m]
% hlim............ max altitude to limit refraction calculations 
% temp............ temperature, matrix containing all levels
% hum............. water vapour, , matrix containing all levels 
% pres............ pressure altitude to limit refraction calculations 
% rwgs............ max altitude to limit refraction calculations 

%% Output data:
% refr_ray........ refractivity for POI [ppm]
% inp_data
%       ref_h..... hydrostatic refractivity [ppm], matrix containing all levels
%       ref_w..... wet refractivity [ppm], matrix containing all levels
%       lat_N..... ellipsoidal latitude of refractivity [degrees]
%       lon_N..... ellipsoidal longitude of refractivity [degrees]
%       alt_N..... ellipsoidal height of refractivity [m]
%       lat_ray... latitude of POI [degrees]
%       lon_ray... latitude of POI [degrees]
%       alt_ray... altitude of POI [m]
%       hlim...... max altitude to limit refraction calculations 
%       temp...... temperature, matrix containing all levels
%       hum....... water vapour, , matrix containing all levels 
%       pres...... max altitude to limit refraction calculations 
%       Nh........ hydrostatic refractivity at ray point
%       Nw........ wet refractivity at ray point
%       N......... total refractivity at ray point
%       T......... temperature at ray point [K]
%       p......... pressue refractivity at ray point [hPa]
%       h......... water vapour refractivity at ray point [kg kg-1]
%       ref0...... mean of total refraction profiles from voxels around ray point  [kg kg-1]

if alt_ray > hlim
    refr_ray = 0;
    return
end
if isempty(inp_data) == 1
    n_lat = length(unique(lat_N));
    n_lon = length(unique(lon_N));
    inp_data.lat = lat_N;
    inp_data.lon = lon_N;
    inp_data.alt = alt_N;
    inp_data.ref = ref_a;
    inp_data.refw = ref_w;
    inp_data.Temp = Temp;
    inp_data.pres= pres;
    inp_data.hum = hum;
    if inp_data.lon(2) >= inp_data.lon(1) 
        inp_data.lon_new = (inp_data.lon(:,1)-inp_data.lon(1,1))';
    else
        inp_data.lon_new = (inp_data.lon(:,1) +360 - inp_data.lon(1,1))';
    end
    inp_data.lon_new = repmat(inp_data.lon_new',[1 n_lat]);
    if min(inp_data.lat(1,:)) < 0
        inp_data.lat_new = abs(inp_data.lat(1,:)-inp_data.lat(1,1));
        inp_data.lat_new = repmat(inp_data.lat_new,[n_lon 1]);
    else
        inp_data.lat_new = inp_data.lat;
    end
    inp_data.hlayer = inp_data.alt;
    inp_data.rwgs= rwgs;
end
% Change ray point position to fit transformed GRID
if lon_ray < inp_data.lon(1,1)
    lon_ray_new = lon_ray +360 - inp_data.lon(1,1);
else
    lon_ray_new = lon_ray - inp_data.lon(1,1);
end 
if min(inp_data.lat(1,:)) < 0
    lat_ray_new = 180 - lat_ray; %180 - lat_ray ;
else
    lat_ray_new =  90 - lat_ray;%90-lat_ray;
end
if lat_ray_new >= min(min(inp_data.lat_new)) && lat_ray_new <= max(max(inp_data.lat_new)) && lon_ray_new >= min(min(inp_data.lon_new)) && lon_ray_new <= max(max(inp_data.lon_new))+5 && alt_ray > min(inp_data.alt) && alt_ray <= max(inp_data.alt)  
    if set.refron % For interpolation of only refractivity
        N = IDW_atom(inp_data.lat*pi()/180, inp_data.lon*pi()/180, inp_data.Temp, inp_data.hum, inp_data.pres,inp_data.ref , inp_data.refw ,inp_data.rwgs,[lat_ray*pi()/180 lon_ray*pi()/180 alt_ray],inp_data.hlayer,set);
        refr_ray = N;
        for i = 1:4
            prof0(:,i) = inp_data.ref(N.coord(i,1),N.coord(i,2),:);
        end
        inp_data.ref0 = mean(prof0,2);
    else % For refractivity from meteo parameters
        N = IDW_atom(inp_data.lat*pi()/180, inp_data.lon*pi()/180, inp_data.Temp, inp_data.hum, inp_data.pres,inp_data.ref ,inp_data.refw , inp_data.rwgs,[lat_ray*pi()/180 lon_ray*pi()/180 alt_ray],inp_data.hlayer,set);
        %% lon_ray_new!!!
        refr_ray.Nh = N.Nh;
        refr_ray.Nw = N.Nw;
        refr_ray.N = N.Nt;
        refr_ray.T = N.T;
        refr_ray.p = N.p;
        refr_ray.h = N.h;
        for i = 1:4
            prof0(:,i) = inp_data.ref(N.coord(i,1),N.coord(i,2),:);
        end
        inp_data.ref0 = mean(prof0,2);
    end
else
    % Find closest grid nodes for each voxel center point
    [~,id_lat] = min(abs(inp_data.lat_new(1,:,1)-lat_ray_new));
    [~,id_lon] = min(abs(inp_data.lon_new(:,1,1)-lon_ray_new));
    inp_data.ref0 = squeeze(mean(mean(ref_a)));
    % Below the given grid
    if alt_ray <= min(min(min(inp_data.alt))) %do 225m
        % Select the lowest available value
        if set.refron
            refr_ray = inp_data.ref(id_lon,id_lat,1); 
        else
            N = IDW_atom(inp_data.lat*pi()/180, inp_data.lon*pi()/180, inp_data.Temp, inp_data.hum, inp_data.pres,inp_data.ref ,inp_data.refw , inp_data.rwgs,[lat_ray*pi()/180 lon_ray_new*pi()/180 min(inp_data.alt)],inp_data.hlayer,set);
            refr_ray.Nh = N.Nh;
            refr_ray.Nw = N.Nw;
            refr_ray.N = N.Nt;
            refr_ray.T = N.T;
            refr_ray.p = N.p;
            refr_ray.h = N.h;
        end
    % Within the height range but outside the grid
    else 
        if set.refron
            refr_ray = interp1(inp_data.alt, squeeze(inp_data.ref(id_lon,id_lat,:)),alt_ray,'linear');
        else
            refr_ray.Nw = 0; %interp1(inp_data.alt, squeeze(inp_data.refw(id_lon,id_lat,1:end-1)),alt_ray,'linear');
            if size(inp_data.alt,1) == size(inp_data.ref,3)
                refr_ray.Nh = interp1(inp_data.alt, squeeze(inp_data.ref(id_lon,id_lat,1:end)),alt_ray,'linear');
            else
                refr_ray.Nh = interp1(inp_data.alt, squeeze(inp_data.ref(id_lon,id_lat,1:end-1)),alt_ray,'linear');
            end
            refr_ray.N = refr_ray.Nh+refr_ray.Nw;
            refr_ray.T = NaN;
            refr_ray.p = NaN;
            refr_ray.h = NaN;
            inp_data.ref0  = inp_data.ref(id_lon,id_lat,1:end);
        end
    end 
end
