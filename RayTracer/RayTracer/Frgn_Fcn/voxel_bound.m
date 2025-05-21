function [lat_v,lon_v,alt_v] = voxel_bound(lat_N,lon_N,alt_N)
%
% This function computes the inner and outer voxel model boundaries from
% given grid points.
%
% Copyright (C) 2019 Gregor Moeller
% All rights reserved.
% Email: gregor.moeller@jpl.nasa.gov

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% Input
% lat_N.....   ellipsoidal latitude of refractivity [degrees]
% lon_N.....   ellipsoidal longitude of refractivity [degrees]
% alt_N.....   ellipsoidal height of refractivity [m]
%
% Output
% lat_v.....   ellipsoidal latitude of voxel boundaries [degrees]
% lon_v.....   ellipsoidal longitude of voxel boundaries [degrees]
% alt_v.....   ellipsoidal height of voxel boundaries [m]
%

% Existing latitudes and longitudes of voxel model
lat_list = unique(lat_N*180/pi());
lon_list = unique(lon_N*180/pi());
alt_list = unique(alt_N);

% Compute voxel model boundaries - inner model
for i = 2:length(lat_list)-1
    lat_v(i-1) = lat_list(i) - (lat_list(3)-lat_list(2))/2;
    lat_v(i)   = lat_list(i) + (lat_list(3)-lat_list(2))/2;    
end
for i = 2:length(lon_list)-1
    lon_v(i-1) = lon_list(i) - (lon_list(3)-lon_list(2))/2;
    lon_v(i)   = lon_list(i) + (lon_list(3)-lon_list(2))/2;
end

if length(lat_list) >= 2
    % Compute voxel model boundaries - outer model
lat_v = [2*lat_list(1)-lat_v(1),lat_v,2*lat_list(end)-lat_v(end)];
else
    lat_v = [];
end

if length(lon_list) >=2
    lon_v = [2*lon_list(1)-lon_v(1),lon_v,2*lon_list(end)-lon_v(end)];
else
    lon_v = [];
end

% Vertical voxel model boundaries
alt_v(1) = alt_list(1) - (alt_list(2)-alt_list(1))/2;
for i = 2:length(alt_list)+1
    alt_v(i) = alt_v(i-1) + 2*(alt_list(i-1)-alt_v(i-1));
end
