function [pos,i_pos] = find_num(lat_p,lon_p,alt_p,lat_v,lon_v,alt_v,set)

% The function find_num checks in which voxel (lat_v,lon_v)
% a point (lat_p,lon_p) is located and delivers coordinates as well as
% id of voxel in 3D voxel model. 

% Input
% lat_p    ... latitude point [degrees]
% lon_p    ... longitude point [degrees]
% alt_p    ... altitude point [m]
% lat_v    ... latidudes voxel model [degrees]
% lon_v    ... longitudes voxel models [degrees]
% alt_p    ... altitude voxel models  [m]

% Output
% pos    ... coordinates voxel boundaries [lat1 lat2 lon1 lon2 alt1 alt2]
% i_pos    ... id of ray point 

% Computde differences
dlat = lat_p-lat_v;
dlon = lon_p-lon_v;
dalt = alt_p-alt_v;

% Get index of phi / lam
n_lat1 = find(dlat>=0,1,'last');
n_lon1 = find(dlon>=0,1,'last');
n_alt1 = find(dalt>=0,1,'last');
n_lat2 = find(dlat<0,1,'first');
n_lon2 = find(dlon<0,1,'first');
n_alt2 = find(dalt<0,1,'first');

% Correct index if outside
if set == 0 || set == 1 || set == 2 || set == 5
    if isempty(n_lat1)
        n_lat1 = 1;
        n_lat2 = 2;
    end
    if isempty(n_lon1)
        n_lon1 = 1;
        n_lon2 = 2;
    end
    if set == 0 || set == 2
        if n_lon1 == length(lon_v)
            n_lon1 = n_lon1 - 1;
            n_lon2 = n_lon1 - 2;
        end
        if n_lat1 == length(lat_v)
            n_lat1 = n_lat1 - 1;
            n_lat2 = n_lat1 - 2;
        end
    else 
        if n_lon1 == length(lon_v)
            n_lon1 = n_lon1;
            n_lon2 = n_lon1 - 1;
        end
        if n_lat1 == length(lat_v)
            n_lat1 = n_lat1;
            n_lat2 = n_lat1 - 1;
        end
    end
elseif set ==  3 || set == 4
    if isempty(n_lat1)
        n_lat1 = NaN;
    end 
    if isempty(n_lat2)
        n_lat2 = NaN;
    end 
    if isempty(n_lon1)
        n_lon1 = NaN;
    end 
    if isempty(n_lon2)
        n_lon2 = NaN;
    end 
    if isempty(n_lat1)
        n_lat1 = NaN;
    end 
    if isempty(n_lat2)
        n_lat2 = NaN;
    end 
end
% Compute voxel coordinates
if set == 0 || set == 4 || set == 5
    pos = [lat_v(n_lat1),lat_v(n_lat2),lon_v(n_lon1),lon_v(n_lon2),alt_v(n_alt1),alt_v(n_alt2)];
elseif set == 1|| set == 2 || set == 3
    pos = [n_lat1,n_lat2,n_lon1,n_lon2,n_alt1,n_alt2];    
end
    
if set == 2
    i_pos = (n_alt1-1)*(length(lon_v)-1)*(length(lat_v)-1) + n_lat1 + (length(lat_v)-1)*(n_lon1-1);
else
    i_pos = (n_alt1-1)*(length(lon_v)-1)*(length(lat_v)-1) + n_lon1 + (length(lon_v)-1)*(n_lat1-1);
end








