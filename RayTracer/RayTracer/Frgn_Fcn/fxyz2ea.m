function [el, az, range, ENU] = fxyz2ea(xyzr, XYZs, solution)
%% ABOUT: 19.08.2016
% xyz2ea calculates elevation angle and azimuth between satellite and 
% receiver with know cartesian coordinates (ECEF)
%% Input data: 
%               xyzr:  	  origin coordinates (i x 3)   [km]
%               XYZs:  	  satellite coordinates (i x 3)[km]
%               solution: e (ellipsoid) or s (sphere)  [string]
% Output data: 
%               az:       azimuth wrt N (clockwise)	   [rad]
%               el:       elevation angle              [rad]
%               range:    slant distance               [km]

%% CONSTANTS
a = 6378.137;               % semi-major axis [km]

% Ellipsoidal or spherical solution
if solution == 'e'
    b = 6356.7523142;       % semi-minor axis [km]
elseif solution == 's' 
    b = a;
end

f = 1/298.257223563;        % inverse flattening
e2 = (a^2-b^2)/a^2;         % power eccentricy (e^2)

% Convinient variables
xr = xyzr(1,1);
yr = xyzr(1,2);
zr = xyzr(1,3);

Xs = XYZs(1,1);
Ys = XYZs(1,2);
Zs = XYZs(1,3);


% Convert origin cartesian coordinates to latitude / longitude
[lat, lon, h] = fcart2geo(xyzr, solution);

% Rotation matrix
R = [-sin(lat).*cos(lon) -sin(lat).*sin(lon) cos(lat)
        -sin(lon)             cos(lon)          0
    cos(lat).*cos(lon)   cos(lat).*sin(lon)   sin(lat)];

% Convert to topocentric coordinates
N = R(1,1)*(Xs-xr) + R(1,2)*(Ys-yr) + R(1,3)*(Zs-zr);
E = R(2,1)*(Xs-xr) + R(2,2)*(Ys-yr) + R(2,3)*(Zs-zr);
U = R(3,1)*(Xs-xr) + R(3,2)*(Ys-yr) + R(3,3)*(Zs-zr);

% Topocentric coordinates
ENU = [E N U];

% Distance between satellite and origin
range = sqrt(N.^2 + E.^2 + U.^2);

% Azimuth [rad]
az = atan2(E, N);

% Elevation angle [rad]
el = atan2(U, sqrt(E.^2 + N.^2));

% Correct negative azimuths
if az < 0
   az = az + 2*pi;
end

end