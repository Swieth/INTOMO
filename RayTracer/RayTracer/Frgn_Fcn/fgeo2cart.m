function [x, y, z] = fgeo2cart(llh, solution)
%% ABOUT: 13.01.2017 VECTORIZED FOR SPEED
% geo2cart converts geocentic coordinates expressed by latitude, longitude
% and altitude to cartesian coordinates (ECEF): x, y, z together with radius.
%
% Input data: 
%               llh:       geocentic latitude          [rad]
%               llh:       geocentic longitude         [rad]
%               llh:       geometric altitude          [km]
%               solution:  e (ellipsoid) or s (sphere) [string]
% Output data: 
%               x:         cartesian x coordinate      [km]
%               y:         cartesian y coordinate      [km]
%               z:         cartesian z coordinate      [km]
%% WGS84 CONSTANTS
a = 6378.137;
b = 6356.7523142;           % semi-minor axis [km]
f = 1/298.257223563;        % inverse flattening
e2 = (a^2-b^2)/a^2;         % power eccentricy (e^2)

% Control input data: should be in 3 columns
assert(size(llh,2) == 3, 'Provide input data columnwise (i x 3)')

% Input coordinates
lat = llh(:,1);
lon = llh(:,2);
alt = llh(:,3);

if solution == 'e'
	%% Ellipsoidal
	Rn = a./sqrt(1 - e2*sin(lat).^2);        %prime vertical
	x = (Rn + alt).*cos(lat).*cos(lon);
	y = (Rn + alt).*cos(lat).*sin(lon);
	z = ((1 - e2).*Rn + alt).*sin(lat);
elseif solution == 's'
	%% Spherical
	x = (a + alt).*cos(lat).*cos(lon);
	y = (a + alt).*cos(lat).*sin(lon);
	z = (a + alt).*sin(lat);
	r = sqrt(x.^2 + y.^2);
    p = sqrt(x.^2 + y.^2 + z.^2);
end
end