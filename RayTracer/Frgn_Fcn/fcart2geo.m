function [lat, lon, h] = fcart2geo(xyz, solution)
%% ABOUT: 13.01.2017 VECTORIZED FOR SPEED
% cart2geo converts cartesian geocentric coordinates: x, y, z to geodetic 
% coordinates expressed by latitude, longitude and height. 
% FEW OPTIONS AVAILABLE
%% Input data: 
%               xyz:  	  cartesian coordinates (i x 3)[km]
%               solution: e (ellipsoid) or s (sphere)  [string]
% Output data: 
%               lat:      geodetic latitude            [rad]
%               lon:      geodetic longitude           [rad]
%               h:        height over ellipsoid/sphere [km]
%% CONSTANTS
a = 6378.137;               % semi-major axis [km]

% Ellipsoidal or spherical solution
if solution == 'e'
    b = 6356.7523142;           % semi-minor axis [km]
elseif solution == 's' 
    b = a;
end

f = 1/298.257223563;        % inverse flattening
e2 = (a^2-b^2)/a^2;         % power eccentricy (e^2)

% Control input data: should be in 3 columns
assert(size(xyz,2) == 3, 'Provide input data columnwise (i x 3)')

% Input coordinates
x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

% Spherical radius
r = sqrt(x.^2 + y.^2);

% Spherical latitude
lats = atan2(z, r);

% Longitude
lon = atan2(y,x);

%% SOLUTION #1
% Ref: S.P. Keeler and Y. Nievergelt, "Computing geodetic coordinates", 
% SIAM Rev. Vol. 40, No. 2, pp. 300-309, June 1998

s = sqrt(r.^2 + z.^2).*(1 - a*((1-e2)./sqrt((1 - e2)*r.^2 + z.^2)));
t0 = 1 + s.*sqrt(1 - (sqrt(e2)*z).^2./(r.^2 + z.^2))/a;
dzeta1 = z.*t0;
xi1 = r.*(t0 - e2);
rho1 = sqrt(xi1.^2 + dzeta1.^2);
c1 = xi1./rho1;
s1 = dzeta1./rho1;
b1 = a./sqrt(1 - (sqrt(e2)*s1).^2);
u1 = b1.*c1;
w1 = b1.*s1.*(1 - e2);

% Output variables #1
lat2 = atan(s1./c1);
h2 = sqrt((r - u1).^2 + (z - w1).^2);

%% SOLUTION #2
% Spherical latitude
lat0 = atan2(z,((1 - e2).*r)); %atan2(z, r)

% Iterate for geodetic latitude
for i = 2:10
    Rn0(:,i) = a./sqrt(1 - e2*sin(lat0(:,i-1)).^2);        %prime vertical
    h0(:,i) = r./cos(lat0(:,i-1)) - Rn0(:,i);
    lat0(:,i) = atan2(z,((1 - e2*Rn0(:,i)./(Rn0(:,i) + h0(:,i))).*r));
    
    % Break loop if precision is sufficient
    if abs(lat0(:,i) - lat0(:,i-1)) < eps;
        break
    end
end

% Output variables #2
lat = lat0(:,end);
Rn = Rn0(:,end);
h = h0(:,end);

% to diverge at poles
[m, ~] = find(lat == pi/2);

% Correct altitudes at the poles
if ~isempty(m);
    fprintf('Diverges at north or south pole \n')
    h(m,:) = z(m,:) + e2*Rn(m,:).*sin(lat(m,:));
end

end