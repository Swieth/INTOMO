function hgeop = fgeom2geop(lat, hgeom)
%% ABOUT: 28.12.2016: vectorized for speed
% fgeom2geop converts geometric to geopotential height expressed in
% geopotential kilometers
% 
% Somigliana's Equation is defined for normal gravity on the ellipsoid  
% while geopotential height is strictly relative to the geoid. 
%
% 1) hgeom wrt ellipsoid -> hgeop wrt geoid
% Geoid undulation need to be substracted from hgeop to compute 
% geopotential height as there is no geopotential wrt ellipsoid.
%
% 2) hgeom wrt geoid -> hgeom wrt geoid
% Somagliana's equation is strictly based on computations wrt ellipsoid, 
% but the error in assuming these apply for conversions wrt geoid is small.
%
% Reference ellipsoid: WGS84
%% Input data: 
%               lat:     geocentric latitude    [rad]
%               hgeom:   geometric height       [km]
% Output data: 
%               hgeop:   geopotential height    [km]
%% CONSTANTS
a = 6378.137;               % semi-major axis [km]
b = 6356.7523142;           % semi-minor axis [km]
e2 = (a^2-b^2)/a^2;         % power eccentricy (e^2)
g0 = 9.80665;               % WMO gravity [m2/s2]
ga = 9.7803253359;          % equatorial gravity [m2/s2]
gr = 0.003449787;           % gravity ratio

% Normal gravity on ellipsoid (Somigliana, 1929)
g = ga*(1 + 1.9*1e-3*sin(lat).^2)./sqrt(1 - e2*sin(lat).^2);

% Effective radius 
Reff = a./(1 + e2/2 + gr - e2*sin(lat).^2);

% Geopotential height
hgeop = (g./g0).*(Reff.*hgeom./(Reff + hgeom));
end