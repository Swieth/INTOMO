function [lat, lon, h] = aer2geodetic2( ...
    az, elev, slantRange, lat0, lon0, h0, spheroid, angleUnit)
%AER2GEODETIC Local spherical AER to geodetic
%
%   [LAT, LON, H] = AER2GEODETIC(AZ, ELEV, SLANTRANGE, LAT0, LON0, H0, ...
%       SPHEROID) transforms point locations in 3-D from local spherical
%   coordinates (azimuth angle, elevation angle, slant range) to geodetic
%   coordinates (LAT, LON, H), given a local coordinate system defined by
%   the geodetic coordinates of its origin (LAT0, LON0, H0).  The geodetic
%   coordinates refer to the reference body specified by the spheroid
%   object, SPHEROID. The slant range and ellipsoidal height H0 must be
%   expressed in the same length unit as the spheroid.  Ellipsoidal height
%   H will be expressed in this unit, also.  The input azimuth and
%   elevation angles, and input and output latitude and longitude angles,
%   are in degrees by default.
%
%   [...] = AER2GEODETIC(..., angleUnit) uses angleUnit, which matches
%   either 'degrees' or 'radians', to specify the units of the azimuth,
%   elevation, latitude, and longitude angles.
%
%   Class support for inputs AZ, ELEV, SLANTRANGE, LAT0, LON0, H0:
%      float: double, single
%
%   See also AER2ECEF, ENU2GEODETIC, GEODETIC2AER, NED2GEODETIC

% Copyright 2012-2020 The MathWorks, Inc.

%#codegen
radii = [6378.137, 6356.752314245];
inDegrees = (nargin < 8 || map.geodesy.isDegree(angleUnit));

[x, y, z] = aer2ecefFormula2( ...
    az, elev, slantRange, lat0, lon0, h0, spheroid, inDegrees);

%if inDegrees
    %[lat, lon, h] = ecef2geodetic2(spheroid, x, y, z);
%else
    %[lat, lon, h] = ecef2geodetic2(spheroid, x, y, z, 'radian');
%end
[lon_ray,lat_ray,alt_ray] = cspice_recgeo([x/1000;y/1000;z/1000],radii(1),(radii(1)-radii(2))/radii(1));
lat = lat_ray*180/pi();
lon = lon_ray*180/pi();
h = alt_ray*1000;

