function [x, y, z] = aer2ecefFormula2( ...
    az, elev, slantRange, lat0, lon0, h0, spheroid, inDegrees)
% Transform position from local spherical (AER) to geocentric (ECEF)
%
%   The outputs are the same as in AER2ECEF. The inputs are also the same,
%   except that the angleUnit string is replaced by logical scalar,
%   inDegrees, which is true when AZ, ELEV, LAT0 and LON0 are in degrees
%   and false when they are in radians.
%
%   See also AER2ECEF, AER2GEODETIC

% Copyright 2012 The MathWorks, Inc.
radii = [6378.137, 6356.752314245];
if inDegrees
    % Sine and cosine function handles, for input in degrees.
    sinfun = @sind;
    cosfun = @cosd;
    [X_T2bis] = cspice_georec( lon0*pi()/180, lat0*pi()/180, h0/1000, radii(1),(radii(1)-radii(2))/radii(1));
    % Origin of the local system in geocentric coordinates.
    %[x0, y0, z0] = spheroid.geodetic2ecef(lat0, lon0, h0);
    x0 = X_T2bis(1)*1000;
    y0 = X_T2bis(2)*1000;
    z0 = X_T2bis(3)*1000;
else
    % Sine and cosine function handles, for input in radians.
    sinfun = @sin;
    cosfun = @cos;
    
    % Origin of the local system in geocentric coordinates.
    [x0, y0, z0] = spheroid.geodetic2ecef(lat0, lon0, h0, 'radian');
    
end

% Transform local spherical AER to Cartesian ENU.
[xEast, yNorth, zUp] = aer2enuFormula2( ...
        az, elev, slantRange, sinfun, cosfun);

% Offset vector from local system origin, rotated from ENU to ECEF.
[dx, dy, dz] = enu2ecefvFormula2( ...
    xEast, yNorth, zUp, lat0, lon0, sinfun, cosfun);

% Origin + offset from origin equals position in ECEF.
x = x0 + dx;
y = y0 + dy;
z = z0 + dz;
