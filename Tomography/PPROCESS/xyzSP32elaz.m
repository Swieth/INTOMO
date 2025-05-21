function [el,az] = xyzSP32elaz(XYZ,xyz)
%          system with origin at X.
%          Both parameters are 3 by 1 vectors.
%          Output: D    vector length in units like the input
%                  Az   azimuth from north positive clockwise, degrees
%                  El   elevation angle, degrees


[lla] = ECEF2LLA(XYZ,1);
lat = lla(1);
lon = lla(2);
h = lla(3);
cl = cos(lon); sl = sin(lon);
cb = cos(lat); sb = sin(lat);
F = [-sl -sb*cl cb*cl;
      cl -sb*sl cb*sl;
       0    cb   sb];
% local_vector = F'*xyz'; %wrong
local_vector = F'*(xyz-XYZ)'; % corrected 30.11.2021
E = local_vector(1);
N = local_vector(2);
U = local_vector(3);
hor_dis = sqrt(E^2+N^2);
if hor_dis < 1.e-20
   az = 0;
   el = pi();
else
   az = atan2(E,N);
   el = atan2(U,hor_dis);
end
if az < 0
   az = az+2*pi();
end

% % 30.11.2021
% % The same results when the following formulas are used
% % https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
% ed = [-sl;cl;0];
% nd = [-cl*sb;-sl*sb;cb];
% ud = [cl*cb;sl*cb;sb];
% ro = (xyz'-XYZ')./norm(xyz'-XYZ');
% el = asin(dot(ud,ro));
% eldeg = el*180/pi();
% azi = atan2(dot(ro,ed),dot(ro,nd));
% azideg = azi*180/pi();










