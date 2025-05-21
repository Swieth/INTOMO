function [tp,ac,tr]=helmertprojective3d(datum1,datum2)

% HERLMERTPROJECTIVE3D    overdetermined cartesian 3D projective transformation to 2D
%
% [param, accur, resid] = helmertprojective3D(datum1,datum2)
%
% Inputs:  datum1  n x 3 - matrix with coordinates in the 3D origin datum (x y z)
%                  datum1 may also be a file name with ASCII data to be processed. No point IDs, only
%                  coordinates as if it was a matrix.
%
%          datum2  n x 2 - matrix with coordinates in the 2D destination datum (x y)
%                  datum2 may also be a file name with ASCII data to be processed. No point IDs, only
%                  coordinates as if it was a matrix.
%                  If either datum1 and/or datum2 are ASCII files, make sure that they are of same
%                  length and contain corresponding points. There is no auto-assignment of points!
%
% Outputs:  param  11 x 1 Parameter set of the 3D projective transformation
%
%           accur  11 x 1 accuracy of the parameters
%
%           resid  n x 2 - matrix with the residuals datum2 - f(datum1,param)
%
% Used to calculate projective transformation parameters when at least 6 identical points in both systems
% are known. An 3D projective transformation is mainly used for image mapping.
% Parameters can be used with d3projectivetrafo.m

% 11/12/21 Peter Wasmeier - Technische Universität München
% p.wasmeier@bv.tum.de
 
%% Argument checking and defaults

% Load input file if specified
if ischar(datum1)
    datum1=load(datum1);
end
if ischar(datum2)
    datum2=load(datum2);
end

if (size(datum1,1)==3)&&(size(datum1,2)~=3)
    datum1=datum1'; 
end
if (size(datum2,1)==2)&&(size(datum2,2)~=2)
    datum2=datum2'; 
end

s1=size(datum1);
s2=size(datum2);
if any(s1~=s2)
    error('The datum sets are not of equal size')
elseif s1(2)~=3)
    error('Datum1 needs to be 3D.')
elseif s2(2)~=2)
    error('Datum2 needs to be 2D.')
elseif any([s1(1) s2(1)]<6)
    error('At least 6 points in each datum are necessary for calculating')
end

%% Projective transformation:
G=zeros(2*s1(1),11);
t=zeros(2*s2(1),1);
for i=1:s1(1)
   G(2*i-1:2*i,:)=[datum1(i,:) 1 0 0 0 0 -datum1(i,:).*datum2(i,1);
                   0 0 0 0 datum1(i,:) 1 -datum1(i,:).*datum2(i,2)];
   t(2*i-1:2*i)=datum2(i,:)';
end
tp=inv(G'*G)*G'*t;
if (size(G,1)>11)
    v=G*tp-t;
    sig0p=sqrt((v'*v)/(size(G,1)-8));
    ac=sqrt(diag(sig0p^2*inv(G'*G)));
else
    ac=zeros(8,1);
end

%% Transformation residuals
idz=zeros(s1);
for i=1:s1(1)
    idz(i,:)=[tp(1)*datum1(i,1)+tp(2)*datum1(i,2)+tp(3)*datum1(i,3)+tp(4) ...
        tp(5)*datum1(i,1)+tp(6)*datum1(i,2)+tp(7)*datum1(i,3)+tp(8)]...
        /(tp(9)*datum1(i,1)+tp(10)*datum1(i,2)+tp(11)*datum1(i,3)+1);
end
tr=datum2-idz;
