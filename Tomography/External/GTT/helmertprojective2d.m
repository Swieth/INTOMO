function [tp,ac,tr]=helmertprojective2d(datum1,datum2)

% HERLMERTPROJECTIVE2D    overdetermined cartesian 2D projective transformation
%
% [param, accur, resid] = helmertprojective2D(datum1,datum2)
%
% Inputs:  datum1  n x 2 - matrix with coordinates in the origin datum (x y)
%                  datum1 may also be a file name with ASCII data to be processed. No point IDs, only
%                  coordinates as if it was a matrix.
%
%          datum2  n x 2 - matrix with coordinates in the destination datum (x y)
%                  datum2 may also be a file name with ASCII data to be processed. No point IDs, only
%                  coordinates as if it was a matrix.
%                  If either datum1 and/or datum2 are ASCII files, make sure that they are of same
%                  length and contain corresponding points. There is no auto-assignment of points!
%
% Outputs:  param  8 x 1 Parameter set of the 2D projective transformation
%
%           accur  8 x 1 accuracy of the parameters
%
%           resid  n x 2 - matrix with the residuals datum2 - f(datum1,param)
%
% Used to calculate projective transformation parameters when at least 4 identical points in both systems
% are known. An projective transformation is mainly used for image rectification.
% Parameters can be used with d2projectivetrafo.m

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

if (size(datum1,1)==2)&&(size(datum1,2)~=2)
    datum1=datum1'; 
end
if (size(datum2,1)==2)&&(size(datum2,2)~=2)
    datum2=datum2'; 
end

s1=size(datum1);
s2=size(datum2);
if any(s1~=s2)
    error('The datum sets are not of equal size')
elseif any([s1(2) s2(2)]~=[2 2])
    error('At least one of the datum sets is not 2D')
elseif any([s1(1) s2(1)]<4)
    error('At least 4 points in each datum are necessary for calculating')
end

%% Projective transformation:
G=zeros(2*s1(1),8);
t=zeros(2*s2(1),1);
for i=1:s1(1)
   G(2*i-1:2*i,:)=[datum1(i,:) 1 0 0 0 -datum1(i,:).*datum2(i,1);
                   0 0 0 datum1(i,:) 1 -datum1(i,:).*datum2(i,2)];
   t(2*i-1:2*i)=datum2(i,:)';
end
tp=inv(G'*G)*G'*t;
if (size(G,1)>8)
    v=G*tp-t;
    sig0p=sqrt((v'*v)/(size(G,1)-8));
    ac=sqrt(diag(sig0p^2*inv(G'*G)));
else
    ac=zeros(8,1);
end

%% Transformation residuals
idz=zeros(s1);
for i=1:s1(1)
    idz(i,:)=[tp(1)*datum1(i,1)+tp(2)*datum1(i,2)+tp(3) tp(4)*datum1(i,1)+tp(5)*datum1(i,2)+tp(6)]/(tp(7)*datum1(i,1)+tp(8)*datum1(i,2)+1);
end
tr=datum2-idz;

