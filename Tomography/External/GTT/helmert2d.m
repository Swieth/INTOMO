function [tp,ac,tr]=helmert2d(datum1,datum2,WithOutScale)

% HERLMERT2D    overdetermined cartesian 2D similarity transformation ("Helmert-Transformation")
%
% [param, accur, resid] = helmert2D(datum1,datum2,DontUseScale)
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
%    DontUseScale  if this is not 0, do not calculate scale factor but set it to the inputted value
%                  Default: 0 (Use scale)
%
% Outputs:  param  4 x 1 Parameter set of the 2D similarity transformation
%                      2 translations (x y) in [Unit of datums]
%                      1 rotation (ez) in [rad]
%                      1 scale factor
%
%           accur  4 x 1 accuracy of the parameters (or 3 x 1 if scale factor is set to be 1)
%
%           resid  n x 2 - matrix with the residuals datum2 - f(datum1,param)
%
% Used to determine transformation parameters e.g. for cadastral purposes (transforming local system
% to 2D mapping system) when at least 2 identical points in both datum systems are known.
% Parameters can be used with d2trafo.m

% 09/12/11 Peter Wasmeier - Technische Universität München
% p.wasmeier@bv.tum.de

%% Argument checking and defaults

if nargin<3 || isempty(WithOutScale)
    WithOutScale=0;
end

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
elseif any([s1(1) s2(1)]<2)
    error('At least 2 points in each datum are necessary for calculating')
end

%% Adjustment

naeh=[0 0 0 1];
if WithOutScale
    naeh(4)=WithOutScale;
end
WertA=[1e-8 1e-8];
zaehl=0;

x0=naeh(1);
y0=naeh(2);
ez=naeh(3);
m=naeh(4);

tp=[x0 y0 ez m];

Qbb=eye(2*s1(1));

while(1)
    A=zeros(2*s1(1),4);
    for i=1:s1(1)
        A(i*2-1,1)=1;
        A(i*2,2)=1;
        A(i*2-1,3)=m*(-sin(ez)*datum1(i,1)+cos(ez)*datum1(i,2));
        A(i*2-1,4)=cos(ez)*datum1(i,1)+sin(ez)*datum1(i,2);
        A(i*2,3)=m*(-cos(ez)*datum1(i,1)-sin(ez)*datum1(i,2));
        A(i*2,4)=-sin(ez)*datum1(i,1)+cos(ez)*datum1(i,2);
 
        w(i*2-1,1)=datum2(i,1)-x0-m*(cos(ez)*datum1(i,1)+sin(ez)*datum1(i,2));
        w(i*2,1)=datum2(i,2)-y0-m*(-sin(ez)*datum1(i,1)+cos(ez)*datum1(i,2));
    end

    if WithOutScale
        A=A(:,1:3);
    end
    
    warning off;
    r=size(A,1)-size(A,2);
    Pbb=inv(Qbb);
    deltax=inv(A'*Pbb*A)*A'*Pbb*w;
    v=A*deltax-w;
    sig0p=sqrt((v'*Pbb*v)/r);
    Qxxda=inv(A'*Pbb*A);
    Kxxda=sig0p^2*Qxxda;
    ac=sqrt(diag(Kxxda));
    warning on;

    testv=sqrt((deltax(1)^2+deltax(2)^2)/2);
    testd=deltax(3);
    zaehl=zaehl+1;
    x0=x0+deltax(1);
    y0=y0+deltax(2);
    ez=ez+deltax(3);
    if ~WithOutScale && (m+deltax(4))>1e-15     % This condition is to prevent numerical problems with m-->0
        m=m+deltax(4);
    end
    tp=[x0 y0 ez m]';
    if abs(testv) < WertA(1) && abs(testd) < WertA(2)
        break;
    elseif zaehl>1000
        warning('Helmert2D:Too_many_iterations','Calculation not converging after 1000 iterations. I am aborting. Results may be inaccurate.')
        break;
    end
end

%% Transformation residuals
idz=zeros(s1);
for i=1:s1(1)
    idz(i,2)=tp(2)+tp(4)*(-sin(tp(3))*datum1(i,1)+cos(tp(3))*datum1(i,2));
    idz(i,1)=tp(1)+tp(4)*(cos(tp(3))*datum1(i,1)+sin(tp(3))*datum1(i,2));
end
tr=datum2-idz;