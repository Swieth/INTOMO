function [ azimuthabgrad ] = azimuth( Pa,Pb )
%UNTITLED Function to calculate azimuth AB in grads
%   PA and PB are points of coordinates (x,y)

xa = Pa(1,1);
ya = Pa(1,2);
xb = Pb(1,1);
yb = Pb(1,2);

dx = xb - xa;
dy = yb - ya;

phi=atan(dy/dx);

azimuthab=NaN;

if dy~=0 && dx>0
    if dy>0
        azimuthab = phi;
    elseif dy<0
        azimuthab = phi + 2*pi();
    end
    
elseif dy~=0 && dx<0
    azimuthab = phi + pi();
    
elseif dy~=0 && dx==0
    if ya<yb
        azimuthab = pi()/2;
    elseif ya>yb
        azimuthab = pi()*(3/2);
    end
elseif dy==0
    if xa<xb
        azimuthab = 0;
    elseif xa>xb
        azimuthab = pi();
    end
end

if azimuthab<0
    azimuthab = azimuthab + 2*pi();
end

azimuthabgrad = azimuthab*200/pi();


end

