%Function to tranform gradients between cartesian and spherical
%coordinates.
%% Input data
% p........... cartesian coordinates of ray point
% grad_g...... refractivity gradient in cartesian coordinates

%% Output data
% grad_r...... refractivity gradient in spherical coordinates 

function grad_r = recgrad(p,grad_p,radii);
    %Earth radius
    R = radii(1);
    %Point coordiantes
    lat = p(1);
    lon = p(2);
    h = p(3);
    %Gradient values in cartesian coordinates
    i_n  = grad_p(1);
    j_n  = grad_p(2);
    k_n = grad_p(3); 
    %Tranformation
    r = sin(pi()-lat)*cos(lon)*(k_n) + sin(pi()-lat)*sin(lon)*i_n + cos(pi()-lat)*j_n;
    o = cos(pi()-lat)*cos(lon)*(k_n) + cos(pi()-lat)*sin(lon)*i_n - sin(pi()-lat)*j_n;
    l = -sin(lon)*(k_n) + cos(lon)*i_n;
    grad_r = [r,o,l];
end
