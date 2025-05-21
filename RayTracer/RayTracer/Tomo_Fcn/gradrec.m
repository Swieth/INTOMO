%Function to tranform gradients between spherical and cartesian
%coordinates.
%% Input data
% p........... spherical coordinates of ray point
% grad_g...... refractivity gradient in spherical coordinates

%% Output data
% grad_r...... refractivity gradient in cartesian coordinates 

function grad_r = gradrec(p,grad_g)
    %Point coordiantes
    lat = p(2);
    lon = p(1);
    h = p(3);
    %Gradient values in spherical coordinates
    lat_n  = grad_g(1);
    lon_n  = grad_g(2);
    h_n = grad_g(3); 
    %a = pi()/2 - lat;
    a = lat;
    b = lon;
    i = sin(a)*cos(b)*(h_n) + cos(a)*cos(b)*lat_n - sin(b)*lon_n;
    j = sin(a)*sin(b)*(h_n) + cos(a)*sin(b)*lat_n + cos(b)*lon_n;
    k = cos(a)*(h_n) - sin(a)*lat_n;
    grad_r = [i,j,k];
end

