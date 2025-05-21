% Function for voxel vertical gradient computation and it's interpolation to
% the ray point
%% Input data
% dh........... height level difference in ray point's voxel [m]
% dh2.......... hegiht level difference in voxel below or above ray point [m]
% dref......... refraction value difference in ray point's voxel
% dref2........ refraction value difference in voxel below or above ray point 
% alt.......... altitude of the ray point
% n_ray........ number of interation
% alt_N........ altitude profile
% m2........... index of height layer below ray point
%% Output data
% dpa   - voxel gradient
% dp    - gradient in raypoint

function [dpa,dp] = dpacalc(dh,dh2,dref,dref2,n_ray,dpa,alt,alt_N,m2)
    hlimit = 1;
    dpa(n_ray-1) =  dref(hlimit:end)./dh(hlimit:end);
    dpa2 =  dref2(hlimit:end)./dh2(hlimit:end);
    ind = 0;
    dh_s = 10;
    if dh_s(end) < abs(dh)
        while dh_s(end) < abs(dh)
            ind = ind + 1;
            dh_s(ind) = dh_s(end) + 1;          
        end
        [~, c3] = find(dpa ~= dpa2);
        try
            c3 = c3(end);
        end
        dp_s = [dpa(c3):((dpa2-dpa(c3))/(ind)):dpa2];          
        if size(dp_s) ~= ind
            dp_s  = [];
            dp_s = [dpa(c3):((dpa2-dpa(c3))/(ind-1)):dpa2]; 
        end
        ind2 = abs((alt*1000-alt_N(m2)));
        if ind2 < min(dh_s)
            ind2 = min(dh_s);
        end
        try
            dp = interp1(dh_s,flip(dp_s) ,ind2);
        catch
            dp = dpa(n_ray-1);
        end
        if isnan(dp) 
            %dp = max(dp_s); Before used for ground based simulations
             dp = dpa(n_ray-1);
        end
     else
        dp = dref(hlimit:end)./dh(hlimit:end);
     end
end