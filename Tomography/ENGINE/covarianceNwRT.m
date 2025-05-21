function Q = covarianceNwRT(Nw_apr,A,num_lat,num_lon,num_lev,BLh_pudel)
%function designed to get first estimate of dynamic disturbance in the
%wet refractivity based on provided Nw data (best working NWP forecasts)
%%%INPUT
%       Nw_apr.........apriori values of refractivities
%       A..............observational matrix with signal deriveratives
%       num_lat........number of voxels in latitudinal direction
%       num_lon........number of voxels in longitudinal direction
%       num_lev........number of voxels in vertical direction
%       BLh_pudel......coordinates of voxel edges/centers
%%%OUTPUT
%       Q..............  the first guess for covariance matrix

%%%%%%%%%
h = BLh_pudel(3,1,:);
h = h(:);
clear wie kol wart
for i = 1: size(Nw_apr,1)
        NWPone =  reshape(Nw_apr(i,:),num_lat*num_lon,num_lev-1);
        NWPext(i,:) =  mean(NWPone);
end

dtt = ones(size(A,2),1);
add_dtt = zeros(size(A,2),1);

modQ = cov(NWPext);
modQ = ((sqrt(diag(modQ))));

if var(modQ) > 0.1  %in case we have real data from NWP
    modQ = repmat(modQ,1,num_lat*num_lon);
    modQ = reshape(modQ',num_lat*num_lon*(num_lev-1),1);
else %in case there are only models available so very limited change over time thus the empirical variance function introduced based on the Trzcina et al.,2016
    % in the folder /TOMO-LAB/PPROCES/Q_poly_fit shows fiting pice-wise
    % function bottom - from 0 to 3875 m linear from 3876 exponential
    % function is hard-coded into this software
    clear modQ
    p1 = [-0.0552 0.6214];
    p2 = [-0.7693 7.7347];
    h = h/1000;
    [lin_wie,~,~] = find(h<=3.875);
    modQ(1,lin_wie) = p1(1,1)*h(lin_wie,1) + p1(1,2);
    clear lin_wie kol_wie wart_wie
    [lin_wie,~,~] = find(h>3.875);
    modQ(1,lin_wie) = p2(1,2)*exp(p2(1,1)*h(lin_wie,1));
    clear lin_wie kol_wie wart_wie
    modQ = repmat(modQ,num_lat*num_lon,1);
    modQ = 9*reshape(modQ,num_lat*num_lon*(num_lev-1),1);
end
if size(modQ,1)<size(A,2)
    %outer model currently over-constrained
    modQ(size(modQ,1)+1:size(A,2),1) = 0.00018;
end
Q = diag((modQ).^2);
