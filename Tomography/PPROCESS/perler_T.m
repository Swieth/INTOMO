function [ Tmatrix,nodes_columns,Nwapr_columns ] = perler_T(levels,num_lev,BLh_pudel,Nw_apr,epoch)
% Perler bilinear-h - create matrix that will be used for the
% interpolation
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       levels.....    altitude layers of tomography model
%       num_lev....    number of voxels in vertical direction
%       BLh_pudel..    coordinates of voxel edges/centers
%       Nw_apr.....    apriori values of refractivities
%       epoch......    processing epoch number
%%%OUTPUT:
%       Tmatrix....    second derivates of refractivity for use in bilin interp
%       nodes_columns  segragated numbers of voxels for use in bilin interp
%       Nw_apr_columns apriori values of wet refractivity sorted by
%       tomography domain columns


dh = diff(levels);
%Set a, b, c, d parameters
a = dh;
b = 6./dh;
for k = 2:size(dh,1)
    c(k,1) = 2*(dh(k-1,1)+dh(k,1));
    d(k,1) = (6/dh(k-1,1))-(6/dh(k,1));
end
% Create C, D matrices
for k = 1:size(levels,1) % For each layer
    if k==1
        C(k,:) = [1 zeros(1,num_lev-1)];
        D(k,:) = zeros(1,num_lev);
    elseif k==num_lev
        C(k,:) = [zeros(1,num_lev-1) 1];
        D(k,:) = zeros(1,num_lev);
    else
        C(k,:) = [zeros(1,k-2) a(k-1,1) c(k,1) a(k,1) zeros(1,num_lev-(k-2)-3)];
        D(k,:) = [zeros(1,k-2) -1*b(k-1,1) d(k,1) b(k,1) zeros(1,num_lev-(k-2)-3)];
    end
end
deriv_matrix = (-1)*C^(-1)*D; %The matrix for computation of the second derivatives of refractivity
deriv_mat_one = [eye(num_lev); deriv_matrix];
clear dh a b c d C D k

% Save numbers of nodes (each vertical column separately)
nodes_columns=1:size(BLh_pudel,3)*size(BLh_pudel,2);
nodes_columns=reshape(nodes_columns,size(BLh_pudel,2),size(BLh_pudel,3))';
% Reshape Nw a priori (each vertical column separately)
Nwapr_columns=reshape(Nw_apr(epoch,:)',size(BLh_pudel,2),size(BLh_pudel,3))';
Tmatrix = zeros(num_lev*8,num_lev*4);
for num = 1:4
    Tmatrix((num_lev*2)*(num-1)+1:num_lev*2*num,(num_lev)*(num-1)+1:num_lev*num)=deriv_mat_one;
end
clear deriv_matrix deriv_mat_one
end

