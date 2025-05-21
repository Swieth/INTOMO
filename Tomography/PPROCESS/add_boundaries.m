function [ nodes_tomo_all,mean_voxels,Nw_apr_nodes_all ] = add_boundaries( BLh_pudel_proj_num,BLh_pudel_proj,nodes_tomo,switches,...
                                                                              BLh2D_pudel_proj,Nw_apr,czas,planes,nodes_columns,Nwapr_columns,...
                                                                              planes_nr,Tmatrix,num_lat,num_lon,num_inner)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% Get numbers of nodes for a priori interpolation from grid
[~,~, ~, ~, ~, mean_voxels] = set_num_pud_nodes(BLh_pudel_proj_num,BLh_pudel_proj,nodes_tomo,switches);

% Calculate ALADIN data in nodes (for a priori and for validation); bilinear/spline interpolation from the regular grid
Np=apriori_irregular(mean_voxels,BLh2D_pudel_proj,Nw_apr,czas,nodes_tomo,planes,nodes_columns,Nwapr_columns,planes_nr,Tmatrix);
% Add bottom layers of nodes
if switches.number_of_layers>0
    [nodes_tomo,Np] = add_regular_bottom( switches, BLh2D_pudel_proj,num_lat,num_lon,Nw_apr,czas,nodes_tomo,Np,num_inner);
end
% Add top layers of nodes
if switches.number_of_layers_top>0
    [nodes_tomo,Np] = add_regular_top( switches,BLh2D_pudel_proj,num_lat,num_lon,Nw_apr,czas,nodes_tomo,Np);
end
Nw_apr_nodes_all=Np;
% Add boundry nodes
if switches.add_boundary_nodes>0
    [nodes_tomo,Nw_apr_nodes_all]=add_regular_boundary(BLh2D_pudel_proj,Nw_apr,czas,nodes_tomo,Nw_apr_nodes_all);
end
nodes_tomo_all=nodes_tomo;
end

