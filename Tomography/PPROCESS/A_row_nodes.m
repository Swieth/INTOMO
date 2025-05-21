function [ SWD_sum_all,A_rows,p_coord_sect,Np_sect,SWD2_p,p_Nw ] = A_row_nodes( switches,nodes_tomo_all,p_coord,p,Nw_apr_nodes_all,mean_voxels,sect,Nw_apr,czas,const,SWD_sum_all,A_rows,p_coord_sect,Np_sect,mult_fact,A_nodes,p_Nw )
% Create a row that will be put into the A matrix
% For the solution with switches.parametrization set to 'nodes'

% Get nearest nodes for interpolation
if strcmp(switches.regular,'no')
    [nearest_nodes_Nw,nearest_nodes_dist,nearest_nodes_nr] = get_nearest_nodes(switches,nodes_tomo_all,p_coord,p,Nw_apr_nodes_all);
elseif strcmp(switches.regular,'yes')
    [nearest_nodes_Nw,nearest_nodes_dist,nearest_nodes_nr,nearest_nodes_Nw_NEW,ommm,nearest_nodes_nr_NEW] = get_nearest_nodes_regular(mean_voxels,sect,nodes_tomo_all,p_coord,p,Nw_apr,czas);
end

% Plot the points on signal and neighbouring nodes
plot_points=0;
if plot_points==1
    if strcmp(switches.regular,'yes')
        scatter3(p_coord(p,1),p_coord(p,2),p_coord(p,3),[],'r','filled');
    elseif strcmp(switches.regular,'no')
        scatter3(p_coord(p,1),p_coord(p,2),p_coord(p,3),[],'r','filled');
    end
    hold on;
    scatter3(nodes_tomo_all(nearest_nodes_nr_NEW,1),nodes_tomo_all(nearest_nodes_nr_NEW,2),nodes_tomo_all(nearest_nodes_nr_NEW,3),[],'b','filled');
end

% Interpolation according to Perler
Np(p,1) = prod(ommm,2)'*nearest_nodes_Nw;
om(p,:) = prod(ommm,2);
om_mult(p,:) = prod(ommm,2).*mult_fact(1,p).*const;

% Save value of Nw in p (for diagnostic)
p_Nw(p,1) = om(p,:)*nearest_nodes_Nw_NEW;

% Calculate the SWD integral without A matrix (for diagnostics)
SWD_sum_all = SWD_sum_all + om_mult(p,:)*nearest_nodes_Nw_NEW;

SWD2_p = om_mult(p,:)*nearest_nodes_Nw_NEW;

% A matrix based on interpolation by Perler
A_row=zeros(1,size(A_nodes,2));
A_row(1,nearest_nodes_nr_NEW)=om_mult(p,:);
A_rows = [A_rows; A_row];

p_coord_sect = [p_coord_sect; p_coord(p,:)];
Np_sect = [Np_sect; Np(p,1)];
end

