function [ nodes_tomo_all,Nw_apr_nodes_all ] = add_regular_boundary( BLh2D_pudel_proj,Nw_apr,czas,nodes_tomo_all,Nw_apr_nodes_all)
% Add boundary voxels to the irregular grid

BLh=BLh2D_pudel_proj';

B = sort(unique(BLh(:,1)));
L = sort(unique(BLh(:,2)));
Bmin_one = B(1,1);
Bmax_one = B(end,1);
Lmin_one = L(1,1);
Lmax_one = L(end,1);

Bmin = find(abs(BLh(:,1)-Bmin_one)<0.00000001);
Bmax = find(abs(BLh(:,1)-Bmax_one)<0.00000001);
Lmin = find(abs(BLh(:,2)-Lmin_one)<0.00000001);
Lmax = find(abs(BLh(:,2)-Lmax_one)<0.00000001);

boundary_nodes = unique([Bmin; Bmax; Lmin; Lmax]);
boundary_BLh = BLh(boundary_nodes,:);
boundary_BLh = boundary_BLh(:,[2 1 3]);
boundary_Nw = Nw_apr(czas,boundary_nodes)';

[~, indA, indB] = intersect(nodes_tomo_all,boundary_BLh,'rows');
boundary_nodes(indB,:)=[];
boundary_BLh(indB,:)=[];
boundary_Nw(indB,:)=[];

nodes_tomo_add_bound = [nodes_tomo_all; boundary_BLh];
Np_add = [Nw_apr_nodes_all; boundary_Nw];
[val ind] = sortrows(nodes_tomo_add_bound,[3 2 1]);

Nw_apr_nodes_all = Np_add(ind);
nodes_tomo_all = val;

end

