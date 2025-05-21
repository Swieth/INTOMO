function [ nodes_tomo_all,Nw_apr_nodes_all ] = add_regular_bottom( switches, BLh2D_pudel_proj,num_lat,num_lon,Nw_apr,czas,nodes_tomo,Np,num_inner)
%Add regular bottom layers to the irregular grid

number_of_layers = switches.number_of_layers;

only_inner = 1;

if only_inner==1
    BLh2D_inner = BLh2D_pudel_proj(:,num_inner);
    Nw_apr_inner = Nw_apr(1,num_inner);
    add_layer = BLh2D_inner(:,1:(num_lat-1)*(num_lon-1)*number_of_layers);
    apr_add = Nw_apr_inner(czas,1:(num_lat-1)*(num_lon-1)*number_of_layers);
else
    
    add_layer = BLh2D_pudel_proj(:,1:(num_lat+1)*(num_lon+1)*number_of_layers);
    apr_add = Nw_apr(czas,1:(num_lat+1)*(num_lon+1)*number_of_layers);
end
nodes_tomo_add_bottom = [nodes_tomo add_layer];
% nodes_tomo = nodes_tomo_add_bottom;
Np_add_bottom = [Np; apr_add'];
Np = Np_add_bottom;
[val ind] = sortrows(nodes_tomo_add_bottom',[3 2 1]);
Nw_apr_nodes_all = Np_add_bottom(ind);
nodes_tomo_all = val(:,[2 1 3]); %First lambda then phi
% clear add_layer apr_add nodes_tomo_add_bottom Np_add_bottom val ind number_of_layers
end

