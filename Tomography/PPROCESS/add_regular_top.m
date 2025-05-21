function [ nodes_tomo_all,Nw_apr_nodes_all ] = add_regular_top( switches,BLh2D_pudel_proj,num_lat,num_lon,Nw_apr,czas,nodes_tomo,Np)
%Add regular top layers to the irregular grid

number_of_layers_top = switches.number_of_layers_top;
nr = size(BLh2D_pudel_proj,2);
add_layer = BLh2D_pudel_proj(:,(nr-(num_lat+1)*(num_lon+1)*number_of_layers_top)+1:end);
apr_add = Nw_apr(czas,(nr-(num_lat+1)*(num_lon+1)*number_of_layers_top)+1:end);
nodes_tomo_add_top = [nodes_tomo add_layer];
Np_add = [Np; apr_add'];
[val ind] = sortrows(nodes_tomo_add_top',[3 2 1]);
Nw_apr_nodes_all = Np_add(ind);
nodes_tomo_all = val(:,[2 1 3]); %First lambda then phi
% clear add_layer apr_add nodes_tomo_add_top Np_add val ind Np number_of_layers
end

