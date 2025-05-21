function [nodes_reg, nodes_reg_coord, nodes_reg_apr] = add_regular_selected(s_list,BLh2D_pudel_proj,Nw_apr,czas,stac,station,...
    BLhall,NAME,BLh_pudel_proj,BLh_pudel_proj_num,z1,L0,num_lat,num_lon,cut_off_angle,...
    switches,path_save)

% Find regular nodes in the close vicinity of the signals' paths

% nodes_reg_all=BLh2D_pudel_proj';
% Nw_apr_reg_all=Nw_apr(czas,:)';

nodes_reg_useful = [];
sat = 0;
b = 0;
% Setting the 0th station loop
for nr = 1:stac
    %setting the satelite loop
    sat_numer = length(station.h(nr).satellite);
    while sat < sat_numer
        sat = sat + 1;
        b = b + 1;
        if isempty(find(s_list==b,1))
            BLhst = station.h(nr).parameters;
            elevation = station.h(nr).satellite(sat).elevation;
            azi = station.h(nr).satellite(sat).azi;
            SWDt = station.h(nr).satellite(sat).SWD;
            
            [~,par_st] = inter_plane_line_LAB(BLhall,NAME,BLhst,elevation,azi,SWDt,BLh_pudel_proj,z1,L0,num_lat-1,num_lon-1,cut_off_angle,switches,path_save,b);
            %         [par_outer,par_st_outer] = inter_plane_line_LAB(BLhall,NAME,BLhst,elevation,azi,SWDt,BLh_outer_proj,z1,L0,3,3,cut_off_angle,switches,path_save,b);
            [~,~,~,~,~,~,mean_voxels] = set_num_pud_fast_notOnlyHor(BLh_pudel_proj_num,BLh_pudel_proj,par_st,switches);
            
            nodes_reg_sl = reshape(mean_voxels,numel(mean_voxels),1);
            nodes_coord_sl = BLh2D_pudel_proj([1 2 3],nodes_reg_sl)';
            plot_signal_nodes = 0;
            if plot_signal_nodes == 1
                scatter3(par_st(:,2),par_st(:,1),par_st(:,3),'filled','r');
                hold on;
                scatter3(nodes_coord_sl(:,2),nodes_coord_sl(:,1),nodes_coord_sl(:,3),'filled','k');
            end
            nodes_reg_useful = [nodes_reg_useful; nodes_reg_sl];
        end
    end
    sat=0;
end

nodes_reg=unique(nodes_reg_useful);

nodes_reg_coord = BLh2D_pudel_proj([1 2 3],nodes_reg)'; %First phi then lambda
% nodes_reg_coord = BLh2D_pudel_proj([2 1 3],nodes_reg)'; %First lambda then phi

nodes_reg_apr = Nw_apr(czas,nodes_reg)';
% scatter3(nodes_reg_coord(:,1),nodes_reg_coord(:,2),nodes_reg_coord(:,3),'filled');


end

