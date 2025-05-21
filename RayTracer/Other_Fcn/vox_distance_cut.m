function [voxEM,distanceEm,i_pos,pos]  = vox_distance_cut(voxEM,i,distanceEm,i_pos,pos)
    [~,id] = min(voxEM.dist2tar(i,1:end-1));
    voxEM.dist2tar(i,id+1:end) = 0;
    voxEM.id_alt(i,id+1:end) = 0;
    voxEM.id_lon(i,id+1:end) = 0;
    voxEM.id_lat(i,id+1:end) = 0;
    distanceEm(i,id+1:end) = 0;
    i_pos(i,id+1:end) = 0;
    pos(i,:,id+1:end) = 0;
    voxEM.missed(i) = i;
end