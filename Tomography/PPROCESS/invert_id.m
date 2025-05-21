function id_par= invert_id(model,indRT2,par)
    radii = [6378.137, 6356.752314245]; 
    dpar = diff(par);
    ps_par = par(1:end-1,1:3)+dpar(:,1:3)./2;
    for i = 1:size(ps_par,1)
       [lon(i),lat(i),h(i)] = cspice_recgeo(ps_par(i,:)',radii(1),(radii(1)-radii(2))/radii(1));
    end
    for i = 1:size(ps_par,1)
        [~,indRT(i)] = find_num(lat(i)*180/pi(),lon(i)*180/pi(),h(i)*1000,model.lat_TOMO,model.lon_TOMO,model.levels_TOMO,2);
    end
   % if sum(indRT2 - indRT) ~= 0
    %    disp('change of index')
    %end
    n_lat = length(model.mid_lat_TOMO);
    n_lon = length(model.mid_lon_TOMO);
    n_alt = length(model.mid_levels_TOMO);
    hor = n_lat*n_lon;
id_par = indRT;
   % for i = 1:size(indRT,2)
     %   altsize = idivide(indRT(i),int32(hor));
    %    i_pos = indRT(i) - altsize*hor;
    %    latsize = idivide(i_pos,int32(n_lat));
    %    lonsize = i_pos - (latsize-1)*n_lon;
    %    i_pos = (altsize)*(n_lon)*(n_lat) + lonsize + (n_lon)*(latsize-1);
        %i_pos = (n_alt1-1)*(length(lon_v)-1)*(length(lat_v)-1) + n_lon1 + (length(lon_v)-1)*(n_lat1-1);
%
   %     id_par(i) = (altsize)*(n_lon)*(n_lat) + latsize + (n_lat)*(lonsize-1);
   % end
end
