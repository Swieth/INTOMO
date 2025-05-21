%%Function for voxel horizontal gradient computation 
%% Input Data
% lat.......... ray point latitude
% lon.......... ray point longitude
% m1........... index of height layer above ray point
% m2........... index of height layer below ray point
% dh3.......... height level difference in ray point's voxel [m]
% dh4.......... hegiht level difference in voxel below or above ray point [m]
% lon_ray_b.... previous ray point longitude 
% lat_ray_b.... previous  ray point latitude 
% set.......... settingd
% inp_data
%% Output
% dh3..........  distance between nodes in lonitude direction [m]
% dh4..........  distance between nodes in latitude direction [m]
% dref.........  refraction value difference in lonitude direction
% dref2........  refraction value difference in latitude direction

function [drefm1,drefm2,dh3,dh4] = dhrefcalc(lat,lon,m1,m2,dh3,dh4,lon_ray_b,lat_ray_b,set,inp_data)
    lat = lat*180/pi();
    lon = lon*180/pi();
    lat_ray_b = lat_ray_b*180/pi();
    lon_ray_b = lon_ray_b*180/pi();
    if lon < inp_data.lon(1,1,1)
        lon_new = lon +360 - inp_data.lon(1,1,1);
    else
        lon_new = lon - inp_data.lon(1,1,1);
    end 
    if min(inp_data.lat(1,:,1)) < 0 && set.refron
         lat_new = 180-lat;
    elseif lat < 90
         lat_new = 90-lat;
    else
        lat_new = lat;
    end
    [m5, ~] = find(inp_data.lon_new(:,1,1) > lon_new);
    [m6, ~] = find(inp_data.lon_new(:,1,1) < lon_new);
    [~, m7] = find(inp_data.lat_new(1,:,1) > lat_new);
    [~, m8] = find(inp_data.lat_new(1,:,1) < lat_new);
    if isempty(m5)
        m5 = find(max(inp_data.lon_new(:,1,1)));
    else
        m5 = m5(1);
    end
    if isempty(m6)
        m6 =find(min(inp_data.lon_new(:,1,1)));
    else
        m6 = m6(end);
    end
    if isempty(m7)
        m7 = find(max(inp_data.lat_new(1,:,1)));
    elseif ~set.ground 
       m7 = m7(1);
    elseif set.ground   
       m7 = m7(end);
    end
    if isempty(m8)
        m8 = find(min(inp_data.lat_new(1,:,1)));
    elseif ~set.ground 
        m8 = m8(end);
    elseif set.ground
        m8 = m8(1);
    end
    if m6 > m5
        dref3 = diff(inp_data.ref(m5:m6,m7,m1));   
        dref4 = diff(inp_data.ref(m5:m6,m7,m2));
        dref5 = diff(inp_data.ref(m5:m6,m8,m1));   
        dref6 = diff(inp_data.ref(m5:m6,m8,m2));
    else
        dref3 = diff(inp_data.ref(m6:m5,m7,m1));   
        dref4 = diff(inp_data.ref(m6:m5,m7,m2));
        dref5 = diff(inp_data.ref(m6:m5,m8,m1));   
        dref6 = diff(inp_data.ref(m6:m5,m8,m2));
    end
    drefm1 = abs(mean(dref3 + dref4+dref5 + dref6));
    if lat_new < 90 && lat_ray_b > lat
      drefm1 = -drefm1;  
    elseif lat_new > 90 && lat_ray_b < lat
      drefm1 = -drefm1;   
    end
    lon1 = inp_data.lon_new(m5,1,1);
    lon2 = inp_data.lon_new(m6,1,1);
    if lon < inp_data.lon(1,1,1)
        lon1= lon1 - 360 + inp_data.lon(1,1,1);
        lon2= lon2 - 360 + inp_data.lon(1,1,1);
    else
        lon1 = lon1 + inp_data.lon(1,1,1);
        lon2 = lon2 + inp_data.lon(1,1,1);
    end
    if lat > 90
        lat = -(lat - 90);
    else
        lat = 90 - lat;
    end
    try
        dh3 = vdist(lat,lon1,lat,lon2);
    end            
    if m8 > m7
        dref3 = diff(inp_data.ref(m5,m7:m8,m1));   
        dref4 = diff(inp_data.ref(m5,m7:m8,m2));
        dref5 = diff(inp_data.ref(m6,m7:m8,m1));   
        dref6 = diff(inp_data.ref(m6,m7:m8,m2));
    else
        dref3 = diff(inp_data.ref(m5,m8:m7,m1));   
        dref4 = diff(inp_data.ref(m5,m8:m7,m2));
        dref5 = diff(inp_data.ref(m6,m8:m7,m1));   
        dref6 = diff(inp_data.ref(m6,m8:m7,m2));
    end
   drefm2 = abs(mean(dref3 + dref4+dref5 + dref6));
    if lon_new < 180 && lon_ray_b > lon
      drefm2 = -drefm2;  
    elseif lon_new > 180 && lon_ray_b < lon
      drefm2 = -drefm2;   
    end
   lat1 = inp_data.lat_new(1,m7,1);
   lat2 = inp_data.lat_new(1,m8,1);
   if min(inp_data.lat(1,:,1)) < 0
        lat1 = 180 - lat1 ;
        lat2 = 180 - lat2 ;
        if lat1>90
            lat1 = -(lat1 - 90);
            lat2 = -(lat2 - 90);
        elseif lat1<90
            lat1 = 90 - lat1;
            lat2 = 90 - lat2;
        end
   end
   try
        dh4 = vdist(lat1,lon,lat2,lon);
   end
end