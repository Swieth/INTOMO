function  unduera = Undulation(lat1,lat2,lon1,lon2,n,path,name)
%% Function to generate undulation file for 3DRT software
% lat1, lat2 ........... min and max latitude of area boundary points
% lon1, lon2 ........... min and max longitude of area boundary points 
% n .......... grid resolutions [degrees]
% path ......... path to save

    dlat = abs(lat1 - lat2);
    dlon = abs(lon1 - lon2);
    for i = 0:dlat/n 
        for j = 0:dlon/n 
            unduera(j+1,i+1) = geoidheight2(max(lat1,lat2)-i*n,min(lon1,lon2)+j*n,'egm2008');
        end
    end
    save(strcat(path,'\',name,'.mat'),'unduera')
end
