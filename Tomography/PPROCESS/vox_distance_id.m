function [X_alt,X_lon,X_lat,distAlt,distLon,distLat] = vox_distance_id(alt_p1,alt_p2,pos_X2,lon_p1,lon_p2,lat_p1,lat_p2,lat_v,lon_v,alt_v,X1,X2,radii)
 %Select layer index, project in altitude direction
                if alt_p1 <= alt_p2
                    idalt = min(pos_X2(5),pos_X2(6));
                else
                    idalt = max(pos_X2(5),pos_X2(6));
                end
                [X_T1] = cspice_georec(lon_p2*pi()/180, lat_p2*pi()/180, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_T2] = cspice_georec(lon_p1*pi()/180, lat_p1*pi()/180, alt_v(idalt)/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_alt,~] = lineIntersect3D([X_T1';X2'],[X_T2';X1']);
                % Calculate distance
                distAlt = vecnorm(X_alt' - X2); 
                 %Select layer index, project in longitude direction
                if lon_p1 <= lon_p2
                    idlon = min(pos_X2(3),pos_X2(4));
                else
                    idlon = max(pos_X2(3),pos_X2(4));
                end
                [X_T1] = cspice_georec(lon_v(idlon)*pi()/180, lat_p2*pi()/180, alt_p2/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_T2] = cspice_georec(lon_v(idlon)*pi()/180, lat_p1*pi()/180, alt_p1/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_lon,~] = lineIntersect3D([X_T1';X2'],[X_T2';X1']);  
                 % Calculate distance
                distLon = vecnorm(X_lon' - X2); 
                %Select layer index, project in latitude direction
                if lat_p1 <= lat_p2
                    idlat = min(pos_X2(1),pos_X2(2));
                else
                    idlat = max(pos_X2(1),pos_X2(2));
                end
                [X_T1] = cspice_georec(lon_p2*pi()/180, lat_v(idlat)*pi()/180, alt_p2/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_T2] = cspice_georec(lon_p1*pi()/180, lat_v(idlat)*pi()/180, alt_p1/1000, radii(1),(radii(1)-radii(2))/radii(1));
                [X_lat,~] = lineIntersect3D([X_T1';X2'],[X_T2';X1']); 
                % Calculate distance
                distLat = vecnorm(X_lat' - X2); 
                % Select coordinates of intersection points and distance by
                % shortest path to projected point
end