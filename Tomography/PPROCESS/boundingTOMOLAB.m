function [X,Y,Z,lat,lon,h,H,NAME]=boundingTOMOLAB(east_limit,west_limit,north_limit,south_limit,Xsta,BLh_ori,NAME)
% Function to cut GNSS station outside the bounding model
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       east_limit.   eastern boundary longitude
%       west_limit.   western boundary longitude
%       north_limit   northern boundary latitude
%       south_limit   southern boundary latitude
%       Xsta          XYZ coordinates of GNSS stations
%       BLh_ori       BLh coordinates of GNSS stations
%%%OUTPUT
%       X..........   X coordinate of GNSS station
%       Y..........   Y coordinate of GNSS station
%       Z..........   Z coordinate of GNSS station
%       lat........   latitude of GNSS station
%       lon........   longitude of GNSS station
%       h..........   altitude of GNSS station
%       H..........   altitude of GNSS station
%       NAME.......   name of GNSS station
    %       number.....   id of GNSS stationX = Xsta(1,:)';
     X = Xsta(1,:)';
     Y = Xsta(2,:)';
     Z = Xsta(3,:)';
     lat = BLh_ori(:,2)*pi()/180;
     lon = BLh_ori(:,3)*pi()/180;
     h = BLh_ori(:,5);
     H = BLh_ori(:,4);
     [wie_east, ~] = find(lon*180/pi()>east_limit);
     [wie_west, ~] = find(lon*180/pi()<west_limit);    
     [wie_north, ~] = find(lat*180/pi()>north_limit);    
     [wie_south, ~] = find(lat*180/pi()<south_limit);
     wie_1 = union(wie_east,wie_west);
     wie_2 = union(wie_south,wie_north);
     wie = union(wie_1,wie_2);
     lon(wie,:) = [];
     lat(wie,:) = [];
     h(wie,:) = [];
     H(wie,:) = [];
     X(wie,:) = [];
     Y(wie,:) = [];
     Z(wie,:) = [];
     NAME(wie,:)=[];