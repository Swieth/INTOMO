function [X,Y,Z,lat,lon,h,H,NAME] =  deleteStat(X,Y,Z,lat,lon,h,H,NAME,stat_range)
% Function to filter out stations to selected spatial density
%---------------------------------------------------------------------------------------------------
%%%INPUT/OUTPUT
%       X..........   X coordinate of GNSS station
%       Y..........   Y coordinate of GNSS station
%       Z..........   Z coordinate of GNSS station
%       lat........   latitude of GNSS station
%       lon........   longitude of GNSS station
%       h..........   altitude of GNSS station
%       H..........   altitude of GNSS station
%       NAME.......   name of GNSS station
%       number.....   id of GNSS station

    nStations = length(X);
    keepStation = true(nStations, 1);  
    coords = [X, Y, Z];
    distances = squareform(pdist(coords));
    for i = 1:nStations
        if keepStation(i)
            closeStations = distances(i, :) < stat_range;
            closeStations(i) = false; 
            keepStation(closeStations) = false;
        end
    end
    X = X(keepStation);
    Y = Y(keepStation);
    Z = Z(keepStation);
    lat = lat(keepStation);
    lon = lon(keepStation);
    h = h(keepStation);
    H = H(keepStation);
    NAME = NAME(keepStation);
end