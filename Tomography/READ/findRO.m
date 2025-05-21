function [coordRO] = findRO(pathRO,obs_set)
%% Function to find RO matching to dates of simulated ground-based tomography
% Input
%   pathRO.........  path to folder with RO.nc files
%   obs_set........  matrix with dates of epochs
% Output
%   coordRO........  structure of RO data matching to obs_set epochs
%          .coordT               : ccordinates of the transmitter (ECEF)
%          .coordR               : ccordinates of the receiver (ECEF)
%          .date                 : date of RO start
%          .exL2                 : RO excess phase on L2
%          .exLC                 : RO excess phase on LC    


% Fing RO files
folderpath = fullfile(pathRO, '**');  
filelist  = dir(folderpath);
filelist = filelist(3:end); %
names       = {filelist.name};
f1 = 1575.42;
f2 = 1227.60;
if ~isempty(names)
    j =  0;
    c = 0;
    % Read data
    for i = 1:size(names,2)
        if contains(names(i),'atmPhs') || contains(names(i),'conPhs') 
            if contains(names(i),'atmPhs') 
                atmPhs=true;
            else
                atmPhs=false;
            end
            j = j + 1;
            k(j) = i;
            spirename = char(names(i));
            year = ncreadatt([pathRO,'\', spirename],'/','year');
            month = ncreadatt([pathRO,'\', spirename],'/','month');
            day = ncreadatt([pathRO,'\', spirename],'/','day');
            hour = ncreadatt([pathRO,'\', spirename],'/','hour');
            minute = ncreadatt([pathRO,'\', spirename],'/','minute');
            second = ncreadatt([pathRO,'\', spirename],'/','second');
            date{j,:} = [year,month,day,hour,minute,round(second)];
            year = num2str(year) ;
            month =num2str(month) ;
            day = num2str(day) ;
            hour = num2str(hour) ;
            dateRO{j,:} = char([year,month,day,hour]);
        end
        %Find corresponding epochs
        for epoch = 1:size(obs_set.observation_set,1)
            dateEP = [num2str(obs_set.observation_set(epoch,3)),num2str(obs_set.observation_set(epoch,7)),num2str(obs_set.observation_set(epoch,8)),num2str(obs_set.observation_set(epoch,9))];
            for i = 1:size(dateRO,1)
                if strcmp(dateRO{i,:} ,dateEP )
                    spirename = char(names(k(i)));
                    if contains(spirename,'atmPhs') 
                        atmPhs=true;
                    else
                        atmPhs=false;
                    end
                    if atmPhs
                        xLeo = ncread([pathRO,'\', spirename],'xLeo');
                        yLeo = ncread([pathRO,'\', spirename],'yLeo');
                        zLeo = ncread([pathRO,'\', spirename],'zLeo');
                        xGps = ncread([pathRO,'\', spirename],'xGps');
                        yGps = ncread([pathRO,'\', spirename],'yGps');
                        zGps = ncread([pathRO,'\', spirename],'zGps');
                    else
                        xLeo = ncread([pathRO,'\', spirename],'xLeoLR');
                        yLeo = ncread([pathRO,'\', spirename],'yLeoLR');
                        zLeo = ncread([pathRO,'\', spirename],'zLeoLR');
                        xGps = ncread([pathRO,'\', spirename],'xGnssLR');
                        yGps = ncread([pathRO,'\', spirename],'yGnssLR');
                        zGps = ncread([pathRO,'\', spirename],'zGnssLR');
                    end
                    exL2 = ncread([pathRO,'\', spirename],'exL2');
                    exL1 = ncread([pathRO,'/', spirename],'exL1');
                    % Calculate excess phase on LC band
                    exLC = exL1+(f2^2/(f1^2-f2^2))*(exL1-exL2);
                    if sum(exLC)<0
                        exLCend = abs(exLC(end));
                        exLC = abs(flip(exLCend + exLC));
                    end
                    coordRO{epoch,i,1,:} = [xLeo, yLeo, zLeo];
                    coordRO{epoch,i,2,:} = [xGps, yGps, zGps];
                    coordRO{epoch,i,3,:} = date{i,:};
                    coordRO{epoch,i,4,:} = exL2;
                    coordRO{epoch,i,5,:} = exLC;
                end
            end
        end
    end
    if ~exist('coordRO','var')
        coordRO = [];
        warning('findRO: No matching radiooccultation found')
    end
end