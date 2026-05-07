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

if ~isempty(names)
    j =  0;
    c = 0;

    % --- Build dateRO once ---
    date = cell(0,1);
    dateRO = cell(0,1);
    k = [];

    for fileIdx = 1:numel(names)
        if contains(names(fileIdx),'atmPhs') || contains(names(fileIdx),'conPhs') 
            if contains(names(fileIdx),'atmPhs') 
                atmPhs=true;
            else
                atmPhs=false;
            end
            j = j + 1;
            k(j) = fileIdx;
            spirename = char(names(fileIdx));
            year = ncreadatt([pathRO,'\', spirename],'/','year');
            month = ncreadatt([pathRO,'\', spirename],'/','month');
            day = ncreadatt([pathRO,'\', spirename],'/','day');
            hour = ncreadatt([pathRO,'\', spirename],'/','hour');
            minute = ncreadatt([pathRO,'\', spirename],'/','minute');
            second = ncreadatt([pathRO,'\', spirename],'/','second');
            date{j,1} = [year,month,day,hour,minute,round(second)];
            dateRO{j,1} = sprintf('%04d%02d%02d%02d',year,month,day,hour);
        end
    end

    % Find corresponding epochs
    for epoch = 1:size(obs_set.observation_set,1)
        dateEP = sprintf('%04d%02d%02d%02d', ...
            obs_set.observation_set(epoch,3), ...
            obs_set.observation_set(epoch,7), ...
            obs_set.observation_set(epoch,8), ...
            obs_set.observation_set(epoch,9));

        for roIdx = 1:numel(dateRO)
            if strcmp(dateRO{roIdx} ,dateEP )
                spirename = char(names(k(roIdx)));
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
                exL1 = ncread([pathRO,'\', spirename],'exL1');
                
                f2 = 1227.60; % default (MHz)

                if strcmp(switches.solution,'REAL')
                    try
                        occfreq2 = ncreadatt([pathRO,'\', spirename], '/', 'occfreq2');
                        f2 = occfreq2 / 1e6; % MHz
                    catch
                        warning('findRO: occfreq2 not found, using default f2=1227.60 MHz');
                    end
                end
                
                % Calculate excess phase on LC band
                exLC = exL1+(f2^2/(f1^2-f2^2))*(exL1-exL2);
                if sum(exLC)<0
                    exLCend = abs(exLC(end));
                    exLC = abs(flip(exLCend + exLC));
                end
                coordRO{epoch,roIdx,1} = [xGps, yGps, zGps];
                coordRO{epoch,roIdx,2} = [xLeo, yLeo, zLeo];
                coordRO{epoch,roIdx,3} = date{roIdx};
                coordRO{epoch,roIdx,4} = exL2;
                coordRO{epoch,roIdx,5} = exLC;
            end
        end
    end
    if ~exist('coordRO','var')
        coordRO = [];
        warning('findRO: No matching radiooccultation found')
    end
end