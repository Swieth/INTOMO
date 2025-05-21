function [ eNwp, erafn ] = fera5Epoch(eGps)
%[YYYY, MM, DD, leap hours, wrf epoch (row: 1 -49), wrf forecast run(0,6,12,18)]
%%  ABOUT
%       fera5Epoch converts orbit epochs from sp3 file to NWP epochs that fit
%       ERA5 model forecast files
%       Example of input datetime epochs:
%       eGps = [YYYY, MM, DD, hh, mm, ss]

%       Example of output WRF epochs:
%       eNwp = [YYYY, MM, DD, ERA5 epoch (0-23), layer index]
%
%       Input data: 
%                   eGps:    YYYY, MM, DD, hh, mm, ss     [number]
%       Output data: 
%               eNwp:    ERA dates                    [number]
%               erafn:   forecast file                [string]
 
hour = eGps(:,4); 
minute = eGps(:,5);
second = eGps(:,6);
% hour and minutes to decimal
hours = double(hour + minute./60 + second./3600);
% round decimal hours to find epoch
hoursRound = round(hours);
% initialize output array with ERA5 epochs
eNwp(:,1:4) = [eGps(:,1:3), hoursRound];
% find hours equal to 24
ind = find(hoursRound == 24);
% update to 1 epoch
eNwp(ind, 4) = 1;
% add one day to dates with 24 epoch
t2 = datenum( [ eGps(ind,1), eGps(ind,2), eGps(ind,3) ] );
t1 = datevec(t2+1);
eNwp(ind, 1:3) = t1(:,1:3);
eNwp(:,5) = eNwp(:,4) + 1;

% ERA filename pattern
erafn = strcat({'ERA5_'}, int2str(eNwp(:,1)), {'-'}, num2str(eNwp(:,2),'%02d'), {'-'}, num2str(eNwp(:,3),'%02d'), {'.nc'});
end

