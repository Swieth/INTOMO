 function [NAME,BLh] = readBLh(filename)
 % Function to read BLH with formatted format
 %%% INPUT
 %      filename...... path to .txt file
 %%%OUTPUT
 %      NAME.......... names of GNSS stations
 %      BLh...........coordinates of GNSS stations
delimiter = '\t';
formatSpec = '%s%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
T = table(dataArray{1:end-1}, 'VariableNames', {'Name','num','lat','lon','h','H'});
 BLh(:,1) = T.num;
 BLh(:,2) = T.lat;
 BLh(:,3) = T.lon;
 BLh(:,4) = T.h;
 BLh(:,5) = T.H;
 NAME = T.Name;
 
 