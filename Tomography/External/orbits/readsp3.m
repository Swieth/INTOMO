function [SP3data,numsat,header] = readsp3_2(filename) 
 
% Read precise satellite orbit from SP3 file 
%  
% %%%%%%%%%% HELP %%%%%%%%%% 
%  
% [SP3data,numsat,header] = ReadSP3(filename) 
%  
% Input Data 
% filename = SP3 filename 
%  
% Output Data 
% SP3data = [GPSweek GPSsec PRN X Y Z clk] 
% numsat = number of satellite in SP3 file 
% header = header of SP3 file 
% *** X Y Z -> km 
%     clk -> microsec 
% 
% Written by  Phakphong Homniam 
% December 21, 2002 
%
% adapted to own needs, 02.2012 by Tomasz Hadas (reading all GNSS)
 
%%%%%%%%%% BEGIN %%%%%%%%%% 
 
% See if this file has already been read and formatted for MATLAB 
% save_file = sprintf('%s.mat', deblank(filename)); 
% if exist(save_file) == 2 
%     file_exists = 1; 
% else 
%     file_exists = 0; 
% end 
% if file_exists == 1 
%     load_string = sprintf('load %s', save_file); 
%     eval(load_string); 
%     fprintf('Data loaded from existing file, %s.\n', save_file); 
%     return 
% end 

gps_week = []; 
gps_sec_start = []; 
ep_interval = []; 
numsat = []; 
gps_sec = []; 
header = []; 
SP3data = []; 
time = []; 
gpstime = []; 

noom = 0;

%header
fid = fopen(filename,'r'); 
header_end=0;
while ~header_end
    line = fgetl(fid);
    if strcmp(line(1:2),'* ')
        header_end=1;
    else
        header = [header;line]; 
    end
end
n_header=size(header,1);
fclose(fid);

%data
fid = fopen(filename,'r');
for i = 1:n_header
    line = fgetl(fid); 
end

nrows = numel(textread(filename,'%1c%*[^\n]'));

textprogressbar(' ')
i=n_header;
while feof(fid) == 0
    %textprogressbar(num2str(i/nrows*100));
    i=i+1;
    line = fgetl(fid);
    if length(line) >= 60
        PRN_str2num(line(2:4));
        add = str2num([num2str(PRN_str2num(line(2:4))) line(5:18) line(19:32) line(33:46) line(47:60)]);
        [gps_week gps_sec] = utc2gps(time,0);
        SP3data = [SP3data;gps_week*7*24*3600+gps_sec add];
        %continue 
    elseif length(line) == 31 
        noom = noom+1; 
        time=str2num(line(4:end));
%         addtime = str2num(line(4:end));
%         time = [time;addtime];
        %continue 
    else 
        %continue 
    end 
end 
fclose(fid);
%textprogressbar('done')
% numsat = str2num(header(3,5:6)); 
% for i = 1:size(time,1) 
%     [gps_week gps_sec] = utc2gps(time(i,:),0); 
%     add_gpstime(1:numsat,1) = gps_week; 
%     add_gpstime(1:numsat,2) = gps_sec; 
%     gpstime = [gpstime;add_gpstime]; 
% end 

%%%%%%%% RESULT %%%%%%%% 
        
%SP3data = [gpstime(:,1)*7*24*3600+gpstime(:,2) SP3data];  
 
% save_string = sprintf('save %s', save_file); 
% eval(save_string); 
 
% %%%%%%%%%% END %%%%%%%%%% 