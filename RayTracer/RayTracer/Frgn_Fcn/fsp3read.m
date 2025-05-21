function [xyzGps, eGps, idGps] = fsp3read(sp3path, sp3file, source)
%% ABOUT: 19.08.2016
% sp3read reads orbit data for given input file. Can work either with URL
% data or local files - defined by source variable.

% Example for orbit source as an URL:
% sp3path = 'ftp://ftp.unibe.ch/aiub/CODE/';
% sp3file = 'COD.EPH_U';
%
% Input data: 
%               sp3path: path to orbit data           [string]
%               sp3file: name of orbit file           [string]
%               source:  'u' = url, 'f' = file        [string]
% Output data: 
%               xyzGps:  satellite coordinates        [number]
%               eGps:    YYYY, MM, DD, hh, mm, ss     [number]
%               idGps:   satellites IDs               [string]
%% READ SP3 DATA

if source == 'u'
    % Read content of source file
    contents = urlread([sp3path,sp3file]);
    % Read data to cell
    sp3r = textscan(contents, '%s', 'Delimiter', '\n');
    % Data array
    sp3 = sp3r{1,1};
elseif source == 'f'
    cd(sp3path);
    % Check if query file exists
    A = exist(sp3file, 'file');
    if A == 2
        fid = fopen(sp3file, 'r');
        % Read data to cell
        sp3r = textscan(fid, '%s', 'Delimiter', '\n');
        % Data array
        sp3 = sp3r{1,1};
        fclose(fid);
    else warning('no such file in directory');
    end
end

% Determine numer of columns in cell
[~, ncols] = cellfun(@size, sp3);

% Find colomuns with coordinates (mP) and epoch (mE): alternative - strfind(sp3, '*')
[mP, ~] = find(ncols == 80 | ncols == 60);
[mE, ~] = find(ncols == 31);

% Remove header in 60char lines
rows = mP > min(mE);
mP = mP(rows);

% Create array with cordinates (P rows)
for i = 1:numel(mP)
    % Satellite numbers: start from value 2 to avoid 'P'
    idGps(i,:) = (sp3{mP(i),1}(1,2:4));
    
    % Satellites coordinates: start from value 5 to idGps
    xyzGps(i,:) = str2num(sp3{mP(i),1}(1,5:46));
end

% Create array with epoch (* rows)
for i = 1:numel(mE)
    % Epoch for satellites: start from value 3 to avoid '*  '  
    E(i,:) = str2num(sp3{mE(i),1}(1,3:31));
end

% Duplicate epoch rows for consistent size of E and xyzGps
for i = 1:size(mE,1)
    % Find xyzGps rows for each epoch
    if i < size(mE,1)
        [m, ~] = find(mP(:,1) > mE(i,1) & mP(:,1) < mE(i+1,1));
    else
        [m, ~] = find(mP(:,1) > mE(i,1));
    end
    % Repeat rows with epoch
    eGps(m,:) = repmat(E(i,:), [numel(m),1]);
end

end