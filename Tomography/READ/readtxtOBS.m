function [ZTDA, MZTDA, DGNA, MDGNA, DGEA, MDGEA,NAMES,test] = readtxtOBS(pathATM,NAME, observation_set)
%% Read ZTD and gradients for observation files
%   INPUT:
%       pathATM...   path to obs files
%       NAME......   GNSS station names
%       observation_set   matrix with dates of epochs
% OUTPUT:
%       ZTDA......   matrix of zenith total delays
%       MZTDA.....   matrix of uncertainties of zenith total delays
%       DGNA......   matrix of N gradients of ZTD
%       MDGNA.....   matrix of uncertainties N gradient of ZTD
%       DGEA......   matrix of E gradients of ZTD
%       MDGEA.....   matrix of uncertainties E gradient of ZTD
%       NAMES.....   matrix of names of GNSS station
%       test......   id of stations for which observations are available

folderpath = fullfile(pathATM, '**');  
filelist  = dir(folderpath);
filelist = filelist(3:end); %
names       = {filelist.name};
c = 1;
s = 0;
A = datetime(observation_set(1,3),observation_set(1,7),observation_set(1,8));
epoch = 0;
for i = 1:size(names,2)
    name = char(names(i));
    year = str2num(name(1:4));
    month = str2num(name(5:6));
    day = str2num(name(7:8));
    B = datetime(year,month,day);
    if A == B
        pathATMdir = [pathATM,name];
        s = i;
        break
    end
    clear month day year name

end
if s ~= 0
    if size(observation_set,1) > 24 && c == 1
        days = idivide(size(observation_set,1),int32(24));
    else
        days = 1;
    end
    for ep = 1:days
        name = char(names(ep));
        pathATMdir = [pathATM,name];
        delimiter = '\t';
        startRow = 2;

        %% Read columns of data as text:
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%s%*s%s%*s%*s%*s%*s%s%s%s%s%s%s%[^\n\r]';

        %% Open the text file.
        fileID = fopen(pathATMdir,'r');

        %% Read columns of data according to the format.
        % This call is based on the structure of the file used to generate this
        % code. If an error occurs for a different file, try regenerating the code
        % from the Import Tool.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

        %% Close the text file.
        fclose(fileID);

        %% Convert the contents of columns containing numeric text to numbers.
        % Replace non-numeric text with NaN.
        raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
        for col=1:length(dataArray)-1
            raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
        end
        numericData = NaN(size(dataArray{1},1),size(dataArray,2));

        for col=[3,4,5,6,7,8]
            % Converts text in the input cell array to numbers. Replaced non-numeric
            % text with NaN.
            rawData = dataArray{col};
            for row=1:size(rawData, 1)
                % Create a regular expression to detect and remove non-numeric prefixes and
                % suffixes.
                regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                try
                    result = regexp(rawData(row), regexstr, 'names');
                    numbers = result.numbers;

                    % Detected commas in non-thousand locations.
                    invalidThousandsSeparator = false;
                    if numbers.contains(',')
                        thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                        if isempty(regexp(numbers, thousandsRegExp, 'once'))
                            numbers = NaN;
                            invalidThousandsSeparator = true;
                        end
                    end
                    % Convert numeric text to numbers.
                    if ~invalidThousandsSeparator
                        numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                        numericData(row, col) = numbers{1};
                        raw{row, col} = numbers{1};
                    end
                catch
                    raw{row, col} = rawData{row};
                end
            end
        end


        %% Split data into numeric and string columns.
        rawNumericColumns = raw(:, [3,4,5,6,7,8]);
        rawStringColumns = string(raw(:, [1,2]));


        %% Replace non-numeric cells with NaN
        R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
        rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

        %% Make sure any text containing <undefined> is properly converted to an <undefined> categorical
        idx = (rawStringColumns(:, 2) == "<undefined>");
        rawStringColumns(idx, 2) = "";

        %% Create output variable
        Table = table;
        Table.ID = rawStringColumns(:, 1);
        Table.Date = categorical(rawStringColumns(:, 2));
        Table.ZTD = cell2mat(rawNumericColumns(:, 1));
        Table.mZTD = cell2mat(rawNumericColumns(:, 2));
        Table.GradN = cell2mat(rawNumericColumns(:, 3));
        Table.GradE = cell2mat(rawNumericColumns(:, 4));
        Table.mGradN = cell2mat(rawNumericColumns(:, 5));
        Table.mGradE = cell2mat(rawNumericColumns(:, 6));
        ind = [];
        for i = 1:size(NAME,1)
            id = find(Table.ID == NAME(i));
            ind = [ind;id];
        end
        Table = Table(ind,:);
        clear id ind
        date = Table.Date;
        [~,idx] = sort(date);
        Table = Table(idx,:);
        clear idx
        if ep == 1
            ZTDA  = [];
            MZTDA = [];
            DGNA = [];
            DGEA = [];
            MDGNA = [];
            MDGEA = [];
        end
        aa = char(Table.Date);
        for k = 1:24
            epoch = epoch + 1;
            epochdate = datetime(observation_set(epoch,3),observation_set(epoch,7),observation_set(epoch,8),observation_set(epoch,9),0,0);
            try
                id = find(datetime(str2num(aa(:,1:4)),str2num(aa(:,6:7)),str2num(aa(:,9:10)),str2num(aa(:,12:13)),str2num(aa(:,15:16)),0) == epochdate);
            catch
                disp(a)
            end
                Tableepoch = Table(id,:);
            [idname,ia,ib] = intersect(Tableepoch.ID,NAME);
            ZTD = nan(1,size(NAME,1));
            MZTD= nan(1,size(NAME,1));
            DGN = nan(1,size(NAME,1));
            DGE = nan(1,size(NAME,1));
            MDGN= nan(1,size(NAME,1));  
            MDGE= nan(1,size(NAME,1));
            ZTD(ib) = Tableepoch.ZTD'; 
            MZTD(ib) = Tableepoch.mZTD';
            DGN(ib) = Tableepoch.GradN';
            DGE(ib) = Tableepoch.GradE';
            MDGN(ib) = Tableepoch.mGradN';
            MDGE(ib) = Tableepoch.mGradE';
            ZTDA = [ZTDA;ZTD];
            MZTDA = [MZTDA;MZTD];
            DGNA = [DGNA;DGN];
            DGEA = [DGEA;DGE];
            MDGNA = [MDGNA;MDGN];
            MDGEA = [MDGEA;MDGE];
            clear Tableepoch id ZTD MZTD DGN DGE MDGN MDGE
        end
        NAMES = unique(Table.ID);
        if length(NAMES)~=length(NAME)
            test = contains(NAME,NAMES);
            ZTDA = ZTDA(:,test); MZTDA= MZTDA(:,test); DGNA = DGNA(:,test); MDGNA = MDGNA(:,test); DGEA = DGEA(:,test); MDGEA = MDGEA(:,test);
            warning('redtextOBS: Missing observational total delays data for selected GNSS stations. Removing them from the processing')
        else
            test = contains(NAME,NAMES);
        end
    end
    clear opts Table
else
    disp('No matching files were found')
    ZTDA = []; MZTDA = []; DGNA = []; MDGNA = []; DGEA = []; MDGEA = [];
end


end