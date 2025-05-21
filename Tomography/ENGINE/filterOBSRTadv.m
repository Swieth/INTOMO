function  [A,SWD,R,elevation,SAT,station_name,id_weights] = filterOBSRTadv(A,SWD,R,elevation,SAT,station_name,xP,SWD_obs,switches)
%% This is experimental filtering function for both ground-based and RO observations based on their geometry.

%%%INPUT
%       A..........   .observational matrix with signal deriveratives
%       SWD........   the values of Slant Delay (SD)
%       R_SWD......   the errors of SD 
%       SAT........   id of the satellites
%       station_name  names of GNSS stations
%       Nw_apr.....   apriori values of refractivities
%       swtiches...   structural matrix of settings of the INTOMO processing
%%%OUTPUT
%       A..........   .observational matrix with signal deriveratives
%       SWD........   the values of Slant Delay (SD)
%       R_SWD......   the errors of SD 
%       SAT........   id of the satellites
%       station_name  names of GNSS stations
%       Nw_apr.....   apriori values of refractivities
%       id_weights.   weights of observations with unsuitable geometry

obs1 = size(SWD_obs,1);
[wi,~,~]=find(SWD(1:obs1)<0);
%% Remove RO crossing only top layers of tomography model
[wi,~,~]=find(SWD(1:obs1)<10000);
SWD(wi,:) = []; A(wi,:) = []; R(wi,:) = []; R(:,wi) = [];  elevation(wi,:) = []; SAT(wi,:) = []; SWD_obs(wi,:) = []; station_name(wi,:) = [];
if strcmp(switches.method,'LSQ') %only for space tomography
    wi = find(elevation~=0);
    SWD(wi,:) = []; A(wi,:) = []; R(wi,wi) = [];  elevation(wi,:) = []; SAT(wi,:) = []; SWD_obs(wi,:) = []; station_name(wi,:) = [];
end
obs = size(SWD,1);
close all
C = 1;
%% Remove outliers
for i = 1:1
    e = (SWD(1:obs)-A(1:obs,:)*xP);
    [wi,~,~] = find(abs(e-mean(e)) > 4* std(e));
    idd = find(elevation~=0);
    wi = unique([wi;idd]);
    SWD(wi,:) = []; A(wi,:) = []; R(wi,:) = []; R(:,wi) = []; elevation(wi,:) = []; SAT(wi,:) = []; SWD_obs(wi,:) = []; obs=obs-size(wi,1); station_name(wi,:)=[]; 
end
%% Remove pararell RO ray paths
A_c = A(1:length(elevation),:);
[numRows, ~] = size(A_c);
indeX = cell(numRows, 1);
nonzero = zeros(numRows, 1);
for i = 1:numRows
    indeX{i} = find(A_c(i, :) ~= 0);
    nonzero(i) = numel(indeX{i});
end 
forremoval = false(numRows, 1);
for i = 1:numRows-1
    if forremoval(i), continue; end
        similarrows = i;
        for j = i+1:numRows
            if forremoval(j), continue; end
            commonindices = intersect(indeX{i}, indeX{j});
            prob_fact = numel(commonindices) / min(nonzero(i), nonzero(j));
            if prob_fact > 0.9
                similarrows = [similarrows, j];
            end
        end
        if numel(similarrows) > 1
            [~, idxMax] = max(nonzero(similarrows));
            similarrows(idxMax) = [];
            forremoval(similarrows) = true;
        end
end
wi = 1:obs1;
wi = wi(forremoval);
SWD(wi,:) = []; A(wi,:) = []; R(wi,:) = []; R(:,wi) = [];  elevation(wi,:) = []; SAT(wi,:) = []; SWD_obs(wi,:) = []; station_name(wi,:) = [];
%% Change the weights of selected Ro ray paths based on their geometry
A_c = A(1:length(elevation),:);
[numRows, ~] = size(A_c);
indeX = cell(numRows, 1);
nonzero = zeros(numRows, 1);
for i = 1:numRows
    indeX{i} = find(A_c(i, :) ~= 0);
    nonzero(i) = numel(indeX{i});
end 
forremoval = false(numRows, 1);
for i = 1:numRows-1
    if forremoval(i), continue; end
        similarrows = i;
        for j = i+1:numRows
            if forremoval(j), continue; end
            commonindices = intersect(indeX{i}, indeX{j});
            prob_fact = numel(commonindices) / min(nonzero(i), nonzero(j));
            if prob_fact > 0.9
                similarrows = [similarrows, j];
            end
        end
        if numel(similarrows) > 1
            [~, idxMax] = max(nonzero(similarrows));
            similarrows(idxMax) = [];
            forremoval(similarrows) = true;
        end
end
wi = 1:obs1;
idweight = wi(forremoval);

