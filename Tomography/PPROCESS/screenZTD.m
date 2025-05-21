function [ZTDA, MZTDA, DGNA, MDGNA, DGEA, MDGEA,model] = screenZTD(ZTDA, MZTDA, DGNA, MDGNA, DGEA, MDGEA,miss_stat,model)
%% Function to look for outliers and NAN values of MZTD
%   INPUT:
%       ZTDA......   matrix of zenith total delays
%       MZTDA.....   matrix of uncertainties of zenith total delays
%       DGNA......   matrix of N gradients of ZTD
%       MDGNA.....   matrix of uncertainties N gradient of ZTD
%       DGEA......   matrix of E gradients of ZTD
%       MDGEA.....   matrix of uncertainties E gradient of ZTD
%       NAMES.....   matrix of names of GNSS station
%       miss_stat.   id of stations for which observations are available
%       model......   parameters of tomography and ray tracing models
% OUTPUT:
%       ZTDA......   matrix of zenith total delays
%       MZTDA.....   matrix of uncertainties of zenith total delays
%       DGNA......   matrix of N gradients of ZTD
%       MDGNA.....   matrix of uncertainties N gradient of ZTD
%       DGEA......   matrix of E gradients of ZTD
%       MDGEA.....   matrix of uncertainties E gradient of ZTD
%       NAMES.....   matrix of names of GNSS station
%       model......   parameters of tomography and ray tracing models

nonan = find(isnan(MZTDA)==0);
meanMZTD = mean(MZTDA(nonan));
a = find(MZTDA>1.5*meanMZTD);
MZTDA(a) = NaN;
ZTDA(a) = NaN;
DGNA(a) = NaN;
DGEA(a) = NaN;
MDGEA(a) = NaN;
model.BLh = model.BLh(miss_stat,:);
model.BLH = model.BLH(miss_stat,:); 
model.BLh_ori = model.BLh_ori(:,miss_stat);
model.NAME = model.NAME(miss_stat,:);
model.X = model.X(miss_stat,:);
model.Y = model.Y(miss_stat,:);
model.Z = model.Z(miss_stat,:);
