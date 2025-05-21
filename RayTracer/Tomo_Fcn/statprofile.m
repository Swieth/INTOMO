function [ZHD, ZWD, Nh, Nw, nodev, nodeh, presh,eh,temph] = statprofile(XLAT, XLONG, pGrid, tempGrid, eGrid, hlayers, llRay, rWGS)

% statprofile calculates zenith delays: hydrostatic and wet  using
% impertoladed grid of meteorological parameters stored in ERA5 
% accurate solution

%% Input data:
%       
%               XLAT:               latitude gird [rad]
%               XLONG:              longitude grid [rad]
%               pGrid:              pressure grid [hPa]
%               tempGrid:           temperature grid [K]
%               eGrid:              water vapour partial pressure grid
%               [hPa]
%               sLWCGrid:           specific liquid water cloud content
%               [kg/kg]
%               sIWCGrid:           specific ice water cloud content
%               [kg/kg]
%               sRWCGrid:           specific rain water content
%               [kg/kg]
%               sLWCGrid:           specific snow water content
%               [kg/kg]
%               hlayer:             grid heights
%               llRay:              station coordinates
%               rWGS:               WGS radii
%
%
%% Output data: 
%               ZHD:                zenith hydrostatic delay          [m]
%               ZWD:                zenith wet delay                  [m]
%               Nh, Nw:             refractivity indicies          [unit]
%               nodev               vertical nodes
%               nodeh               horizontal nodes
%               presh               atmopshperic preassure at point [hPa]
%               eh                  water vapour preassure at point [kg/kg-1]
%               temph               temperature at point            [K]


%% DEFINE DATA LIMITS FOR GIVEN RAYPATH
wrfbt = length(hlayers);

% Angular distance (precise)
AD = sqrt((llRay(1,1) - XLAT).^2 + (llRay(1,2) - XLONG).^2);

% Find one lowest value in each column and sort array in ascending order
sdmin = sort(min(AD), 'ascend');
for i = 1:2
    % Find column number for 2 lowest values
    [~, col] = find(AD == sdmin(1,i));
    % Error found: can result in 2 equal values but only 1 can be used
    if numel(col) > 1
        col = col(1,1);
    end
    % Give value and row for found column
    [val, row] = sort(AD(:,col), 'ascend');
    for j = 1:2
        % Create cell with location and value for wrf grid
        node{i,j} = [row(j,1); col; val(j,1)];
    end
end

% Array with 4 nearest nodes: wrf row (1), wrf column (2), distance to ray (3)
nodeh = cell2mat(reshape(node,1,4));

% Weights: inverse squared distance
wNode = nodeh(3,:).^(-2);
[ ePrf, tempPrf, pPrf] = deal(zeros(wrfbt, 4));
Rwgs = zeros(1,4);
%% METEOROLOGICAL PARAMETERS PROFILES
for i=1:4   
    % pressure
    pPrf(:,i) = squeeze(pGrid(nodeh(1,i), nodeh(2,i), 1:wrfbt));

    % water vapor pressure
    ePrf(:,i) = squeeze(eGrid(nodeh(1,i), nodeh(2,i), 1:wrfbt));

    % temperature
    tempPrf(:,i) = squeeze(tempGrid(nodeh(1,i), nodeh(2,i), 1:wrfbt));

    % WGS radius
    Rwgs(:,i) =  rWGS(nodeh(1,i), nodeh(2,i));
end

% find vertical layers above station
geomRay = llRay(1,3);
hpre = hlayers(hlayers>geomRay);

dhSurf = deal(zeros(length(hpre)+1,1));

% interpolate refractivity at station position
% Find layer above and below current ray altitude
dgeom = geomRay - hlayers;
% Layer below
valb = find(dgeom >= 0);
% Layer above
vala = find(dgeom < 0 );
     
nodev = [max(valb); min(vala)];

% Increment between ray altitude and first vertical model layer (ground)
dhSurf(1,:) = geomRay;
dhSurf(2:end,:) = hpre';

%% ERA METEOROLOGICAL PARAMETERS FOR NODES
[ presi, tempi, ei] = deal(zeros(size(nodev,1),4));
% Loop through nodes to read data for 2 nearest levels and 4 nodes
for k = 1:4
    % DIM: west - east (x)/ south - north (y)/ bottom - top / time
    for i = 1:size(nodev,1)
        % pressure [hPa]
        presi(i,k) = pPrf(nodev(i),k);
        
        % temperature
        tempi(i,k) = tempPrf(nodev(i),k);
        
        % water vapor pressure
        ei(i,k) = ePrf(nodev(i),k);
    end
end    

[eh, temph, presh] = deal(zeros(length(dhSurf),1));
%% interpolate parameters at station height
% water vapor pressure
Ce = (hlayers(nodev(2)) - hlayers(nodev(1)))*1000./log(ei(2,:)./ei(1,:));
% Can results in NaN
eh1 = ei(1,:).*exp((geomRay - hlayers(nodev(1)))*1000./Ce);
eh2 = ei(2,:).*exp((geomRay - hlayers(nodev(2)))*1000./Ce);
eh(1,:) = sum((eh1 + eh2)./2.*wNode)./sum(wNode);

% temperature
Ctemp = (tempi(2,:) - tempi(1,:))./(hlayers(nodev(2)) - hlayers(nodev(1)));
% Linear interpolation of temperature
temph1 = tempi(1,:) + (geomRay - hlayers(nodev(1))).*Ctemp;
temph2 = tempi(2,:) + (geomRay - hlayers(nodev(2))).*Ctemp;
temph(1,:) = sum((temph1 + temph2)./2.*wNode)./sum(wNode);

presh(1,:) = sum(mean(presi).*wNode)./sum(wNode); 

% calculate parameters for rest layers
presh(2:end,:) = sum(pPrf(vala',:).*repmat(wNode,length(vala'),1),2)./sum(wNode);
eh(2:end,:) = sum(ePrf(vala',:).*repmat(wNode,length(vala'),1),2)./sum(wNode);
temph(2:end,:) = sum(tempPrf(vala',:).*repmat(wNode,length(vala'),1),2)./sum(wNode);

% calculate refractivites
[ Nh, Nw] = refcalc(presh, temph, eh, 'c');

% total refractivity index
Nt = Nh + Nw; %+ Nhm;
% refractivities
% wet
nw = Nw.*10^-6+1;
% hydrostatic
nh = Nh.*10^-6+1;
% total
nt = Nt.*10^-6+1;

% Mean refractive indeces
nhm = (nh(1:end-1,1) + nh(2:end,1))/2; 
nwm = (nw(1:end-1,1) + nw(2:end,1))/2;  
ntm = (nt(1:end-1,1) + nt(2:end,1))/2; 

% Incremental heights for profiles integration
dhinc = dhSurf(2:end)-dhSurf(1:end-1);

% Partial zenith delays
dZHD = (nhm - 1).*1e6.*dhinc;
dZWD = (nwm - 1).*1e6.*dhinc;
dZTD = (ntm - 1).*1e6.*dhinc;

ZHD = sum(dZHD)/1e6*1000;
ZWD = sum(dZWD)/1e6*1000;
ZTD = sum(dZTD)/1e6*1000;


end