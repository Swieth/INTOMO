function [SP3time,SP3X,SP3Y,SP3Z,SP3clk] = LG_interp(nXYZ,nCLK,igsTime,igsSP3X,igsSP3Y,igsSP3Z,igsSP3clk,delt)

% Lagrange interpolation satellite position
% 
% %%%%%%%%%% HELP %%%%%%%%%%
% 
% [SP3time,SP3X,SP3Y,SP3Z] = LG_interp(nXYZ,nCLK,igsTime,igsSP3X,igsSP3Y,igsSP3Z,igsSP3clk,delt)
% 
% Input Data
% nXYZ = number of points used for the Lagrange interpolation [1 x 1]
% igsTime = time corresponding to satellite position of IGS data [numep x 1], sec
% igsSP3X,igsSP3Y,igsSP3Z = satellite position of IGS data [numep x 1], metres
% igsSP3clk = satellite clock correction of IGS data [numep x 1], metres
% delt = time interval that wants to compute [1 x 1], second
% 
% Output Data
% SP3time = time corresponding to satellite position
% SP3X,SP3Y,SP3Z = satellite position after interpolate [numep x 1], metres
% 
% comment: igsTime is same GPSweek 
%
% Written by Phakphong Homniam
% April 22, 2003

%%%%%%%%%% BEGIN %%%%%%%%%%

delt;
NAJ = 999999.999999;

SP3time = [igsTime(1):delt:igsTime(end)];
numep = length(SP3time);

SP3X(1:numep) = NAJ;
SP3Y(1:numep) = NAJ;
SP3Z(1:numep) = NAJ;
SP3clk(1:numep) = NAJ;

for i = 1:numep
%     fprintf('epoch no: %6u\n',i);
    timeX = SP3time(i);
    [x1,x2,x3,x4] = XYZ(timeX,igsTime,igsSP3X,nXYZ,NAJ);
    SP3X(i) = LG(nXYZ,x1,x2,NAJ);
    [y1,y2,y3,y4] = XYZ(timeX,igsTime,igsSP3Y,nXYZ,NAJ);
    SP3Y(i) = LG(nXYZ,y1,y2,NAJ);
    [z1,z2,z3,z4] = XYZ(timeX,igsTime,igsSP3Z,nXYZ,NAJ);
    SP3Z(i) = LG(nXYZ,z1,z2,NAJ);  
    [c1,c2,c3,c4] = XCLK(timeX,igsTime,igsSP3clk,nCLK,NAJ);
    SP3clk(i) = LG(nCLK,c1,c2,NAJ);
end

%%%%%%% RESULT %%%%%%%

SP3time = SP3time(:);
SP3X = SP3X(:);
SP3Y = SP3Y(:);
SP3Z = SP3Z(:);
SP3clk = SP3clk(:);

%%%%%%%%%% END %%%%%%%%%%