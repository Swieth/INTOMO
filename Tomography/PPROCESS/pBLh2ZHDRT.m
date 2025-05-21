function ZHD = pBLh2ZHDRT(model,station)
% Function to generate Zenith Hydrostatic values based on apriori model.
%%%INPUT
%       model......   parameters of tomography and ray tracing models
%       station....   structural matrix containing all informations regarding GNSS observations
%%%OUTPUT
%       ZHD....... Zenith Hydrostatic Delay Values

RTstat = model.BLh;
RTName = model.NAME;
[~,ia,~] = intersect(RTName, station);
BLh = RTstat(ia,:);
RTstat = RTstat(ia,2:4);
set.refron = false;
[Xq,Yq,Zq] = meshgrid(model.lat_TOMO_RT,model.lon_TOMO_RT,model.levels_TOMO_RT);
p=[];
for epoch = 1:size(model.temp,1)
    for i = 1:size(RTstat,1)
        llRay = RTstat(i,:);
        p(epoch,i) = interp3(Xq,Yq,Zq,model.pres{epoch},llRay(1),llRay(2),llRay(3)./1000,'spline');
    end
end

ZHD = (0.0022767*p)./(1-0.00266*cos(2*((repmat(BLh(:,2),1,size(p,1))/180*pi()))')-0.00000028*(repmat(BLh(:,4),1,size(p,1)))');
