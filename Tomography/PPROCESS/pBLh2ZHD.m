function ZHD = pBLh2ZHD(p,BLh)

ZHD = (0.002267*p)./(1-0.00266*cos(2*((repmat(BLh(:,2),1,size(p,1))/180*pi()))')-0.00000028*(repmat(BLh(:,4),1,size(p,1)))');
