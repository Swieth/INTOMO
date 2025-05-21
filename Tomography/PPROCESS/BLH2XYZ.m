function [X,Y,Z] = BLH2XYZ (B,L,H)
%function to calculate X Y Z coordinates from B, L and H.
%__________INPUT___________
% B - latitude [rad] 
% L - longitude [rad]
% H - height [m]
%__________OUTPUT___________
% X - X coordinate in ECEF [m]
% Y - Y coordinate in ECEF [m]

% Z - Z coordinate in ECEF [m]
 a=6378137.000;
 b=6356752.314;
 eks=sqrt((a^2-b^2)/a^2);
for i=1:length(B)
    R(i,:)=a/sqrt(1-(eks^2)*sin(B(i,:))^2);
    X(i,:)=(R(i,:)+H(i,:))*cos(B(i,:))*cos(L(i,:));
    Y(i,:)=(R(i,:)+H(i,:))*cos(B(i,:))*sin(L(i,:));
    Z(i,:)=(R(i,:)*(1-eks^2)+H(i,:))*sin(B(i,:));
end
