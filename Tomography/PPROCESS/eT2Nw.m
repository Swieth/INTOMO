function  [Nw,Nw3D] = eT2Nw(e, T)
% The function has been developed to calculate the wet refractivity in a
% Voxel
%____________input________________
% e - the water vapour partial pressure (mbar)
% T - the temperature [K]
%__________output_________________
% Nw - the wet refractivity [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k2 = 72;
k3 = 370100;
[wi ko] = size(e);
clear ko
for i = 1 : wi
Zv = 1 + e(i,1)*(1 + 3.7*10 ^-4*e(i,1))*(-2.37321*10^-3 + 2.23366*T(i,1)^-1 - 710.792*T(i,1)^-2 + 7.76147*10^4*T(i,1)^-3);
Nw(i,1) = ( k2 * e(i,1) / T(i,1) + k3 * e(i,1) / T(i,1)^2 )*Zv ;
%Nw1(i,1) = ( k2 * e(i,1) / T(i,1) + k3 * e(i,1) / T(i,1)^2 ); 
clear Zv
Nw3D = [];
end
end