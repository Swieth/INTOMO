function q = e2q(e,T)
%Function designed to compute  water vapour absolut humidity from partial pressure and
%temperature
%_______________INPUT_______________
%%e - water vapour partial pressure at specified location in hPa
%T - temperature at specific location in K
%_______________OUTPUT______________
%q - absolute humidity in g / m-3 amount of water vapour in air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = 2.16679 * (e*100)./T;

 
 