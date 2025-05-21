function sh = ep2sh(e,p)
%%%function designed to comvert from water vapour partial pressure and
%%%total pressure to specific humidyt
%%%%%INPU%%%%%%%%%%%%%%%%%%%%
% e - water vapour partial pressire
% p - pressure
%%%%OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sh - specific humidity 

sh = 0.622 * e./ (p - 0.378 .* e);
