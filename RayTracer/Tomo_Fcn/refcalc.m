
function [Nh, Nw] = refcalc(pres, temp, e, solution)
% refcalc calculates refractivity in a
% function of pressure, temperature and water vapor pressure 
%% Input data: 
% 	pres.........  pressure                              [hPa]
% 	temp.........  temperature                           [K]
%   e............  water vapor pressure                  [hPa]
%   expressed as (mass of condensate / mass of moist air)
%   solution:       a, b, c                               [string]
%% Output data: 
%   Nh..........   hydrostatic refractivity    [-]
%   Nw..........   wet refractivity            [-]


%% CONSTANTS
Rd = 287.058;       % gas constant dry air [J/K/kg]
Rw = 461.525;       % gas constant wet air [J/K/kg]
% Bevis (1994) 
k1 = 77.689;
k2 = 71.295;
k3  = 375463;  

k2p = k2 - k1*Rd/Rw; 


%% CALCULATE REFRACTIVITY 

% Results can affect ZHD in +/- 1cm
if solution == 'a'
    % - use total pressure and k2' coefficient: biggest Nh
    Nh = k1*(pres)./temp;
    Nw = k2p*e./temp +  k3*e./temp.^2;

elseif solution == 'b' 
    % Dry and vapor part (Smith and Weintraub, 1953; Thayer et al., 1974)
    % - use dry pressure and k2 coefficient: smallest Nh, biggest Nw
    Nh = k1*(pres - e)./temp;
    Nw = k2*e./temp +  k3*e./temp.^2;

elseif solution == 'c' 
    % Hydrostatic and nonhydrostatic part (Davis et al., 1985)
    % Air densities
    Dd = (pres - e)./Rd./temp;
    Dw = e./Rw./temp;
    D = Dd + Dw;

    % moderate Nh
    Nh = k1*Rd.*D; %Rd = Ru/Md
    Nw = k2p*e./temp + k3*e./temp.^2;
end

end