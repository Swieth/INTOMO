function [PRN_num] = PRN_str2num(PRN_str)
%PRN_str2num converts string PRN to numeric PRN(non-vectorized version)
% 
%  [author]   : Tomasz Hadaœ, 03.04.2013
%  [argin]    : PRN_str - string containing PRN eg. R10, G07 etc.
%  [argout]   : PRN_num - numeric represantation of PRN_str eg. 60, 7
%
% [changelog] :

% constants
g_add = 0;
r_add = 50;
e_add = 100;

sat_num=str2double(PRN_str(2:3));
switch PRN_str(1)
    case 'G'
        PRN_num=g_add+sat_num;
    case 'R'
        PRN_num=r_add+sat_num;
    case 'E'
        PRN_num=e_add+sat_num;
    otherwise
end

end