function [A,R_SWD,elevat,stat_name,dSWD,SAT,SWD] = obsErrorFilter(A,R_SWD,elevat,stat_name,dSWD,SAT,SWD)
% Function to find outliers in observations based on observation error    
%---------------------------------------------------------------------------------------------------    
%%%INPUT/OUTPUT
%       A..........   the distances of each signal in tomography voxel
%       R_SWD......   the errors of the salnt delays
%       elevat.....   elevation angle to satellite
%       stat_name..   names of GNSS stations
%       dSWD.......   difference between observed and simulated slant delay
%       SAT........   id of GPS satellite
%       SWD........   the values of slant delays

    iddelt = [];
    % Find observations based on elevation level
    for i = 5:5:85
        idel1 = find(elevat>i);
        idel2 = find(elevat<i+5);
        idel = intersect(idel1,idel2);
        R_SWDel  = R_SWD(idel);
        iddelt1 = find(R_SWDel>(mean(R_SWDel)+ 4*std(R_SWDel))); 
        iddelt = [iddelt idel(iddelt1)'];
    end
    % Delete Nan values from dSWD
    if any(isnan(dSWD))
        id = isnan(dSWD);
        dSWD(id) = [];
    end
    % Delete outliers
    A(iddelt,:) = [];
    R_SWD(iddelt) = [];
    elevat(iddelt) = [];
    stat_name(iddelt) = [];
    dSWD(iddelt) = [];
    SAT(iddelt) = [];
    SWD(iddelt) = [];
    if ~isempty(iddelt)
        disp(['Observation outliers found, ',num2str(length(iddelt)),' observation removed from processing.'])
    end
end