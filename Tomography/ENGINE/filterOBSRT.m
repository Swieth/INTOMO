function  [A,SWD,R_SWD,elevation,SAT,station_name,id_del] = filterOBSRT(A,SWD,R_SWD,elevation,SAT,station_name,Nw_apr,SWD_obs,switches)
%Function to remove outlayers based on standard deviation of parameter e
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       A..........   .observational matrix with signal deriveratives
%       SWD........   the values of Slant Delay (SD)
%       R_SWD......   the errors of SD 
%       SAT........   id of the satellites
%       station_name  names of GNSS stations
%       Nw_apr.....   apriori values of refractivities
%       swtiches...   structural matrix of settings of the INTOMO processing
%%%OUTPUT
%       A..........   .observational matrix with signal deriveratives
%       SWD........   the values of Slant Delay (SD)
%       R_SWD......   the errors of SD 
%       SAT........   id of the satellites
%       station_name  names of GNSS stations
%       Nw_apr.....   apriori values of refractivities
%       id_del.....   id of observations deleted from the processing

obs1 = size(SWD_obs,1);
[wi,~,~]=find(SWD(1:obs1)<0);
SWD(wi,:) = []; A(wi,:) = []; R_SWD(wi,:) = [];  elevation(wi,:) = []; SAT(wi,:) = []; SWD_obs(wi,:) = []; station_name(wi,:) = [];
obs = size(SWD,1);
close all
C = 1;
id_del = [];
for i = 1:1
    e = (SWD(1:obs)-A(1:obs,:)*Nw_apr);
    [wi,~,~] = find(abs(e-mean(e)) > 3* std(e));
    id_del = [id_del;wi];
    SWD(wi,:) = []; A(wi,:) = []; R_SWD(wi,:) = [];R_SWD(:,wi) = []; elevation(wi,:) = []; SAT(wi,:) = []; SWD_obs(wi,:) = []; obs=obs-size(wi,1); station_name(wi,:)=[]; 
    if strcmp(switches.method,'LSQ')
        R_SWD = R_SWD./5; 
    end
end