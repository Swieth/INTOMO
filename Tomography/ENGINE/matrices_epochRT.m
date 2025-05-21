function [A_k,SWD_k,R_SWD_k,elevat_k,SAT_k,stat_name_k] = matrices_epochRT(A,A_apriori,A_constr,SWD_apriori,SWD_constr,P_apriori,P_constr,SWD,R_SWD,elevation,SAT,station_name,epoch)
% Function to stack matrices
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       A..........   observational matrix with signal deriveratives
%       A_apriori..   the A apriori matrix with 0 or 1 values
%       A_constr...   constraints given to A matrix
%       SWD_apriori   the apriori values of Slant Delay (SD)
%       SWD_constr.   constraints given to SD 
%       P_apriori..   the apriori weights for the tomography grid
%       P_constr...   constraints given to P 
%       SWD........   the values of slant delays
%       elevat.....   elevation angle to satellite
%       SAT........   id of GPS satellite
%       station_name  names of GNSS stations
%%%INPUT
%       A_k........   stacked A matrices
%       SWD_k......   stacked SWD vector
%       R_SWD_k....   stacked R_SWD vector
%       elevat_k...   stacked elevation angle vector
%       SAT_k......   stacked id of GPS vector
%       stat_name_k   stacked names of GNSS vector


A_k = [A(epoch).A; A_apriori(epoch).A_apriori; A_constr(epoch).A_constr];
A_k = sparse(A_k);
SWD_k = [SWD(epoch).SWD; SWD_apriori(epoch).SWD_apriori; SWD_constr(epoch).SWD_constr];
SWD_k = sparse(double(SWD_k));
R_SWD_k = [R_SWD(epoch).R_SWD; P_apriori(epoch).P_apriori; P_constr(epoch).P_constr];
if any(isnan(R_SWD_k))
   if size(R_SWD_k,1) > size(A_k,1)
       id =  isnan(R_SWD_k);
       R_SWD_k(id) = [];
       R_SWD_k = R_SWD_k(1:size(A_k,1));
   else
       id =  isnan(R_SWD_k);
       R_SWD_k(id) = 0;
       R_SWD_k = R_SWD_k(1:size(A_k,1));
   end
   warning('Matrices Epoch: Unequal sizes of observation error matrix R and derivatatives matrix A_k. Cutting observations')
end
R_SWD_k = double(R_SWD_k);
R_SWD_k = sparse(R_SWD_k);
elevat_k = elevation(epoch).elevation;
SAT_k = SAT(epoch).SAT;
stat_name_k = station_name(epoch).station_name;