function   [A_apriori, SWD_apriori, P_apriori] = apConstRT(num_Nw,Nw_apr,varargin)
%Setup apriori inofrmation in form of a matrices
%%%INPUT
%   num_Nw....   number with number of voxels in particular time having observations
%   Nw_apr....   apriori values of refractivities
%%%OUTPUT
%   A_apriori.   the matrix containing 1 where there is aprioiri
%   observation or 0 where ther is no observation apriori
%   SWD_apriori  the matrix of all apriori values (Nw) wet refractivities
%   P_apriori..  the matrix containig weights for aprioir observations
if isempty(varargin)==0
    P_val = varargin{1};
else
    P_val = 0.1;
end
SWD_apriori = Nw_apr(~isnan(Nw_apr))';
A_apriori(1:size(SWD_apriori,1),1:size(Nw_apr,2)) = 0;
for zz = 1:size(num_Nw,1)
    A_apriori(zz,num_Nw(zz)) = 1;
    P_apriori(zz,1)= P_val;
end
if isempty(SWD_apriori)==1
    P_apriori = [];
end
clear zz

