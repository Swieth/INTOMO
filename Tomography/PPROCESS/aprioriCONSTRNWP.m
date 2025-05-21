function [num_Nw, Nw_obs] = aprioriCONSTRNWP(start_time_n,num_Nw_n,time_c,num_Nw_c,time_stA,num_Nw_stA,NwA_NWP)
%function build to derive information about the aprioiri values
%%%INPUT
%       start_time_n.... all time epochs when radiosonding is available
%       num_Nw_n........ all numbers of radiosounding crossed voxels
%       time_c.......... all time epochs when COSMIC is availible
%       num_Nw_c........ all numbers of COSMIC crossed voxels
%       time_stA........ all times epochs when Ground data are avlilabe 
%       num_Nw_stA...... all numbers of Groud crossed voxles
%       NwA_NWP......... TOMO model with reference refractivity from NWP\
%%%OUTPUT
%       num_Nw..........the structural matrix with num_Nw(i).numbers where i is time and
%       numbers is equal to voxel numbers
%       Nw_obs.......... set of observations from NwA_NWP in the spaces of all availible
%       data.

RADIO_ap = [start_time_n'; num_Nw_n];
COSMIC_ap = [time_c'; num_Nw_c];
GROUND_ap = [time_stA'; num_Nw_stA];

[wie,kie] = find(sum(GROUND_ap')==0);
GROUND_ap(kie,:) = [];
clear wie kie
[wie,kol] = size(NwA_NWP);
Nw_obs(1:wie,1:kol) = NaN;

for i = 1:size(NwA_NWP,1)
    num_Nw(i).number= [];
end
RADIO_ap = unique(RADIO_ap','rows')';
COSMIC_ap = unique(COSMIC_ap','rows')';
GROUND_ap = unique(GROUND_ap','rows')';

for i = 1 : size(RADIO_ap,2)
     if isempty(num_Nw(RADIO_ap(1,i)).number) ==1
        num_Nw(RADIO_ap(1,i)).number =  RADIO_ap(2:end,i);
     else
        num_Nw(RADIO_ap(1,i)).number = union(num_Nw(RADIO_ap(1,i)).number,RADIO_ap(2:end,i));
     end
     Nw_obs(RADIO_ap(1,i), RADIO_ap(2:end,i)) = NwA_NWP(RADIO_ap(1,i), RADIO_ap(2:end,i));
end

for i = 1 : size(COSMIC_ap,2)
   if isempty(num_Nw(COSMIC_ap(1,i)).number) ==1
      num_Nw(COSMIC_ap(1,i)).number =  COSMIC_ap(2:end,i);
   else
      num_Nw(COSMIC_ap(1,i)).number = union(num_Nw(COSMIC_ap(1,i)).number,COSMIC_ap(2:end,i));
   end
      Nw_obs(COSMIC_ap(1,i), COSMIC_ap(2:end,i)) = NwA_NWP(COSMIC_ap(1,i), COSMIC_ap(2:end,i));
end

for i = 1 : size(GROUND_ap,2)
    
    if  GROUND_ap(1,i)<=size(NwA_NWP,1);
    to_unite = GROUND_ap(2:end,i);
    [wie,kol,tt] = find(to_unite>0);
    to_unite = unique(to_unite(wie,1));
    clear wie kol tt
    if isempty(num_Nw(GROUND_ap(1,i)).number) ==1
      %num_Nw(GROUND_ap(1,i)).number =  GROUND_ap(2:end,i);
      num_Nw(GROUND_ap(1,i)).number =  to_unite;
    else
      
      %num_Nw(GROUND_ap(1,i)).number = union(num_Nw(GROUND_ap(1,i)).number,GROUND_ap(GROUND_ap(2:end,i)~=0),i));
      num_Nw(GROUND_ap(1,i)).number = union(num_Nw(GROUND_ap(1,i)).number,to_unite);
    end
      NwA_PACK = NwA_NWP(GROUND_ap(1,i),to_unite);
      [wart,wie,kol] = find(abs(NwA_PACK)>300);
      NwA_PACK(:,wie) = [];
      to_unite(wie,:) = [];
      num_Nw(GROUND_ap(1,i)).number(wie,:) = [];
       Nw_obs(GROUND_ap(1,i), to_unite) = NwA_PACK;
      clear to_unite wart wie kol
    end
end


