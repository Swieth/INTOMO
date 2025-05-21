function [E,E3D] = distr_e_gpt2(BLh_pudel_num,dmj,z1)
%function designed to generate the temperature over model domain
%___________________input_________________
% dmj - day of the year
% BLh_pude_num - the centers of voxels coordinates
% z1 - number of UTM projection
% __________________output_________________
% E, E3D - the water vapour partial pressure distribuition inside model [mbar]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wi2,ki2,wa2] = size(BLh_pudel_num);
BLh_pudel_num_2D = reshape(BLh_pudel_num,[wi2 ki2*wa2]);
B = BLh_pudel_num_2D(1,:);
L = BLh_pudel_num_2D(2,:);
h = BLh_pudel_num_2D(3,:);
B = B';
L = L';
h = h';
if isstruct(z1)==0 & ischar(z1)==0 & license('test','map_toolbox')==1
  [ellipsoid,estr] = utmgeoid(z1);
   utmstruct = defaultm('utm'); 
   utmstruct.zone = z1; 
   utmstruct.geoid = ellipsoid; 
   utmstruct = defaultm(utmstruct); 
  [B,L] = minvtran(utmstruct,L,B);
elseif isstruct(z1)== 1 & license('test','map_toolbox')==1
    utmstruct = z1;
   [B,L] = minvtran(utmstruct,L,B);

elseif license('test','map_toolbox')==0
    if strcmpi(z1,'utm') == 1
    PRO =  tm2ell([L B],'utm');
      L = PRO(:,1);
      B = PRO(:,2);
      z1 = 'utm';
    else
      PRO =  tm2ell([L B],'pl92');
      L = PRO(:,1);
      B = PRO(:,2);
      z1 = 'pl92';
       
    end
end
   B = B/180*pi();
   L = L/180*pi();
%for i = 1 : length(B)
   [p,T,dT,E,ah,aw,undu] = gpt2(dmj,B,L,h,size(B,1),1);
   clear p T dT ah aw undu
%end
clear T P TM 
E = E';
E3D = reshape(E,[1 ki2 wa2]);
