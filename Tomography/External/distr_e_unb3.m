function [E,E3D] = distr_e_unb3( BLh_pudel_num,undu,jd,z1 )
%function designed to generate the water vapor pressure over model domain
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
h = BLh_pudel_num_2D(3,:)*1000;
B = B';
L = L';
h = h';
H = h - undu';
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

doy=jd2doy(jd);

for i = 1:size(doy,1)
    for j = 1:size(B,1)
        [~,~,E(j,i),~]=UNB3MM(B(j),H(j),doy(i));
    end
end

E3D = reshape(E,[1 ki2 wa2]);

end

