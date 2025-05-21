function [E,P] = distr_e_unb3RT( BLh_pudel_num,undu,jd)
%function designed to generate the water vapor pressure over model domain
%___________________input_________________
% dmj - day of the year
% BLh_pude_num - the centers of voxels coordinates
% z1 - number of UTM projection
% __________________output_________________
% E, E3D - the water vapour partial pressure distribuition inside model [mbar]
% P - pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[wi2,ki2,wa2] = size(BLh_pudel_num);
BLh_pudel_num_2D = reshape(BLh_pudel_num,[wi2 ki2*wa2]);
B = BLh_pudel_num_2D(1,:);
L = BLh_pudel_num_2D(2,:);
h = BLh_pudel_num_2D(3,:);
B = B';
L = L';
h = h';
H = h - undu;

B = B/180*pi();
L = L/180*pi();

doy=jd2doy(jd);

for i = 1:size(doy,1)
    for j = 1:size(B,1)
        [~,P,E(j,i),~]=UNB3MM(B(j),H(j),doy(i));
    end
end

E3D = reshape(E,[1 ki2 wa2]);

end

