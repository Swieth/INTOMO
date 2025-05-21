function [T,T3D,undu] = distr_T_gpt2RT(BLh_pudel_num,dmj)
%function designed to generate the temperature over model domain
%___________________input_________________
% dmj - day of the year
% BLh_pudel - the centers of voxels coordinates
% z1 - number of UTM projection
% __________________output_________________
% T, T3D - the temperature distribuition inside model [K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[wi2,ki2,wa2] = size(BLh_pudel_num);
BLh_pudel_num_2D = reshape(BLh_pudel_num,[wi2 ki2*wa2]);
B = BLh_pudel_num_2D(1,:);
L = BLh_pudel_num_2D(2,:);
h = BLh_pudel_num_2D(3,:);
B = B';
L = L';
h = h';



% load('/home/estera/Testowy/model.mat');
% in=model.num_inner;
% B=B(in);
% L=L(in);
% h=h(in);
% fid=fopen('/home/estera/Testowy/WRFcoord_inner','w');
% for i = 1:size(B,1)
%     fprintf(fid,[sprintf('%03.8f', B(i,1)) '	' sprintf('%03.8f', L(i,1)) '	' sprintf('%2.4f', h(i,1))]);
%     fprintf(fid,'\n');
% end
% fclose(fid);

   B = B/180*pi();
   L = L/180*pi();
%for i = 1 : length(B)
   [p,t,dT,E,ah,aw,undu] = gpt2(dmj,B,L,h,size(B,1),0);
   clear p T dT ah aw
%end
clear p N 
T = t + 273.15;
T = T';
T3D = reshape(T,[1 ki2 wa2]);
