function [ p_coord,const] = get_p_coord( par_st,sect,s_sect_nr,ss_sect_nr,azi,elevation )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Length of the segment
delta_s = par_st(sect+1,4)-par_st(sect,4); %length of the ray in section
ds=delta_s/(s_sect_nr*ss_sect_nr); %length of the ray in piece [subsubsection] (1/16 of section)
const=2*ds/1000/45; %/1000 due to km not m
% Get coordinates of the points along the slant path
p_coord=NaN(s_sect_nr*ss_sect_nr+1,3);
p_coord(1,1)=par_st(sect,1);
p_coord(1,2)=par_st(sect,2);
p_coord(1,3)=par_st(sect,3);
%Mozna to przekszta³ciæ pozniej na liczenie na podstawieunit vector przy RO
% Loop for each piece [subsubsection] (coordinates)
for p=2:(s_sect_nr*ss_sect_nr)
    p_coord(p,1) = p_coord(p-1,1) + cos(azi*pi()/180)*ds*cos(elevation*pi()/180);
    p_coord(p,2) = p_coord(p-1,2) + sin(azi*pi()/180)*ds*cos(elevation*pi()/180);
    p_coord(p,3) = round(p_coord(p-1,3) + sin(elevation*pi()/180)*ds,10);
end
p_coord(s_sect_nr*ss_sect_nr+1,1)=par_st(sect+1,1);
p_coord(s_sect_nr*ss_sect_nr+1,2)=par_st(sect+1,2);
p_coord(s_sect_nr*ss_sect_nr+1,3)=round(par_st(sect+1,3),10);
p_coord=p_coord(:,[2 1 3]); % First lambda then phi
% [num_pud,Dsr, ~, node_dist, node_nr, mean_voxels] = set_num_pud_nodes(BLh_pudel_proj_num,BLh_pudel_proj,p_coord,switches);

end

