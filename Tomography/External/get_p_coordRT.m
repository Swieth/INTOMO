function [ p_coord,p_coordBLH,const] = get_p_coordRT( par_st,sect,s_sect_nr,ss_sect_nr,azi,elev )
% Function to calculate the coordinates of inner points (subpoints) between simulated
% signal ray points.
%%%INPUT
%    p_coord........coordiantes of the subpoints in XYZ
%    p_coordBLH.....coordiantes of the subpoints in BLh
%    const..........integral constant
%%%OUTPUT
%    par_st.........coordiantes of the ray points
%    sect...........number of integral sections
%    s_sect_nr......number of integral subsections
%    ss_sect_nr.....number of integral subsubsections
%    azi............azimuth of the satellite
%    elev...........elevation angles of ray traced signals

radii = [6378.137, 6356.752314245]; 
wgs84 = wgs84Ellipsoid;
delta_s = vecnorm(par_st(sect,1:3)-par_st(sect+1,1:3));
ds=delta_s/(s_sect_nr*ss_sect_nr); %length of the ray in piece [subsubsection] (1/16 of section)
const=2*ds/45; %/1000 due to km not m
% Get coordinates of the points along the slant path
p_coord=NaN(s_sect_nr*ss_sect_nr+1,3);
p_coord(1,1)=par_st(sect,1);
p_coord(1,2)=par_st(sect,2);
p_coord(1,3)=par_st(sect,3);
% Loop for each piece [subsubsection] (coordinates)
if elev == 0 && azi == 0
    eX   = (par_st(sect,1:3)-par_st(sect+1,1:3))./(norm(par_st(sect,1:3)-par_st(sect+1,1:3))); 
    for p=2:(s_sect_nr*ss_sect_nr)
        p_coord(p,:) = p_coord(p-1,:) + eX* ds;
    end
else
    eX   = (par_st(sect+1,1:3)-par_st(sect,1:3))'./(norm(par_st(sect+1,1:3)-par_st(sect,1:3))); 
    for p=2:(s_sect_nr*ss_sect_nr)
        p_coord(p,:) = p_coord(p-1,:) + eX'* ds;
    end
end
p_coord(s_sect_nr*ss_sect_nr+1,1)=par_st(sect+1,1);
p_coord(s_sect_nr*ss_sect_nr+1,2)=par_st(sect+1,2);
p_coord(s_sect_nr*ss_sect_nr+1,3)=round(par_st(sect+1,3),10);
for i = 1:size(p_coord,1)
    [p_coordBLH(i,1),p_coordBLH(i,2),p_coordBLH(i,3)] = cspice_recgeo(p_coord(i,:)'  ,radii(1),(radii(1)-radii(2))/radii(1));
end
p_coord=p_coord(:,[2 1 3]); % First lambda then phi
p_coordBLH=p_coordBLH(:,[2 1 3]);

end

