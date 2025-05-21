function BLHstruc = pudel2(model,latT,lonT,levT,switches)
% Function to generate ray tracing and tomography model domain coordinates
% in 2D
%%% INPUT
%           levT......  altitudes of tomography model voxel centers
%           latT......  latitudes of tomography model voxel centers
%           lonT......  longitudes of tomography model voxel centers
%%% OUTPUT
%          BLHstruc...   2D coordinates of tomography model nodes

B_T = latT;
L_T = lonT;
h_T_o = levT;
B_T_o = repmat(B_T,1,size(lonT,2));
L_T_o = repmat(L_T,size(latT,2),1);
L_T_o = reshape(L_T_o,(size(latT,2))*(size(lonT,2)),1);
L_T_o = L_T_o';
B_T_o = repmat(B_T_o,[1 1 size(levT,2)]);
L_T_o = repmat(L_T_o,[1 1 size(levT,2)]);
h_T_o= repmat(h_T_o,size(B_T_o,2),1);   
for w = 1:size(levT,2)
    h_T_new(1,:,w) = h_T_o(:,w);        
end
clear h_T_o, 
h_T_o = h_T_new;
clear h_T_new
BLh_pudel_rad = [B_T_o; L_T_o; h_T_o];
clearvars -except BLh_pudel_rad model latT lonT levT switches
B_T = latT(1:end-1) + diff(latT)/2;
L_T = lonT(1:end-1) + diff(lonT)/2;
h_T_o = levT(1:end-1) + diff(levT)/2;
B_T_o = repmat(B_T,1,size(L_T,2));
L_T_o = repmat(L_T,size(B_T,2),1);
L_T_o = reshape(L_T_o,(size(B_T,2))*(size(L_T,2)),1);
L_T_o = L_T_o';
B_T_o = repmat(B_T_o,[1 1 size(h_T_o,2)]);
L_T_o = repmat(L_T_o,[1 1 size(h_T_o,2)]);
h_T_o = repmat(h_T_o,size(B_T_o,2),1);   

for w = 1:size(h_T_o,2)
    h_T_new(1,:,w) = h_T_o(:,w);        
end
clear h_T_o, 
h_T_o = h_T_new;
clear h_T_new
BLh_pudel_rad_num = [B_T_o; L_T_o; h_T_o];
   
for i = 1:size(BLh_pudel_rad_num,2)
    X_T_num(:,i) = cspice_georec(BLh_pudel_rad_num(2,i,1)*pi()/180, BLh_pudel_rad_num(1,i,1)*pi()/180, BLh_pudel_rad_num(3,i,1)/1000, model.radii(1),(model.radii(1)-model.radii(2))/model.radii(1));
end
for i = 1:size(BLh_pudel_rad,2)
    X_T(:,i) = cspice_georec(BLh_pudel_rad(2,i,1)*pi()/180, BLh_pudel_rad(1,i,1)*pi()/180, BLh_pudel_rad(3,i,1)/1000, model.radii(1),(model.radii(1)-model.radii(2))/model.radii(1));
end
BLh_pudel_num = repmat(X_T_num,1,1,size(BLh_pudel_rad_num,3));
BLh_pudel = repmat(X_T,1,1,size(BLh_pudel_rad,3));

if strcmp(switches.parametrization,'constant')
    removal_inner = 1:size(BLh_pudel_num,2);
    num_outer = 1:size(BLh_pudel_num,2)*size(BLh_pudel_num,3);
    removal_inner = reshape(removal_inner,size(latT,2)-1,size(lonT,2)-1);
    removal_inner = flipud(removal_inner);
    removal_inner(:,1) = [];
    removal_inner(:,size(removal_inner,2)) = [];
    removal_inner(1,:) = [];
    removal_inner(size(removal_inner,1),:) = [];
    removal_inner = flipud(removal_inner);
    removal_inner = reshape(removal_inner,(size(latT,2)-3)*(size(lonT,2)-3),1);
    removal_inner = removal_inner';
    num_inner = removal_inner;
    add_inner = 0:(size(BLh_pudel_num,3)-1);
    add_inner = add_inner';
    add_inner = repmat(add_inner,1,size(num_inner,2));
    num_inner = repmat(num_inner, size(BLh_pudel_num,3),1);
    num_inner = num_inner + add_inner*(size(latT,2)-1)*(size(lonT,2)-1);
    num_inner = reshape(num_inner',(size(latT,2)-3)*(size(lonT,2)-3)*size(BLh_pudel_num,3),1);
    num_outer(:,num_inner) =[];
else
    removal_inner = 1:size(BLh_pudel,2);
    num_outer = 1:size(BLh_pudel,2)*size(BLh_pudel,3);
    removal_inner = reshape(removal_inner,size(latT,2),size(lonT,2));
    removal_inner = flipud(removal_inner);
    removal_inner(:,1) = [];
    removal_inner(:,size(removal_inner,2)) = [];
    removal_inner(1,:) = [];
    removal_inner(size(removal_inner,1),:) = [];
    removal_inner = flipud(removal_inner);
    removal_inner = reshape(removal_inner,(size(latT,2)-2)*(size(lonT,2)-2),1);
    removal_inner = removal_inner';
    num_inner = removal_inner;
    add_inner = 0:(size(BLh_pudel,3)-1);
    add_inner = add_inner';
    add_inner = repmat(add_inner,1,size(num_inner,2));
    num_inner = repmat(num_inner, size(BLh_pudel,3),1);
    num_inner = num_inner + add_inner*(size(latT,2))*(size(lonT,2));
    num_inner = reshape(num_inner',(size(latT,2)-2)*(size(lonT,2)-2)*size(BLh_pudel,3),1);
    num_outer(:,num_inner) =[];
end

for i = 1:size(model.BLh,1)
    Xsta(:,i) = cspice_georec(model.BLh(i,3)*pi()/180, model.BLh(i,2)*pi()/180, model.BLh(i,4)/1000, model.radii(1),(model.radii(1)-model.radii(2))/model.radii(1));
end


BLHstruc.BLh_pudel = BLh_pudel;
BLHstruc.BLh_pudel_num = BLh_pudel_num;
BLHstruc.BLh_pudel_rad = BLh_pudel_rad;
BLHstruc.BLh_pudel_num_rad = BLh_pudel_rad_num;
BLHstruc.BLh_outer = [];
BLHstruc.BLh_outer_num = [];
BLHstruc.BLh_outer_rad = [];
BLHstruc.BLh_outer_num_rad = [];
BLHstruc.Bsta = Xsta(1,:);
BLHstruc.Lsta = Xsta(2,:);
BLHstruc.num_outer = num_outer;
BLHstruc.num_inner = num_inner;    
    
    
    

end