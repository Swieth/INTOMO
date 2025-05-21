function Pzero = PzerosetupRT(A,BLh_pudel_proj_num,num_lat,num_lon)
%function to set initial estimation error
%The assigned error is related to the findings published in Advances in Space Research and
%Atmospheric Research

H = BLh_pudel_proj_num(3,1,:);
H = H(:);
num_lev = size(H,1);
[num,wi,wa] = find(H>7000);
%change due to research done by E. Trzcina 23.05.2016
%mNw(1:num(1)-1,1) = 1.5*5;
%mNw(num,1) = 3*5;
% mNw(1:num(1)-1,1) =4*5;
% mNw(num,1) = 1*5;
mNw(1:num(1)-1,1) =1;
mNw(num,1) = 1;
mNw = repmat(mNw,1,num_lat*num_lon);
mNw = reshape(mNw',num_lat*num_lon*(num_lev),1);
if size(mNw,1)<size(A,2)
    mNw(size(mNw,1)+1:size(A,2),1)=3*5;
end
Pzero = diag(mNw.^2);