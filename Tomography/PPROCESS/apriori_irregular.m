function [ Np ] = apriori_irregular( mean_voxels, BLh2D_pudel_proj,Nw_apr,czas,nodes_tomo,planes,nodes_columns,Nwapr_columns,planes_nr,Tmatrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for point_nr = 1:size(mean_voxels,1)
    nodes=mean_voxels(point_nr,:);
    nodes_coord = BLh2D_pudel_proj([2 1 3],nodes)'; %First lambda then phi
    nodes_val = Nw_apr(czas,nodes)';
    p_coord=[nodes_tomo(2,point_nr)' nodes_tomo(1,point_nr)' nodes_tomo(3,point_nr)']; %First lambda then phi
    
    %Sort nodes according to coordinates (i, j, k)
    [nodes_coord,n_order]=sortrows(nodes_coord,[3,2,1]);
    nodes_val=nodes_val(n_order,:);
    nodes_sorted = nodes(n_order);
    
    % Get heights of nodes
    nodes_bottom = nodes_sorted(:,1:4);
    [row,~]=find(nodes_columns==nodes_bottom(1,1));
    hk = planes(row,1); % height k
    hkp1 = planes(row+1,1); % height k+1
    dhk = hkp1-hk; % height difference dhk
    
    % Get Nw apriori values and numbers of nodes used in interpolation / derivatives calculation
    Nij=[];
    nodes_num = [];
    for n = 1:size(nodes_bottom,2)
        [row,col]=find(nodes_columns==nodes_bottom(1,n));
        Nij = [Nij; Nwapr_columns(:,col)];
        nodes_num = [nodes_num nodes_columns(:,col)'];
    end
    nodes_numbers = nodes_num;
    
    % Create matrix for vertical interpolation
    p=1;
    h = p_coord(p,3);
    omegavert_ad=zeros(planes_nr*2,1);
    omegavert_ad(row,1) = 1-((h-hk)/dhk);
    omegavert_ad(row+1,1) = (h-hk)/dhk;
    omegavert_ad(planes_nr+row) = (h-hk)^2/2 - dhk*(h-hk)/3 - (h-hk)^3/(6*dhk);
    omegavert_ad(planes_nr+row+1) = (h-hk)^3/(6*dhk) - (dhk*(h-hk))/6;
    omega_vert = zeros(4,planes_nr*8);
    for num = 1:4
        omega_vert(num,planes_nr*2*(num-1)+1:planes_nr*2*num) = omegavert_ad';
    end
    
    % Create matrix for horizontal interpolation
    nodes_coord_interp = [nodes_coord(1:4,1:2) repmat(h,4,1)];
    nodes_coord_interp2 = repmat(nodes_coord_interp,2,1); %for comfortable extraction of coordinates for the last nodes
    for n = 1:size(nodes_bottom,2)
        omega(n,1)=abs(nodes_coord_interp2(n+1,1)-p_coord(p,1))/abs(nodes_coord_interp2(n+1,1)-nodes_coord_interp2(n,1));
        omega(n,2)=abs(nodes_coord_interp2(n+2,2)-p_coord(p,2))/abs(nodes_coord_interp2(n+2,2)-nodes_coord_interp2(n,2));
    end
    % Apply both interpolations together (vertical and horizontal)
    Np_temp = omega'*omega_vert*Tmatrix*Nij;
    Np(point_nr,1)=sum(Np_temp)/4;
    %consider: Np(point_nr,1) = prod(omega,2)'*omega_vert*Tmatrix*Nij
end

%             % Plot apriori - scatter with colors
%             % regular
%             figure(1);
%             x = BLh2D_pudel_proj(1,:);
%             y = BLh2D_pudel_proj(2,:);
%             z = BLh2D_pudel_proj(3,:);
%             v = Nw_apr(1,:);
%
%             colormap(jet)
%             scatter3(x,y,z,[],v,'filled');
%             colorbar
%             caxis([0, 60])
%
%             % irregular
%             figure(2)
%             x = nodes_tomo(1,:);
%             y = nodes_tomo(2,:);
%             z = nodes_tomo(3,:);
%             v = Np(:,1);
%
%             colormap(jet)
%             scatter3(x,y,z,[],v,'filled');
%             colorbar
%             caxis([0, 60])
end

