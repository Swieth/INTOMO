function [ SWD_sum_all, A_rows, SWD2_p,Np,nodesC ] = A_row_bilin_regRT(model,p,Nw_apr,epoch,nodes_columns,levels,p_coordBLH,num_lev,Tmatrix,mult_fact,const,SWD_sum_all,A_rows,A_nodes,Np)
% Function to calculate the deriveratives of spline/bilinear functions
% based on tomgroaphy domain nodes
%%%INPUT
%    model......... parameters of tomography and ray tracing models
%    p.............id of subpoint
%    Nw_apr........apriori values of refractivities
%    epoch.........processing epoch number
%    nodes_columns. segragated id's of voxels for use in bilin interp
%    levels........altitude layers of tomography model
%    p_coordBLH.....coordinates of the subpoints in BLh
%    num_lev.......id numbers of voxels/nodes in inner model
%    Tmatrix.......second derivates of refractivity for use in bilin interp
%    mult_fact.....multiplication factors
%    const..........integral constant
%    SWD_sum_all... SWD integral based on a priori data
%    A_rows........id of ray path row in A matrix
%    A_nodes.......nodes surrounding the ray point
%    Np............apriori refractivity test value 
%%%OUTPUT
%    SWD_sum_all....SWD integral based on a priori data
%    A_rows.........id of ray path row in A matrix
%    SWD2_p.........subsection deriverative of refractivity
%    Np.............apriori refractivity test value 
%    nodesC.........nodes surrounding the subpoint
    %---------------------------------------------------------------------------------------------------
    X = model.Xbil; 
    Y = model.Ybil; 
    %% Find id of voxel nodes 
    AD = sqrt((p_coordBLH(p,1)*180/pi() - X').^2 + (p_coordBLH(p,2)*180/pi() - Y').^2);
    sdmin = sort(min(AD), 'ascend');
    for i = 1:2
        [~, col] = find(AD == sdmin(1,i));
        if numel(col) > 1
            col = col(1,1);
        end
        [~, row] = sort(AD(:,col), 'ascend');
        for j = 1:2
            node{i,j} = [row(j,1); col];
        end
    end
    id1 = find(model.bilmeshLatLon(1,:) == model.lat_TOMO(node{1,2}(1,1)));
    id2 = find(model.bilmeshLatLon(2,:) == model.lon_TOMO(node{1,2}(2,1)));
    x1 = intersect(id1,id2);
    id1 = find(model.bilmeshLatLon(1,:) == model.lat_TOMO(node{1,1}(1,1)));
    id2 = find(model.bilmeshLatLon(2,:) == model.lon_TOMO(node{1,1}(2,1)));
    x2 = intersect(id1,id2);
    id1 = find(model.bilmeshLatLon(1,:) == model.lat_TOMO(node{2,1}(1,1)));
    id2 = find(model.bilmeshLatLon(2,:) == model.lon_TOMO(node{2,1}(2,1)));
    x3 = intersect(id1,id2);
    id1 = find(model.bilmeshLatLon(1,:) == model.lat_TOMO(node{2,2}(1,1)));
    id2 = find(model.bilmeshLatLon(2,:) == model.lon_TOMO(node{2,2}(2,1)));
    x4 = intersect(id1,id2);
    dhz = find((p_coordBLH(p,3)*1000 - model.levels_TOMO)<0,1);
    %% Get coordinates of nodes
    if p_coordBLH(p,3)*1000 > model.levels_TOMO(end)
       dhz = size(model.levels_TOMO,2);
    end
    nodesA = [x1(dhz-1) x2(dhz-1) x3(dhz-1) x4(dhz-1) ];
    nodes = [nodesA nodesA + model.num_lat_TOMO*model.num_lon_TOMO];
    nodes = sort(nodes);
    nodesC(i,:) = nodes;
    nodes_coord = model.bilmeshLatLon([2 1 3],nodes)'; %First lambda then phi
    % Sort nodes according to coordinates (i, j, k)
    [nodes_coord,~]=sortrows(nodes_coord,[3,2,1]);
    % Get heights of the nodes
    nodes_bottom = sort([node{1,1}(2,1) + (size(model.lon_TOMO,2))*(node{1,1}(1,1)-1); node{2,1}(2,1) + (size(model.lon_TOMO,2))*(node{2,1}(1,1)-1);node{1,2}(2,1) + (size(model.lon_TOMO,2))*(node{1,2}(1,1)-1);node{2,2}(2,1) + (size(model.lon_TOMO,2))*(node{2,2}(1,1)-1)])';
    %[row,~]=find(nodes_columns==nodes_bottom(1,1));
    hk = levels(dhz-1,1); 
    hkp1 = levels(dhz,1); 
    dhk = hkp1-hk; 
    %% Calculate Nw apriori values and numbers of nodes used in interpolation / derivatives calculation
    nodes_num = [];
    for n = 1:size(nodes_bottom,2)
        [~,col]=find(nodes_columns==nodes_bottom(1,n));
        nodes_num = [nodes_num nodes_columns(:,col)'];
    end
    % Create matrix for vertical interpolation
    h = p_coordBLH(p,3)*1000; %h = 93;
    omegavert_ad=zeros(num_lev*2,1);
    omegavert_ad(dhz-1,1) = 1-((h-hk)/dhk);
    omegavert_ad(dhz,1) = (h-hk)/dhk;
    omegavert_ad(num_lev+dhz-1) = (h-hk)^2/2 - dhk*(h-hk)/3 - (h-hk)^3/(6*dhk);
    omegavert_ad(num_lev+dhz) = (h-hk)^3/(6*dhk) - (dhk*(h-hk))/6;
    omega_vert = zeros(4,num_lev*8);
    for num = 1:4
        omega_vert(num,num_lev*2*(num-1)+1:num_lev*2*num) = omegavert_ad';
    end
    % Horizonthal interpolation according to Perler
    nodes_coord_interp = [nodes_coord(1:4,1:3)];
    ncoorfud = flipud(nodes_coord_interp(:,1:2));
    % Find weights for nodes
    Bp = flip(p_coordBLH(p,1:2)*180/pi());
    for n = 1:size(ncoorfud,1)    
        for c = 1:2
            ommm(n,c) = abs(ncoorfud(n,c)-Bp(1,c))/abs(ncoorfud(n,c)-nodes_coord_interp(n,c));
        end
    end
    % Apply both interpolations together (vertical and horizontal)
    Np(p,1)=(prod(ommm,2)'*omega_vert*Tmatrix)*Nw_apr(epoch,nodes_num)'; % For diagnostics - Nw values based on a priori
    om_mult(p,:) = (prod(ommm,2)'*omega_vert*Tmatrix).*mult_fact(1,p).*const; % For A matrix construction
    % Calculate the SWD integral without A matrix, based on a priori data
    SWD_sum_all = SWD_sum_all + om_mult(p,:)*Nw_apr(epoch,nodes_num)';
    SWD2_p = om_mult(p,:)*Nw_apr(epoch,nodes_num)'; % For diagnostics - SWD value corresponding to this fragment of the signal
    %% Create a row that will be put into A matrix
    A_row=zeros(1,size(A_nodes,2));
    A_row(1,nodes_num)=om_mult(p,:);
    A_rows = [A_rows; A_row];
end

