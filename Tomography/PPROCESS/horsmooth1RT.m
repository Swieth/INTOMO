function [A_constr, SWD_constr, P_constr] = horsmooth1RT(not_cro,BLh_pudel_proj_num,num_lat, num_lon,constraints_r)
%_________________________________________
%function used to obtain horizonthal smoothing equations
%________________input__________________
%not_cro - set of not crossed voxel here set to all
%BLh_pudel_proj_num - set of projected coordinates of the central point voxel
%constraints_r - distance for constraints
%num_lat - number of voxels in north direction
%num_lon - number of voxels in est direction
%______________output_______________
%A_constr - the matrix containing weights for each voxel surrounding
%SWD_constr - the matrix of zeros
%P_constr - the matrix containig weights for constraints
%_____________________________________
if isempty(not_cro) == 1;
    W = [];
else
    [wi1,ki1] = size(not_cro);
    [wi2, ki2, wa2] = size (BLh_pudel_proj_num);
    BLh_pudel_proj_num_2D = reshape(BLh_pudel_proj_num,[wi2 ki2*wa2]);
    clear wi2 ki2 ki1 wa1
    for p = 1 : wi1
        dowag = not_cro(p,1);
        warstwa = floor((dowag-0.0001)/(num_lon*num_lat)) + 1;
        BLh_pudel_proj_num_1D = BLh_pudel_proj_num(:,:,warstwa);
        BLh_pudel_proj_num_1D = BLh_pudel_proj_num_1D';
        if dowag > num_lat*num_lon
            dowag = dowag - (warstwa-1)*num_lat*num_lon;
            if dowag == 0
                dowag = 1;
            end
        end
        num_pud_ide(:,1) = 1:num_lon*num_lat ;
        num_pud_ide(dowag,:) =[];
        point = BLh_pudel_proj_num_1D(dowag,:);
        BLh_pudel_proj_num_1D(dowag,:)=[];
        point = repmat(point,num_lon*num_lat-1,1);
        to_dist =[BLh_pudel_proj_num_1D point];
        to_dist(:,end+1) = ((to_dist(:,1) - to_dist(:,4)).^2 + (to_dist(:,2) - to_dist(:,5)).^2).^(0.5);
        [to_dist(:,end+1),wier] =sort(to_dist(:,7),1);
        clear wier
        % apply the change of correlation %
        
        if constraints_r >= 1 & constraints_r <= to_dist(1,7)
            wR = find(to_dist(:,7)>to_dist(1,7));
            to_dist(wR,7) = 0;
            clear wR
        elseif constraints_r > to_dist(1,7)
            wR = find(to_dist(:,7)>constraints_r);
            to_dist(wR,7) = 0;
            clear wR
            %        if warstwa < 7
            %        wR = find(to_dist(:,7)>60000);
            %        to_dist(wR,7) = 0;
            %        clear wR
        end
        Waga(dowag,:) = -1;
        
        to_dist(:,7) = to_dist(:,7)/sum(to_dist(:,7));
        Waga(num_pud_ide,:) =to_dist(:,7);
        zer = zeros((num_lat*num_lon),1);
        if warstwa > 1
            Przed = repmat(zer,1,warstwa - 1);
            Za = repmat(zer,1,wa2 - warstwa);
        else
            Przed = [];
            Za = repmat(zer,1,wa2 - warstwa);
        end
        Wexp = [Przed Waga Za];
        Wexp = reshape(Wexp,1,num_lon*num_lat*wa2);
        Wend(p,:)=Wexp;
        %Wend(not_cro(p,1),:)=Wexp;
        A_constr = Wend;
        SWD_constr(p,1) = 0;
        P_constr(p,1) = 1;
        clear dowag warstwa BLh_pudel_proj_num_1D num_pud_ide point to_dist Waga Wexp Za Przed zer
        
    end
    % [wi,ki] = size(W);
    %   if wi < ile *ile* wa2
    %       W(ile*ile* wa2,ile*ile * wa2) = 0;
    %   end
end