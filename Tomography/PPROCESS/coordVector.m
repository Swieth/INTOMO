function coordV = coordVector(coordV,i_pos)
% Function to find local direction vectors of simulated ray path (in each traversed voxel).
%%%INPUT
%      coordV.... coordinates of the vectors is XYZ directions
%%%OUTPUT
%      coordV.... coordinates of the vectors is XYZ directions
%      i_pos..... id number of ray points inside tomography model domain

    id_coord = 1:size(coordV,1);
    coordV = squeeze(coordV(:,:,1));
    id = find(isnan(coordV(:,1)));
    id_coord(id) = [];
    for i = 1:size(id_coord,2)-1
        len_coord = length(id_coord(i)+1:id_coord(i+1));
        for k = id_coord(i)+1:id_coord(i+1)-1
            coordV(k,:) = coordV(k-1,:) + (coordV(id_coord(i+1),:) - coordV(id_coord(i),:))./len_coord;
        end
    end
    coordV = diff(coordV(id_coord,:));
    coordV(size(coordV,1)+1,:) = [NaN,NaN,NaN];
    id = find(isnan(coordV(:,1)));
    coordV(id,:) = repmat(coordV(min(id)-1,:),length(id),1); 
end