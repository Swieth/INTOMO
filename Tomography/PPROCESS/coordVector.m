function dirVec = coordVector(coordIn, ~)
% COORDVECTOR  Compute local direction vectors of a simulated ray path
%              (one vector per traversed voxel).
%
% INPUT
%   coordIn - [nPts x 3 x 2] XYZ coordinates of voxel entry/exit points
%             (3-D array as produced by the ray tracer; NaN rows mark gaps)
%   ~       - (unused) i_pos voxel index array passed by callers;
%             retained in signature for backward compatibility with
%             groundRT.m and spaceRT.m
%
% OUTPUT
%   dirVec  - [nVoxels x 3] direction vectors for each traversed voxel

    % Flatten to a 2-D [nPts x 3] array using the first slice
    coordIn = squeeze(coordIn(:,:,1));

    % Identify rows that carry valid (non-NaN) coordinates
    id_coord = find(~isnan(coordIn(:,1)))';

    % Compute difference vectors between successive valid boundary points
    dirVec = diff(coordIn(id_coord,:));   % [(numel(id_coord)-1) x 3]

    % Append a NaN row as sentinel at the end (preserves legacy behaviour)
    dirVec(end+1, :) = NaN;

    % Replicate the last valid direction into any trailing NaN rows
    id_nan = find(isnan(dirVec(:,1)));
    if ~isempty(id_nan)
        firstNaN = min(id_nan);
        if firstNaN < 2
            error('coordVector: no valid direction row precedes the first NaN row; cannot fill direction vectors');
        end
        dirVec(id_nan, :) = repmat(dirVec(firstNaN - 1, :), numel(id_nan), 1);
    end
end
