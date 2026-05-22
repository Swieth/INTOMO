function dirVec = coordVector(coordIn, i_pos)
% COORDVECTOR  Compute chord direction vectors for a simulated ray path.
%
%   For each unique voxel V traversed by the ray, the chord is:
%       exit_V - entry_V
%   where entry_V is coordIn(i_first,:,1) and exit_V is coordIn(i_last,:,2),
%   with i_first and i_last being the first and last ray-point indices whose
%   voxel id equals V.
%
%   This correctly uses both boundary-coordinate sheets of coordIn (as
%   produced by vox_distance.m) and handles multiple ray points per voxel.
%
% INPUT
%   coordIn - [nPts x 3 x 2] XYZ voxel-boundary intersection coordinates.
%             Sheet 1 (:,:,1): entry-side boundary points.
%             Sheet 2 (:,:,2): exit-side boundary points.
%             NaN marks a step where the ray did not cross a boundary.
%   i_pos   - [nPts x 1] or [1 x nPts] voxel index for each ray point
%             (as returned by interateModelsA via find_num).
%
% OUTPUT
%   dirVec  - [nVoxels x 3] chord vectors, one row per element of
%             unique(i_pos) in ascending voxel-id order.
%             Callers must index Avec with unique(i_pos), not ind_sum(1,:).

    entry_pts = squeeze(coordIn(:,:,1));   % [nPts x 3]
    exit_pts  = squeeze(coordIn(:,:,2));   % [nPts x 3]

    i_pos = i_pos(:);
    assert(numel(i_pos) == size(entry_pts, 1), ...
        'coordVector: numel(i_pos) (%d) must equal size(coordIn,1) (%d).', ...
        numel(i_pos), size(entry_pts, 1));

    [unique_vox, ~, ic] = unique(i_pos);   % ascending order
    n_vox  = numel(unique_vox);
    dirVec = nan(n_vox, 3);

    for k = 1:n_vox
        idx = find(ic == k);               % indices of ray points in voxel k
        e   = entry_pts(idx(1),   :);      % entry boundary at first ray point
        x   = exit_pts (idx(end), :);      % exit  boundary at last  ray point
        if any(isnan(e)) || any(isnan(x))
            error('coordVector: missing entry/exit boundary for voxel id %d (row first=%d, last=%d).', ...
                  unique_vox(k), idx(1), idx(end));
        end
        dirVec(k, :) = x - e;
    end
end
