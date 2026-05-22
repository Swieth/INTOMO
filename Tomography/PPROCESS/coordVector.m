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
%   Edge cases handled gracefully (zero chord, warning):
%     - Ray enters the tomography domain already inside a voxel: sheet 1
%       at row 1 is NaN because no voxel-face crossing was recorded by
%       vox_distance.  The forward scan finds the first available entry; if
%       none exists (entire traversal started at domain boundary) a zero
%       chord is assigned.
%     - Ray exits the tomography domain while still inside a voxel: sheet 2
%       at the last row is NaN for the same reason; the backward scan is used.
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
        % Fallback: scan for first non-NaN entry and last non-NaN exit within
        % this voxel's rows. This handles the domain-entry edge case where the
        % ray enters the tomography domain already inside a voxel so vox_distance
        % never records a sheet-1 entry for those initial rows.
        if any(isnan(e))
            entry_rows   = entry_pts(idx, :);
            first_valid  = find(~any(isnan(entry_rows), 2), 1, 'first');
            if ~isempty(first_valid)
                e = entry_rows(first_valid, :);
            end
        end
        if any(isnan(x))
            exit_rows   = exit_pts(idx, :);
            last_valid  = find(~any(isnan(exit_rows), 2), 1, 'last');
            if ~isempty(last_valid)
                x = exit_rows(last_valid, :);
            end
        end

        if any(isnan(e)) || any(isnan(x))
            warning('coordVector: no resolvable entry/exit boundary for voxel id %d (n_pts=%d); zero chord assigned.', ...
                    unique_vox(k), numel(idx));
            dirVec(k, :) = 0;
            continue
        end
        dirVec(k, :) = x - e;
    end
end
