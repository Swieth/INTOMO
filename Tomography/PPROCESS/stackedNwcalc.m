function [ZWDROstack,idout] = stackedNwcalc(Nwstack,height)
%% Function to calculate Zenith Delays from refractivities in the intersection point between ray path and voxel boundaries
    idout = [];
    nH = numel(height);

    for i = 1:numel(Nwstack)
        Nw = Nwstack{i};
        if size(Nw,1) < 2
            idout(i) = true;
        end
        
        for j = 1:size(Nw,1)
            try
                % Extract layer-midpoint refractivity profile (one value per layer)
                prof = Nw(j,:)';
                % Profile must have exactly nH elements (one per layer midpoint)
                if numel(prof) ~= nH
                    % Degenerate ray: profile/height dimension mismatch
                    ZWDRO(j) = NaN;
                    ZWDROstack(i,j) = ZWDRO(j);
                    continue
                end
                ZWDRO(j) = excessphasezenith(prof,height); 
            catch ME
                warning('stackedNwcalc: %s', ME.message);
                ZWDRO(j) = NaN;
            end
            ZWDROstack(i,j) = ZWDRO(j);
         end
    end
    
    id = [1:size(Nwstack,2)];
    try
        idout = id(idout);
    end
    
    ZWDROstack(ZWDROstack==0) = NaN;
end