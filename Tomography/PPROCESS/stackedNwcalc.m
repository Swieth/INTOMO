function [ZWDROstack,idout] = stackedNwcalc(Nwstack,height)
%% Function to calculate Zenith Delays from refractivities in the intersection point between ray path and voxel boundaries
    idout = [];
    for i = 1:size(Nwstack,2)
        Nw = Nwstack{1,i};
        if size(Nw,1) < 2
            idout(i) = true;
        end
         for j = 1:size(Nw,1)
            try
                ZWDRO(j) = excessphasezenith(Nw(j,:)',height); 
            catch
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