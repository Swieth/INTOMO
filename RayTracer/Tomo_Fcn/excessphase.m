%% Function to calculate excess phase based od ray path points coordinates.
%% Input
% set........... settings
% ray
%% Output
% dLs.......... excess phase [m]
% ray
%      nstepb....... bent segment lenght [km]
%      stepb........ geometric segment lenght [km]
%      step......... straight segment lenght [km]


function [ray,dLs] = excessphase(ray,set)
    if set.refopt == 1
        refr = ray.refr;
    elseif set.refopt == 2
        refr = ray.refw'; 
    elseif set.refopt == 3
        refr = ray.refh';
    end
    X_ray_b = ray.X_ray_b(:,1:size(refr,1));
    X_ray = ray.X_ray(:,1:size(refr,1));
    alt_ray_b = ray.de_alt_ray_fin(1:size(refr,1));
    if set.ground && ~set.stop
        id = size(refr,2);
    elseif ~set.ground && ~set.stop
        [~,id] = min(alt_ray_b);
    end
    if isempty(refr) == 0  &&  isempty(find(refr >= 1)) == 0 
        id_dL = find(refr >= 1);
        for i = 1:size(id_dL,1)-1
            bent_l(i) = (refr(id_dL(i))+refr(id_dL(i+1)))/2*norm(X_ray_b(:,id_dL(i+1))-X_ray_b(:,id_dL(i))); %Bent line length
            straight_l(i) = norm(X_ray(:,id_dL(i+1))-X_ray(:,id_dL(i))); %Straight line length
            geom_l(i) = norm(X_ray_b(:,id_dL(i+1))-X_ray_b(:,id_dL(i))); %Geometrical line length
        end
        dL = (sum(bent_l))*1000 - sum(straight_l)*1000;
        ray.nstepb = bent_l;
        ray.step = straight_l;
        ray.stepb = geom_l;
        dLs = dL;                   
        alt = ray.de_alt_ray_fin(:,end);
        %fprintf('%9.4f %9.6f\n',alt(id),dL)
        ray.altend = alt(id);
    else
        ray.dL = 0;
        ray.refr = 0;
        ray.X_ray_b = 0; 
        dLs = 0;
    end
end