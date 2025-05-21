function mean_voxels = m_voxelsRT(model,par,num_pud,switches)
% to trzeba przepiac tak zeby szukal id wezlow na podstawie id woksela,
% teraz sa problemy a beda tez przy ro pewnie

 radii = [6378.137, 6356.752314245]; 
[wi,ki,wa] = size(par);
clear ki wa
if wi > 1
    hor = model.num_lat_TOMO*model.num_lon_TOMO;
    horind1 = 1:hor-model.num_lon_TOMO+1;
    horind2 = horind1(2:end-1);
    horindx = sort([horind1';horind2']);
    horindx = reshape(horindx,2,size(horindx,1)/2)';
    id = model.num_lon_TOMO-1:model.num_lon_TOMO-1: hor-model.num_lon_TOMO;
    horindx(id,:) = [];  
    horindy = horindx + model.num_lat_TOMO;
    horindh = [horindx+hor,horindy+hor];
    horin = [horindx,horindy, horindh];
    num_vox = [];
    c = diff(num_pud);
    for i = 1:size(model.mid_levels_TOMO,1)
        num_vox = [num_vox;horin+(i-1)*hor];
    end
    if max(num_pud) > size(num_vox,1) 
        if min(num_pud) > size(horin,1)
            mean_voxels = num_vox(num_pud-210,:);
        else
            disp('incorrect coordinates')
        end
    else
        mean_voxels = num_vox(num_pud,:);
    end

 
end
