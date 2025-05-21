function [altstack,Nwstack,zwderr] = intersetionNw(MI,radii,model,epoch,switches)
%% Function to calculate refractivites in the hoorizontal coordinates of the intersection point between ray path and tomography vertical layer
    MI = [squeeze(MI(:,:,1));squeeze(MI(:,:,2))];
    id = isnan(MI(:,1));
    MI(id,:) = [];
    errfunc = [0.219740324581235 0.275687659153536 0.298534562541324 0.142097536297436 0.118997178178812 0.134314930804440 0.118454306153678 0.0695602766473414 0.0728226412061322 0.0523875402957172 0.0532846858017076 0.0356764999929888 0.0284876860767259 0.0166562060524089 0.0103045162085157 0.00703281625698092 0.00331718955353701 0.000414634522075665 0.00128079170310092 0.000205125164377180 0.000138772925997198 0.000472264287511310 0.000453395445643634 0.000521857805132862 0.000323401175654430 0.000268492344632992];
    erralt = [2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.25 8.75 9.25 9.75 10.25 10.75 11.25 11.75 12.25 12.75 13.25 13.75 14.25 14.75];
    if size(MI,1) > 0
        for j = 1:size(MI,1)
            [lon_ray(j),lat_ray(j),alt_ray(j)] = cspice_recgeo(MI(j,:)'  ,radii(1),(radii(1)-radii(2))/radii(1));
            if  alt_ray > 0
                for z = 1:size(model.levels_TOMO,2)
                    N = IDW_atom(model.LAT, model.LON, model.temp{epoch,1}, model.wvpr{epoch,1}, model.pres{epoch,1},model.refr{epoch,1},model.refrRT{epoch,1} ,model.rWGS,[(pi()/2-lat_ray(j)) lon_ray(j), model.levels_TOMO(z)],model.levels_TOMO_RT.*1000',switches);
                    Nw(j,z) = N.Nw;
                end 
                zwderr(1,j) = interp1(erralt,errfunc,alt_ray(j));
            else
                for z = 1:size(model.levels_TOMO,2)
                    Nw(j,z) = NaN;
                end  
                zwderr(1,j) = NaN;
            end
        end
        try
            Nwstack = Nw;
        catch
            Nwstack = []; 
        end
        try
           altstack = min(alt_ray);
        catch
           altstack = NaN;
        end
    else
        altstack = NaN;
        Nwstack = NaN;
	zwderr = NaN;
    end
end

  