%% Function to save variables
function ray = savevar(lat_ray_b,lon_ray_b,alt_ray_b,lat_ray,lon_ray,alt_ray,refr,X_ray,X_ray_b,diff_dist,grad_n,grad_ng,ts,ds,ray,n_iter,refrh,refrw,Tempprof,presprof,humprof,id)
    if ds > 100
        ray.de_lat_ray_fin{:,n_iter-1} = lat_ray_b';
        ray.de_lon_ray_fin{:,n_iter-1} = lon_ray_b';
        ray.de_alt_ray_fin{:,n_iter-1} = alt_ray_b';
        ray.de_lat_ray_fins{:,n_iter-1} = lat_ray';
        ray.de_lon_ray_fins{:,n_iter-1} = lon_ray';
        ray.de_alt_ray_fins{:,n_iter-1} = alt_ray';
        ray.refr{:,n_iter-1} = refr';
        ray.X_ray{:,:,n_iter-1} = X_ray;
        ray.X_ray_b{:,:,n_iter-1} = X_ray_b;
        ray.diff_dist(1,n_iter-1) = diff_dist;
        ray.grad_n = grad_n;
        ray.grad_ng = grad_ng;
        ray.t = ts;
        ray.refh = refrh;
        ray.refw = refrw;
        ray.Temp = Tempprof;
        ray.pres = presprof;
        ray.hum = humprof;
        ray.id = id;
    else
        ray.de_lat_ray_fin = lat_ray_b';
        ray.de_lon_ray_fin = lon_ray_b';
        ray.de_alt_ray_fin = alt_ray_b';
        ray.de_lat_ray_fins = lat_ray';
        ray.de_lon_ray_fins = lon_ray';
        ray.de_alt_ray_fins = alt_ray';
        ray.refr = refr';
        ray.X_ray = X_ray;
        ray.X_ray_b = X_ray_b;
        ray.diff_dist(1,n_iter-1) = diff_dist;
        ray.grad_n = grad_n;
        ray.grad_ng = grad_ng;
        ray.t = ts;
        ray.refh = refrh;
        ray.refw = refrw;
        ray.Temp = Tempprof;
        ray.pres = presprof;
        ray.hum = humprof;
        ray.id = id;
end