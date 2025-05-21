%Funcion for CoordFigure calc
function[height,ray2,sums,kh,kh2] = vdistsum(lat,lon,alt,refr,k,hprog)
j = 0;
try
    for i = 1:k:size(lat,2)-k
        j = j + 1;
        s(j) = vdist(lat(i)*180/pi,lon(i)*180/pi(),lat(i+k)*180/pi,lon(i+k)*180/pi());
        if isnan(s(j))
            s(j) = s(j-1);
        end
    end
end
sums = cumsum(s);
%ray = refr(1:k:j*k);
ray = refr(1:k:size(refr,2));
%ray = refr; %refr(1:k:j*k);
%height = alt(1:k:j*k);
height = alt(1:k:size(refr,2));
kh = min(find(height<hprog));
kh2 = min(find(height(kh+1:end)>hprog))+kh - 1;
ray2 = ray;
try
    for i = 1:size(ray,2)
        if ray(i) < 1
           ray2(i) = 0;
        else
           ray2(i) = (ray(i)-1)*10^6;
        end

    end
end
end