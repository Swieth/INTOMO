function es=t2es(t)
for i = 1: size(t,1)
    for j = 1:size(t,2)
        if t(i,1)>273.16
          es(i,j)=6.11*exp((-2.500*10^6/461.525)*(1/t(i,j)-1/273.16));
         else
          es(i,j)=6.11*exp((-2.834*10^6/461.525)*(1/t(i,j)-1/273.16));
        end
    end
end