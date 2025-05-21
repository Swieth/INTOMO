function  [K,Cov_v,SWD_k,A_k,R_SWD_k] = localtestRT(K,Cov_v,SWD_k,A_k,R_SWD_k,xPminus,loc_value,path_save,epoch,n_obs)
%simple function to remove outlayers

v = SWD_k(1:n_obs)-A_k(1:n_obs,:)*xPminus;
mv = diag(sqrt(Cov_v(1:n_obs,1:n_obs)));
loc = abs(v)./mv;
outl = find(loc>loc_value);
K(:,outl)=[];
keep = find(loc<=loc_value);
Cov_v=Cov_v(keep,keep);
SWD_k(outl)=[];
A_k(outl,:)=[];
R_SWD_k(outl)=[];

close all
fig=figure(1);
plot(SWD_k-A_k*xPminus);
saveas(fig,[path_save '/afterlocaltest_epoch' num2str(epoch,'%02d') '.jpg']);
close(fig);