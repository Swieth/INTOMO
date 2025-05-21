function [K,Wd,thin] = getGainR(Pminus,A,R_SWD,SWD,xPminus)
%Robust estimation of gain matrix
%%%INPUT
%       K..........Kalman gain matrix
%       Wd.........modified Kalman weights 
%       thin.......structure to save Kalman filtering parameters
%%%OUTPUT
%       A..........observational matrix with signal deriveratives
%       R_SWD......the errors of the slant delays
%       Pminus.....Kalman filtering predicted covariance matrix of estimation uncertainties
%       xPminus....Kalman filtering predicted covariance matrix of
%       estimation uncertainties from previous epoch

%R = diag(R);
sig = 1; %the apriori sigma for observations
c = 2;
W = (1./(R_SWD.^2));
try
    e = SWD - A*double(xPminus);
catch
    warning('getGainR: Kalman processing failed')
end
clear U S;
for i = 1:size(W,1)
   if sqrt(W(i,i))*abs(e(i,1)) > sig*c
        W(i,1) = c*sig*sqrt(W(i,1))/abs(e(i,1));
    end
end
Wd = (1./sqrt(W));
Anew = (A*Pminus*A' + sig^2*(diag(W.^(-1))));   
[ile_w, ile_kol] = size(Anew);
if issparse(Anew) == 1
    [U,S,V] = svds(Anew,ile_w);
else
    try
        [U,S,V] = svd(Anew,'econ');
    catch
        [U,S,V] = svdecon(Anew);
        warning('getGainR: SVD inversion failed. Using economic SVD function')
    end
end
proc1 = log(diag(S));
[bzd,bzd2 ,bzd3] = find(proc1>(-30));
proc2 = diff(proc1(bzd),2);
theta = atan(proc2);
clear bzd*
[sorted, theta_ind] = sortrows(abs(theta),-1);
k = 1;
%to make sure singular values from the beginnig are not selected
num = theta_ind(1);
while (theta_ind(k)<0.1*size(Anew,1) & k<size(theta_ind,1))
    k = k + 1;
    num = (theta_ind(k));
end
SD = diag(S);
if num(1,1)>1
    con = SD(1,:)./SD((num(1,1)-1),:);
else
    con = NaN;
end
if con > 300
    con_test = SD(1,:)./SD(1:(num(1,1)-1),:);
    [row,col,num] = find(con_test<300);
    clear col num
    num = row(end,1);
    con = con_test(num,1);
end
SD = SD(1:(num(1,1)-1),:);
SD = 1./SD;
SD = diag(SD);
VD = V';
VD = VD(1:(num(1,1)-1),:);
VD = VD';
Splus = [SD zeros((num(1,1)-1),ile_kol-(num(1,1)-1)); zeros(ile_w-(num(1,1)-1),ile_kol)];
Ainv = V*Splus'*U';

K = Pminus*A'*Ainv;
thin.thind = sort(theta_ind(1:num-1));
thin.theta = theta;
thin.Anew = diag(Anew);