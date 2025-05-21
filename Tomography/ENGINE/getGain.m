function [K,con,Ainv, Anew,thin] = getGain(Pminus,A,R_SWD)
%Function of estimation of gain matrix
%%%INPUT
%       K..........Kalman gain matrix
%       con........condition number
%       Ainv.......inverted A matrix
%       Anew.......weighted Kalman A matrix
%       thin.......structure to save Kalman filtering parameters
%%%OUTPUT
%       A..........observational matrix with signal deriveratives
%       R_SWD......the errors of the slant delays
%       Pminus.....Kalman filtering predicted covariance matrix of estimation uncertainties

Anew = (A*Pminus*A' + R_SWD.^2);
[ile_w, ile_kol] = size(Anew);
if issparse(Anew) == 1
    [U,S,V] = svds(Anew,ile_w);
else
    try
        [U,S,V] = svd(Anew,'econ');
    catch
        [U,S,V] = svdecon(Anew);
        warning('getGain: SVD inversion failed. Using economic SVD function')
    end
end

proc1 = log(diag(S));

[bzd,~ ,~] = find(proc1>(-30));
proc2 = diff(proc1(bzd),2);
theta = atan(proc2);
clear bzd*
[~, theta_ind] = sortrows(abs(theta),-1);
k = 1;
%to make sure singular values from the beginnig are not selected
num = theta_ind(1);
while (theta_ind(k)<0.1*size(A,1) & k<size(theta_ind,1))
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
    [row,~,~] = find(con_test<300);
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