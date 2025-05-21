function [wie] = Adecor(A)
%producing number of rows/columns from SWD and A R that has to be removed
%Filtering A matrix to get into processing only independent columns based on the
%assumption that:
%   1. There is an inflection point in singular values
%   2. Condition number is limited to 300 -> multiplication of Covariance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ile_w, ile_kol] = size(A);
if issparse(A) == 1
    [U,S,V] = svd(full(A));
else
    [U,S,V] = svd(A);
end
proc1 = log(diag(S));

[bzd,bzd2 ,bzd3] = find(proc1>(-30));
proc2 = diff(proc1(bzd),2);
theta = atan(proc2);
clear bzd*
[sorted, theta_ind] = sortrows(abs(theta),-1);

num = (theta_ind(1));
k = 1;
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
con_test = SD(1,:)./SD(:,1);
if con > 300
    [row,col,num] = find(con_test<300);
    clear col num
    num = row(end,1);
    con = con_test(num,1);
end

SD = SD(1:(num(1,1)-1),:);
Sminus = SD;
Sminus = [diag(Sminus) zeros((num(1,1)-1),ile_kol-(num(1,1)-1)); zeros(ile_w-(num(1,1)-1),ile_kol)];
Adec = U*Sminus*V';
wie = find(mean((A-Adec)')>2* std(mean((A-Adec)')));