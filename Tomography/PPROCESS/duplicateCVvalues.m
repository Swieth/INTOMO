function [var_obs,var_obsR] = duplicateCVvalues(sZ,sR,cZWD,cR_ZWD,A_stC,ZWD)  
    var_obs = [];
    var_obsR = [];
    % covariance
    try
        fZWD = zeros(sZ-sR,sZ-sR);
        fZWDR = zeros(sZ-sR,sZ-sR);
        h = size(cZWD,1);
        c = 1;
        for i = 1:size(ZWD,2)
            A_st = A_stC{i,:,:};
            vara = repmat(cZWD(i,i),size(A_st,1),size(A_st,1));
            varaR = repmat(cR_ZWD(i,i),size(A_st,1),size(A_st,1));
            c1 = c + size(A_st,1);
            fZWD(c:c1-1,c:c1-1) = vara;
            fZWDR(c:c1-1,c:c1-1) = varaR;
            b = c1;
            for j = i+1:h
                A_st2 = A_stC{j,:,:};
                covar = repmat(cZWD(i,j),size(A_st,1),size(A_st2,1));
                covarR = repmat(cR_ZWD(i,j),size(A_st,1),size(A_st2,1));
                b1 = b + size(A_st2,1);
                fZWD(c:c1-1,b:b1-1) = covar;
                fZWD(b:b1-1,c:c1-1) = covar';
                fZWDR(c:c1-1,b:b1-1) = covarR;
                fZWDR(b:b1-1,c:c1-1) = covarR';
                b = b1;
            end
            c = c1;
         end
         var_obs = fZWD;
         var_obsR = fZWDR;
    catch
      error('Covariance matrix calculation failed')  
    end