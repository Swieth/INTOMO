function output = TomoLSQ(A,A_apriori,A_constr,SWD_apriori,SWD_constr,P_apriori,P_constr,SWD,R_SWD,elevat,SAT,station_name,Nw_apr,solution_epochs, model, path_save, switches)
% Function to calculate tomography with Least Square method.
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       A..........   the distances of each signal in tomography voxel
%       A_apriori..   the A apriori matrix with 0 or 1 values
%       A_constr...   constraints given to A matrix
%       SWD_apriori   the apriori values of Slant Delay (SD)
%       SWD_constr.   constraints given to SD 
%       P_apriori..   the apriori weights for the tomography grid
%       P_constr...   constraints given to P 
%       SWD........   the values of slant delays
%       elevat.....   elevation angle to satellite
%       SAT........   id of GPS satellite
%       station_name  names of GNSS stations
%       Nw_apr.....   apriori filed of refractivity
%       path_save..   path to save matrices
%       swtiches...   settings of the processing
%%%OUTPUT:
%       output
%             xP......  estiamted with tomgoraphy processing refractivity values
%             mxP.....  error of estimated refractivity values
%             dxP.....  difference between estimated and apriori values

%% Find id of nodes in columns 
NwDET = Nw_apr.Nw_apr_num;
NwEra = Nw_apr.aprioriEra;
if switches.cutModel 
    if strcmp(switches.parametrization,'constant')
        num_lat = length(model.mid_lat_TOMO);
        num_lon = length(model.mid_lon_TOMO);
        num_lvl = length(model.levels_TOMO)-1;
    elseif strcmp(switches.parametrization,'nodes')
        num_lat = length(model.lat_TOMO);
        num_lon = length(model.lon_TOMO);
        num_lvl = length(model.levels_TOMO);
    end
    idvox = [];
    range_intervals = 1:num_lat*num_lon:size(A(1).A,2);
    num_ranges = length(range_intervals) - 1;
    for i = 1:num_lat*num_lon
        intervals_hor(i,:) = i:num_lat*num_lon:size(A(1).A,2);
    end
    intervals_hor = reshape(intervals_hor,num_lat,num_lon,num_lvl);
end
num_latO = num_lat;
num_lonO = num_lon;
num_lvlO = num_lvl;
switches.roonly = true;
for epoch = 1 : size(solution_epochs,2)
    isro = true;
    %% Find intersected voxels and layers
    if switches.cutModel
         elev = elevat(epoch).elevation;
         if switches.roonly
             id = find(elev~=0);
             A(epoch).A(id,:) = [];
             SWD(epoch).SWD(id) = [];
             R_SWD(epoch).R_SWD(id) = [];
             elevat(epoch).elevation(id) = [];
             elev(id) = [];
         end
         int_vox = find(sum(A(epoch).A(1:length(elev),:),1)~=0);
         counts_ind = zeros(1, num_ranges);
         for i = 1:num_ranges
            counts_ind(i) = sum(unique(int_vox) >= range_intervals(i) & unique(int_vox) < range_intervals(i+1));
         end
         counts_ind = find(counts_ind>0);  % find(counts_ind>5); 
         for x = 1:size(intervals_hor,1)
            if any(intersect(int_vox,squeeze(intervals_hor(x,:,:)))) % length(intersect(int_vox,squeeze(intervals_hor(x,:,:))))>5
                counts_horX(x) = true;
            else
                counts_horX(x) = false;
            end
         end
         for y = 1:size(intervals_hor,2)
            if any(intersect(int_vox,squeeze(intervals_hor(:,y,:)))) % length(intersect(int_vox,squeeze(intervals_hor(:,y,:))))>5
              counts_horY(y) = true;
            else
              counts_horY(y) = false;
            end
         end
         verLeft = [range_intervals(min(counts_ind)):size(A(epoch).A,2)];
         horLeft = intervals_hor(counts_horX,counts_horY,:);
         horLeft = reshape(horLeft,size(horLeft,1)*size(horLeft,2)*size(horLeft,3),1);
         indLeft = intersect(verLeft,horLeft);
         indLeftTrue = indLeft;
         if ~isempty(indLeft)
             num_lat = sum(counts_horX); 
             num_lon = sum(counts_horY);
             num_lvl = length(model.levels_TOMO) -  min(counts_ind)+1;
             %% Redefine apriori matrices
             BLHstruc = pudel2(model,model.lat_TOMO(counts_horX),model.lon_TOMO(counts_horY),model.levels_TOMO(min(counts_ind):end),switches);
             %[~,BLHstruc.num_inner] = intersect(BLHstruc.num_inner,indLeft);
             [num_Nw, Nw_obs] = aprioriCONSTRNWP([],[],[],[],(1:size(NwDET,1))',(repmat(BLHstruc.num_inner',size(NwDET,1),1))',NwDET(:,indLeft));
             [A_aprioriE, SWD_aprioriE, ~] = apConstRT(num_Nw(epoch).number,Nw_obs(epoch,:),1);
             Acut = A(epoch).A;
             Acut = Acut(:,indLeft);
             A(epoch).A = Acut;
             SWD_apriori(epoch).SWD_apriori = SWD_aprioriE;
             A_apriori(epoch).A_apriori = A_aprioriE;
             values.Nw_DETER = NwDET(:,indLeft);
             if epoch == 1
                range_intervals = [range_intervals, size(NwDET,2)];
             end
             for lvl = 1:length(range_intervals)-1
                Pmean1(lvl,1) = abs(mean(NwDET(epoch,range_intervals(lvl):range_intervals(lvl+1)) - NwEra(epoch,range_intervals(lvl):range_intervals(lvl+1))));
                if lvl>5
                    Pmean1(lvl,1) = Pmean1(lvl,1)./1000;    
                end
             end
             Pmean2 = repmat(Pmean1',num_latO*num_lonO,1);
             Pmean = reshape(Pmean2,num_latO*num_lonO*(num_lvlO),1);
             P_apriori(epoch).P_apriori = Pmean(num_Nw(epoch).number)./100;
             clear Pmean
         else
             disp(['No RO profiles on epoch ', num2str(epoch)])
             isro = false;
         end
         %indLeft = int_vox;  % only intersected voxels how to fit apriori
         %data
    else
         indLeft = 1:size(A(epoch).A,2);
         indLeftTrue = indLeft;
    end
    if isro
        %% Stack matrices
        [A_k,SWD_k,R_SWD_k,elevation_k,SAT_k,station_name_k] = matrices_epochRT(A,A_apriori,A_constr,SWD_apriori,SWD_constr,P_apriori,P_constr,SWD,R_SWD,elevat,SAT,station_name,epoch);
        indLeft = 1:size(A_k,2);
        if strcmp(switches.regular,'yes')
            R_SWD_k = diag(R_SWD_k);
            [A_k,SWD_k,R_SWD_k,elev,~,~,idweight] = filterOBSRTRO(A_k,SWD_k,R_SWD_k,elevation_k,SAT_k,station_name_k,(NwEra(solution_epochs(epoch),indLeftTrue))',SWD(epoch).SWD,switches);
        end   
        if ~isempty(A_k)
            %Tikhonov Damped Least Squares (Levenberg-Marquardt)
            %% Calculate Singluar Value Decomposition    
            Pdist = 1./abs(SWD_k(1:length(elev)) - A_k(1:length(elev),indLeft)*NwEra(epoch, indLeftTrue)')*1000;
            Pdist(idweight) = Pdist(idweight)./1000;
            R_SWD_k = diag(R_SWD_k);
            P=diag([Pdist;R_SWD_k(length(Pdist)+1:end,1)]); 
            Pall(epoch).P = P;
            A_k = A_k(:,indLeft);
            [ile_w, ~] = size(A_k);
            if issparse(A_k) == 1
                [U,S,V] = svds(A_k,ile_w);
            else
                [U,S,V] = svd(A_k,'econ');
            end
            %% Find condition number
            cond_A = max(diag(S)) / max(min(diag(S)), 1e-12)
            disp(epoch)
            %if epoch == 1
                if  cond_A < 1e3 
                    lambda = 1e-6; % Low damping
                elseif cond_A >= 1e3 && cond_A < 1e6
                    lambda = 1e-3; % Medium damping
                else
                    lambda = 1e-1; % High damping
                end
            %end
            %D = diag(diag(A_k' * A_k));
            Alfa = diag(1./(diag(S^2)+lambda));
            D = V*Alfa*V';
            %% Calculate estimated refractivity
            Anew = A_k'*P*A_k+lambda*D;
            b = A_k'*P*SWD_k+lambda*D* NwDET(epoch, indLeftTrue)';
            Xtemp = Anew \ b;
            xP(epoch,:) = NwDET(epoch, :);
            xP(epoch,indLeftTrue) = Xtemp;
            xP(epoch,isnan(xP(epoch,:))) = NwDET(epoch, isnan(xP(epoch,:)));
            %% Adaptive Lambda Update Based on Residuals
            residual = norm(A_k*xP(epoch,indLeft)'-SWD_k);
            if exist('residual_old', 'var') && residual > residual_old
                lambda = lambda * 10; % Increase lambda if error increases
            else
                lambda = lambda * 0.1; % Reduce lambda if error decreases
            end
             residual_old = residual; % Store residual for next iteration
	     dxP(epoch).dxP=xP(epoch,:)-NwDET(epoch,:);
        else
            xP(epoch).xP=NwDET(epoch, :)';
            dxP(epoch).dxP=[];
        end
    else
        xP(epoch,:) = NwDET(epoch, :)';
        dxP(epoch).dxP=[];
    end
end
%% Save values
if strcmp(switches.regular,'yes')
    output.xP = full(xP);
end
output.mxP = NaN(size(output.xP,1),size(output.xP,2));
output.dxP = dxP;
output.P = Pall;
output.S = Sc;