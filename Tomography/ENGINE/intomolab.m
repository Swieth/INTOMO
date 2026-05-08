function [output] = intomolab(station,model,values,pathTOMO,path_save,observation_set,switches)
%The newest version of tomography software INTOMO v. 1.0
%%%%%INPUT%%%%%%%%%%%%%%%%%

%   station.......... structural matrix containing all observations as follows                                     
%   model............ structural matrix containing all informations regarding tomography model structures
%   values........... structural matrix containing inner/outer model data 
%        .Nw_apr                                       : apriori wet refractivity
%        .Nw_out                                       : apriori wet refractivity for outer model
%        .num_Nw_out                                   : number referring to the location of additional wet refractivity for outer model voxel
%        .Nw_obs_out                                   : additional wet refractivity for outer model voxel
%        .num_Nw                                       : number referring to the location of additional wet refractivity for inner model voxel
%        .Nw_obs                                       : additional wet refractivity for inner model voxel
%       - apriori data and additional observations
%        .WV_apr                                       : apriori water vapour
%        .WV_out                                       : aprioriwater vapour for outer model
%        .num_WV_out                                   : number referring to the location of additional water vapour for outer model voxel
%        .WV_obs_out                                   : additional water vapour for outer model voxel
%        .num_WV                                       : number referring to the location of additional /water vapour for inner model voxel
%        .WV_obs                                       : additional water vapour for inner model voxel
%   path_TOMO - location of work directory
%   path_save - location for observation matrices storage files
%   observation_set - matrix with dates of epochs
%   switches -  structural matrix of settings of the INTOMO processing
%%%%%%OUTPUT%%%%%%%%%%%%%%%%%
%   output - output variables from tomography processing
%        .xP                                       : estimated values of refractivities
%        .mxP                                      : errors of refractivities estimation

%% Set initial parameters
if isempty(switches.observations{1,1})  && isempty(switches.observations{1,2})
    sprintf('%s','User did not make a choice of observation type, DEFAULT: SWD \n');
    switches.observations{1,1} = 'SWD';
end
if strcmpi(switches.observations{1,1},'SWD')
    switches.run_SWD = true;
else
    switches.run_SWD = false;
end
% Set apriori parameters
if isempty(switches.apriori{1,1}) && isempty(switches.apriori{1,2})
    sprintf('%s','User did not make a choice of apriori information, checking values variable \n');
    if isempty(values.num_Nw)
        switches.apriori{1,1} = 'INNER';
        sprintf('%s','INNER apriori was detected and set \n');
    end
    if isempty(values.num_Nw_out)
        switches.apriori{1,2} = 'OUTER';
        sprintf('%s','OUTER apriori was detected and set \n');
    end    
end
if strcmpi(switches.apriori{1,1},'INNER')
    switches.run_INNER = true;
else
    switches.run_INNER = false;
end
if strcmpi(switches.apriori{1,2},'OUTER')
    switches.run_OUTER = true;
else
    switches.run_OUTER = false;
end
% Set contraints
if isempty(switches.constraints{1,1})  && isempty(switches.constraints{1,2}) 
    sprintf('%s','User did not make a choice of constraints information, DEFAULT: HORIZONTAL NO, VERTICAL NO  \n');
    switches.constraints{1,1} = 'NO';
    switches.constraints{1,2} = 'NO';
end

if strcmpi(switches.constraints{1,1},'HORIZONTAL')
    switches.run_HORIZONTAL = true ;
else
    switches.run_HORIZONTAL = false; 
end

if strcmpi(switches.constraints{1,2},'VERTICAL') 
    switches.run_VERTICAL = true;
else
    switches.run_VERTICAL = false;
end
% Set stacking parameters
if isempty(switches.stacking{1,1})
    sprintf('%s','User did not make a choice of use of stacking, DEFAULT: NO  \n');
    switches.constraints{1,1} = 'NO';
end
if strcmpi(switches.stacking{1,1},'NO')
    switches.run_STACKING = 1;
    solution_epochs = 1:size(observation_set,1);
elseif isnumeric(switches.stacking{1,1})
    switches.run_STACKING = switches.stacking{1,1};
    solution_epochs = 1:switches.run_STACKING:size(observation_set,1);
end
% Set filter parameters
if isempty(switches.filter{1,1})
    sprintf('%s','User did not make a choice of use of FILTER, DEFAULT: KALMAN  \n');
    switches.constraints{1,1} = 'KALMAN';
end
if strcmpi(switches.filter{1,1},'KALMAN')
    switches.run_KALMAN = 0;
elseif strcmpi(switches.filter{1,1},'ROBUST') %not tested
    switches.run_KALMAN = 1;
end
% Set solution paramaters
if isempty(switches.solution{1,1})
    sprintf('%s','User did not make a choice of use of type of solution, DEFAULT: REAL  \n');
    switches.solution{1,1} = 'REAL';
end
% Set decorelation parameters
if isempty(switches.decorelation{1,1}) %not tested
    sprintf('%s','User did not make a choice of decorrelation of matrix A, DEFAULT: NO  \n');
    switches.decorelation{1,1} = 'NO';
end
if strcmpi(switches.decorelation{1,1},'YES') 
    switches.run_DECOR = 1;
else
    switches.run_DECOR = 0;
end
%% Unpack variables from structures
BLh_pudel_proj = model.BLh_pudel_rad;
BLh_pudel_proj_num = model.BLh_pudel_num_rad;
BLh_outer_proj_num = model.BLh_outer_num_rad;   
num_lat = model.num_lat_TOMO-1;
num_lon = model.num_lon_TOMO-1;
rem_SAT = values.rem_SAT;
rem_REC = values.rem_REC;
if ~strcmp(rem_REC,'') %not tested
    b = 0;
    for j = 1:size(rem_REC,2)
        [wie,~,~] = find(strcmp(model.NAME,rem_REC(1,j))==1);
        if isempty(wie) == 0
            b = b + 1;
            rem_REC_NUM(b,1) = wie;
        end
        clear wie kol num
    end
    clear b i rem_REC
    rem_REC = rem_REC_NUM;
    clear rem_REC_NUM;
else
    rem_REC = [];
end
if switches.run_SWD
     if strcmp(switches.parametrization,'constant')
        Nw_apr = values.Nw_apr;
        num_Nw = values.num_Nw_num;
        Nw_obs = values.Nw_obs_num;
        if strcmp(switches.apriori{1,2},'OUTER')
            Nw_out = values.Nw_out_num;
            num_Nw_out = values.num_Nw_out_num;
            Nw_obs_out = values.Nw_obs_out_num;
        end
     else 
         Nw_apr = values.Nw_apr;
         num_Nw = values.num_Nw;
         Nw_obs = values.Nw_obs;
         if strcmp(switches.apriori{1,2},'OUTER')
            Nw_out = values.Nw_out;
            num_Nw_out = values.num_Nw_out;
            Nw_obs_out = values.Nw_obs_out;
         end 
     end
end
if any(any(isnan(Nw_apr)))
    [num_Nw_out, Nw_obs_out] = aprioriCONSTRNWP([],[],[],[],(1:size(observation_set,1))',(repmat(model.num_outer,size(observation_set,1),1))',Nw_apr); 
    [num_Nw, Nw_obs] = aprioriCONSTRNWP([],[],[],[],(1:size(observation_set,1))',(repmat(model.num_inner',size(observation_set,1),1))',Nw_apr);        
end
if isa(Nw_apr,'single')
    Nw_apr = double(Nw_apr);
end
[~,~,ile_wys] = size(BLh_pudel_proj_num);
%% Calculate basic matrices for tomograpy processing
for epoch = 1:size(observation_set,1)
    disp(['Epoch: ' num2str(epoch) ' / ' num2str(size(observation_set,1))]);
    [A(epoch).A,SWD(epoch).SWD,SIWV(epoch).SIWV,R_SWD(epoch).R_SWD,R_SIWV(epoch).R_SIWV,dt(epoch).not_cro,...
        elevation(epoch).elevation,SAT(epoch).SAT,station_name(epoch).station_name,coord(epoch).coord,dSWD(epoch).dSWD] = ...
        initialOBS(station(epoch),epoch,rem_SAT,rem_REC, pathTOMO, switches, Nw_apr,model);   

    %% Check the consistency of initialOBS matrices with 'station' structure
    try 
        if strcmp(switches.solution,'REAL')
             st = unique(station_name(epoch).station_name);
             % Check whether the same GNSS stations are used for processing
             for stat = 1:size(station(epoch).h,2)
                namestat{stat,1} =  char(station(epoch).h(stat).nazwa);
             end
             if any(elevation(epoch).elevation==0)
                namestat{stat+1,1} = 'RO';
             end
             [~,indst] = intersect(st,namestat);
             if length(indst) ~= length(st)
                 warning('intomolab: Cutting GNSS observation, incosistent station.mat file')
                 id_statf = [];
                for id_del = 1:size(indst,1) 
                    id_stat = find(strcmp(namestat(id_del),station_name(epoch).station_name));
                    id_statf = [id_statf; id_stat];
                end
                % Copy RO rows to new matrices
                id_RO = find(elevation(epoch).elevation==0);
                A(epoch).A = [A(epoch).A(id_statf,:);A(epoch).A(id_RO,:)];
                SWD(epoch).SWD = [SWD(epoch).SWD(id_statf);SWD(epoch).SWD(id_RO)];
                R_SWD(epoch).R_SWD = [R_SWD(epoch).R_SWD(id_statf);R_SWD(epoch).R_SWD(id_RO)];
                elevation(epoch).elevation = [elevation(epoch).elevation(id_statf);elevation(epoch).elevation(id_RO)];
                SAT(epoch).SAT = [SAT(epoch).SAT(id_statf);SAT(epoch).SAT(id_RO)];
                station_name(epoch).station_name = [station_name(epoch).station_name(id_statf);station_name(epoch).station_name(id_RO)];
                % Delete NaN values from dSWD
                id = isnan(dSWD(epoch).dSWD);
                dRO = dSWD(epoch).dSWD(id_RO);
                dSWD(epoch).dSWD(id) = [];
                id = isnan(dRO);
                dRO(id) = nanmean(dRO);
                dSWD(epoch).dSWD = [dSWD(epoch).dSWD(id_statf),dRO];
             end
        end
    end
    %% Calculate apriori matrices
    if switches.run_SWD
        if  switches.run_INNER % INNNER
            A_apriori_t = [];
            SWD_apriori_t = [];
            P_apriori_t = [];
        end
        if  switches.run_OUTER %OUTER
            A_apriori_out = [];
            SWD_apriori_out = [];
            P_apriori_out = [];
        end        
        if strcmp(switches.regular,'yes')
            if  switches.run_INNER
                [A_apriori_t, SWD_apriori_t, ~] = apConstRT(num_Nw(epoch).number,Nw_obs(epoch,:),1);
                if strcmp(switches.parametrization,'constant')
                    P_apriori_t = covarianceAprRT(Nw_apr,A(epoch).A,num_lat,num_lon,size(BLh_pudel_proj,3),BLh_pudel_proj_num,num_Nw(epoch).number,values,epoch,switches);
                elseif strcmp(switches.parametrization,'nodes') || strcmp(switches.parametrization,'bilinear-h')
                    P_apriori_t = covarianceAprRT(Nw_apr,A(epoch).A,num_lat+1,num_lon+1,size(BLh_pudel_proj,3)+1,BLh_pudel_proj,num_Nw(epoch).number,values,epoch,switches);
                end
            end
            if switches.run_OUTER
                [A_apriori_out, SWD_apriori_out, ~] = apConstRT(num_Nw_out(epoch).number,Nw_obs_out(epoch,:),1);
                if strcmp(switches.parametrization,'constant')
                    P_apriori_out = covarianceAprRT(Nw_apr,A(epoch).A,num_lat,num_lon,size(BLh_pudel_proj,3),BLh_pudel_proj_num,num_Nw_out(epoch).number,values,epoch,switches);
                elseif strcmp(switches.parametrization,'nodes') || strcmp(switches.parametrization,'bilinear-h')
                    P_apriori_out = covarianceAprRT(Nw_apr,A(epoch).A,num_lat+1,num_lon+1,size(BLh_pudel_proj,3)+1,BLh_pudel_proj,num_Nw_out(epoch).number,values,epoch,switches);
                end
                if ~isempty(BLh_outer_proj_num)
                    [~,~,ile_wys] = size(BLh_pudel_proj_num);
                    clear ile ko
                    A_apriori_out(:,5:3*3:ile_wys*3*3) = [];
                    A_apriori_out(5:3*3:ile_wys*3*3,:) =[];
                    SWD_apriori_out(5:3*3:ile_wys*3*3,:) = [];
                    P_apriori_out(5:3*3:ile_wys*3*3,:) = [];
                end
            end
        end
        %Vertical/horizontal contraints
        A_consth = [];
        SWD_consth = [];
        P_consth = [];
        A_constv = [];
        SWD_constv = [];
        P_constv = [];
        A_constr(epoch).A_constr = [A_consth; A_constv];
        SWD_constr(epoch).SWD_constr = [SWD_consth; SWD_constv];
        P_constr(epoch).P_constr = [P_consth; P_constv];        
        %% Merge simulated with apriori data
        if  ~isempty(BLh_outer_proj_num) 
            if switches.run_INNER && switches.run_OUTER
                A_apriori(epoch).A_apriori = blkdiag(A_apriori_t,A_apriori_out);
            elseif switches.run_INNER && ~switches.run_OUTER
                fill_outer = zeros(size(A_apriori_t,1),size(BLh_outer_proj_num,2)*size(BLh_outer_proj_num,3));
                fill_outer(:,5:3*3:size(BLh_outer_proj_num,3)*3*3) = [];
                A_apriori(epoch).A_apriori = [A_apriori_t fill_outer];
                clear fill_outer
            elseif ~switches.run_INNER && switches.run_OUTER
                fill_inner = zeros(size(A_apriori_out,1),size(BLh_pudel_proj_num,2)*size(BLh_pudel_proj_num,3));
                A_apriori(epoch).A_apriori = [fill_inner A_apriori_out];
                clear fill_outer
            elseif ~switches.run_INNER && ~switches.run_OUTER
                A_apriori(epoch).A_apriori = blkdiag(A_apriori_t,A_apriori_out);
            end
        else
            A_apriori(epoch).A_apriori = [A_apriori_t; A_apriori_out];
        end
        SWD_apriori(epoch).SWD_apriori = [SWD_apriori_t; SWD_apriori_out];
        P_apriori(epoch).P_apriori = [P_apriori_t; P_apriori_out];
        
        clear A_apriori_t SWD_apriori_t P_apriori_t A_apriori_out SWD_apriori_out P_apriori_out wi ki wa
        clear A_consth A_constv SWD_consth SWD_constv P_consth P_constv
    end
end
%% Start tomography processing
if switches.run_SWD
    % Run LSQ solution
    if strcmp(switches.method,'LSQ') 
        output = TomoLSQ(A,A_apriori,A_constr,SWD_apriori,SWD_constr,P_apriori,P_constr,SWD,R_SWD,elevation,SAT,station_name,Nw_apr, pathTOMO,switches);        
    else
        for epoch = 1 : size(solution_epochs,2)
            %% Initialize Kalman Filter
            if epoch == 1
                %Find outliers in observations
                [A_k,SWD_k,R_SWD_k,elevation_k,SAT_k,station_name_k] = matrices_epochRT(A,A_apriori,A_constr,SWD_apriori,SWD_constr,P_apriori,P_constr,SWD,R_SWD,elevation,SAT,station_name,epoch);
                [P2,id_del]  = weightObs(A(epoch).A,elevation(epoch).elevation,R_SWD_k(1:size(elevation(epoch).elevation,1)),SWD(epoch).SWD, station_name(epoch).station_name,epoch ,model,station,(Nw_apr(solution_epochs(epoch),:))',coord(epoch).coord,dSWD(epoch).dSWD,switches);  
                if ~isempty(P2)
                    A_k(id_del,:) = [];
                    R_SWD_k(id_del,:) = [];
                    SWD_k(id_del,:) = [];
                    elevation_k(id_del,:) = [];
                    R_SWD_ktemp(1:size(P2,1),1:size(P2,1)) = double(P2);
                    R_SWD_ktemp(size(P2,1)+1:size(R_SWD_k,1),size(P2,1)+1:size(R_SWD_k,1)) = full(diag(R_SWD_k(size(double(P2),1)+1:end)));
                    R_SWD_k = R_SWD_ktemp;
                    %checking for outlayers in the data
                    [A_k,SWD_k,R_SWD_k,~,~,~,id_weights] = filterOBSRT(A_k,SWD_k,R_SWD_k,elevation_k,SAT_k,station_name_k,(Nw_apr(solution_epochs(epoch),:))',SWD(epoch).SWD,switches);
                    %R_SWD_k(id_weights,id_weights) = R_SWD_k(id_weights,id_weights)*1000;
                    clear Nw_out2 R_SWD_ktemp;
                    if isempty(BLh_outer_proj_num)==0 %separate model inner and outer
                        Nw_ap = Nw_out(1,:);
                        Nw_ap(:,5:3*3:ile_wys*3*3) = [];
                        xPminus = [Nw_apr(1,:) Nw_ap];
                        clear Nw_ap
                    else   
                        xPminus = [Nw_apr(1,:)];
                    end
                    xPminus = xPminus';
                    n_c = size(A_apriori(epoch).A_apriori,1);
                    n_obs = size(SWD_k,1)-n_c;
                    %state predictor 
                    Phi(1:size(A_k,2)) = 1;
                    Phi = diag(Phi);
                    %state predicted disturbance noise 
                    if strcmp(switches.parametrization,'constant')
                        Q = covarianceNwRT(Nw_apr,A_k,num_lat,num_lon,size(BLh_pudel_proj,3),BLh_pudel_proj_num);
                    elseif strcmp(switches.parametrization,'bilinear-h')
                        Q = covarianceNwRT(Nw_apr,A_k,num_lat+1,num_lon+1,size(BLh_pudel_proj,3)+1,BLh_pudel_proj);
                    end
                    if strcmp(switches.parametrization,'constant')
                        Pzero = covarianceAprRT(Nw_apr,A(epoch).A,num_lat,num_lon,size(BLh_pudel_proj,3),BLh_pudel_proj_num,1:size(A(epoch).A,2),values,epoch,switches);
                    elseif strcmp(switches.parametrization,'bilinear-h')
                        Pzero = covarianceAprRT(Nw_apr,A(epoch).A,num_lat+1,num_lon+1,size(BLh_pudel_proj,3)+1,BLh_pudel_proj,1:size(A(epoch).A,2),values,epoch,switches);
                    end
                    Pzero=diag(Pzero.^2);
                    Pminus = Phi * Pzero * Phi+ Q;
                    if switches.run_KALMAN==0 || switches.run_KALMAN==1
                        id=isnan(diag(R_SWD_k));
                        R_SWD_k(id,id)=0.1;
                        [K, con(epoch,1),Y(epoch).Y,~,thin] = getGain(Pminus,A_k,R_SWD_k);
                    end
                    thin.id_del = id_del;
                    %Corrected state estimate
                    xPplus = xPminus + K * (SWD_k - A_k*xPminus);
                    %Corrected state covariance
                    Pplus = Pminus - K *  A_k * Pminus;
                    Pplusall(epoch).Pplus = Pplus;
                    thin.Pplus = diag(Pplus);
                    thin.Pminus = diag(Pminus);
                    %thin.K = K;
                    if ~isvector(P2)
                       thin.P2 = diag(P2);
                       thin.R_SWD_k = diag(R_SWD_k);
                    else
                       thin.P2 = P2;
                       thin.R_SWD_k = R_SWD_k;
                    end
                    thin.Q = diag(Q);
                    %% Save Kalman parameters to output file
                    try
                        save([path_save,switches.project_name,'/KAL/kalman_mat_',num2str(epoch),'.mat'],'thin')
                    catch
                        warning('intomolab: Unable to save kalman filtering errors. Missing path')
                    end
                    clear A_ll A_c SWD_ll Cov_ll P_ll Nw_0 Cov_c P_c res_SWD res_Nw Cov_ll_adj r_ll Cov_c_adj r_c P_ll P_c R_ll R_c R_SWD_k
                    output.n_c(epoch) = n_c;
                    output.n_obs(epoch) = n_obs;
                    xPplus = full(xPplus);
                    xPminus = full(xPminus);
                    clear A_k SWD_k R_SWD_k
                    mxP(epoch,:) = sqrt(diag(Pplus));
                    xP(epoch,:) = xPplus;
                else
                    error(['intomolab: Empty matrix P2. Processing failed on epoch ', num2str(epoch)])
                end
            else
                %% Normal cycle for Kalman processing
                A_t = [];
                SWD_t = [];
                R_SWD_t = [];
                %%%%%%Stacking of multiple epochs
                for j = 0:switches.run_STACKING-1
                    [A_t,SWD_t,R_SWD_t,elevation_t,SAT_t,station_name_t] = matrices_epochRT(A,A_apriori,A_constr,SWD_apriori,SWD_constr,P_apriori,P_constr,SWD,R_SWD,elevation,SAT,station_name,solution_epochs(epoch)-j);
                    clear A_k SWD_k R_SWD_k
                    n_c =  size(SWD_apriori(solution_epochs(epoch)-j).SWD_apriori,1);
                end
                [P2,id_del]  = weightObs(A(epoch).A,elevation(epoch).elevation,R_SWD_t(1:size(elevation(epoch).elevation,1)),SWD(epoch).SWD, station_name(epoch).station_name,epoch ,model,station,(Nw_apr(solution_epochs(epoch),:))',coord(epoch).coord,dSWD(epoch).dSWD,switches);  
                A_t(id_del,:) = [];
                R_SWD_t(id_del,:) = [];
                SWD_t(id_del,:) = [];
                elevation_t(id_del,:) = [];
                try
                    R_SWD_ttemp(1:size(P2,1),1:size(P2,1)) = double(P2);
                    R_SWD_ttemp(size(P2,1)+1:size(R_SWD_t,1),size(P2,1)+1:size(R_SWD_t,1)) = diag(R_SWD_t(size(P2,1)+1:end));
                    R_SWD_t = R_SWD_ttemp;
                    %checking for outlayers in the data
                    [A_t,SWD_t,R_SWD_t,elevation_t,SAT_t,station_name_t,id_weights] = filterOBSRT(A_t,SWD_t,R_SWD_t,elevation_t,SAT_t,station_name_t, (Nw_apr(solution_epochs(epoch),:))', SWD(solution_epochs(epoch)).SWD,switches);
                    %R_SWD_k(id_weights,id_weights) = R_SWD_k(id_weights,id_weights)*1000;
                    clear Nw_out2 R_SWD_ttemp
                    n_obs = size(SWD_t,1)-n_c;
                    if switches.run_DECOR == 1 %decorelation of A and removal of SWD observations [1 -on, 0 -off]
                        remA = Adecor(A_t);
                        A_t(remA,:) = [];
                        SWD_t(remA,:) = [];
                        R_SWD_t(remA,:) = [];
                        clear remA;
                    end
                    %% Kalman state prediction
                    Phi = (Nw_apr(epoch,:)')./(Nw_apr(epoch-1,:)');
                    Phi(find(Phi==Inf))=1;
                    Phi(find(isnan(Phi)))=1;
                    Phi = diag(Phi);
                    %state prediction
                    try
                        xPminus =Phi * xPplus;
                    catch
                        xPminus =Phi * xPplus';    
                    end
                    if strcmp(switches.parametrization,'constant')
                        Q = covarianceNwRT(Nw_apr,A_t,num_lat,num_lon,size(BLh_pudel_proj,3),BLh_pudel_proj_num);
                    elseif strcmp(switches.parametrization,'nodes') || strcmp(switches.parametrization,'bilinear-h')
                        Q = covarianceNwRT(Nw_apr,A_t,num_lat+1,num_lon+1,size(BLh_pudel_proj,3)+1,BLh_pudel_proj);
                    end
                    Qall(epoch).Q=Q;
                    %covariance state predictor
                    Pminus = Phi * Pplus * Phi' + Q;
                    if switches.run_KALMAN ==1 %switching on or of robust filtering (DOWNWEIGHTING) [1 -on, 0 -off]
                        id=isnan(diag(R_SWD_t));
                        R_SWD_t(id,id)=1;
                        [K,~,thin] = getGainR(Pminus,A_t,R_SWD_t,SWD_t,xPminus);
                        thin.id_del = id_del;
                        xPminusit = xPminus;
                        Pminusit = Pminus;
                        m_alter = 0;
                        while m_alter<1 % usually 2 repeats of robust filter is enought
                            m_alter = m_alter + 1;
                            xPminusit = xPminusit + K*(SWD_t - A_t*xPminusit);
                            try
                                [K, ~] = getGainR(Pminus,A_t,R_SWD_t,SWD_t,xPminusit);
                            catch
                                warning('intomolab: Kalman inversion failed in robust filtering')
                            end
                        end
                        clear xPminusit
                     elseif switches.run_KALMAN ==0
                        [K, con(epoch,1),Y(epoch).Y,~] = getGain(Pminus,A_t,R_SWD_t);
                     end
                     %% Kalman state correction
                     %Corrected state estimate
                     xPplus = xPminus + K * (SWD_t - A_t*xPminus);
                     %Corrected state covariance
                     Pplus = Pminus - K *  A_t * Pminus;
                     Pplusall(epoch).Pplus = Pplus;
                     thin.Pplus = diag(Pplus);
                     thin.Pminus = diag(Pminus);
                     if ~isvector(P2)
                        thin.P2 = diag(P2);
                        thin.R_SWD_k = diag(R_SWD_t);
                     else
                        thin.P2 = P2;
                        thin.R_SWD_k = R_SWD_t;
                     end
                     %thin.K = K;
                     thin.Q = diag(Q);
                    %% Save Kalman parameters to output file
                    try
                        save([path_save,switches.project_name,'/KAL/kalman_mat_',num2str(epoch),'.mat'],'thin')
                    catch
                        warning('intomolab: Unable to save kalman filtering errors. Missing path')
                    end
                    output.n_obs(epoch) = n_obs;
                    output.n_c(epoch) = n_c;
                    CovX = Pplus;
                    mxP(epoch,:) = sqrt(diag(Pplus));
                    xP(epoch,:) = xPplus;
                catch
                     error(['intomolab: Epoch processing failed on epoch ',num2str(epoch)])
                end
            end
        end
        CovxP = [];
        output.xP = xP;
        output.mxP = mxP;
        clear xPlus xPminus Pminus Q Pzero Phi
    end
end
