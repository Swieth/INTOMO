function [C,idcol]  = weightObs(A,elev,R_SWD,SWD_obs,station_name,epoch,model,station,Nw_apr,coord,dSWD,switches)
% Function to calculate weights of observations based on Zenith Wet Delays
% and LSQ outliers detection
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       A..........   observational matrix with signal deriveratives
%       R_SWD......   the errors of the salnt delays
%       SWD_obs....   the values of slant delays
%       station_name  names of GNSS stations
%       epoch......   processing epoch number
%       model......   parameters of tomography and ray tracing models
%       station....   structural matrix containing all informations regarding GNSS observations
%       Nw_apr.....   apriori values of refractivities
%       coord......   coord. of interesetcions between RO and voxel model
%       dSWD.......   difference between observed and simulated slant delay
%       swtiches...   settings of the processing
%%%OUTPUT
%       C..........   the weighting matrix for observations
%       idcol......   id numbers of observations deleted for further processing


%% Initial values
radii = model.radii; 
idoutRO = [];
icut = 0;
%Save original variables
A_b = A;
SWD_obs_b = SWD_obs;
% Chack apriori Nw matrix
if any(isnan(Nw_apr))
    error(['weightObs: There are NaN values in apriori Nw at epoch '],num2str(epoch))
end
%% Recalculate observed to simulated SWD for real data
if strcmp(switches.solution,'REAL')
    nsim = 1000;
    if ~size(dSWD,2)== 1
        dSWD = dSWD';
    end
    dSWD(isnan(dSWD)) = [];
    try
        SWD_obs = SWD_obs-dSWD(1:size(SWD_obs,1))* nsim;
    catch
        id = find(elev==0);
        l_dSWD = length(dSWD);
        sSWD_obs = length(SWD_obs)-length(id);
        dSWD = [dSWD(1:sSWD_obs,1);zeros(length(id),1)+0.001];
        SWD_obs = SWD_obs-dSWD(1:size(SWD_obs,1))* nsim;
        if size(A,1) ~= l_dSWD+length(id)
            warning(['weightObs: Incorrect dSWD. Check the amtrix data. Results might be corrupted on epoch ',num2str(epoch)])
        end
    end
else
     nsim = 1;
end
%% Calculate ZWD covariance matrix  
% Find RO observations
un_name = unique(station_name);
namestat =  {station(epoch).h.nazwa}';
if any(strcmp(un_name,'RO'))
    un_name(strcmp(un_name,'RO')) = [];
    idLenGr = min(find(elev == 0)) - 1;
else
    idLenGr = size(A,1);
end
if size(R_SWD,2)>1
    R_SWD = R_SWD';
end
% Find ZWD for each station
for i = 1:size(un_name,1)
    try
        stat = station(epoch).h(find(strcmp(namestat,un_name(i)))).satellite;
    catch
        warning('weightObs: Incorrect reference station file. Check the reference station.mat or amtrix data')
    end
    id = find(strcmp(station_name,un_name(i)));
    elev_st = elev(id);
    A_st = A(id,:);
    SWD_st = SWD_obs(id);
    R_SWD_st = R_SWD(id); 
    % Find VMF
    elev_check = cell2mat({stat.elevation})';
    vmfh = cell2mat({stat.vmf1h})';
    vmfw = cell2mat({stat.vmf1w})';      
    elev_check(isnan(elev_check)) = [];
    vmfh(isnan(vmfh)) = [];
    vmfw(isnan(vmfw)) = [];
    % Find common observed and ray traced rays
    [id_st,id_check] = ismember(round(elev_st,2),round(elev_check,2));
    id_temp = 1:size(round(elev_st,2),1);
    id_st = id_temp(id_st);
    elev_st = elev_st(id_st);
    vmfh = vmfh(id_check);
    vmfw = vmfw(id_check);
    SWD_st = SWD_st(id_st);
    % Temporalily save data for validation
    A_stC{i,:,:} = A_st(id_st,:);
    SWD_stC{i,:} = SWD_st;
    R_SWD_stC{i,:} = R_SWD_st(id_st);
    if switches.totalN
        ZWD{i,:} = SWD_st./((vmfw+vmfh)./2);
    else
        ZWD{i,:} = SWD_st./vmfw;
    end
    R_ZWD{i,:} = R_SWD_st(sort(id_st))*1000;
    clear idstat stat id elevcheck idelev vmfh vmfw
end  
% Look for missing ZWD values
[ZWDnummiss,~] = cellfun(@size,ZWD);
id_missZWD = find(ZWDnummiss<max(ZWDnummiss));
% Delete station from processing if only one signal left
if any(ZWDnummiss<2)
    [ZWDnummiss,id_missZWD,lenStatout,A,SWD_obs,R_SWD,un_name,ZWD,A_stC,SWD_stC,~,R_ZWD] = deleteStatZWD(id_missZWD,ZWDnummiss,A,SWD_obs,R_SWD,un_name,ZWD,A_stC,SWD_stC,R_SWD_stC,R_ZWD,station_name,switches);
else
    lenStatout = 0;
end
% Fill missing signals with NaN values
if ~isempty(id_missZWD)
    for i = 1:size(id_missZWD,1)
        ZWD{id_missZWD(i),1}(ZWDnummiss(id_missZWD(i))+1:max(ZWDnummiss),1) = NaN; 
        R_ZWD{id_missZWD(i),1}(ZWDnummiss(id_missZWD(i))+1:max(ZWDnummiss),1) = NaN; 
    end
end
% Calculate coviariance matrix
ZWD = reshape(cell2mat(ZWD),max(ZWDnummiss),size(un_name,1));
cR_ZWD = full(reshape(cell2mat(R_ZWD),max(ZWDnummiss),size(un_name,1)));
cZWD = cov(ZWD,"omitrows");
cR_ZWD = cov(cR_ZWD,"omitrows");
% Copy covariance values to observation matrix size
[var_obs,var_obsR] = duplicateCVvalues(idLenGr,lenStatout,cZWD,cR_ZWD,A_stC,ZWD); 
if strcmp(switches.solution,'REAL')
    SWD_obs = SWD_obs_b;
    var_obs = sqrt(var_obs.^2+var_obsR.^2);
end
% Check the quality of the processing
if sum(ZWDnummiss) ~= size(var_obs,1)
    warning('weightObs:  Ground-based part of weighting matrix do not match in size to observations. Results might be corrupted') 
end
if any(isnan(var_obs(1,:)))
    warning('weightObs: Ground-based part of weighting matrix contains NaN values. Results might be corrupted') 
end
%% Calculate RO covariance matrix
% Calculate zenith wet delays in RO intresection points
try
    % Variance only solution
    if  ~isempty(coord) 
        for i = 1:size(coord,1)
            % Find RO-voxel intersection points
            OU = coord(i).OU;
            IN = coord(i).IN;
            if ~isempty(OU) && size(OU,3) == 1
                OU = repmat(OU,1,1,2);
                OU(1,:,2)= [NaN,NaN,NaN];
            end
            if size(OU,3) > 2
               OU = OU(:,:,1:2);
            end
            MI = [IN;OU];
            % Calculate refractivities in point's lat, lon, alt
            if ~isempty(MI)
                [altstack(i),Nwstack{i},zwderr{i}] = intersetionNw(MI,radii,model,epoch,switches);
            end
        end
        % Calculate Zenith Wet Delay
        height = diff(model.levels_TOMO)./1000;
        % Calculate Zenith Wet Delay
        try
            [ZWDROstack,~] = stackedNwcalc(Nwstack,height);
        catch
            warning(['weightObs: ROZWD calculation failed at epoch ',num2str(epoch)])
        end
        id_temp = 1:size(altstack,2);
        idoutRO = unique([find(altstack>model.levels_TOMO(end-1)),id_temp(isnan(altstack)),find(altstack<=0)]);
        ZWDROstack(idoutRO,:) = [];
        altstack(idoutRO) = [];
        zwderr(idoutRO) = [];
        A(idLenGr+idoutRO,:) = [];
        SWD_obs(idLenGr+idoutRO,:) = [];
        R_SWD(idLenGr+idoutRO,:) = [];
        idoutRO = idoutRO+size(var_obs,1);
        if ~isempty(ZWDROstack)
            [~,szwd] = cellfun(@size,zwderr);
            iderrzwd = find(szwd<max(szwd));
            for i = 1:size(iderrzwd,2)
                zwderr{1,iderrzwd(i),1}(szwd(iderrzwd (i))+1:max(szwd)) = NaN; 
            end
            zwderr = reshape(cell2mat(zwderr),max(szwd),size(zwderr,2))';
            ZWDROstackAlt = [altstack',ZWDROstack];
            id = find(ZWDROstackAlt(:,1)==0);
            ZWDROstackAlt(id,:) = [];
            A(idLenGr+id,:) = [];
            SWD_obs(idLenGr+id,:) = [];
            R_SWD(idLenGr+id,:) = [];
            for i = round(model.levels_TOMO(end-1)./4):-1:0
                id_lev = intersect(find(i*4<ZWDROstackAlt(:,1)),find(ZWDROstackAlt(:,1)<(i+1)*4));
                ZWDROlvl = ZWDROstackAlt(id_lev,:);
                zwderrlvl = zwderr(id_lev,:);
                if size(ZWDROlvl,1) > 1
                    zwderrlvl = zwderrlvl(:,~all(isnan(zwderrlvl)));
                    ZWDROprof = cov(ZWDROlvl(:,2:end)',"omitrows");
                    ZWDROproferr = cov(zwderrlvl(:,2:end)',"omitrows");
                    if sum(sum(ZWDROprof)) ~= 0
                        var_obs(size(var_obs,1)+1:size(var_obs,1)+size(ZWDROprof,1),size(var_obs,1)+1:size(var_obs,1)+size(ZWDROprof,1)) = sqrt(ZWDROprof.^2+ZWDROproferr.^2);
                    end
                end
            end
        end
    end
    if size(var_obs,1) ~= size(A,1) 
        %In case of only one observation in ZWD group with lowest altitude
        icut = size(A,1)-size(var_obs,1); 
        SWD_obs = SWD_obs(1:size(var_obs,1));
        R_SWD= R_SWD(1:size(var_obs,1));
        A = A(1:size(var_obs,1),:);
        if length(icut) > 1
            warning(['weightObs: Incorrect weighting matrix size. Check the amtrix data. Results might be corrupted on epoch ',num2str(epoch)]);
        end
    end
catch
    warning('weightObs: Failed to calulcate RO covariance matrix (weightObs)')
end   

 %% LSQ for outliers
 % Check matrix dimensions
n = 3;
Aout = 1:size(A_b,1);
Aout(idoutRO) = [];
if icut > 0
    Aout(end-icut+1:end) = [];
end
% Calculate inital weighting matrix
C = pinv(var_obs); 
C(isinf(C)|isnan(C)) = 0;
% Delete outliers
if ~isempty(SWD_obs)
    for i = 1:1
        x = pinv(A'*C*A)*A'*C*SWD_obs; 
        r = A*x-SWD_obs;
        try
            idGr =  find(abs(r(1:idLenGr))>n*std(r(1:idLenGr)));
        catch
            idGr =  find(abs(r(1:size(var_obs,1)))>n*std(r(1:size(var_obs,1))));
        end
        idRO = [];
        % Outliers check done separately for RO
        if size(r,1) > idLenGr
            idRO =  find(abs(r(idLenGr:end))>n*std(r(idLenGr:end)));
            idRO = idRO+idLenGr-1;
        end
        id = unique([idGr;idRO]);
        r(id) = [];
        C(id,:) = [];
        C(:,id) = [];
        A(id,:) = [];
        SWD_obs(id,:) = [];
        R_SWD(id,:) = [];
        var_obs(id,:) = [];
        Aout(id) = [];
    end
    it = 0;
    % Recalculate weighting matrix to m0
    while it < 3
        it = it + 1;
        try
            x = pinv(A'*C*A)*A'*C*SWD_obs;  
            r = A*x-SWD_obs;
            var_est = sqrt(abs((r'*diag(diag(C))*r)/(size(A,1)-size(A,2))));
            if and(var_est <1.1,var_est > 0.9)
                break
            else
                C =  1./var_est * C;
            end
        catch
            warning('weightObs: Failed to invert matrix. Results might be corrupted')
        end
    end
else
    C = [];
end
idcol = 1:size(A_b,1);
idcol(Aout) = [];
C = C;
end
