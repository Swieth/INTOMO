function [A,SWD,SIWV,R_SWD,R_SIWV,not_cro,elev,SAT,station_name,coord,dSWD] = initialOBS(station,epoch,rem_SAT,rem_REC,path_save, switches, Nw_apr,model)
%Function used to build basic matrices using information from
%station and satellite coordinates in integrated and non intregrated mode
%---------------------------------------------------------------------------------------------------
%%%INPUT
%       station....   structural matrix containing all informations regarding tomography model 
%       epoch......   processing epoch number
%       rem_SAT....   id of satellites to delete from apriori data
%       rem_REC....   id of receivers to delete from apriori data
%       path_save..   path to save the output matrices
%       swtiches...   structural matrix of settings of the INTOMO processing
%       Nw_apr.....   apriori values of refractivities
%       model......   parameters of tomography and ray tracing models
%%%OUTPUT:
%       A..........   observational matrix with signal deriveratives
%       SWD........   the values of Slant Delay (SD) or excess phase (EP)
%       SIWV.......   the values of Slant Integrated Water Vapour
%       R_SWD......   the errors of SD or EP
%       R_SIWV.....   the errors of Slant Integrated Water Vapour
%       not_cro....   voxels not crossed by ray path
%       elev.......   elevation angles of ray traced signals
%       SAT........   id of the satellites
%       station_name  names of GNSS stations
%       coord......   coord. of interesetcions between RO and voxel model
%       dSWD.......   difference between observed and simulated slant delay

%% Check if file exists
save_filename = [path_save switches.amtrix '/amtrix_' num2str(epoch) '.mat'];
if ~exist(save_filename,'file')
    %% Create empty A matrix for GNSS and A_RO for RO
    assum_sat = 32;
    stac = size(station.h,2);
    if strcmp(switches.regular,'yes') && strcmp(switches.parametrization,'constant')
        A((model.num_lat_TOMO-1)*(model.num_lon_TOMO-1)*(model.num_levels_TOMO-1),assum_sat*stac) = 0;
        Avec((model.num_lat_TOMO-1)*(model.num_lon_TOMO-1)*(model.num_levels_TOMO-1),assum_sat*stac,3) = 0;
        A_RO((model.num_lat_TOMO-1)*(model.num_lon_TOMO-1)*(model.num_levels_TOMO-1),1000) = 0;  
        AvecRO((model.num_lat_TOMO-1)*(model.num_lon_TOMO-1)*(model.num_levels_TOMO-1),1000,3) = 0; 
        exPh = zeros(1000,1);
    elseif strcmp(switches.regular,'yes') && strcmp(switches.parametrization,'bilinear-h')
        A_nodes(assum_sat*stac,model.num_lat_TOMO*model.num_lon_TOMO*model.num_levels_TOMO) = 0;
        Avec=[];
        A_RO(1000,model.num_lat_TOMO*model.num_lon_TOMO*model.num_levels_TOMO) = 0;
        exPh = zeros(1000,1);
        AvecRO = [];
    end
    %% Create matrices used for the interpolation in bilinear approach
    if strcmp(switches.parametrization,'bilinear-h') 
       % Get heights of the tomographic layers
        planes = model.levels_TOMO';
        planes_nr = size(planes,1);
        interp_h='Perler'; %other options may be added
        if strcmp(interp_h,'Perler')
            [Tmatrix,nodes_columns,~]  = perler_T( planes,planes_nr,model.BLh_pudel_rad,Nw_apr,epoch);
        end
    end
    %% Create empty vectors and matrices to fill in later
    SWD(assum_sat*stac,:) = 0;
    SIWV(assum_sat*stac,:) = 0;
    R_SIWV(assum_sat*stac,:) = 0;
    R_SWD(assum_sat*stac,:) = 0;
    SAT(assum_sat*stac,:) = 0;
    REC(assum_sat*stac,:) = 0;
    elev(assum_sat*stac,:) = 0;
    station_name = cell(assum_sat*stac,1);
    %% Calculate excess phase and fill A_RO matrix
    if strcmp(switches.integrated,'yes') %%yes
        if strcmp(switches.parametrization,'constant')
           [A_RO, exPh,~,R_SWD,dexPh,coord,AvecRO] = spaceRT(model,A_RO,AvecRO,station,epoch,0,0,0,0,0,switches); 
        else 
           [A_RO, exPh,~,R_SWD,dexPh,coord,AvecRO] = spaceRT(model,A_RO,AvecRO,station,epoch,nodes_columns,planes,planes_nr,Tmatrix,Nw_apr,switches); 
        end
    else
        A_RO = [];
        exPh = [];
        RelDistRO = [];
        R_SWD = [];
        dexPh = [];
        coord = [];
    end  
   %% Calculate slant delays and fill A matrix
   sat = 0;
   b = 0;
   try
       for nr = 1:stac
            sat_number = length(station.h(nr).satellite);
            while sat < sat_number
                sat = sat + 1;
                b = b + 1;
                disp(['Slant number: ' num2str(b,'%03d') '/' num2str(sat_number*stac,'%4d') '      ' datestr(datetime)  ' Station: ' num2str(nr) ' Satellite: ' num2str(sat) ]);
                BLhst = station.h(nr).parameters;
                X_st = cspice_georec(BLhst(3)*cspice_rpd,BLhst(2)*cspice_rpd,BLhst(4)/1000,model.radii(1),(model.radii(1)-model.radii(2))/model.radii(1));
                elevat = station.h(nr).satellite(sat).elevation;
                satcoord = station.h(nr).satellite(sat).coord;
                % Simulate data
                if strcmp(switches.parametrization,'constant')
                    [A, SWD(b,1),SIWV(b,1),REC(b,1),elev(b,1),SAT(b,1),station_name(b),R_SIWV(b,1),R_SWD(b,1),SWD_nodes_integ,RelDist(b,1),dSWD(b,1),Avec] = groundRT(A,Avec,model,station,epoch,X_st,satcoord,elevat,b,nr,sat,0,0,0,0,0,switches); 
                elseif strcmp(switches.parametrization,'bilinear-h')
                    [A_nodes, SWD(b,1),SIWV(b,1),REC(b,1),elev(b,1),SAT(b,1),station_name(b),R_SIWV(b,1),R_SWD(b,1),SWD_nodes_integ,RelDist(b,1),dSWD(b,1),Avec] = groundRT(A_nodes,Avec,model,station,epoch,X_st,satcoord,elevat,b,nr,sat,nodes_columns,planes,planes_nr,Tmatrix,Nw_apr,switches); 
                end  
            end
            sat = 0;
       end
   catch
       warning('initialOBS:rtFailed', 'initial OBS: Failed to calculate RT for station number: %s', num2str(nr))
   end    
   if strcmp(switches.parametrization,'bilinear-h')
        A = A_nodes';
        A_RO = A_RO';
   end
    %% Join A and A_RO matrices
     elev(size(A,2)+1:size(A,2)+size(exPh,1),1) = 0;
     SAT(size(A,2)+1:size(A,2)+size(exPh,1),1) = 100; 
     REC(size(A,2)+1:size(A,2)+size(exPh,1),1) = 100;
     R_SWD(size(A,2)+1:size(A,2)+size(exPh,1),1) = 0;
     R_SIWC(size(A,2)+1:size(A,2)+size(exPh,1),1) = 0;
     RelDist(size(A,2)+1:size(A,2)+size(exPh,1),1) = 0;
     station_name(size(A,2)+1:size(A,2)+size(exPh,1),1) = {'RO'};
     SWD = [SWD;exPh]*1000;
     A = [A,A_RO];
     Avec = [Avec,AvecRO];
     dSWD = [dSWD;dexPh];
     A = A(:,1:size(elev,1));
     if strcmp(switches.parametrization,'constant')
        Avec = Avec(:,1:size(elev,1),:);
     end         
    %% Filter out zero and NaN values from matrices
    if switches.run_SWD
        [wi1,~,~] = find(SWD==0);
        [wi2,~,~] = find(isnan(SWD)==1);
    else
        [wi1,~,~] = find(SIWV==0);
        [wi2,~,~] = find(isnan(SIWV)==1);
    end
    wi=union(wi1,wi2);
    SWD(wi,:) = [];
    SIWV(wi,:) =[];
    A(:,wi) =[];
    if strcmp(switches.parametrization,'constant')
        Avec(:,wi,:) =[];
    end
    SAT(wi,:) = [];
    REC(wi,:) = [];
    elev(wi,:) = [];
    dSWD(wi,:) = [];
    RelDist(wi,:) = [];
    station_name(wi,:) = [];
    R_SWD(wi,:) =[];
    R_SIWV(wi,:) =[];
    clear wi wi1 wi2
    A = A';
    [~,not_cro] = find(sum(A)==0);
    if strcmp(switches.parametrization,'constant')
        Avec = permute(Avec,[2,1,3]);
    end
    not_cro = not_cro';
    clear wi ki wa
    %% Remove unwanted observations from receivers or satellites 
    wi1 = [];
    if ~isempty(rem_SAT)
        for i = 1:size(rem_SAT,1)
            [wik,~,~] = find((SAT)==rem_SAT(i,1));
            wi1 = [wi1; wik];
            clear wik ki wa
        end
    end
    wi2 = [];
    if ~isempty(rem_REC)
        for i = 1:size(rem_REC,1)
            [wik,~,~] = find((REC)==rem_REC(i,1));
            wi2 =[wi2; wik];
            clear wik ki wa
        end
    end
    wi = union(wi1,wi2);
    SWD(wi,:) = [];
    SIWV(wi,:) = [];
    A(:,wi) =[];
    if strcmp(switches.parametrization,'constant')
        Avec(:,wi,:) =[];
    end
    SAT(wi,:) = [];
    REC(wi,:) = [];
    elev(wi,:) = [];
    station_name(wi,:) = [];
    dSWD(wi,:) = [];
    RelDist(wi,:) = [];
    R_SWD(wi,:) =[];
    R_SIWV(wi,:) = [];
    clear wi ki wa wi1 wi2
    SWD_nodes_integ = A*Nw_apr(epoch,:)';
    %% Save matrices
    save(save_filename,'A','Avec','SWD','SIWV','R_SWD','R_SIWV','not_cro','elev','SAT','SWD_nodes_integ','station_name','RelDist','dSWD','coord','-v7.3');
else
    %% Load matrices
    load(save_filename);
    if ~exist('coord','var')
        coord = [];
    end
end

