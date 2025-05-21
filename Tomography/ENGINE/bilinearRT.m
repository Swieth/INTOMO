function [A,SWD_nodes_integ,nodesD] = bilinearRT(model,A,elevat,azi,nodes_columns,levels,num_lev,Tmatrix,Nw_apr,epoch,ray_num,coord,switches)                    
% Bilinear/spline interpolation function accroding to 
% Perler, Donat, Alain Geiger, and Fabian Hurter. "4D GPS water vapor tomography: 
% new parameterized approaches." Journal of Geodesy 85 (2011): 539-550; and 
% Trzcina, Estera, Witold Rohm, and Kamil Smolak. "Parameterisation of the GNSS 
% troposphere tomography domain with optimisation of the nodes’ distribution." 
% Journal of Geodesy 97.1 (2023): 2.
%%%INPUT
%       model......   parameters of tomography and ray tracing models
%       A..........   observational matrix with signal deriveratives
%       elevat.....   elevation angle of satellite
%       azi........   azimuth of the satellite
%       nodes_columns segragated id's of voxels for use in bilin interp
%       levels.....   altitude layers of tomography model
%       numa_lev..    number of voxels in vertical direction
%       Tmatrix....   second derivates of refractivity for use in bilin interp
%       Nw_apr.....   apriori values of refractivities
%       epoch......   number of processing epoch
%       ray_num....   number of row in A matrix
%       coord......   coord. of interesetcions between RO and voxel model
%       swtiches...   structural matrix of settings of the INTOMO processing
%%%OUTPUT:
%       A..........   parameters of tomography and ray tracing models
%       SWD_nodes_integ    SWD integral based on a priori data
%       nodesD....    id of nodes assigned to the voxels intersected by signal
%---------------------------------------------------------------------------------------------------

    if ~isnan(coord) % If the ray intersects with the domain's grid
        if strcmp(switches.parametrization,'nodes') || strcmp(switches.parametrization,'bilinear-h')           
            %% Loop for each section (each layer; layer is a space between two consecutive horizontal planes of the domain)
            p_coord_sect = [];
            Np_sect = [];
            A_rows=[];
            p_Nw = [];
            % Number of subsections and subsubsections
            s_sect_nr=4;
            ss_sect_nr=4;
            mult_fact = [7 32 12 32 14 32 12 32 14 32 12 32 14 32 12 32 7];
            SWD_sum_all=0;
            for sect=1:size(coord,1)-1
                %% Get coordinates of the consecutive points
                [p_coord,p_coordBLH,const] = get_p_coordRT(coord,sect,s_sect_nr,ss_sect_nr,azi,elevat);
                clear Np
                %% Calculate integrals of interpolation
                for p = 1:size(p_coord,1)
                    [SWD_sum_all,A_rows,SWD_p(p,1),p_Nw,nodesC] = A_row_bilin_regRT(model,p,Nw_apr,epoch,nodes_columns,levels,p_coordBLH,num_lev,Tmatrix,mult_fact,const,SWD_sum_all,A_rows,A,p_Nw);
                    nodesDa(p,:) = nodesC(end,:);
                end % end of all p_coord in section (end of section)
                nodesD(sect).nodes = nodesDa;
                SWD_sect(sect,1) = sum(SWD_p);
                ray(ray_num).Nw(sect).Nw = p_Nw;
                ray(ray_num).p_coord(sect).coord = p_coord;
                len(ray_num).par = coord;
            end % end of all sections (end of one slant observation)
            SWD_nodes_integ = SWD_sum_all;
            %% Assign integrals to A matrix 
            A_row_true = sum(A_rows);
            A(ray_num,:) = A_row_true;
        end
        clear SWD_sect SWD2_sect
    end
                  