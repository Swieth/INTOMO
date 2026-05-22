function Q = covarianceAprRT(Nw_apr,A,num_lat,num_lon,levels,BLh_pudel,num_Nw,values,epoch,switches)

%function designed to get first estimate of dynamic disturbance in the
%wet refractivity based on provided Nw data (best working NWP forecasts)
%%%INPUT
%       Nw_apr.........apriori values of refractivities
%       A..............observational matrix with signal deriveratives
%       num_lat........number of voxels in latitudinal direction
%       num_lon........number of voxels in longitudinal direction
%       levels.........altitude layers of tomography model
%       BLh_pudel......coordinates of voxel edges/centers
%       num_Nw.........id numbers of voxels/nodes in inner model
%       values.........structural matrix containing inner/outer model data
%       epoch..........processing epoch number
%       switches.......structural matrix of settings of the INTOMO processing
%%%OUTPUT
%       Q..............  the first guess for covariance matrix
%---------------------------------------------------------------------------------------------------


if exist('switches','var') && strcmp(switches.regular,'no')
    Nw_NWP_irr = Nw_apr;
    try
        if strcmp(switches.apriori,'WRF')
            Nw_apr = values.Nw_WRF;
        elseif strcmp(switches.apriori,'DETER')
            Nw_apr = values.Nw_DETER;
        end
    end
end
if exist('nodes_BLh','var')
    if ~isempty(nodes_BLh)
        h = nodes_BLh(:,3);
    else
        h = [];
    end
else
    h = BLh_pudel(3,1,:);
    if strcmp(switches.integrated,'yes')
        h = h(:);
    end
end

clear kol wart
for i = 1: size(Nw_apr,1)
    NWPone =  reshape(Nw_apr(i,:),num_lat*num_lon,levels-1);
    NWPext(i,:) =  mean(NWPone);
end
modQ = cov(NWPext);
modQ = sqrt(diag(modQ));

if var(modQ) > 0.1  %in case we have real data from NWP
    modQ = repmat(modQ,1,num_lat*num_lon);
    modQ = reshape(modQ',num_lat*num_lon*(levels-1),1);
else
    % If GPT2 data used as aproiri, use statistics calculated based on
    % a comparison between GPT2 and WRF (lowcost case)
    deterconst = 1;
    if strcmp(switches.aprModel,'DETER')
        if isfield(values,'Nw_WRF') && deterconst==0
            NwWRF = reshape(values.Nw_WRF(epoch,:),num_lat*num_lon,levels-1);
            NwGPT = reshape(values.Nw_DETER(epoch,:),num_lat*num_lon,levels-1);
            dNw = NwGPT-NwWRF;
            Nwb = mean(dNw);
            Nws = std(dNw);
            modQ = sqrt(Nwb.^2+Nws.^2)';
        else
            href = [150, 275, 400, 575, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500,...
                2750, 3250, 3750, 4500, 5250, 6000, 6750, 7500, 8250, 10000, 11750];
            hlow = find(href<3000);
            modQref = [13.29, 12.35, 12.07, 11.63, 11.55, 11.03, 10.20, 8.98, 7.91, 6.92, 5.98, 5.17...
                4.52, 3.43, 2.79, 2.18, 1.75, 1.47, 1.19, 0.98, 0.77, 0.44, 0.22];
            modQ = interp1(href,modQref,h,'spline');
            plot_check = 0;
            if plot_check==1
                plot(modQref,href,'-o');
                hold on
                plot(modQ,h,'or');
            end
        end
        if ~exist('switches','var') || strcmp(switches.regular,'yes')
            modQl = repmat(modQ',num_lat*num_lon,1);
            modQ = reshape(modQl,num_lat*num_lon*(levels-1),1);
        end
        if switches.totalN
             range_intervals = 1:num_lat*num_lon:size(A,2);
             range_intervals = [range_intervals, size(A,2)];
             for lvl = 1:length(range_intervals)-1
                Pmean1(lvl,1) = abs(mean(values.Nw_apr_num(epoch,range_intervals(lvl):range_intervals(lvl+1)) - values.aprioriEra(epoch,range_intervals(lvl):range_intervals(lvl+1))));
                if lvl>5
                   % Pmean1(lvl,1) = Pmean1(lvl,1)./1000;    
                end
             end
             Pmean2 = repmat(Pmean1',num_lat*num_lon,1);
             modQ = reshape(Pmean2,num_lat*num_lon*(levels-1),1);
        end
    else
        %in case there are only models available so very limited change over time thus the empirical variance function introduced based on the Trzcina et al.,2016
        % in the folder /TOMO-LAB/PPROCES/Q_poly_fit shows fiting pice-wise
        % function bottom - from 0 to 3875 m linear from 3876 exponential
        % function is hard-coded into this software
        clear modQ
        p1 = [-0.0552 0.6214];
        p2 = [-0.7693 7.7347];
        h = h/1000;
        [lin_wie,kol_wie,~] = find(h<=3.875);
        modQ(1,lin_wie) = p1(1,1)*h(lin_wie,1) + p1(1,2);
        clear lin_wie kol_wie wart_wie
        [lin_wie,kol_wie,wart_wie] = find(h>3.875);
        modQ(1,lin_wie) = p2(1,2)*exp(p2(1,1)*h(lin_wie,1));
        clear lin_wie kol_wie wart_wie
        modQ = repmat(modQ,num_lat*num_lon,1);
        modQ = 9*reshape(modQ,num_lat*num_lon*(levels-1),1); % ???
        modQ = modQ*2;
    end
end
if size(modQ,1)<size(A,2)
    %outer model currently over-constrained
    modQ(size(modQ,1)+1:size(A,2),1) = 0.00018*10;
end
Q = modQ(num_Nw)*10;
