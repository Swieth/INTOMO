function N = IDW_atom( LAT, LON, Tlayer,Hlayer,Player,Nlayer,Nwlayer ,rWGS, llRay,hlayers,set)
    % IDW_atom provides IDW interpolation on ray point coordiantes to 
    % calculate refraction values from temperature, pressure and water
    % vapour pressure
    
    %% Input data
    % LAT.............. ellipsoidal latitude of refractivity [radians]
    % LON.............. ellipsoidal longitude of refractivity [radians] 
    % Tlayer........... temperature, matrix containing all levels
    % Hlayer........... water vapour pressure, matrix containing all levels
    % Player........... pressure, matrix containing all levels
    % rWGS............. Earth radius matrix
    % llRay............ ray point coordinates
    % hlayers.......... vertical layers heights
    % set.............. settings 
    %% Output data
    % N
    %       Nw......... wet refractivity at ray point
    %       Nh......... hydrostatic refractivity at ray point
    %       Nt......... total refractivity at ray point
    %       T.......... temperature  at ray point [K]
    %       p.......... pressure  at ray point [hPa]
    %       h.......... total refractivity at ray point [kg kg-1]
    
    % Constants   
    Rd = 287.058;       % gas constant dry air [J/K/kg]
    Md = 28.9644;       % molar mass dry air [kg/mol]
    Mw = 18.0152;       % molar mass wet air [kg/mol]
    g0 = 9.80665;       % WMO gravity (optional: 9.7803253359)
    % Find distances to nodes
    AD = sqrt((pi()/2 - llRay(1,1) - LAT).^2 + (llRay(1,2) - LON).^2);
    % Find one lowest value in each column and sort array in ascending order
    sdmin = sort(min(AD), 'ascend');
    if set.refron
          k = 1;
        for i = 1:2
            % Find column number for 2 lowest values
            [~, col] = find(AD == sdmin(1,i));
            % Error found: can result in 2 equal values but only 1 can be used
            if numel(col) > 1
                col = col(1,1);
            end
            % Give value and row for found column
            [val, row] = sort(AD(:,col), 'ascend');
            for j = 1:2
                % Create cell with location and value for wrf grid
                node{i,j} = [row(j,1); col; val(j,1)];
                Nprf(:,k) = squeeze(Nlayer(row(j,1),col,:));
                Wprf(:,k) = squeeze(Nwlayer(row(j,1),col,:));
                N.coord(k,:) = [row(j,1),col];
                k = k + 1;  
            end
        end
        % Array with 4 nearest nodes: wrf row (1), wrf column (2), distance to ray (3)
        nodeh2 = cell2mat(reshape(node,1,4));
        nodeh2(3,1);
        % Weights: inverse squared distance
        wNode = nodeh2(3,:).^(-2);   
        geomRay = llRay(1,3);  
        % interpolate refractivity at station position
        % Find layer above and below current ray altitude
        dgeom = geomRay - hlayers;
        % Layer below
        valb = find(dgeom >= 0);
        % Layer above
        vala = find(dgeom < 0 ); 
        nodev = [max(valb); min(vala)];
        %% ERA REFARCTION VALUES FOR NODES
        % Loop through nodes to read data for 2 nearest levels and 4 nodes
        for k = 1:4
            % DIM: west - east (x)/ south - north (y)/ bottom - top / time
            for i = 1:size(nodev,1)
                % pressure [hPa]
                Nprf2(i,k) = Nprf(nodev(i),k);  
                Wprf2(i,k) = Wprf(nodev(i),k); 
            end
        end   
        try
            N.Nt = sum((Nprf2(1,:)+Nprf2(2,:))./2.*wNode)./sum(wNode);  
        catch
            N.Nt = sum((Nprf2(1,:).*wNode))./sum(wNode);  
        end
        try
            N.Nw = sum((Wprf2(1,:)+Wprf2(2,:))./2.*wNode)./sum(wNode); 
       catch
            N.Nw = sum((Wprf2(1,:).*wNode))./sum(wNode);  
        end

    else
        k = 1;
        for i = 1:2
            % Find column number for 2 lowest values
            [~, col] = find(AD == sdmin(1,i));
            % Error found: can result in 2 equal values but only 1 can be used
            if numel(col) > 1
                col = col(1,1);
            end
            % Give value and row for found column
            [val, row] = sort(AD(:,col), 'ascend');
            for j = 1:2
                % Create cell with location and value for wrf grid
                node{i,j} = [row(j,1); col; val(j,1)];
                Tprf(:,k) = squeeze(Tlayer(row(j,1),col,:));
                Hprf(:,k) = squeeze(Hlayer(row(j,1),col,:));
                Pprf(:,k) = squeeze(Player(row(j,1),col,:));
                rWGprf(:,k) = squeeze(rWGS(row(j,1),col));
                N.coord(k,:) = [row(j,1),col];
                k = k + 1;  
            end
        end

        % Array with 4 nearest nodes: wrf row (1), wrf column (2), distance to ray (3)
        nodeh2 = cell2mat(reshape(node,1,4));
        nodeh2(3,1);
        % Weights: inverse squared distance
        wNode = nodeh2(3,:).^(-2);   
        geomRay = llRay(1,3)+20;  
        % interpolate refractivity at station position
        % Find layer above and below current ray altitude
        dgeom = geomRay - hlayers;
        % Layer below
        valb = find(dgeom >= 0);
        % Layer above
        vala = find(dgeom < 0 ); 
        nodev = [max(valb); min(vala)];
        %% ERA METEOROLOGICAL PARAMETERS FOR NODES
        % Loop through nodes to read data for 2 nearest levels and 4 nodes
        for k = 1:4
            % DIM: west - east (x)/ south - north (y)/ bottom - top / time
            for i = 1:size(nodev,1)
                % pressure [hPa]
                Tprf2(i,k) = Tprf(nodev(i),k); 
                Hprf2(i,k) = Hprf(nodev(i),k);  
                Pprf2(i,k) = Pprf(nodev(i),k);  
            end
        end    
        if geomRay < 181
            %% interpolate parameters at station height
            % water vapor pressure
            Ce = (hlayers(nodev(2)) - hlayers(nodev(1)))./log(Hprf2(2,:)./Hprf2(1,:));
            % Can results in NaN
            eh1 = Hprf2(1,:).*exp((geomRay - hlayers(nodev(1)))./Ce);
            eh2 = Hprf2(2,:).*exp((geomRay - hlayers(nodev(2)))./Ce);
            eh(1,:) = sum((eh1 + eh2)./2.*wNode)./sum(wNode);
            if isnan(eh) && sum(Hprf2(2,:)) == 0 || sum(Hprf2(1,:)) == 0 
                eh = 0;
            end
            if isnan(eh)
                disp(eh)
            end
            % temperature
            Ctemp = (Tprf2(2,:) - Tprf2(1,:))./(hlayers(nodev(2)) - hlayers(nodev(1)));
            % Linear interpolation of temperature
            temph1 = Tprf2(1,:) + (geomRay - hlayers(nodev(1))).*Ctemp;
            temph2 = Tprf2(2,:) + (geomRay - hlayers(nodev(2))).*Ctemp;
            temph(1,:) = sum((temph1 + temph2)./2.*wNode)./sum(wNode);
            % pressure

            for i = 1:size(nodev,1)
                % Acceleration due to gravity (Kraus, 2004): Die atmosphere der Erde
                gm(i,:) = g0*(1 - 0.0026373*cos(2*llRay(1,2)*pi()/180) + 5.9*1e-6*cos(2*llRay(1,2)*pi()/180)).^2*(1./(1 + (geomRay/1000)./rWGprf).^2);

                % Virtual temperature [K]
                Tv(i,:) = Tprf2(i,:).*Pprf2(i,:)./(Pprf2(i,:) - (1 - Mw/Md)*Hprf2(i,:));

                % Extrapolate pressure from grid's layer
                presh(i,:) = Pprf2(i,:).*exp(-gm(i,:).*(geomRay - hlayers(nodev(i)))./(Rd.*Tv(i,:)));
            end
            presh = sum(mean(presh).*wNode)./sum(wNode);
        else
            try
                presh = sum((Pprf2(1,:)+Pprf2(2,:))./2.*wNode)./sum(wNode);
                temph = sum((Tprf2(1,:)+Tprf2(2,:))./2.*wNode)./sum(wNode);
                eh = sum((Hprf2(1,:)+Hprf2(2,:))./2.*wNode)./sum(wNode);
            catch
                presh = sum((Pprf2(1,:)+Pprf2(1,:))./2.*wNode)./sum(wNode);
                temph = sum((Tprf2(1,:)+Tprf2(1,:))./2.*wNode)./sum(wNode);
                eh = sum((Hprf2(1,:)+Hprf2(1,:))./2.*wNode)./sum(wNode);
            end
        end   
        [ Nh, Nw] = refcalc(presh, temph, eh, 'c'); 
        N.Nh = Nh;
        N.Nw = Nw;
        N.Nt = Nh+Nw;
        N.T = temph;
        N.p = presh;
        N.h = eh;
    end
end