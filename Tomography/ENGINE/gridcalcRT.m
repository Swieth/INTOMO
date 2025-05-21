function [ pGrid, eGrid, tempGrid, radWGS] = gridcalcRT( XLAT,XLON, undugpt, pres, temp, specHum, geop, layers,lat,lon,struc,model)

    %% ABOUT
    %       function for whole ERA5 grid meteorological parameters calculation
    %       US standard atmosphere is applied above 20 km of altitude to provide
    %       temperature. Pressure is extrapolated to from model top layer. 
    %
    %% Input 
    %
    %           XLAT        -   Latitude array [rad]
    %           undugpt     -   Undulation grid calculated in advance
    %           pres        -   Presure grid in [hPa]
    %           temp        -   Temperature grid in [K]
    %           specHum     -   Specific humidity grid in [kg/kg]
    %           geop        -   geopotential heights
    %           layers      -   spacing, layers at which interpolate parameters
    %           lat         -   latitiudes of processed domain
    %           lon         -   longotudes of processed domain    
    %           struc       -   variable to modify undulation grid
    %           model       -   parameters of tomography and ray tracing models
    %% Output
    %           pGrid      -       pressure grid
    %           eGrid      -       water vapour grid
    %           tempGrid   -       temperature grid

    %% METEOROLOGICAL CONSTANTS
    Rd = 287.058;       % gas constant dry air [J/K/kg]
    Md = 28.9644;       % molar mass dry air [kg/mol]
    Mw = 18.0152;       % molar mass wet air [kg/mol]
    g0 = 9.80665;       % WMO gravity (optional: 9.7803253359)

    %% ELLIPSOID PARAMETERS
    a = model.radii(1);               % semi-major axis [km]
    b = model.radii(2);           % semi-minor axis [km]
    e2 = (a^2-b^2)/a^2;         % power eccentricy (e^2)
    if size(XLAT,2) ~= size(lat,2)
        [~,id] = intersect(round(XLAT(1,:)*180/pi(),2),lat);
        pres =  flip(pres(:,id,:),2);
        specHum =  flip(specHum(:,id,:),2);
        temp =  flip(temp(:,id,:),2);
        geop = flip(geop(:,id,:),2);
        XLAT = flip(XLAT(:,id),2);
    end
    if size(XLAT,1) ~= size(lon,2)
        [~,id] = intersect(round(XLON(:,1)*180/pi(),2),lon);
        geop = geop(id,:,:);    
        pres =  pres(id,:,:);
        specHum =  specHum(id,:,:);
        temp =  temp(id,:,:);
        XLAT = XLAT(id,:);
    end
    if ~isempty(struc)
        if size(struc.LAT,2) ~= size(lat,2)
            [~,id] = intersect(round(struc.LAT(1,:)*180/pi(),2),lat);
            undugpt =  flip(undugpt(:,id,:),2);
        end
        if size(struc.LAT,1) ~= size(lon,2)
            [~,id] = intersect(round(struc.LON(:,1)*180/pi(),2),lon);
            undugpt = undugpt(id,:,:);    
        end
    end

    % convert to ellipsoidal heights
    undu = repmat(undugpt,1,1,size(geop,3))./1000;
    XLAT2 = repmat(XLAT,1,1,size(geop,3));
    geomPrf = fgeop2geom(XLAT2, geop) + undu;


    % Water vapor partial pressure profile [hPa]
    eg = specHum.*pres./(Mw/Md + (1 - Mw/Md).*specHum);  % from specific humidity

    % find maximum maximum height - above this layer all nodes have standard
    % atmosphere
    lastH = max(max(geomPrf(:,:,end)));
    % find heights above maximum model height
    satmH = layers(layers>lastH);
    % find heights below maximum model height
    intH = layers(layers<=lastH);
    % 3D array with all heights
    satmHa = zeros(size(XLAT,1),size(XLAT,2),length(layers));
    for ih=1:length(layers)
        satmHa(:,:,ih) = repmat(layers(ih),size(XLAT,1),size(XLAT,2),1);
    end

    % Convert to geopotential height [km]
    XLATa = repmat(XLAT,1,1,size(satmH,2));
    geopUS76 = fgeom2geop(XLATa, satmHa(:,:,length(intH)+1:end));

    % US standard atmosphere (US76) for temperature only (geopotential height)
    [~, ~, tempRay, ~, ~, ~] = stdatmo(geopUS76*1000);

    % Earth radius
    Rns = a^2*b^2.*(a^2.*cos(XLAT).^2 + b^2.*sin(XLAT).^2).^(-3/2);
    Rew = a./((1-e2.*sin(XLAT).^2).^(1/2));

    % latitude array for acceleration
    XLATacc = repmat(XLAT,1,1,size(layers,2));
    Rwgs = repmat(sqrt(Rns.*Rew),1,1,size(layers,2));

    % Acceleration due to gravity (Kraus, 2004): Die atmosphere der Erde
    gm = g0.*(1 - 0.0026373.*cos(2.*XLATacc) + 5.9*1e-6.*cos(2.*XLATacc)).^2.*(1./(1 + satmHa./Rwgs).^2);

    % Extrapolate pressure from model's top
    presh = repmat(pres(:,:,end),1,1,size(satmH,2)).*exp(-gm(:,:,length(intH)+1:end).*(satmHa(:,:,length(intH)+1:end) - repmat(geomPrf(:,:,end),1,1,size(satmH,2))).*1000./(Rd.*repmat(temp(:,:,end),1,1,size(satmH,2))));
    presh(presh<0) = 0;

    % specifiv cloud water and ice content, water vapor partial pressure above model = 0 
    %[ eAb, sLWCAb, sIWCAb, sRWCAb, sSWCAb] =  deal(zeros(size(XLAT,1),size(XLAT,2),size(satmH,2)));
     eAb =  deal(zeros(size(XLAT,1),size(XLAT,2),size(satmH,2)));

    % initialize zeros array for parameters below max model height
    %[ tempW, sLWCW, sIWCW, sRWCW, sSWCW, presW, eW] = deal(zeros(size(XLAT,1),size(XLAT,2),size(intH,2)));

    [ tempW, presW, eW] = deal(zeros(size(XLAT,1),size(XLAT,2),size(intH,2)));

    % interpolate parameters - each profile individually
    for i=1:size(XLAT,1)
        for j=1:size(XLAT,2)
            clearvars tempbt tempab Hab
            [tbl, tempab] = deal([]);
            % profiles of meteorological paramateres for given  lon/lat node
            hPrf = squeeze(geomPrf(i,j,:));
            tdsf(i,j) = hPrf(1);
            tPrf = squeeze(temp(i,j,:));
            pPrf = squeeze(pres(i,j,:));
            ePrf = squeeze(eg(i,j,:));
            gPrf = squeeze(gm(i,j,:));

            % heights inside ERA-I layers
            Hbet = intH(intH>=hPrf(1,:) & intH<=hPrf(end,:));
            % below
            Hbl = intH(intH<hPrf(1,:));
            % above
            hr = find(intH>hPrf(end,:));
            % find layers where want to interpolate below model
            indB = find(layers<=hPrf(1,:));
            % if there are RT heights below model
            if ~isempty(indB)
                % Virtual temperature [K]
                Tv = tPrf(1).*pPrf(1)./(pPrf(1) - (1 - Mw/Md)*ePrf(1));
                % Extrapolate pressure from model's 1st layer
                presW(i,j,indB) = pPrf(1).*exp(-gPrf(indB).*(layers(indB)' - hPrf(1))*1000./(Rd.*Tv));
                % Keep water vapor pressure unchanged;
                eW(i,j,indB) = repmat(ePrf(1),length(indB),1);
                % interpolate temperature below
                tbl = tPrf(1,:) + (hPrf(1,:)-Hbl).*6.5;

            end

            % RT heights inside model
            % linear interpolation of temperature
            tempbt = interp1(hPrf,tPrf,Hbet,'linear');
            Cea = (hPrf(2:end) - hPrf(1:end-1))*1000./log(ePrf(2:end)./ePrf(1:end-1));
            % Virtual temperature [K]
            Tv = tPrf.*pPrf./(pPrf - (1 - Mw/Md)*ePrf);
            for hw=1:length(hPrf)-1
                indW = find(layers>hPrf(hw) & layers<=hPrf(hw+1));
                % Extrapolate pressure from model's top
                presw1 = pPrf(hw,:).*exp(-gPrf(indW).*(layers(indW)' - hPrf(hw,:))*1000./(Rd.*Tv(hw,:)));
                presw2 = pPrf(hw+1,:).*exp(-gPrf(indW).*(layers(indW)' - hPrf(hw+1,:))*1000./(Rd.*Tv(hw+1,:)));
                presw1(presw1<0) = 0;
                presw2(presw2<0) = 0;
                presw1(isnan(presw1)) = 0;
                presw2(isnan(presw2)) = 0;
                presW(i,j,indW) = (presw1 + presw2)./2;

                % Can results in NaN
                eh1 = ePrf(hw,:).*exp((layers(indW)' - hPrf(hw,:))*1000./Cea(hw));
                eh2 = ePrf(hw+1,:).*exp((layers(indW)' - hPrf(hw+1,:))*1000./Cea(hw));
                eh1(eh1<0) = 0;
                eh2(eh2<0) = 0;

                 if ~isreal((eh1 + eh2)./2)
                     eW(i,j,indW) = 0;
                else
                    eW(i,j,indW) = (eh1 + eh2)./2;
                 end
            end

            % RT heights above top model's layer
            if ~isnan(hr)
                Hab = intH(hr);
                % Convert to geopotential height layers above [km]
                XLATa = repmat(XLAT(i,j), 1, size(Hab,2));
                geopAbove = fgeom2geop(XLATa, Hab);
                % US standard atmosphere (US76) for temperature above only (geopotential height)
                [~, ~, tempab, ~, ~, ~] = stdatmo(geopAbove*1000);
                % extrapolate pressure above
                presW(i, j, hr) = repmat(pPrf(end,:),length(hr),1).*exp(-gPrf(hr,:).*(Hab' - repmat(hPrf(end,:),length(hr),1)).*1000./(Rd.*repmat(tPrf(end,:),length(hr),1)));
                % calculate temperature corrections between last model layer
                % and the first RT height above the model
                % calculate temperature correction
                Htemp = hPrf(end,:);
                geopTemp = fgeom2geop(XLAT(i,j), Htemp);
                [~, ~, tempTemp, ~, ~, ~] = stdatmo(geopTemp*1000);
                corr = tPrf(end) - tempTemp;
                tempab = tempab + corr;
            end  
             tempW(i,j,:) = [tbl'; tempbt'; tempab'];     
        end
    end

    for w=1:size(eW,3)
        % find nan values
        [row, col, ~] = find(isnan(reshape(eW(:,:,w),size(XLAT,1),size(XLAT,2))));
        if ~isempty(row)
            for l=1:length(row)
                rmin = row - 1;
                rmax = row + 1;
                cmin = col - 1;
                cmax = col + 1;
                if row == 1
                   rmin = 2; 
                end
                if row == size(eW,1)
                    rmax = row - 1;
                end
                if col == size(eW,2)
                    cmax = col - 1;
                end
                if col == 1
                    cmin = 2;
                end
                dsr = [eW(rmin,col,w);eW(rmax,col,w);eW(row,cmin,w);eW(row,cmax,w)];
                dsr(isnan(dsr)) = [];
                eW(row,col,w) = mean(dsr);
            end
        end
    end
    presW(presW<0) = 0;
    % create additional layer for satellite param
    zeroGrid = zeros(size(XLAT,1), size(XLAT,2), 1);
    onesGrid = ones(size(XLAT,1), size(XLAT,2), 1);
    tempGrid = cat(3, tempW, tempRay, onesGrid);
    eGrid = cat(3, eW, eAb, zeroGrid);
    pGrid = cat(3, presW, presh, zeroGrid);
    radWGS = sqrt(Rns.*Rew);


end

