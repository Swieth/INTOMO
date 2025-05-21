function [station] = construct_station_LAB(BLh,BLH,ZTD,M_ZTD,ZHD,NAME,PRN,SP3X,SP3Y,SP3Z,observation_set,cut_off_angle,ro,switches,varargin)
% Function to convert the specific ZTD data into the station structural matrix an input for tomography model.
%-------------------------------------------------------------------------------------
%%% OUTPUT
% station(t).                                          : t-number of epochs
%            h(nr).                                    : nr-number of stations
%                  ZTD                                 : initial ZTD, gradients and conversion factor if present [ZTD GN GE Q]
%                  parameters                          : id and coordinates of GNSS receiver [id lat lon h H]
%                  name                                : station 4 letters abbrev
%                  satellite(sat)                      : sat - satellite number currently 1:32
%                                .PRN                  : PRN number
%                                .elevation            : satellite elevation [deg]
%                                .SWD                  : slant wet delay
%                                .M_SWD                : SWD error
%                                .STD                  : slant total delay
%                                .M_STD                : STD error
%                                .SIVW                 : slant integrated water vapour
%                                .M_SIVW               : SIWV error
%                                .vmf1w                : Vienna Mapping Function wet
%                                .vmf1h                : Vienna Mapping Function hydrtostatic
%                                .coord                : satellite position (ECEF)
%            ro.                                       : radio occultation (RO) data
%                                .coordT               : ccordinates of the transmitter (ECEF)
%                                .coordR               : ccordinates of the receiver (ECEF)
%                                .date                 : date of RO start
%                                .exL2                 : RO excess phase on L2
%                                .exLC                 : RO excess phase on LC                             

% The are number of data that can be stored in this matrix but :
% 1. In case of non integrated solution or lack of RO data, the RO
% structure part remains empty
% 2. SIWV/SWD struture parts may remain empty depending on choosen solution
% 3. exLC and exL2 values are optional however the information from them
% may decrease computation time

idx = [];
if strcmp(switches.solution,'REAL')
    if isempty(varargin)==0 & size(varargin,2)==1
        Q = varargin{1};
    elseif isempty(varargin)==0 & size(varargin,2)==4
        DNA = varargin{1};
        DEA = varargin{2};
        MDNA = varargin{3};
        MDEA = varargin{4};

    elseif isempty(varargin)==0 & size(varargin,2)==6
        DNA = varargin{1};
        DEA = varargin{2};
        MDNA = varargin{3};
        MDEA = varargin{4};
        GdN = varargin{5};
        GdE = varargin{6};
    elseif isempty(varargin)==0 & size(varargin,2)==7
        DNA = varargin{1};
        DEA = varargin{2};
        MDNA = varargin{3};
        MDEA = varargin{4};
        GdN = varargin{5};
        GdE = varargin{6};
        Q = varargin{7};
    end  
end

for t = 1:size(SP3X,1) 
    if ~isempty(ro)
        empt = ~cellfun(@isempty,ro);
    end
    try
        for i = 1:size(ro,2)
            if empt(t,i,1)
                station(t).ro(i).cordT = ro{t,i,1};
                station(t).ro(i).cordR = ro{t,i,2};
                station(t).ro(i).date = ro{t,i,3};
                station(t).ro(i).exL2 = ro{t,i,4};
                if strcmp(switches.solution,'REAL')
                    station(t).ro(i).exLC = ro{t,i,5};
                end
            end
        end
    end
    for nr = 1:size(BLh,1)
        station(t).h(nr).parameters = [];
        station(t).h(nr).nazwa = [];
        for sat = 1:max(max(PRN))
            if strcmp(switches.solution,'REAL')
                station(t).h(nr).satellite(sat).PRN = sat;
                station(t).h(nr).satellite(sat).elevation = NaN;
                station(t).h(nr).satellite(sat).azi = NaN;
                station(t).h(nr).satellite(sat).SWD = NaN;
                station(t).h(nr).satellite(sat).M_SWD = NaN;
                station(t).h(nr).satellite(sat).STD = NaN;
                station(t).h(nr).satellite(sat).M_STD = NaN;
                station(t).h(nr).satellite(sat).SIWV = NaN;
                station(t).h(nr).satellite(sat).M_SIWV = NaN;
                station(t).h(nr).satellite(sat).vmf1h = NaN;
                station(t).h(nr).satellite(sat).vmf1w = NaN;
                station(t).h(nr).satellite(sat).coord = [NaN NaN NaN];
            else
                station(t).h(nr).satellite(sat).PRN = sat;
                station(t).h(nr).satellite(sat).elevation = NaN;
                station(t).h(nr).satellite(sat).azi = NaN;
                station(t).h(nr).satellite(sat).vmf1h = NaN;
                station(t).h(nr).satellite(sat).vmf1w = NaN;
                station(t).h(nr).satellite(sat).coord = [NaN NaN NaN];
            end
        end
    end
    clear nr sat
    for s = 1 : length(BLh)
        if strcmp(switches.solution,'REAL')
            if size(varargin,2)==1
                station(t).h(s).ZTD(1,4) = Q(t,s);
            end
            if isnan(ZTD(t,s)) ==0 & ZTD(t,s)>0
                ZTDkop = ZTD(t,s)/1000;
            else
                ZTDkop = NaN;
            end
            station(t).h(s).ZTD(1,1) = ZTDkop;
            if size(varargin,2)==4
                station(t).h(s).ZTD(1,2) = DNA(t,s)/1000;
                station(t).h(s).ZTD(1,3) = DEA(t,s)/1000;
            end
            if size(varargin,2)==6
                station(t).h(s).ZTD(1,2) = DNA(t,s)/1000;
                station(t).h(s).ZTD(1,3) = DEA(t,s)/1000;
            end

            if size(varargin,2)==7
                station(t).h(s).ZTD(1,2) = DNA(t,s)/1000;
                station(t).h(s).ZTD(1,3) = DEA(t,s)/1000;
                station(t).h(s).ZTD(1,4) = Q(t,s);
            end
            clear  ZTDkop
            if any(isnan(station(t).h(s).ZTD))
                idx = [idx s];
            end
        end
        count = 0;
        for nr = 1 :size(PRN,2)
            doy_num = floor(observation_set(1,2) - floor(observation_set(1,2)))+1;
            sat = PRN(doy_num,nr);
            [X,Y,Z] = BLH2XYZ ( BLh(s,2)/180*pi(),BLh(s,3)/180*pi(), BLh(s,4));
            % Note that X Y Z (receiver) are in meters, SP3X SP3Y SP3Z (satellite) are in kilometers
            [elev , azi] = xyzSP32elaz([X Y Z],[SP3X(t,nr) SP3Y(t,nr) SP3Z(t,nr)].*1000);
            clear X Y Z
            elev = elev/pi()*180;
            azi = azi/pi()*180;
            coord = [SP3X(t,nr) SP3Y(t,nr) SP3Z(t,nr)];
            [w,~,~] = find(elev(:,1)>cut_off_angle);
            if isempty(w) == 0
                count = count +1;
                elev = elev(w,:);
                azi = azi(w,:);
                %VMF1 with GPT2 function used to calculate mapping fum
                mjdday = doy2jd(observation_set(t,3),floor(observation_set(t,2)));
                mjdday = jd2mjd(mjdday);
                if t == 1 && s ==1
                    mjdday_old = 0;
                    
                end
                if mjdday > mjdday_old
                    clear ah aw
                    [~,~,~,~,ah,aw,~] = gpt2 (mjdday,BLh(s,2)*pi()/180,BLh(s,3)*pi()/180,BLh(s,4),1,0);
                    mjdday_old = mjdday;
                end
                [vmf1h,vmf1w] = vmf1_ht (ah,aw,mjdday,BLh(s,2)*pi()/180,BLh(s,4),(pi()/2 - elev*pi()/180));  
            if strcmp(switches.solution,'REAL')
                if isempty(varargin) == 1
                    STD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + vmf1h*ZHD(t,s);
                    SWD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s));
                    M_SWD = vmf1w*M_ZTD(t,s)/1000;
                    M_STD = M_SWD+0.2*M_SWD;
                    SIWV = NaN;
                    M_SIWV = NaN;
                end
                if  isempty(varargin)==0 & size(varargin,2)==1
                    STD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + vmf1h*ZHD(t,s);
                    SWD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s));
                    M_SWD = vmf1w*M_ZTD(t,s)/1000;
                    M_STD = M_SWD+0.2*M_SWD;
                    SIWV = SWD./Q(t,s);
                    M_SIWV = M_SWD*0.15;
                end
                if  isempty(varargin)==0 & size(varargin,2)==4
                    STD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + vmf1h*ZHD(t,s) + cot(elev*pi()/180)*((DNA(t,s)/1000)*cos(azi*pi()/180)+(DEA(t,s)/1000)*sin(azi*pi()/180));
                    SWD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + cot(elev*pi()/180)*((DNA(t,s)/1000)*cos(azi*pi()/180)+(DEA(t,s)/1000)*sin(azi*pi()/180));
                    M_SWD = vmf1w*M_ZTD(t,s)/1000;
                    M_STD = M_SWD+0.2*M_SWD;
                    SIWV = NaN;
                    M_SIWV = NaN;
                end
                if  isempty(varargin)==0 & size(varargin,2)==6
                    STD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + vmf1h*ZHD(t,s) + cot(elev*pi()/180)*((DNA(t,s)/1000)*cos(azi*pi()/180)+(DEA(t,s)/1000)*sin(azi*pi()/180));
                    SWD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + cot(elev*pi()/180)*((DNA(t,s)/1000-GdN(t,s)/1000)*cos(azi*pi()/180)+(DEA(t,s)/1000 - GdE(t,s)/1000)*sin(azi*pi()/180));
                    M_SWD = vmf1w*M_ZTD(t,s)/1000;
                    M_STD = M_SWD+0.2*M_SWD;
                    SIWV = NaN;
                    M_SIWV = NaN;
                end
                if  isempty(varargin)==0 & size(varargin,2)==7
                    STD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + vmf1h*ZHD(t,s) + cot(elev*pi()/180)*((DNA(t,s)/1000)*cos(azi*pi()/180)+(DEA(t,s)/1000)*sin(azi*pi()/180));
                    SWD = vmf1w*(ZTD(t,s)/1000 - ZHD(t,s)) + cot(elev*pi()/180)*((DNA(t,s)/1000-GdN(t,s)/1000)*cos(azi*pi()/180)+(DEA(t,s)/1000 - GdE(t,s)/1000)*sin(azi*pi()/180));
                    SIWV = SWD./Q(t,s);
                    M_SWD = vmf1w*M_ZTD(t,s)/1000;
                    M_STD = M_SWD+0.2*M_SWD;
                    M_SIWV = M_SWD*0.15; %Very simplified version of error propagation!!!!
                end
            end
            if strcmp(switches.solution,'REAL')
                station(t).h(s).satellite(sat).SWD = SWD;
                station(t).h(s).satellite(sat).M_SWD = M_SWD;
                station(t).h(s).satellite(sat).STD = STD;
                station(t).h(s).satellite(sat).M_STD = M_STD;
                station(t).h(s).satellite(sat).SIWV= SIWV;
                station(t).h(s).satellite(sat).M_SIWV= M_SIWV;
                clear SWD M_SWD STD M_STD SIWV M_SIWV
            end
            else
                elev = NaN;
                azi = NaN;
                vmf1h = NaN;
                vmf1w = NaN;  
                coord = [NaN NaN NaN];
            end
            station(t).h(s).parameters = [BLh(s,1:3) BLh(s,4) BLH(s,4)];
            station(t).h(s).nazwa = NAME{s};
            station(t).h(s).satellite(sat).PRN = sat;
            station(t).h(s).satellite(sat).elevation = elev;
            station(t).h(s).satellite(sat).azi = azi;
            station(t).h(s).satellite(sat).vmf1h = vmf1h;
            station(t).h(s).satellite(sat).vmf1w = vmf1w;
            station(t).h(s).satellite(sat).coord = coord;
            clear  elevation azi m vmf1w vmf1w 
        end
    end
    if ~isempty(idx)
        station(t).h(idx) = [];
    end
    idx = [];
end


