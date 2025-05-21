function [obs_set] = time_listing(switches,observation_start,observation_end, est_interval_TOMO,obs_interval_SP3,obs_interval_ZTD,obs_interval_NWP,int_interval_METEO)
%Function to create matrix with dates of epochs
%-------------------------------------------------------------------------------------------
%%%INPUT
%       switches.......... structural matrix of settings of the INTOMO processing
%       observation_start.  time of initial epoch of the processing
%       observation_start.  time of final epoch of the processing
%       est_interval_TOMO. interval of TOMO processing [s]
%       obs_interval_SP3..  interval of SP3 data [s]
%       obs_interval_ZTD..  interval of ZTD data [s]
%       obs_interval_NWP..  interval of NWP data [s]
%       int_interval_METEO interval of meteorological data [s]
%%%OUTPUT
%       obs_set......    matrix with dates of epochs
%           .obervations_set...... tomography processing epochs
%           .obervations_set_SP3.. SP3 processing epochs
%           .interpolation_set.... ray tracing processing epochs
%               .all above in the format [time_since:DOY:YEAR:GPS_WEEK:seconds_since:DOW:MM:DD:HH:MM]

%%%%%%TOMO observation/TRO files%%%%%%%%
start_time = observation_start(1,4)/24 + observation_start(1,5)/60/24 + observation_start(1,6)/60/60/24 ;
jd_start = cal2jd(observation_start(1,1),observation_start(1,2),observation_start(1,3)+start_time);

end_time = observation_end(1,4)/24 + observation_end(1,5)/60/24 + observation_end(1,6)/60/60/24 ;
jd_end = cal2jd(observation_end(1,1),observation_end(1,2),observation_end(1,3)+end_time);

observation_set = jd_start:est_interval_TOMO/60/60/24:jd_end;

observation_set = observation_set';
[observation_set(:,2), observation_set(:,3)] = jd2doy(observation_set(:,1));
[observation_set(:,4), observation_set(:,5)] = jd2gps(observation_set(:,1));
[observation_set(:,6)] = jd2dow(observation_set(:,1)); 
observation_set(:,6) = observation_set(:,6)-1;% The - 1 due to the script producing SUN=1 SAT=7 GPS SUN=0 SAT=6
[~, observation_set(:,7), dayT] = jd2cal(observation_set(:,1));
clear yearT;
observation_set(:,8) = floor(dayT);
[hour,minutes,~] = cal2time(dayT);
observation_set(:,9) = hour;
observation_set(:,10) = minutes;

observation_set(1,11) = 1;

clear dayT hour minutes seconds file_mark
obs_set.observation_set = observation_set;

%%%%%%%%WRF%%%%%%%%%
if strcmp(switches.aprModel,'WRF') == 1
    observation_set_NWP = jd_start:obs_interval_NWP/60/60/24:jd_end;
    observation_set_NWP = observation_set_NWP';
    [observation_set_NWP(:,2), observation_set_NWP(:,3)] = jd2doy(observation_set_NWP(:,1));
    [observation_set_NWP(:,4), observation_set_NWP(:,5)] = jd2gps(observation_set_NWP(:,1));
    [observation_set_NWP(:,6)] = jd2dow(observation_set_NWP(:,1)); 
    observation_set_NWP(:,6) = observation_set_NWP(:,6)-1;% The - 1 due to the script producing SUN=1 SAT=7 GPS SUN=0 SAT=6
    [~, observation_set_NWP(:,7), dayT] = jd2cal(observation_set_NWP(:,1));
    clear yearT;
    observation_set_NWP(:,8) = floor(dayT);
    [hour,minutes,~] = cal2time(dayT);
    observation_set_NWP(:,9) = hour;
    observation_set_NWP(:,10) = minutes;
    file_mark = find(((hour == 0)|(hour == 6)|(hour == 12)|(hour == 18))&(minutes==0));
    observation_set_NWP(file_mark,11) = 1;
    clear dayT hour minutes seconds file_mark
    obs_set.observation_set_NWP = observation_set_NWP;
end

%%%%%%%%SP3 Files%%%%%%%%%
observation_set_SP3 = jd_start:obs_interval_SP3/60/60/24:jd_end;
observation_set_SP3 = observation_set_SP3';
[observation_set_SP3(:,2), observation_set_SP3(:,3)] = jd2doy(observation_set_SP3(:,1));
[observation_set_SP3(:,4), observation_set_SP3(:,5)] = jd2gps(observation_set_SP3(:,1));
[observation_set_SP3(:,6)] = jd2dow(observation_set_SP3(:,1)); 
observation_set_SP3(:,6) = observation_set_SP3(:,6)-1;% The - 1 due to the script producing SUN=1 SAT=7 GPS SUN=0 SAT=6
[~, observation_set_SP3(:,7), dayT] = jd2cal(observation_set_SP3(:,1));
clear yearT;
observation_set_SP3(:,8) = floor(dayT);
[hour,minutes,~] = cal2time(dayT);
observation_set_SP3(:,9) = hour;
observation_set_SP3(:,10) = minutes;
file_mark = find((hour == 0)&(minutes==0));
observation_set_SP3(file_mark,11) = 1;
clear dayT hour minutes seconds file_mark
obs_set.observation_set_SP3 = observation_set_SP3;



%%%%% Time interpolation limits %%%%%%%%%%%%%%

interpolation_set = jd_start:int_interval_METEO/60/60/24:jd_end;
interpolation_set = interpolation_set';
[interpolation_set(:,2), interpolation_set(:,3)] = jd2doy(interpolation_set(:,1));
[interpolation_set(:,4), interpolation_set(:,5)] = jd2gps(interpolation_set(:,1));
[interpolation_set(:,6)] = jd2dow(interpolation_set(:,1)); 
interpolation_set(:,6) = interpolation_set(:,6)-1;% The - 1 due to the script producing SUN=1 SAT=7 GPS SUN=0 SAT=6
[~, interpolation_set(:,7), dayT] = jd2cal(interpolation_set(:,1));
clear yearT;
interpolation_set(:,8) = floor(dayT);
[hour,minutes,~] = cal2time(dayT);
interpolation_set(:,9) = hour;
interpolation_set(:,10) = minutes;
clear dayT hour minutes seconds 
obs_set.interpolation_set = interpolation_set;
%%%%%%%%%%%%%%%%%%%%%%%%%




