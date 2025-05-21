function SP3data = readSP3dat(obs_set,pathORB);
    disp('SP3 processing');
    SP3data = [];
    gps_week_day = unique([obs_set.observation_set_SP3(:,4) obs_set.observation_set_SP3(:,6)], 'rows');
    if obs_set.observation_set_SP3(end,9)==23
       [gps_week_next, gps_s_next] = jd2gps(obs_set.observation_set_SP3(end,1)+1);
       gps_week_day = [gps_week_day;[gps_week_next floor(gps_s_next/(3600*24))]];
    end
    if obs_set.observation_set_SP3(end,9)==0 || obs_set.observation_set_SP3(end,9)==1
        [gps_week_prev, gps_s_prev] = jd2gps(obs_set.observation_set_SP3(end,1)-1);
        gps_week_day = [[gps_week_prev floor(gps_s_prev/(3600*24))];gps_week_day];
    end
    for j = 1: size(gps_week_day,1)
        filename{1} = ['igs' num2str(gps_week_day(j,1)) num2str(gps_week_day(j,2)) '.sp3'];
        filename{2} = ['igr' num2str(gps_week_day(j,1)) num2str(gps_week_day(j,2)) '.sp3'];
        filename{3} = ['igu' num2str(gps_week_day(j,1)) num2str(gps_week_day(j,2)) '_18.sp3'];
        filename{4} = ['igu' num2str(gps_week_day(j,1)) num2str(gps_week_day(j,2)) '_12.sp3'];
        filename{5} = ['igu' num2str(gps_week_day(j,1)) num2str(gps_week_day(j,2)) '_06.sp3'];
        filename{6} = ['igu' num2str(gps_week_day(j,1)) num2str(gps_week_day(j,2)) '_00.sp3'];
        for f = 1:size(filename,2)
           if exist([pathORB filename{f}],'file')==2
               file_ORB = dir([pathORB filename{f}]);
               break
           end
        end
        if exist('file_ORB','var')==1
           [tempSP3data,~,~] = readsp3([pathORB file_ORB.name]);
           SP3data = [SP3data; tempSP3data];
        end
        clear tempSP3data numsat header file_ORB filename1 filename2 filename3 filename4 filename5 filename6
    end   
end