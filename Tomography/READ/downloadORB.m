function downloadORB(pathORB,observation_set_SP3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function developed to download sp3 files from ftp server igs.bkg.bund.de
%%%%%%%%%%%%%%%%%%%%%%%input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pathORB - download path
%observation_set - the matrix containing time

% gps time
gps_week_day = unique([observation_set_SP3(:,4) observation_set_SP3(:,6)], 'rows');
if observation_set_SP3(end,9)==23
    [gps_week_next, gps_s_next] = jd2gps(observation_set_SP3(end,1)+1);
    gps_week_day = [gps_week_day;[gps_week_next floor(gps_s_next/(3600*24))]];
end
if observation_set_SP3(end,9)<6
    [gps_week_prev, gps_s_prev] = jd2gps(observation_set_SP3(end,1)-1);
    gps_week_day = [gps_week_day;[gps_week_prev floor(gps_s_prev/(3600*24))]];
end

[gps_week_prev, gps_s_prev] = jd2gps(observation_set_SP3(end,1)-1);
gps_week_day_prev = [gps_week_prev floor(gps_s_prev/(3600*24))];

%download sp3 files
for i = 1:size(gps_week_day,1)
gps_names{1} = ['igs' num2str(gps_week_day(i,1)) num2str(gps_week_day(i,2)) '.sp3'];
gps_names{2} = ['igr' num2str(gps_week_day(i,1)) num2str(gps_week_day(i,2)) '.sp3'];
gps_names{3} = ['igu' num2str(gps_week_day(i,1)) num2str(gps_week_day(i,2)) '_18.sp3'];
gps_names{4} = ['igu' num2str(gps_week_day(i,1)) num2str(gps_week_day(i,2)) '_12.sp3'];
gps_names{5} = ['igu' num2str(gps_week_day(i,1)) num2str(gps_week_day(i,2)) '_06.sp3'];
gps_names{6} = ['igu' num2str(gps_week_day(i,1)) num2str(gps_week_day(i,2)) '_00.sp3'];
gps_names{7} = ['igu' num2str(gps_week_day_prev(1,1)) num2str(gps_week_day_prev(1,2)) '_18.sp3'];
gps_names{8} = ['igu' num2str(gps_week_day_prev(1,1)) num2str(gps_week_day_prev(1,2)) '_12.sp3'];

for j = 1:size(gps_names,2)
    if exist([pathORB gps_names{j}],'file')==0
        p = ftp('igs.bkg.bund.de');
        if isempty(dir(p,['/IGS/products/orbits/' num2str(gps_week_day(i,1)) '/' gps_names{j} '.Z']))==0
            cd(p, ['/IGS/products/orbits/' num2str(gps_week_day(i,1)) '/']);
            mget(p,[gps_names{j} '.Z'],pathORB);
            system(['uncompress ' pathORB gps_names{j} '.Z']);
            break
        end
        close(p);
    else
        break
    end
end
end

end

