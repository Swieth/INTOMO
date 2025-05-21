function  [ZWDnummiss,id_missZWD,sR,A,SWD_obs,R_SWD,un_name,ZWD,A_stC,SWD_stC,R_SWD_stC,R_ZWD] = deleteStatZWD(id_missZWD,ZWDnummiss,A,SWD_obs,R_SWD,un_name,ZWD,A_stC,SWD_stC,R_SWD_stC,R_ZWD,name,switches)
      try
        idcolR = [];
        idstatdel = find(ZWDnummiss<2);
        namedelstat = un_name(idstatdel);
        num = 1:length(name);
        for i = 1:size(namedelstat,1)
            id_missZWD(find(id_missZWD == idstatdel(i))) = [];
            idcolR = [idcolR,num(strcmp(namedelstat(i),name))];
            id_missZWD(id_missZWD>idstatdel(i)) = id_missZWD(id_missZWD>idstatdel(i))-1;
        end
        sR = length(idcolR);
        A(idcolR,:) = [];
        SWD_obs(idcolR) = [];
        R_SWD(idcolR) = [];
        un_name(idstatdel) = [];
        ZWD(idstatdel,:)= [];
        A_stC(idstatdel) =[];
        SWD_stC(idstatdel) =[];
        R_SWD_stC(idstatdel) =[];
        ZWDnummiss(idstatdel) =[];
        R_ZWD(idstatdel) =[];
        if strcmp(switches.solution,'REAL')
            SWD_obs_b(idcolR) = [];
        end
    catch
        warning('Deleting GNSS stations from weigting process failed (weightObs)')
    end
        idcolR = sort(idcolR);
        disp(['Cutting ',num2str(length(namedelstat)),' GNSS stations due to the impossiblity of covariance function estimation'])
end