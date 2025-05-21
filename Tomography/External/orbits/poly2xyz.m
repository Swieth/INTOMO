%input: epoch_sat, curr_prns, arc_inx
% output: Xst, Yst, Zst, vXst, vYst, vZst
temp.emptyD=NaN(n_systems*50,1);
Xst=temp.emptyD;
Yst=temp.emptyD;
Zst=temp.emptyD;
vXst=temp.emptyD;
vYst=temp.emptyD;
vZst=temp.emptyD;
temp.dtime=0.1;

for i=1:length(curr_prns)
    p=curr_prns(i);
    p_epoch=epoch_sat(p,1)*7*24*3600+epoch_sat(p,2);
    if sum(isnan(SP3.poly(p,arc_inx).X))+sum(isnan(SP3.poly(p,arc_inx).Y))+sum(isnan(SP3.poly(p,arc_inx).Z))==0
        %X and vX
        [Xst(p),temp.delta]=polyval(SP3.poly(p,arc_inx).X,p_epoch,SP3.poly(p,arc_inx).X_s,SP3.poly(p,arc_inx).X_mu);
        [temp.X,temp.delta]=polyval(SP3.poly(p,arc_inx).X,p_epoch-temp.dtime,SP3.poly(p,arc_inx).X_s,SP3.poly(p,arc_inx).X_mu);
        vXst(p)=(Xst(p)-temp.X)/temp.dtime;
        %Y and vY    
        [Yst(p),temp.delta]=polyval(SP3.poly(p,arc_inx).Y,p_epoch,SP3.poly(p,arc_inx).Y_s,SP3.poly(p,arc_inx).Y_mu);
        [temp.Y,temp.delta]=polyval(SP3.poly(p,arc_inx).Y,p_epoch-temp.dtime,SP3.poly(p,arc_inx).Y_s,SP3.poly(p,arc_inx).Y_mu);
        vYst(p)=(Yst(p)-temp.Y)/temp.dtime;
        %Z and vZ
        [Zst(p),temp.delta]=polyval(SP3.poly(p,arc_inx).Z,p_epoch,SP3.poly(p,arc_inx).Z_s,SP3.poly(p,arc_inx).Z_mu);
        [temp.Z,temp.delta]=polyval(SP3.poly(p,arc_inx).Z,p_epoch-temp.dtime,SP3.poly(p,arc_inx).Z_s,SP3.poly(p,arc_inx).Z_mu);
        vZst(p)=(Zst(p)-temp.Z)/temp.dtime;
    else
        Xst(p)=NaN; vXst(p)=NaN;
        Yst(p)=NaN; vYst(p)=NaN;
        Zst(p)=NaN; vZst(p)=NaN;
    end
end

% scale from km to m
Xst=Xst*1000;
Yst=Yst*1000;
Zst=Zst*1000;
vXst=vXst*1000;
vYst=vYst*1000;
vZst=vZst*1000;

% apply Sat. Antena Offsets
if ~exist('curr_date','var')
    curr_date=RNX.epoch';
end
[mjd full_mjd frac_mjd]=date2mjd(curr_date); % modified julian date for the current date
sid_time=mjd2sdt(mjd,curr_date); % siderial time
temp.xsun=suncrd(mjd,sid_time); % sun coordinates
for i=1:length(ico)
    temp.prn=ico(i);
    temp.XYZ=(santoff([Xst(temp.prn) Yst(temp.prn) Zst(temp.prn)]',SatOffsets(temp.prn,:),temp.xsun))';
    Xst(temp.prn)=temp.XYZ(1);
    Yst(temp.prn)=temp.XYZ(2);
    Zst(temp.prn)=temp.XYZ(3);
end

clear temp