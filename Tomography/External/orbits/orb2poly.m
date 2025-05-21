function [poly,fr,to]=orb2poly(data,deg,fitlen,arclen)
%orb2poly fits polynomials into SP3 data (separately for X,Y,Z,dTs)
%  [author]   : Tomasz Hadaœ, 18.01.2014
%  [argin]    : data [n,6] - SP3 data [epoch prn X[m] Y[m] Z[m] dTs[s*10^6]]
%                            (epoch is GPSweek*7*24*3600+GPSsec)
%               deg - polynomial degree, warning - interpolation in orb2poly will not work properly, if deg<16
%               fitlen - length of data used for polynomial fit [s]
%               arclen - max length of the arc [s] (only this part will be used for interpolation)
%  [argout]   : poly - array of structures (poly(prn) - indicates PRN):
%                             .from - epoch from which arc is valid
%                             .to   - epoch to which arc is valid
%                             .X [deg+1,1] - coefficients of polynomial: a(1)*x^n+a(2)*x^(n-1)+...+a(n)*x+a(n+1)
%                             .X_s and X_m - additional polynomials parameters from polyfit that remove bad scaling a
%                             .Y [deg+1,1] -  --||--
%                             .Z [deg+1,1] -  --||--
%                             .dT [deg+1,1] -  --||--
%               fr - [n_arc,1] with time indicating arc beginings
%               to - [n_arc,1] with time indicating arc ends
% WARNING - use [val delta]=polyval(.X,time,.X_s,.X_mu) to get value 'val' at epoch 'time'
% [changelog] 

if arclen>fitlen
    fitlen=arclen;
end

%ref_epoch=(max(data(:,1))-min(data(:,1)))/2; %reference epoch of data
%data(:,1)=data(:,1)-ref_epoch; % for better numerical performance of polyfit function (read WARNING above)

epoch0=min(data(:,1)); %first epoch of new data (should be 0)
epoch1=max(data(:,1)); %last epoch of new data
prns=unique(data(:,2));

%prealocate poly as array of structure
max_prns=fix((50+max(prns))/50)*50; 
n_arc=ceil((epoch1-epoch0)/arclen); %number of arcs
s=struct('from',NaN,'to',NaN,'X',NaN(deg+1,1),'X_s',NaN(deg+1,1),'X_mu',NaN(deg+1,1),'Y',NaN(deg+1,1),'Y_s',NaN(deg+1,1),'Y_mu',NaN(deg+1,1),'Z',NaN(deg+1,1),'Z_s',NaN(deg+1,1),'Z_mu',NaN(deg+1,1),'dT',NaN(deg+1,1),'dT_s',NaN(deg+1,1),'dT_mu',NaN(deg+1,1));
poly=repmat(s,max_prns,n_arc);

arc_to=epoch1;
% epoch1-epoch0
% n_arc
fr=NaN(n_arc,1);
to=NaN(n_arc,1);

for a=n_arc:-1:1 %start with the last arc - this is due to near real-time purpose
    arc_fr=arc_to-arclen+0.5;%
    fit_to=arc_to+(fitlen-arclen)/2;
    fit_fr=arc_fr-(fitlen-arclen)/2;
    % num2str([min(data(:,1)) fit_fr arc_fr arc_to fit_to max(data(:,1))])
    if fit_to>epoch1
        fit_to=epoch1;
    end
    if fit_fr<epoch0
        fit_fr=epoch0;
    end
    fdata=data((and(data(:,1)>=fit_fr,data(:,1)<=fit_to)),:);
    for i=1:length(prns)
        if i==1
            fr(a)=arc_fr;
            to(a)=arc_to;
        end
        p=prns(i);
        pfdata=fdata(fdata(:,2)==p,:);
        poly(p,a).to=arc_to;
        poly(p,a).from=arc_fr;
        if size(pfdata,1)>=deg
            [poly(p,a).X poly(p,a).X_s poly(p,a).X_mu]=polyfit(pfdata(:,1),pfdata(:,3),deg);
            [poly(p,a).Y poly(p,a).Y_s poly(p,a).Y_mu]=polyfit(pfdata(:,1),pfdata(:,4),deg);
            [poly(p,a).Z poly(p,a).Z_s poly(p,a).Z_mu]=polyfit(pfdata(:,1),pfdata(:,5),deg);
            [poly(p,a).dT poly(p,a).dT_s poly(p,a).dT_mu]=polyfit(pfdata(:,1),pfdata(:,6),deg);
        else
            % zostaj¹ wszêdzie NaN
        end
    end
    arc_to=arc_fr-0.5;
end

end %function