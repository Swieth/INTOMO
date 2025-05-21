function [yr,doy,sec] = jd2gps2(jd)
%function developed to extend the capability of geodetic package to be able
%change jd to year, day of the year (doy), seconds of the day.
%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%
%jd - joulian date
%%%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%
%yr - year
%doy - day of year
%sec - seconds of the day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[doy,yr] = jd2doy(jd);
doy = floor(doy);
[yr2,mm,day] = jd2cal(jd);
[hh,min,ss] = cal2time(day);
sec = ss + min*60 + hh*3600;
