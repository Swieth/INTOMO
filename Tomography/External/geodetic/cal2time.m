function [hh,mm,sec] = cal2time(day)
%function developed to extend the capability of geodetic package to be able
%change day (given as day with decimal parts) into hh mm and sec
%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%
%day - decimal number showing number of calendar day
%%%%%%%%%%%%%%%%%%OUTPUT%%%%%%%%%%%%%%%%
%hh - hours
%mm - minutes
%sec - seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seconds = round((day - floor(day))*24*60*60);
hh = floor(seconds/60/60);
mm = floor(seconds/60) - hh*60;
sec = (seconds - mm*60 - hh*60*60);