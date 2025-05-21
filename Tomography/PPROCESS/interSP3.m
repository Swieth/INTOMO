function [SP3X,SP3Y,SP3Z,PRN,observation_set] = interSP3(SP3data,observation_set)
%Function developed to interpolate orbits to the time specified in the
%observation
%%%INPUT
%   SP3data......... orbit data from sp3 files
%   observation_set. matrix with dates of epochs
%%%OUTPUT
%   SP3X............ X coordinate of GNSS satellite
%   SP3Y............ Y coordinate of GNSS satellite
%   SP3Z............ Z coordinate of GNSS satellite
%   PRN............. PRN of GNSS satellite
%   observation_set. matrix with dates of epochs

%%%%%%%%checking the time in SP3
startSP3 = min(SP3data(:,1));
endSP3 = max(SP3data(:,1));
%%%%%%%%checking the time in observations
startOBS = observation_set(1,4)*7*24*3600 + observation_set(1,5);
endOBS = observation_set(end,4)*7*24*3600 + observation_set(end,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%checking time overlap%%%%%%%%%%%%%%%
if startOBS > endSP3 
    sprintf(['Check your orbits! Starting time for observation is: ' num2str(startOBS) ' and end time for orbits is: ' num2str(endSP3) '\n'])
    sprintf(['this makes ' num2str((endSP3 - startOBS)/3600) 'hours difference'])
    SP3X = [];
    SP3Y = [];
    SP3Z = [];
    PRN = [];
    observation_set = [];
    return;
    
    
elseif startSP3 > endOBS
    sprintf(['Check your orbits! End time for observation is: ' num2str(endOBS) ' and start time for orbits is: ' num2str(startSP3) '\n'])
    sprintf(['this makes ' num2str((startSP3 - endOBS)/3600) 'hours difference'])
    SP3X = [];
    SP3Y = [];
    SP3Z = [];
    PRN = [];
    observation_set = [];   
    return;

elseif startOBS >= startSP3 && endSP3 < endOBS
    endOBS = endSP3;
    [wie,kol,wart] = find((observation_set(:,4)*7*24*3600 + observation_set(:,5))>endOBS);
    observation_set(wie,:) = [];
    clear wie kol wart
elseif startSP3 >= startOBS && endOBS < endSP3 
    startOBS = startSP3;
    [wie,kol,wart] = find((observation_set(:,4)*7*24*3600 + observation_set(:,5))<startOBS);
    observation_set(wie,:) = [];
    clear wie kol wart
    
elseif startSP3>startOBS && endOBS > endSP3
    startOBS = startSP3;
    endOBS = endSP3;     
    [wie,kol,wart] = find((observation_set(:,4)*7*24*3600 + observation_set(:,5))>endOBS);
    observation_set(wie,:) = [];
    clear wie kol wart
    
    [wie,kol,wart] = find((observation_set(:,4)*7*24*3600 + observation_set(:,5))<startOBS);
    observation_set(wie,:) = [];
    clear wie kol wart
        
end
    p = unique(SP3data(:,2));
    maxp = max(p);
    selectSP3sat = SP3data(find(SP3data(:,2)==p(1)==1),:); 
    resolutionTIME = median(diff(selectSP3sat(:,1)));
    clear selectSP3sat
    
    SP3X(1:size(observation_set,1),1:maxp) = NaN;
    SP3Y(1:size(observation_set,1),1:maxp) = NaN;
    SP3Z(1:size(observation_set,1),1:maxp) = NaN;
    PRN(1,1:maxp) = 1:maxp;
    PRN = repmat(PRN,size(observation_set,1),1);

[poly,fr,to]=orb2poly(SP3data,16,resolutionTIME*12*3,resolutionTIME*12);
% figure(1)
%    hold on
for j = 1:size(p,1)
%for j = 1
   for i = 1:size(observation_set,1) %wspolrzedne co 15 minut 
    
    timeCURRENT = observation_set(i,4)*7*24*3600 + observation_set(i,5);
    test = [fr - timeCURRENT to - timeCURRENT];
    w = find(test(:,1)< 0.5 & test(:,2)>=0.5);
    if  sum(isnan(poly(p(j),w).X))==0
       [SP3X(i,p(j)),temp.delta]=polyval(poly(p(j),w).X,timeCURRENT,poly(p(j),w).X_s,poly(p(j),w).X_mu);
       [SP3Y(i,p(j)),temp.delta]=polyval(poly(p(j),w).Y,timeCURRENT,poly(p(j),w).Y_s,poly(p(j),w).Y_mu);
       [SP3Z(i,p(j)),temp.delta]=polyval(poly(p(j),w).Z,timeCURRENT,poly(p(j),w).Z_s,poly(p(j),w).Z_mu);
        PRN(i,p(j))=p(j);
    end
%    scatter3(SP3X(i,1),SP3Y(i,1),SP3Z(i,1),'or')
%    text(SP3X(i,1),SP3Y(i,1),SP3Z(i,1),num2str(i))
   clear test w;
   end
end