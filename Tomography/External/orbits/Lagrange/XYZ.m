function [stdX,XtractPoArr,TimeArr,index] = XYZ(timeX,TimeArr,PoArr,d,NA)

% Returns Standardzied Time & extracts PoArr for Lagrange Interpolation (LG)
% and time at the start point LG
% Called by SatPo
% timeX = GPS time epoch at which you need to interpolate
% TimeArr = GPS time array corresponding to array of PoArr
% PoArr = GPS coordinate X, Y or Z corresponding to TimeArr
% d = number of points used for the Lagrange interpolation
%
% Written by  Phakphong Homniam
% September 13, 2002
% Original Mathcad source code by Boonsap Witchayangkoon, 2000

rowTimeArr = size(TimeArr,1);
rowPoArr = size(PoArr,1);
if rowTimeArr ~= rowPoArr
    fprintf('Err(TimeArr, PoArr)\n');
    return
elseif d > rowTimeArr
    fprintf('Not enough data using d\n');
    return
end

% timeInterval = frequent GPS time to provide Satellite Position
timeInterval = TimeArr(3,1)-TimeArr(2,1);
PositionTimeArr = ceil((timeX-TimeArr(1,1))/timeInterval)+1;
if PoArr(PositionTimeArr,1) == NA
    stdX = PositionTimeArr;
    XtractPoArr = 0;
    TimeArr = 0;
    index = 0;
    return
end
if PositionTimeArr < ceil(d/2)
    for i = 1:d
        XtractPoArr(i,1) = PoArr(i,1);
    end
    stdX = ((timeX-TimeArr(1,1))/timeInterval)+1;
    XtractPoArr;
    TimeArr = TimeArr(1,1);
    index = 1;
    return
else
    NearEnd = rowTimeArr-PositionTimeArr;
    if NearEnd < floor(d/2)
        StartPt = rowTimeArr-d+1;
        for i = 1:d
            XtractPoArr(i,1) = PoArr(StartPt+i-1,1);
        end
    else
        for i = 1:d
            % Extract Satellite Position Array
            ptr = PositionTimeArr-floor(d/2)+i-1;
            XtractPoArr(i,1) = PoArr(ptr,1);
        end
        StartPt = PositionTimeArr-floor(d/2);
    end
    stdX = ((timeX-TimeArr(StartPt,1))/timeInterval)+1;
    XtractPoArr;
    TimeArr = TimeArr(StartPt,1);
    index = 2;
    return
end

%%%%%%%%%END%%%%%%%%%
        
