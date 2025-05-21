function [stdX,XtractClkArr,index,ClkArr] = XCLK(timeX,TimeArr,ClkArr,d,NA)

% Returns Standardized Time & extracts ClkArr for Lagrange Interpolation(LG)
% Called by SVCLK
% timeX = GPS time that wants to compute Sat Clk Correction
% TimeArr = GPS time corresponding to Sat positions in PoArr
% ClkArr = GPS Clock Correction corresponding to time in TimeArr
% d = number of points used for the Lagrange interpolation
%
% Written by  Phakphong Homniam
% September 16, 2002
% Original Mathcad source code by Boonsap Witchayangkoon, 2000

rowTimeArr = size(TimeArr,1);
rowClkArr = size(ClkArr,1);
if rowTimeArr ~= rowClkArr
    fprintf('TimeArr not equal ClkArr\n');
    return
elseif d > rowTimeArr
    fprintf('Not enough data for using d\n');
    return
end

% timeInterval = frequent GPS time to provide Sat Position
timeInterval = TimeArr(3,1)-TimeArr(2,1);
PositionTimeArr = ceil((timeX-TimeArr(1,1))/timeInterval)+1;
if ClkArr(PositionTimeArr,1) == NA
    stdX = PositionTimeArr;
    XtractClkArr = 0;
    index = 0;
    ClkArr = NA;
    return
end
if PositionTimeArr < ceil(d/2)
    for i = 1:d
        XtractClkArr(i,1)= ClkArr(i,1);
    end
    stdX = ((timeX-TimeArr(1,1))/timeInterval)+1;
    XtractClkArr;
    index = 1;
    ClkArr = ClkArr(PositionTimeArr,1);
    return
else
    NearEnd = rowTimeArr-PositionTimeArr;
    if NearEnd < floor(d/2)
        StartPt = rowTimeArr-d+1;
        for i = 1:d
            XtractClkArr(i,1) = ClkArr(StartPt+i-1,1);
        end
    else
        for i = 1:d
            % Extract Satellite Position Array
            ptr = PositionTimeArr-floor(d/2)+i-1;
            XtractClkArr(i,1) = ClkArr(ptr,1);
        end
        StartPt = PositionTimeArr-floor(d/2);
    end
    stdX = ((timeX-TimeArr(StartPt,1))/timeInterval)+1;
    XtractClkArr;
    index = 2;
    ClkArr = ClkArr(PositionTimeArr,1);
    return
end

%%%%%%%%%%END%%%%%%%%%%
            
