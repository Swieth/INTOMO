function [idlon,idlat,idalt,switches] = select_id(switches,pos_X1,idlon,idlat,idalt,b)
    if b == 1
        if switches.direction(1) 
            idlon = max(pos_X1(3),pos_X1(4));
        else
            idlon = min(pos_X1(3),pos_X1(4));
        end
        if switches.direction(2)
            idlat = max(pos_X1(1),pos_X1(2));
        else
            idlat = min(pos_X1(1),pos_X1(2));
        end
        if switches.direction(3)
            idalt = max(pos_X1(5),pos_X1(6));
        else
            idalt = min(pos_X1(5),pos_X1(6));
        end
    else
        if switches.direction(1) && switches.change(1)
            idlon = idlon + 1;
        elseif switches.change(1)
            idlon = idlon - 1;
        end
        if switches.direction(2) && switches.change(2)
            idlat = idlat + 1;
        elseif switches.change(2)
            idlat = idlat - 1;
        end
        if switches.direction(3) && switches.change(3)
            idalt = idalt + 1;
        elseif switches.change(3)
            idalt = idalt - 1;
        end
        switches.change = [false;false;false];
    end
end