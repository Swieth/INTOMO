% Function to obtain mean voxel difference of refraction value and height
% between nodes
%% Input data
% m1........... index of height layer above ray point
% m2........... index of height layer below ray point
% alt_N........ altitude profile
% inp_data
%% Output data
% dh........... height level difference in ray point's voxel [m]
% dh2.......... hegiht level difference in voxel below or above ray point [m]
% dref......... refraction value difference in ray point's voxel 
% dref2........ refraction value difference in voxel below or above ray point  m1/m2 - upper and lower index of grid (alt_N) layers
function [dref,dref2,dh,dh2,m1,m2] = dp_alt(m1,m2,alt_N,inp_data)
    try
        m1 = m1(1);
    end
    if isempty(m2)
        m2 = 1;
    else
        m2 = m2(end);
    end
    if m2 == size(alt_N,1) 
        m1 = m2-1;
    elseif m1 == 1
        m2 = 2;
    end
    if m1 < m2
        dref = diff(inp_data.ref0(m1:m2)); 
        dh = diff(alt_N(m1:m2)); 
        if m1 == 1
            m3 = 1;
            m4 = 2;
        else
            m3 = m1 - 1;
            m4 = m2 - 1;
        end
        dref2 = diff(inp_data.ref0(m3:m4)); 
        dh2 = diff(alt_N(m3:m4)); 
    else 
        try
            dref = diff(inp_data.ref0(m2:m1)); 
        catch
            dref = 0;
        end
        dh = diff(alt_N(m2:m1));
        if m2 == 1
            m3 = 2;
            m4 = 1;
        else
            m3 = m1 - 1;
            m4 = m2 - 1;
        end      
        dref2 = diff(inp_data.ref0(m4:m3)); 
        dh2 = diff(alt_N(m4:m3));
    end
end