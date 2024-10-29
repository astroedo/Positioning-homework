function [R1, R3_w, R3_W] = ORStoITRF(i, w, W)
% Rotation matrices from ORS to ITRF
    R1 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
    R3_w = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];
    R3_W = [cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];
end
