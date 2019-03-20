function rotz = rotz(degrees)

radians = deg2rad(degrees);
rotz = [cos(radians)    -sin(radians)   0; 
        sin(radians)    cos(radians)    0;
        0               0               1];
end