function roty = roty(degrees)

radians = deg2rad(degrees);
roty = [cos(radians)    0   sin(radians);
        0               1   0;
        -sin(radians)   0   cos(radians)];

end