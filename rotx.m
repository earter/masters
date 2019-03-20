function rotx = rotx(degrees)

radians = deg2rad(degrees);
rotx = [1               0             0; 
        0               cos(radians)  -sin(radians); 
        0               sin(radians)  cos(radians)];

end