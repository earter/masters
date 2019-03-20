% % rosetta SC
sc_mass = 1230;        %[kg]
sc_z = 2.8;      % dimensions [m]
sc_y = 2;
sc_x = 2.1;
I_z = 1/12*sc_mass*(sc_y^2+sc_x^2); % moments of inertia [kg*m^2]
I_x = 1/12*sc_mass*(sc_z^2+sc_y^2);
I_y = 1/12*sc_mass*(sc_z^2+sc_x^2);
SC_LOS = [0; 0; 1];

% % orbita: kolowa 30km: GMP Global Mapping Phase
r = 30000; % [m]
T_all = 14*24*3600; % [s]
T = T_all/2;
l = 2*pi*r;

% Keplarian elements
eccentricity = 0;
epoch_time = 0;
inclination = pi/180 * 30;
right_ascension = 0;
arg_perigee = 0; % zalezne od czasu symulacji, perigee zgodne z osią X asteroidy

G = 6.672e-11;
M = 9.9828e12;
v_kepler = G*M/r;
