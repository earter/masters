close all; clear all; clc
% % best?
ptCloud = pcread('Churyumov-Gerasimenko SPC 2017 - 96k poly.ply');
% figure; pcshow(ptCloud);

% % too round
% ptCloud = pcread('Ceres OpNav5 24k poly.ply');
% figure; pcshow(ptCloud);

% % too dense
% ptCloud = pcread('Churyumov-Gerasimenko SPC 2017 - 199k poly.ply');
% figure; pcshow(ptCloud);

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


% simulation time: polowka orbity
sim_time = [0:1:T-1] ;
aa=1/8*length(sim_time);
% sim_time = [sim_time(1):sim_time(aa)];
% linear velocity in orbit
v_lin = l/T_all; % [m/s]
vel_ang = v_lin / r;

% CB - celestial body is in [0,0,0] of frame
CB_loc = [0;0;0];
% SC location at time=0
sc_loc_0 = [sqrt(2)/2*r; sqrt(2)/2*r;0];

% figure; scatter3([CB_loc(1) sc_loc_0(1)], [CB_loc(2) sc_loc_0(2)], [CB_loc(3) sc_loc_0(3)]);

d_t = 1; % time resolution (step)
d_sc = v_lin * d_t / sqrt(2); 



% sc_loc_1 = sc_loc_0;
% for i=1:length(sim_time)
%     sc_loc(:,i) = sc_loc_1 + [d_sc; d_sc; 0];
%     
%     sc_loc_1 = sc_loc(:,i);
% end
% figure;plot(sc_loc(1,:), sc_loc(2,:))

d_k = d_t*2*pi/T; % delta argument (parameter)
k_1 = -d_k; % so k(1)=0

rotx = [1 0 0; 
        0 cos(inclination) -sin(inclination); 
        0 sin(inclination) cos(inclination)];


sc_loc = zeros(3, length(sim_time));

for i=1:length(sim_time)
    k = k_1 + d_k;
    sc_loc(1,i) = r*cos(k);
    sc_loc(2,i) = r*sin(k);
    sc_loc(3,i) = 0;
    
    sc_loc(:,i) = rotx*sc_loc(:,i);
    
    k_1 = k;
end

% figure; plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  hold on; axis equal;
% scatter3(sc_loc(1,end), sc_loc(2,end), sc_loc(3,end));
% scatter3(0,0,0,'ok');
% legend("orbit","SC"); grid on;

% % SC LOS
% ROT_CB_SC = R(Z,arg) * R(Y,0) * R(X,30deg)
arg = 1/8*pi;
arg = 30000;

% rotx = [1 0 0; 
%         0 cos(xx) -sin(xx); 
%         0 sin(xx) cos(xx)];
yy = arg;
roty = [cos(yy) 0 sin(yy);
        0 1 0;
        -sin(yy) 0 cos(yy)];
    zz = 90* pi/180;
rotz = [cos(zz) -sin(zz) 0; 
        sin(zz) cos(zz) 0;
        0 0 1;];
yy = -90 * pi/180;    
      roty2 = [cos(yy) 0 sin(yy);
        0 1 0;
        -sin(yy) 0 cos(yy)];  
    
ROT_CB_SC = roty * rotz * roty2;

CB_LOS = ROT_CB_SC * SC_LOS;

scale=10000;

arg = 30000;
figure; 
plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  
hold on; axis equal;
plot3([sc_loc(1,arg) sc_loc(1,arg)+CB_LOS(1)*scale], [sc_loc(2,arg) sc_loc(2,arg)+CB_LOS(2)*scale], [sc_loc(3,arg) sc_loc(3,arg)+CB_LOS(3)*scale]);
scatter3(sc_loc(1,arg), sc_loc(2,arg), sc_loc(3,arg));
scatter3(0,0,0,'ok');
legend("orbit", "LOS","SC","CB"); grid on;

zz = ptCloud.Location';
cnt=0;
for i=1:length(zz)
    if mod(i,100)==0
        cnt = cnt+1;
        aa(:,cnt)=zz(:,i);
    end
end

% figure; 
% scatter3(1000*aa(1,:),1000*aa(2,:),1000*aa(3,:),'.'); 
% hold on; axis auto;
% plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));
% zlim([-2000 2000]);

% figure; scatter3(aa(1,:),aa(2,:),aa(3,:),'b.') % rozciagniete jakies
% legend("hh"); axis equal