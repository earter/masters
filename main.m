close all; clear all;
% % best?
ptCloud = pcread('Churyumov-Gerasimenko SPC 2017 - 96k poly.ply');
pcshow(ptCloud);

% % round
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

% % orbita: kolowa 30km
r = 30000; % [m]
T_all = 14*24*3600; % [s]
T = T_all/2;
l = 2*pi*r;

% simulation time: polowka orbity
sim_time = [0:1:T-1] ;

% linear velocity in orbit
v_lin = l/T_all; % [m/s]
vel_ang = v_lin / r;

% CB - celestial body is in [0,0,0] of frame
CB_loc = [0;0;0];
% SC location at time=0
sc_loc_0 = [sqrt(2)/2*r; sqrt(2)/2*r;0];

% figure; scatter3([CB_loc(1) sc_loc_0(1)], [CB_loc(2) sc_loc_0(2)], [CB_loc(3) sc_loc_0(3)]);

d_t = 1;
d_sc = v_lin * d_t / sqrt(2); 

sc_loc = zeros(3, length(sim_time));

% sc_loc_1 = sc_loc_0;
% for i=1:length(sim_time)
%     sc_loc(:,i) = sc_loc_1 + [d_sc; d_sc; 0];
%     
%     sc_loc_1 = sc_loc(:,i);
% end
% figure;plot(sc_loc(1,:), sc_loc(2,:))

d_k = d_t *2*pi/T;
k_1 = -d_k; % so k(1)=0

for i=1:length(sim_time)
    k = k_1 + d_k;
    sc_loc(1,i) = r*cos(k);
    sc_loc(2,i) = r*sin(k);

    k_1 = k;
end
figure; plot(sc_loc(1,:), sc_loc(2,:)); axis equal; hold on;
scatter(sc_loc(1,end), sc_loc(2,end));
legend("orbit","SC");

zz = ptCloud.Location';
cnt=0;
for i=1:length(zz)
    if mod(i,100)==0
        cnt = cnt+1;
        aa(:,cnt)=zz(:,i);
    end
end

figure; plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  hold on;
scatter3(aa(1,:),aa(2,:),aa(3,:)); hold off;


figure; scatter3(aa(1,:),aa(2,:),aa(3,:),'b.') % rozciagniete jakies
legend("hh")