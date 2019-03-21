close all; clear all; clc
run parameters.m;

% ptCloud = pcread('Churyumov-Gerasimenko SPC 2017 - 96k poly.ply');
% figure; pcshow(ptCloud); xlabel('x axis'); ylabel('y axis');
% ptCloud = (rotx(30) * ptCloud.Location')';
% figure; pcshow(ptCloud); xlabel('x axis'); ylabel('y axis');

% simulation time: polowka orbity
sim_time = [0:1:T-1] ;
% aa=1/8*length(sim_time);          % skrócić Tsim 
% sim_time = [sim_time(1):sim_time(aa)];

% velocity in orbit
v_lin = l/T_all; % [m/s]
vel_ang = v_lin / r;

% CB - celestial body% figure; 
% plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  
% hold on; axis equal;
% plot3([sc_loc(1,arg) sc_loc(1,arg)+CB_LOS(1)*scale], [sc_loc(2,arg) sc_loc(2,arg)+CB_LOS(2)*scale], [sc_loc(3,arg) sc_loc(3,arg)+CB_LOS(3)*scale]);
% scatter3(sc_loc(1,arg), sc_loc(2,arg), sc_loc(3,arg));
% scatter3(0,0,0,'ok');
% legend("orbit", "LOS","SC","CB"); grid on;
% xlabel('x axis'); ylabel('y axis')
% %CB  is in [0,0,0] of CB (global) frame
CB_loc = [0;0;0];
% SC location at time=0
% sc_loc_0 = [sqrt(2)/2*r; sqrt(2)/2*r; 0];
sc_loc_0 = [r; 0; 0];
% figure; 
% plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  
% hold on; axis equal;
% plot3([sc_loc(1,arg) sc_loc(1,arg)+CB_LOS(1)*scale], [sc_loc(2,arg) sc_loc(2,arg)+CB_LOS(2)*scale], [sc_loc(3,arg) sc_loc(3,arg)+CB_LOS(3)*scale]);
% scatter3(sc_loc(1,arg), sc_loc(2,arg), sc_loc(3,arg));
% scatter3(0,0,0,'ok');
% legend("orbit", "LOS","SC","CB"); grid on;
% xlabel('x axis'); ylabel('y axis');
% figure; scatter3([CB_loc(1) sc_loc_0(1)], [CB_loc(2) sc_loc_0(2)], [CB_loc(3) sc_loc_0(3)]);

dt = 1;                        % time resolution (step) [s]
d_sc = v_lin * dt / sqrt(2);   % delta linear SC position
d_k = dt*2*pi/T;               % delta angle argument (parameter)
k_1 = -d_k;                    % so k(1)=0

sc_loc = zeros(3, length(sim_time));
% sc_loc(:,1) = sc_loc_0;

for i=1:length(sim_time)
    if i==30000
        i;
    end
    k = k_1 + d_k;
    sc_loc(1,i) = r*cos(k);
    sc_loc(2,i) = r*sin(k);
    sc_loc(3,i) = 0;
   
    
    sc_loc(:,i) = rotx(30)*sc_loc(:,i);
    
    k_1 = k;
end

% figure; plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  hold on; axis equal;
% scatter3(sc_loc(1,end), sc_loc(2,end), sc_loc(3,end));
% scatter3(0,0,0,'ok');
% legend("orbit","SC"); grid on; xlabel('x axis'); ylabel('y axis')

%% SC LOS

ROT_CB_SC = rotx(90) * roty(-90) * roty(rad2deg(0.2917)) * rotx(asind(sc_loc(3,30000)/r));
CB_LOS = ROT_CB_SC * SC_LOS;

scale = r;

arg = 30000;
% arg = 1;
figure; 
plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));  
hold on; axis equal;
plot3([sc_loc(1,arg) sc_loc(1,arg)+CB_LOS(1)*scale], [sc_loc(2,arg) sc_loc(2,arg)+CB_LOS(2)*scale], [sc_loc(3,arg) sc_loc(3,arg)+CB_LOS(3)*scale]);
scatter3(sc_loc(1,arg), sc_loc(2,arg), sc_loc(3,arg));
scatter3(0,0,0,'ok');
legend('orbit', 'LOS','SC','CB'); grid on;
xlabel('x axis'); ylabel('y axis')

% zz = ptCloud.Location';
% cnt=0;
% for i=1:length(zz)
%     if mod(i,100)==0
%         cnt = cnt+1;
%         aa(:,cnt)=zz(:,i);
%     end
% end

% figure; 
% scatter3(1000*aa(1,:),1000*aa(2,:),1000*aa(3,:),'.'); 
% hold on; axis auto;
% plot3(sc_loc(1,:), sc_loc(2,:), sc_loc(3,:));
% zlim([-2000 2000]);

% figure; scatter3(aa(1,:),aa(2,:),aa(3,:),'b.') % rozciagniete jakies
% legend("hh"); axis equal