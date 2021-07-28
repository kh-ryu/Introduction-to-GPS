clear all;
close all;
clc;

load(['EXP_20200922_1_data.mat']) % Group 1 data
%load(['EXP_20191016_G2_data.mat']) % Group 2 data

%% time setting
ini_epoch = 10*60;
n_epoch = size(true.GPSTime,1); % num of GPS time steps
maskangle = 15;                 
nt = 32;

%% Visible satellite selection
PRN = [1:32];
vis_sat = zeros(1,32);

for ii = 1:32
    
    % mask angle 15 deg for both user and reference station
    temp_el1 = find(ref.El(:,ii) < maskangle| user.El(:,ii) < maskangle); 
    temp_el2 = find(ref.pr1(:,ii)<10 |user.pr1(:,ii)<10);
    if isempty(temp_el1)
        if isempty(temp_el2)
            vis_sat(ii)=1;
        end
    end
end
SV_vis = find(vis_sat==1);
n_vis = length(SV_vis);

%% Standalone 
user.pos_standalone = zeros(n_epoch,3);
user.pos_standalone_error = zeros(n_epoch,3);
user.pos_standalone_enu = zeros(n_epoch,3);
R_ur = zeros(n_epoch,3);

VDOP = zeros(n_epoch,1);
HDOP = zeros(n_epoch,1);
PDOP = zeros(n_epoch,1);
TDOP = zeros(n_epoch,1);
GDOP = zeros(n_epoch,1);
 
 for ti = 1: n_epoch
     dx = 100;
     
     % Step 1.
     R_user = zeros(3,1).'; % Initially in center of the earth
     B_ur = 0;  % What is B
     x_old = [R_user.'; B_ur];
     
     while (dx > 10^(-4)) 
         H_standalone = zeros(n_vis,4);
         z_standalone =  zeros(n_vis,1);

         for jj = 1: n_vis
             
             % Step 2. e_hat
             R_su = user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)) - R_user;
             e_hat = R_su/norm(R_su);
             
             % Step 3. 
             % Step 3.1. H matrix 
             H_standalone(jj,1:3) = e_hat;
             H_standalone(jj,4) = -1;

             % Step 3.2. measurement vector
             z_standalone(jj,1) = e_hat*user.svpos(ti,3*SV_vis(jj)-2:3*SV_vis(jj)).'-(user.pr1(ti,SV_vis(jj))+user.b(ti,SV_vis(jj))*c);
         end
         
         % Step 3.3. user position calculation
         x = pinv(H_standalone)*z_standalone;
         R_user = x(1:3).';
         B_ur = x(4);   
         
         % Step 4. new position update
         dx = norm(x-x_old);
         x_old = x;
         
     end
     
     % Calculated result
     user.pos_standalone(ti,:) = R_user;
     user.pos_standalone_enu(ti,:) = ( Rtran*(user.pos_standalone(ti,:) - ref_xyz).').';
     error.standalone_enu(ti,:) = user.pos_standalone_enu(ti,:) - true.enu(ti,:);
     
     % Calculate DOP -> Convert to ENU!
     e_enu = zeros(3, n_vis);    % Assembly of e-values in enu
     
     for i = 1:n_vis
         e_enu(:, i) = H_standalone(i, 1:3)';
         e_enu(:, i) = Rtran*e_enu(:, i);
         length = sqrt(sum(e_enu(:,i).*e_enu(:,i)));
        
         if length > 1
             e_enu(:,i) = e_enu(:,i)/length;
         end
     end
     
     H = horzcat(e_enu', -ones(n_vis,1));
     
     DOP = inv(H'*H);
     HDOP(ti) = sqrt(DOP(1,1) + DOP(2,2));
     VDOP(ti) = sqrt(DOP(3,3));
     PDOP(ti) = sqrt(DOP(1,1) + DOP(2,2) + DOP(3,3));
     TDOP(ti) = sqrt(DOP(4,4));
     GDOP(ti) = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3)+DOP(4,4));
 end

% DGPS from pseudorange -> additional point
%% Measured position vs real position
figure(1)
plot3(user.pos_standalone_enu(1:end,1),user.pos_standalone_enu(1:end,2),user.pos_standalone_enu(1:end,3),'k.');
hold on
plot3(true.enu(:,1), true.enu(:,2), true.enu(:,3),'r.')

title('Standalone Trajectory')
xlabel('East(meter)')
ylabel('North(meter)')
zlabel('Altitude(meter)')
legend('Standalone', 'True')
grid on

%% Plot on East-North, Up
figure(2)
plot(user.pos_standalone_enu(1:n_epoch,1),user.pos_standalone_enu(1:n_epoch,2), 'x');
hold on
plot(true.enu(:,1), true.enu(:,2))
title('East-North Trajectory')
xlabel('East(meter)')
ylabel('North(meter)')
legend('Standalone', 'True')
grid on

figure(3)
plot(true.GPSTime(1:n_epoch), user.pos_standalone_enu(1:n_epoch,3),'x')
hold on
plot(true.GPSTime(1:n_epoch), true.enu(1:n_epoch,3))
title('Time-Vertical Trajectory')
xlabel('Time(second)')
ylabel('Vertical(meter)')
legend('DGPS', 'True')
grid on

%% Plot Dop
toffset = true.GPSTime(1);
t = true.GPSTime(1:n_epoch)-toffset;

figure(4)
plot(t, VDOP)
title("VDOP vs time")
xlabel("time(sec)")
ylabel("VDOP")
grid on

figure(5)
plot(t, HDOP)
title("HDOP vs time")
xlabel("time(sec)")
ylabel("HDOP")
grid on

figure(6)
plot(t, PDOP)
title("PDOP vs time")
xlabel("time(sec)")
ylabel("PDOP")
grid on

figure(11)
plot(t, TDOP)
title("TDOP vs time")
xlabel("time(sec)")
ylabel("TDOP")
grid on

figure(12)
plot(t, GDOP)
title("GDOP vs time")
xlabel("time(sec)")
ylabel("GDOP")
grid on

%% Tendency of error
figure(7)
plot(t, error.standalone_enu(:, 3))
xlabel("time(sec)")
ylabel("Error(m)")
title('Vertical Error')
grid on

figure(8)
plot(t, sqrt(error.standalone_enu(:, 1).^2 + error.standalone_enu(:, 2).^2))
xlabel("time(sec)")
ylabel("Error(m)")
title("Horizontal Error")
grid on

figure(9)
plot(t, error.standalone_enu(:, 1))
xlabel("time(sec)")
ylabel("Error(m)")
title("East Error")
grid on

figure(10)
plot(t, error.standalone_enu(:, 2))
xlabel("time(sec)")
ylabel("Error(m)")
title("North Error")
grid on
%% RMS error calculation
VRMS = 0;
HRMS = 0;
PRMS = 0;

VRMS = rms(error.standalone_enu(:,3));
HRMS = rms(sqrt(error.standalone_enu(:, 1).^2 + error.standalone_enu(:, 2).^2));
PRMS = rms(sqrt(error.standalone_enu(:, 1).^2 + error.standalone_enu(:, 2).^2 + error.standalone_enu(:,3).^2));

