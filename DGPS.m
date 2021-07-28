clear all;
close all;
clc;

load(['EXP_20201006_1_data.mat'])
%load(['EXP_20190925_G2_data.mat'])

%% time & basic setting
ini_epoch = 10*60; 
n_epoch = size(true.GPSTime,1);
maskangle = 15;                 
nt = 32;

%% Visible satellite selection
PRN = [1:32];
vis_sat = zeros(1,32);

for ii = 1:32
    % mask angle보다 elevation angle이 높고, 측정치가 있는 위성 찾기
    temp_el1 = find(ref.El(:,ii) < maskangle| user.El(:,ii) < maskangle); 
    temp_el2 = find(ref.pr1(:,ii)<10 | user.pr1(:,ii)<10);
    if isempty(temp_el1)
        if isempty(temp_el2)
            vis_sat(ii)=1;
        end
    end
end
SV_vis = find(vis_sat==1);
n_vis = length(SV_vis);

%% DGPS
 user.pos_DGPS = zeros(n_epoch,3);
 user.pos_DGPS = zeros(n_epoch,3);
 user.pos_DGPS_enu = zeros(n_epoch,3);
 
 R_ur = zeros(n_epoch,3);
 del_rj = zeros(n_epoch, n_vis);
 d_rj = zeros(n_epoch, n_vis);
 
 HDOP = zeros(n_epoch,1);
 VDOP = zeros(n_epoch,1);
 PDOP = zeros(n_epoch,1);
 TDOP = zeros(n_epoch,1);
 GDOP = zeros(n_epoch,1);
 
 % Position calculation loop
 for ti = 1: n_epoch
     dx = 100;
     
     R_user = zeros(3,1).'; % Assume position is origin
     B_user = 0;              % Negelct time error
     x_old = [R_user';B_user];
     
     while(dx > 10^-4)
         H_DGPS = zeros(n_vis,4);
         H_enu = zeros(n_vis,4);
         z_DGPS = zeros(n_vis,1);
         
         % Calculation for H, Z matrix
         for jj = 1: n_vis
             % Visible satellite data index
             SV_z = SV_vis(jj)*3;
             SV_x = SV_z - 2;
             
             % Calculate e_hat
             R_su = user.svpos(ti, SV_x:SV_z) - R_user;
             e_hat = R_su/norm(R_su);
             
             % Calculate H matrix
             H_DGPS(jj,1:3) = e_hat;
             H_DGPS(jj,4) = -1;
             
             % Calculate correction factor
             d_rj(ti, jj) = norm(ref.svpos(ti, SV_x:SV_z) - ref_xyz);
             del_rj = ref.pr1(ti,SV_vis(jj)) - d_rj(ti, jj);

             % Build measurement error matrix with correction factor
             pr_rj = user.pr1(ti, SV_vis(jj)) - del_rj;
             z_DGPS(jj,1) = e_hat*user.svpos(ti, SV_x:SV_z)' - pr_rj;
         end
         
         % Calculate user position
         x = pinv(H_DGPS)*z_DGPS;
         R_user = x(1:3)';
         B_user = x(4);
         
         % Calculate DOP -> Convert to ENU!
         e_enu = zeros(3, n_vis);    % Assembly of e-values in enu
         
         % Convert H_DOP into H_enu
         for i = 1:n_vis
             e_enu(:, i) = H_DGPS(i, 1:3)';
             e_enu(:, i) = Rtran*e_enu(:, i);
             length = sqrt(sum(e_enu(:,i).*e_enu(:,i)));
        
            if length > 1
                e_enu(:,i) = e_enu(:,i)/length;
            end
          end
     
         H_enu = horzcat(e_enu', -ones(n_vis,1));

         DOP = inv(H_enu'*H_enu);
         HDOP(ti) = sqrt(DOP(1,1) + DOP(2,2));
         VDOP(ti) = sqrt(DOP(3,3));
         PDOP(ti) = sqrt(DOP(1,1) + DOP(2,2) + DOP(3,3));
         
         % Calculated updated amount
         dx = norm(x - x_old);
         x_old = x;
     end
     
     user.pos_DGPS(ti,:) = R_user;
     user.pos_DGPS_enu(ti,:) = (Rtran*(user.pos_DGPS(ti,:)-ref_xyz).').';
     error.DGPS(ti,:) = user.pos_DGPS(ti,:) - true.xyz(ti,:);
     error.DGPS_enu(ti,:) = (Rtran*(error.DGPS(ti,:)'))';
     
     DOP = inv((H_enu)'*(H_enu));
     HDOP(ti) = sqrt(DOP(1,1)+DOP(2,2));
     VDOP(ti) = sqrt(DOP(3,3));
     PDOP(ti) = sqrt(DOP(1,1) + DOP(2,2) + DOP(3,3));
     TDOP(ti) = sqrt(DOP(4,4));
     GDOP(ti) = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3)+DOP(4,4));
 end
 
 %% Measured position vs real position
 figure(1)
 plot3(user.pos_DGPS_enu(1:end,1),user.pos_DGPS_enu(1:end,2),user.pos_DGPS_enu(1:end,3),'k.');
 hold on
 plot3(true.enu(:,1), true.enu(:,2), true.enu(:,3),'r.')

 title('DGPS Trajectory')
 xlabel('East(meter)')
 ylabel('North(meter)')
 zlabel('Altitude(meter)')
 legend('DGPS', 'True')
 grid on
  
 %% Plot
figure(2)
plot(user.pos_DGPS_enu(1:n_epoch,1),user.pos_DGPS_enu(1:n_epoch,2), 'x');
hold on
plot(true.enu(:,1), true.enu(:,2))
title('East-North Trajectory')
xlabel('East(meter)')
ylabel('North(meter)')
legend('DGPS', 'True')
grid on

figure(3)
plot(true.GPSTime(1:n_epoch), user.pos_DGPS_enu(1:n_epoch,3),'x')
hold on
plot(true.GPSTime(1:n_epoch), true.enu(1:n_epoch,3))
title('Time-Vertical Trajectory')
xlabel('Time(second)')
ylabel('Vertical(meter)')
legend('DGPS', 'True')

%% Plot Dop

toffset = true.GPSTime(1);
t = true.GPSTime(1:n_epoch)-toffset;
%{
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
%}
%% Tendency of error
figure(7)
plot(t, error.DGPS_enu(:, 3))
xlabel("time(sec)")
ylabel("Error(m)")
title('Vertical Error')
grid on

figure(8)
plot(t, sqrt(error.DGPS_enu(:, 1).^2 + error.DGPS_enu(:, 2).^2))
xlabel("time(sec)")
ylabel("Error(m)")
title("Horizontal Error")
grid on

figure(9)
plot(t, error.DGPS_enu(:, 1))
xlabel("time(sec)")
ylabel("Error(m)")
title("East Error")
grid on

figure(10)
plot(t, error.DGPS_enu(:, 2))
xlabel("time(sec)")
ylabel("Error(m)")
title("North Error")
grid on

%% RMS error calculation
VRMS = 0;
HRMS = 0;
PRMS = 0;

VRMS = rms(error.DGPS_enu(:,3))
HRMS = rms(sqrt(error.DGPS_enu(:, 1).^2 + error.DGPS_enu(:, 2).^2))
PRMS = rms(sqrt(error.DGPS_enu(:, 1).^2 + error.DGPS_enu(:, 2).^2 + error.DGPS_enu(:,3).^2))
