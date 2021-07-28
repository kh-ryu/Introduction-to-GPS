clear all;
close all;
clc;

load(['EXP_20201006_1_data.mat'])

%% time setting
ini_epoch = 5*60;              % Data aquisition interval for determining integer ambiguity (5 mins)
n_epoch = size(user.GPSTime,1); % Overall data length
maskangle = 15;                 
nt = 32;

%% Visible satellite selection & find maximum elevation satllite
PRN = [1:32];
vis_sat = zeros(1,32);
ref_sat = 1; % index of reference(maximum elevation) satellite. set 1 as default.
El = 0; % mean elevation of visible satellite. set 0 as default.
El_max = 0; % to save maximum elevation

for ii = 1:32
    % find satellites that have larger elevation than mask angle with effectvie data
    temp_el1 = find(ref.El(:,ii) < maskangle| user.El(:,ii) < maskangle); 
    temp_el2 = find(ref.pr1(:,ii)<10 |user.pr1(:,ii)<10);
    if isempty(temp_el1)
        if isempty(temp_el2)
            vis_sat(ii)=1;
        end
    end
    
    % Calculate mean elevation
    El = mean(ref.El(:,ii));
    % Multiply vis_sat value - El will be 0 if not visible
    El = El*vis_sat(ii);
    
    % Compare with previous maximum value & update
    if El_max <= El
       ref_sat = ii;
       El_max = El;
    end
    
end
vis_sat(ref_sat) = 0; % Remove maximum elevation satellite
SV_vis = find(vis_sat==1);
n_vis = length(SV_vis);

%% Values for CDGPS Calculation
 phi_j_ru = zeros(n_epoch, n_vis);  % phase difference for satellite other than reference satellite
 phi_k_ru = zeros(n_epoch, 1);  % phase diffrence with reference satellite
 phi_jk_ru = zeros(n_epoch, n_vis); % double difference
 
 d_j_ru = zeros(n_epoch, n_vis); % distance between satellite other than reference satellite
 d_k_ru = zeros(n_epoch, 1); % distance with reference satellite
 d_jk_ru = zeros(n_epoch, n_vis); % double difference
 
 N_jk_ru = zeros(n_epoch, n_vis); % Integer ambiguity. Should round it
 
 %% N_jk_ru(Integer ambiguity) calculation
 
 for ti = 1:n_epoch
     % user position fixed -> user_xyz0
     
     phi_k_ru(ti) = ref.cp1(ti, ref_sat) - user.cp1(ti, ref_sat);   % k = reference satellite
     
     % svpos has vertical stack of x y z coordinates
     refsat_z = ref_sat*3;
     refsat_x = refsat_z - 2;
     d_k_ru(ti) = norm(ref.svpos(ti, refsat_x:refsat_z) - ref_xyz) - norm(user.svpos(ti, refsat_x:refsat_z) - user_xyz0); % difference between distance at (r, u)
     
     for jj = 1:n_vis
         phi_j_ru(ti,jj) = ref.cp1(ti, SV_vis(jj)) - user.cp1(ti, SV_vis(jj));
         phi_jk_ru(ti,jj) = phi_j_ru(ti, jj) - phi_k_ru(ti);
         
         SV_z = SV_vis(jj)*3;
         SV_x = SV_z - 2;
         d_j_ru(ti,jj) = norm(ref.svpos(ti, SV_x:SV_z) - ref_xyz) - norm(ref.svpos(ti, SV_x:SV_z) - user_xyz0);
         d_jk_ru(ti,jj) = d_j_ru(ti,jj) - d_k_ru(ti);
         N_jk_ru(ti,jj) = (phi_jk_ru(ti,jj) - d_jk_ru(ti,jj))/lambda_L1;
     end
 end
 
 % N value is obtained by taking average of calculated Ns 
 N = zeros(1, n_vis);
 
 for jj = 1:n_vis
     N(jj) = round(mean(N_jk_ru(1:ini_epoch, jj)));
 end
 
 N
%% Position calculation
 H = zeros(n_vis, 3);
 Z = zeros(n_vis, 1);
 
 user.pos_CDGPS = zeros(n_epoch, 3);
 user.pos_CDGPS_error = zeros(n_epoch, 3);
 user.pos_CDGPS_enu = zeros(n_epoch, 3);
 
 % Find position by iteration
 for ti = 1: n_epoch
    dx = 100;
     
    % Step 1. Suppose the user postion as O, then iterate
    R_u = zeros(3,1).';
    R_j = zeros(3,1).';
    R_k = zeros(3,1).';
    x_old = R_u.';
    
    % Step 2. Find position by iteration.
    while (dx > 10^-4)
        H = zeros(n_vis,3);
        Z = zeros(n_vis,1);
        
        R_uk = user.svpos(ti, refsat_x:refsat_z) - R_u;
        e_hat_uk = R_uk/norm(R_uk);
        
        R_k = ref.svpos(ti, refsat_x:refsat_z);
        R_rk = R_k - ref_xyz;
        e_hat_rk = R_rk/norm(R_rk);
        
        % Step 2.1. H, Z matrix calculation
        for jj = 1: n_vis
             % index of each satellite
             SV_z = SV_vis(jj)*3;
             SV_x = SV_z - 2;
             
             % Calculate e hat vector
             R_uj = user.svpos(ti, SV_x:SV_z) - R_u;
             e_hat_uj = R_uj/norm(R_uj);
             
             R_j = ref.svpos(ti, SV_x:SV_z);
             R_rj = R_j - ref_xyz;
             e_hat_rj = R_rj/norm(R_rj);

             % Calculate H matrix entries
             H(jj,1:3) = e_hat_uj - e_hat_uk;
             
             % Calculate Z matrix entries
             Z(jj,1) = phi_jk_ru(ti,jj) - N(jj)*lambda_L1 - ((e_hat_rj - e_hat_uj)*R_j.' - (e_hat_rk - e_hat_uk)*R_k.' - (e_hat_rj - e_hat_rk)*ref_xyz.');
         end
         
        % Step 2.2. Calculate user position
        x = pinv(H)*Z;
        R_u = x(1:3).';
         
        % Step 2.3. Compare with previous position and update
        dx = norm(x - x_old);
        x_old = x;
         
    end
    
    % Step 3. Transform into enu coordinate
    user.pos_CDGPS(ti,:) = x.';
    user.pos_CDGPS_enu(ti,:) = ( Rtran*(user.pos_CDGPS(ti,:) - ref_xyz).' ).'; 
    if ti == 1
        user.pos_CDGPS_error(ti,:) = 0;
    else
        user.pos_CDGPS_error(ti,:) = user.pos_CDGPS_enu(ti,:) - true.enu(ti - 1,:);
    end
    
    %DOP
    H_enu = horzcat((Rtran*H')', -ones(n_vis,1));
    DOP = inv((H_enu)'*(H_enu));
    HDOP(ti) = sqrt(DOP(1,1)+DOP(2,2));
    VDOP(ti) = sqrt(DOP(3,3));
    PDOP(ti) = sqrt(DOP(1,1) + DOP(2,2) + DOP(3,3));
    TDOP(ti) = sqrt(DOP(4,4));
    GDOP(ti) = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3)+DOP(4,4));
 end
 
 %% Measured position vs real position
 figure(1)
 plot3(user.pos_CDGPS_enu(ini_epoch+1:end,1),user.pos_CDGPS_enu(ini_epoch+1:end,2),user.pos_CDGPS_enu(ini_epoch+1:end,3),'k.');
 hold on
 plot3(true.enu(:,1), true.enu(:,2), true.enu(:,3),'r.')

 title('CDGPS Trajectory')
 xlabel('East(meter)')
 ylabel('North(meter)')
 zlabel('Altitude(meter)')
 legend('CDGPS', 'True')
 grid on
  
 
 %% Plot result
figure(2)
plot(user.pos_CDGPS_enu((ini_epoch+1):end, 1), user.pos_CDGPS_enu((ini_epoch+1):end, 2), 'o');
hold on
plot(true.enu((ini_epoch+1):end, 1), true.enu((ini_epoch+1):end, 2));
title('East-North Trajectory')
xlabel('East(meter)')
ylabel('North(meter)')
legend('CDGPS', 'True position')
grid on

figure(3)
plot(0:1:(n_epoch-ini_epoch-1), user.pos_CDGPS_enu((ini_epoch+1):end, 3), 'x');
hold on
plot(0:1:(n_epoch-ini_epoch-1), true.enu(ini_epoch:end, 3));
title('Up Trajectory')
xlabel('time(sec)')
ylabel('Up(meter)')
legend('CDGPS', 'True position')
grid on

%% Tendency of error
figure(4)
plot(ini_epoch+1:n_epoch, user.pos_CDGPS_error(ini_epoch+1:end, 3))
xlabel("time(sec)")
ylabel("Error(m)")
title('Vertical Error')
grid on

figure(5)
plot(ini_epoch+1:n_epoch, sqrt(user.pos_CDGPS_error(ini_epoch+1:end, 1).^2 + user.pos_CDGPS_error(ini_epoch+1:end, 2).^2))
xlabel("time(sec)")
ylabel("Error(m)")
title("Horizontal Error")
grid on

figure(6)
plot(ini_epoch+1:n_epoch, user.pos_CDGPS_error(ini_epoch+1:end, 1))
xlabel("time(sec)")
ylabel("Error(m)")
title("East Error")
grid on

figure(7)
plot(ini_epoch+1:n_epoch, user.pos_CDGPS_error(ini_epoch+1:end, 2))
xlabel("time(sec)")
ylabel("Error(m)")
title("North Error")
grid on

%% RMS error calculation
VRMS = 0;
HRMS = 0;
PRMS = 0;

VRMS = rms(user.pos_CDGPS_error(:,3))
HRMS = rms(sqrt(user.pos_CDGPS_error(:, 1).^2 + user.pos_CDGPS_error(:, 2).^2))
PRMS = rms(sqrt(user.pos_CDGPS_error(:, 1).^2 + user.pos_CDGPS_error(:, 2).^2 + user.pos_CDGPS_error(:,3).^2))

