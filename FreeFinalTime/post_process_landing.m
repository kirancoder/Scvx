% post_process_landing.m
% Analyzes and exports individual, full-size plots for 6-DoF Rocket Landing

clc; clear; close all;

% --- Setup Output Directory ---
output_dir = 'outputs';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
fprintf('Saving individual plots to: %s\n', output_dir);

% --- 1. Load Data ---
% Load MAT file (for parameters and 3D history)
load('LandingTrajectoryData.mat', 'all_X', 'all_U', 'm', 'K');
num_iters = length(all_X);

% FIX: Restore 'm' to physical units (kg, N) so plot limits match CSV data!
m.redimensionalize();

% Load CSV file (for 2D plotting)
csv_filename = fullfile(output_dir, 'Final_Converged_Trajectory.csv');
if ~isfile(csv_filename)
    error('CSV file not found! Ensure Final_Converged_Trajectory.csv is in the outputs folder.');
end
csv_data = readtable(csv_filename);
fprintf('Successfully loaded CSV data.\n');

% Extract CSV columns
t_vec = csv_data.Time_s;
K_csv = length(t_vec);
mass = csv_data.Mass_kg;
thrust_mag = csv_data.Thrust_Magnitude_N;

% --- 2. Calculate Constraints ---
gimbal_angle_deg = zeros(K_csv, 1);
tilt_angle_deg = zeros(K_csv, 1);
glideslope_margin = zeros(K_csv, 1);
omega_mag_deg = zeros(K_csv, 1);
vel_norm = zeros(K_csv, 1);

for k = 1:K_csv
    U_vec = [csv_data.Control_Ux_N(k); csv_data.Control_Uy_N(k); csv_data.Control_Uz_N(k)];
    if thrust_mag(k) > 1e-6
        gimbal_angle_deg(k) = rad2deg(acos(U_vec(3) / thrust_mag(k)));
    end
    
    q = [csv_data.Quat_w(k); csv_data.Quat_x(k); csv_data.Quat_y(k); csv_data.Quat_z(k)];
    CBI = m.dir_cosine(q);
    body_z_in_inertial = CBI' * [0; 0; 1];
    tilt_angle_deg(k) = rad2deg(acos(max(min(body_z_in_inertial(3), 1), -1)));
    
    r_xy_norm = norm([csv_data.Pos_X_East_m(k); csv_data.Pos_Y_North_m(k)]);
    glideslope_margin(k) = csv_data.Pos_Z_Up_m(k) - (r_xy_norm / m.tan_gamma_gs);
    
    omega_vec = [csv_data.Omega_X_radps(k); csv_data.Omega_Y_radps(k); csv_data.Omega_Z_radps(k)];
    omega_mag_deg(k) = rad2deg(norm(omega_vec));
    
    vel_vec = [csv_data.Vel_X_mps(k); csv_data.Vel_Y_mps(k); csv_data.Vel_Z_mps(k)];
    vel_norm(k) = norm(vel_vec);
end

% Keep windows hidden while generating to prevent screen clutter
fig_vis = 'off';

% --- Plot 1: Thrust ---
fig1 = figure('Name', 'Thrust', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
hold on; grid on;
yline(m.T_max / 1000, 'r--', 'Upper Limit (T_{max})', 'LineWidth', 2);
yline(m.T_min / 1000, 'r--', 'Lower Limit (T_{min})', 'LineWidth', 2);
plot(t_vec, thrust_mag / 1000, 'b', 'LineWidth', 4); % Extra thick line
title('Thrust Magnitude Profile');
xlabel('Time [s]'); ylabel('Thrust [kN]');
% Pad Y-limits by 50 kN so the data line never touches the edge
ylim([(m.T_min/1000) - 50, (m.T_max/1000) + 50]);
exportgraphics(fig1, fullfile(output_dir, '01_Thrust_Profile.png'), 'Resolution', 300);

% --- Plot 2: Gimbal Angle ---
fig2 = figure('Name', 'Gimbal', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
hold on; grid on;
yline(m.max_gimbal, 'r--', 'Upper Limit (\delta_{max})', 'LineWidth', 2);
yline(0, 'r--', 'Lower Limit (0)', 'LineWidth', 2);
plot(t_vec, gimbal_angle_deg, 'b', 'LineWidth', 3);
title('Gimbal Deflection Angle');
xlabel('Time [s]'); ylabel('Angle [deg]');
ylim([-2, m.max_gimbal + 2]);
exportgraphics(fig2, fullfile(output_dir, '02_Gimbal_Angle.png'), 'Resolution', 300);

% --- Plot 3: Tilt Angle ---
fig3 = figure('Name', 'Tilt', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
hold on; grid on;
yline(m.max_angle, 'r--', 'Upper Limit (\theta_{max})', 'LineWidth', 2);
yline(0, 'r--', 'Lower Limit (0)', 'LineWidth', 2);
plot(t_vec, tilt_angle_deg, 'b', 'LineWidth', 3);
title('Vehicle Tilt Angle');
xlabel('Time [s]'); ylabel('Angle [deg]');
ylim([-5, m.max_angle + 10]);
exportgraphics(fig3, fullfile(output_dir, '03_Tilt_Angle.png'), 'Resolution', 300);

% --- Plot 4: Glideslope Margin ---
fig4 = figure('Name', 'Glideslope', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
hold on; grid on;
yline(0, 'r--', 'Lower Limit (Cone Edge)', 'LineWidth', 2);
plot(t_vec, glideslope_margin, 'b', 'LineWidth', 3);
title('Glideslope Margin (Positive = Inside Cone)');
xlabel('Time [s]'); ylabel('Margin [m]');
g_min = min(glideslope_margin); g_max = max(glideslope_margin);
ylim([min(g_min - 20, -10), g_max + 20]);
exportgraphics(fig4, fullfile(output_dir, '04_Glideslope_Margin.png'), 'Resolution', 300);

% --- Plot 5: Angular Velocity ---
fig5 = figure('Name', 'Omega', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
hold on; grid on;
yline(rad2deg(m.w_B_max), 'r--', 'Upper Limit (\omega_{max})', 'LineWidth', 2);
yline(0, 'r--', 'Lower Limit (0)', 'LineWidth', 2);
plot(t_vec, omega_mag_deg, 'b', 'LineWidth', 3);
title('Angular Velocity Magnitude');
xlabel('Time [s]'); ylabel('\omega [deg/s]');
ylim([-5, rad2deg(m.w_B_max) + 10]);
exportgraphics(fig5, fullfile(output_dir, '05_Angular_Velocity.png'), 'Resolution', 300);

% --- Plot 6: Mass ---
fig6 = figure('Name', 'Mass', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
hold on; grid on;
yline(m.m_wet, 'r--', 'Upper Limit (Wet Mass)', 'LineWidth', 2);
yline(m.m_dry, 'r--', 'Lower Limit (Dry Mass)', 'LineWidth', 2);
plot(t_vec, mass, 'b', 'LineWidth', 4); % Extra thick line
title('Vehicle Mass Depletion');
xlabel('Time [s]'); ylabel('Mass [kg]');
% Pad Y-limits by 1000 kg so flat lines don't get squashed
ylim([m.m_dry - 1000, m.m_wet + 1000]);
exportgraphics(fig6, fullfile(output_dir, '06_Mass_Profile.png'), 'Resolution', 300);

% --- Plot 7: Position ---
fig7 = figure('Name', 'Position', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
plot(t_vec, csv_data.Pos_X_East_m, 'r', 'LineWidth', 2); hold on; grid on;
plot(t_vec, csv_data.Pos_Y_North_m, 'g', 'LineWidth', 2);
plot(t_vec, csv_data.Pos_Z_Up_m, 'b', 'LineWidth', 2);
title('Position Profiles');
xlabel('Time [s]'); ylabel('Position [m]');
legend('X (East)', 'Y (North)', 'Z (Altitude)', 'Location', 'best');
exportgraphics(fig7, fullfile(output_dir, '07_Position.png'), 'Resolution', 300);

% --- Plot 8: Velocity ---
fig8 = figure('Name', 'Velocity', 'Position', [100, 100, 800, 500], 'Visible', fig_vis);
plot(t_vec, csv_data.Vel_X_mps, 'r--', 'LineWidth', 2); hold on; grid on;
plot(t_vec, csv_data.Vel_Y_mps, 'g--', 'LineWidth', 2);
plot(t_vec, csv_data.Vel_Z_mps, 'b--', 'LineWidth', 2);
plot(t_vec, vel_norm, 'k', 'LineWidth', 3);
title('Velocity Profiles');
xlabel('Time [s]'); ylabel('Velocity [m/s]');
legend('V_x (East)', 'V_y (North)', 'V_z (Up)', '||V|| (Magnitude)', 'Location', 'best');
exportgraphics(fig8, fullfile(output_dir, '08_Velocity.png'), 'Resolution', 300);

% --- Plot 9: 3D Trajectory ---
fig9 = figure('Name', '3D Trajectory', 'Position', [100, 100, 900, 700], 'Visible', fig_vis);
start_iter = max(1, num_iters - 10);
colors = jet(num_iters - start_iter + 1);
hold on; grid on; view(3);
xlabel('X (East) [m]'); ylabel('Y (North) [m]'); zlabel('Z (Up) [m]');
title('3D Landing Trajectory Progression');

% Plot target pad
th = linspace(0, 2*pi, 50);
fill3(20*cos(th), 20*sin(th), zeros(1,50), [0.8 0.8 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.5);

c_idx = 1;
for i = start_iter:num_iters
    X_iter = all_X{i};
    if i == num_iters
        plot3(X_iter(2, :), X_iter(3, :), X_iter(4, :), 'Color', 'k', 'LineWidth', 3, 'DisplayName', 'Final Converged');
    else
        plot3(X_iter(2, :), X_iter(3, :), X_iter(4, :), 'Color', colors(c_idx, :), 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', sprintf('Iter %d', i-1));
    end
    c_idx = c_idx + 1;
end
legend('Target Pad', 'Location', 'best');
axis equal;
exportgraphics(fig9, fullfile(output_dir, '09_Trajectory_3D.png'), 'Resolution', 300);

fprintf('\nDone! Check the "%s" folder for 9 separate high-resolution plots.\n', output_dir);