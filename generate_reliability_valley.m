%% generate_reliability_valley.m
% Generates the two-panel Reliability Valley figure for main.tex
% Uses calibrated baseline parameters directly (no Dynare needed)
%
% Output: reliability_valley.png

clear; clc;

%% Calibrated baseline parameters
p_psi      = 2.0;
p_sigma_ren = 0.12;
p_mu       = 0.16;
p_rho_flex = 0.40;
p_delta_g  = 0.0125;
p_delta_b  = 0.030;
p_theta_ren = 0.30;    % baseline penetration

% Steady-state capital stocks (from steadystate file)
K_b_ss  = 0.19;
K_g_ss  = 1.14;
A_bat_ss = 1.0;

% phi_int calibrated so u = 0.97 at steady state
u_ss = 0.97;
Vol_ren_bar = p_theta_ren * (0.045/p_theta_ren) * p_sigma_ren;  % = 0.0054
s_f = (p_rho_flex - 1) / p_rho_flex;
F_ss = (p_mu * (A_bat_ss * K_b_ss)^s_f + (1-p_mu) * K_g_ss^s_f)^(1/s_f);
phi_int_val = -p_psi * F_ss / (Vol_ren_bar * log(1 - u_ss));

fprintf('Calibrated phi_int = %.4f\n', phi_int_val);
fprintf('Steady-state F = %.4f\n', F_ss);

%% Transition path parameters
T_trans = 100;           % 100 quarters = 25 years
t_trans = (1:T_trans)';

% S-curve: theta_ren from 15% to 50%
theta_start = 0.15;
theta_end   = 0.50;
k_steep = 0.08;          % steepness
t_mid   = 35;            % front-loaded (midpoint at Q35)
theta_path = theta_start + (theta_end - theta_start) ./ (1 + exp(-k_steep * (t_trans - t_mid)));

%% Simulate baseline transition
% K_ren_base calibrated at theta_start=0.15 (matching run_dynare.m convention)
% This gives Vol_ren = theta_start * K_ren_base * sigma_ren = 0.15*0.30*0.12 = 0.0054
K_ren_base = 0.045 / theta_start;   % = 0.30
grid_growth_base = 0.0135;
bat_growth_base  = 0.020;

K_b_t   = K_b_ss;
K_g_t   = K_g_ss;
A_bat_t = A_bat_ss;

u_path = zeros(T_trans, 1);

for qq = 1:T_trans
    K_ren_t   = K_ren_base * (theta_path(qq) / theta_start);
    Vol_ren_t = theta_path(qq) * K_ren_t * p_sigma_ren;

    F_t = (p_mu * (A_bat_t * K_b_t)^s_f + (1-p_mu) * K_g_t^s_f)^(1/s_f);
    z_t = p_psi * F_t / (phi_int_val * Vol_ren_t);
    u_t = 1 - exp(-z_t);
    u_path(qq) = u_t;

    % Reliability gap drives investment acceleration
    u_gap = max(0, 0.97 - u_t) / 0.97;

    grid_growth = grid_growth_base * (1 + 2*u_gap);
    I_grid_t = (p_delta_g + grid_growth) * K_g_t;
    K_g_t = (1 - p_delta_g) * K_g_t + I_grid_t;

    bat_growth = bat_growth_base * (1 + 3*u_gap);
    I_bat_t = (p_delta_b + bat_growth) * K_b_t;
    K_b_t = (1 - p_delta_b) * K_b_t + I_bat_t;

    A_bat_t = A_bat_t * (1 + 0.10 * u_gap);
end

u_pct = u_path * 100;

%% Valley characteristics
[min_u_val, min_u_idx] = min(u_pct);
below_97 = u_pct < 97;
first_below = find(below_97, 1, 'first');
last_below  = find(below_97, 1, 'last');
valley_dur  = sum(below_97);

fprintf('\nReliability Valley results:\n');
fprintf('  Nadir:    %.2f%% at Q%d (year %.1f)\n', min_u_val, min_u_idx, min_u_idx/4);
fprintf('  Duration: %d quarters (%.1f years) below 97%%\n', valley_dur, valley_dur/4);
fprintf('  First below target: Q%d\n', first_below);
fprintf('  Recovery: Q%d\n', last_below);

%% Two-panel figure
fig = figure('Position', [100 100 1100 440], 'Color', 'w');

%--- Left panel: S-curve renewable penetration ---
ax1 = subplot(1, 2, 1);
hold on;

% Shade the front-loaded acceleration region (before midpoint)
x_shade = [0, t_mid/4, t_mid/4, 0];
y_shade = [10, 10, 55, 55];
fill(x_shade, y_shade, [0.9 1.0 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

plot(t_trans/4, theta_path*100, 'Color', [0.13 0.55 0.13], 'LineWidth', 2.8);

% PDP8 target
yline(50, '--', 'LineWidth', 1.6, 'Color', [0.8 0 0]);
text(0.8, 51.2, 'PDP8 target (50%)', 'Color', [0.8 0 0], 'FontSize', 8.5);

% Baseline penetration
yline(15, ':', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
text(0.8, 13.5, 'Baseline (15%)', 'Color', [0.5 0.5 0.5], 'FontSize', 8.5);

% Inflection point annotation
[~, mid_q] = min(abs(theta_path - (theta_start + theta_end)/2));
plot(mid_q/4, theta_path(mid_q)*100, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', [0.13 0.55 0.13]);
text(mid_q/4 + 1.0, theta_path(mid_q)*100 - 2.5, ...
    sprintf('Inflection\n(Q%d)', mid_q), 'FontSize', 8, 'Color', 'k');

% Annotation: "Front-loaded"
text(3, 37, {'Front-loaded';'acceleration'}, 'FontSize', 8.5, ...
    'Color', [0.0 0.4 0.0], 'FontAngle', 'italic');

hold off;
grid on; box on;
xlabel('Years from transition start', 'FontSize', 10);
ylabel('Renewable penetration (theta_{ren}, %)', 'FontSize', 10);
title('(a)\ VRE Deployment: S-Curve Path', 'FontSize', 11, 'FontWeight', 'bold', ...
    'Interpreter', 'none');
xlim([0 T_trans/4]);
ylim([10 55]);
set(ax1, 'FontSize', 9);

%--- Right panel: Reliability valley ---
ax2 = subplot(1, 2, 2);
hold on;

% Red shaded valley region
if ~isempty(first_below) && ~isempty(last_below)
    x_vfill = [(first_below-1)/4, last_below/4, last_below/4, (first_below-1)/4];
    y_vfill = [min_u_val - 2, min_u_val - 2, 97, 97];
    fill(x_vfill, y_vfill, [1.0 0.85 0.80], 'EdgeColor', 'none', 'FaceAlpha', 0.65);
end

% Reliability path
plot(t_trans/4, u_pct, 'b-', 'LineWidth', 2.8, 'DisplayName', 'Grid reliability');

% Target line
yline(97, '--', 'LineWidth', 1.6, 'Color', [0.8 0 0]);
text(0.8, 97.5, 'Target (97%)', 'Color', [0.8 0 0], 'FontSize', 8.5);

% Nadir marker
plot(min_u_idx/4, min_u_val, 'v', 'MarkerSize', 10, ...
    'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0.6 0 0]);
text(min_u_idx/4 + 0.8, min_u_val - 0.5, ...
    sprintf('Nadir: Q%d\n(u = %.1f%%)', min_u_idx, min_u_val), ...
    'FontSize', 8.5, 'Color', [0.7 0 0]);

% Valley label
if ~isempty(first_below) && ~isempty(last_below)
    valley_center = (first_below + last_below) / 2;
    text(valley_center/4, min_u_val + 1.2, sprintf('Reliability Valley\n(%d quarters)', valley_dur), ...
        'FontSize', 9, 'Color', [0.75 0 0], 'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');
end

% Recovery annotation
if ~isempty(last_below)
    plot([last_below/4, last_below/4], [min_u_val - 2, 97.5], ...
        'k:', 'LineWidth', 1.0);
    text(last_below/4 + 0.5, 97.9, 'Recovery', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end

hold off;
grid on; box on;
xlabel('Years from transition start', 'FontSize', 10);
ylabel('Grid reliability (u, %)', 'FontSize', 10);
title('(b)\ The Reliability Valley: Baseline Transition', 'FontSize', 11, ...
    'FontWeight', 'bold', 'Interpreter', 'none');
xlim([0 T_trans/4]);
ylim([min_u_val - 3, 100]);
set(ax2, 'FontSize', 9);

%% Save
print('reliability_valley.png', '-dpng', '-r300');
fprintf('\nFigure saved: reliability_valley.png\n');
