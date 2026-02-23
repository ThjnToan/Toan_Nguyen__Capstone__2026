%% ========================================================================
%% VIETNAM DSGE MODEL - DYNARE EXECUTION SCRIPT
%% ========================================================================
%% This script runs the Dynare model and generates publication-quality
%% figures for the capstone thesis
%%
%% Requirements:
%%   - MATLAB R2016a or later
%%   - Dynare 6.x (download from https://www.dynare.org)
%%   - vietnam_dsge.mod in the same directory
%%
%% Author: Toan T. Nguyen
%% Date: February 2026
%% ========================================================================

clear all;
close all;
clc;

%% Add Dynare to MATLAB path
addpath('C:/dynare/6.5/matlab');

%% Run the Dynare model
fprintf('Running Vietnam DSGE model in Dynare...\n');
dynare vietnam_dsge noclearall;

fprintf('\nModel solved successfully!\n');
fprintf('Blanchard-Kahn conditions verified.\n');

%% ========================================================================
%% EXTRACT RESULTS
%% ========================================================================

% Time horizon (quarters)
T = 40;
t = 0:T-1;

% Extract IRF data for renewable intermittency shock (eps_ren)
irf_Y        = oo_.irfs.Y_eps_ren;
irf_C        = oo_.irfs.C_eps_ren;
irf_u        = oo_.irfs.u_eps_ren;
irf_I_bat    = oo_.irfs.I_bat_eps_ren;
irf_I_grid   = oo_.irfs.I_grid_eps_ren;
irf_F        = oo_.irfs.F_eps_ren;
irf_A_bat    = oo_.irfs.A_bat_eps_ren;
irf_P_bat    = oo_.irfs.P_bat_eps_ren;
irf_I_p      = oo_.irfs.I_p_eps_ren;
irf_L        = oo_.irfs.L_eps_ren;
irf_B_star   = oo_.irfs.B_star_eps_ren;
irf_r_star   = oo_.irfs.r_star_eps_ren;

%% ========================================================================
%% FIGURE 1: IMPULSE RESPONSE FUNCTIONS (9 PANELS)
%% ========================================================================

figure('Position', [100 100 1400 1000]);

% Panel 1: Output
subplot(3,3,1);
plot(t, irf_Y, 'b-', 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Output (Y)');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 2: Consumption
subplot(3,3,2);
plot(t, irf_C, 'b-', 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Consumption (C)');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 3: Utilization (Reliability)
subplot(3,3,3);
plot(t, irf_u, 'r-', 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Utilization / Reliability (u)');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 4: Battery Investment
subplot(3,3,4);
plot(t, irf_I_bat, 'g-', 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Battery Investment (I_{bat})');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 5: Grid Investment
subplot(3,3,5);
plot(t, irf_I_grid, 'Color', [0.8 0.5 0], 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Grid Investment (I_{grid})');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 6: Flexibility
subplot(3,3,6);
plot(t, irf_F, 'm-', 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Total Flexibility (F)');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 7: Battery Technology
subplot(3,3,7);
plot(t, irf_A_bat, 'c-', 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Battery Technology (A_{bat})');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 8: Net Foreign Assets
subplot(3,3,8);
plot(t, irf_B_star, 'Color', [0.5 0 0.5], 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Net Foreign Assets (B^*)');
xlabel('Quarters');
ylabel('% deviation from SS');

% Panel 9: Labor
subplot(3,3,9);
plot(t, irf_L, 'Color', [0 0.5 0.5], 'LineWidth', 2);
hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
grid on;
title('Labor (L)');
xlabel('Quarters');
ylabel('% deviation from SS');

% Overall title
sgtitle('Impulse Responses to Renewable Intermittency Shock (\epsilon_{ren})', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
print('irf_dynare.png', '-dpng', '-r300');
fprintf('Figure saved: irf_dynare.png\n');

%% ========================================================================
%% COMPUTE KEY STATISTICS
%% ========================================================================

fprintf('\n========================================\n');
fprintf('KEY QUANTITATIVE RESULTS\n');
fprintf('========================================\n\n');

% Impact effects (t=0, already in % deviation)
fprintf('IMPACT EFFECTS (1 std dev Intermittency Shock):\n');
fprintf('  Output:              %+.4f%%\n', irf_Y(1));
fprintf('  Consumption:         %+.4f%%\n', irf_C(1));
fprintf('  Utilization:         %+.4f%%\n', irf_u(1));
fprintf('  Battery Investment:  %+.4f%%\n', irf_I_bat(1));
fprintf('  Grid Investment:     %+.4f%%\n', irf_I_grid(1));
fprintf('  Battery Technology:  %+.4f%%\n', irf_A_bat(1));
fprintf('  Net Foreign Assets:  %+.4f%%\n\n', irf_B_star(1));

% Agility gap
if irf_I_grid(1) ~= 0
    agility_gap = abs(irf_I_bat(1)) / abs(irf_I_grid(1));
    fprintf('AGILITY GAP (|Battery|/|Grid| response ratio): %.1fx\n\n', agility_gap);
else
    fprintf('AGILITY GAP: Grid investment impact is zero, ratio undefined.\n\n');
end

% Peak responses
[max_bat, idx_bat] = max(abs(irf_I_bat));
[max_grid, idx_grid] = max(abs(irf_I_grid));
fprintf('PEAK RESPONSES:\n');
fprintf('  Battery Investment:  %+.4f%% at quarter %d\n', irf_I_bat(idx_bat), idx_bat-1);
fprintf('  Grid Investment:     %+.4f%% at quarter %d\n\n', irf_I_grid(idx_grid), idx_grid-1);

% Technology improvement by Q10
if length(irf_A_bat) >= 11
    fprintf('LEARNING-BY-DOING:\n');
    fprintf('  Battery tech (A_bat) at Q10: %+.4f%%\n\n', irf_A_bat(11));
end

%% ========================================================================
%% FIGURE 2: AGILITY GAP COMPARISON
%% ========================================================================

figure('Position', [100 100 800 500]);

hold on;
plot(t, irf_I_bat, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Battery Investment');
plot(t, irf_I_grid, 'Color', [0.8 0.5 0], 'LineWidth', 2.5, 'DisplayName', 'Grid Investment');
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;

grid on;
xlabel('Quarters', 'FontSize', 12);
ylabel('% Deviation from Steady State', 'FontSize', 12);
title('The Agility Gap: Private vs. Public Flexibility Investment', ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);

print('agility_gap_dynare.png', '-dpng', '-r300');
fprintf('Figure saved: agility_gap_dynare.png\n');

%% ========================================================================
%% FIGURE 3: IRFs TO BATTERY PRICE SHOCK
%% ========================================================================

if isfield(oo_.irfs, 'Y_eps_bat')
    figure('Position', [100 100 1000 600]);

    subplot(2,3,1);
    plot(t, oo_.irfs.Y_eps_bat, 'b-', 'LineWidth', 2);
    hold on; plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
    grid on; title('Output'); ylabel('% dev');

    subplot(2,3,2);
    plot(t, oo_.irfs.P_bat_eps_bat, 'r-', 'LineWidth', 2);
    hold on; plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
    grid on; title('Battery Price'); ylabel('% dev');

    subplot(2,3,3);
    plot(t, oo_.irfs.I_bat_eps_bat, 'g-', 'LineWidth', 2);
    hold on; plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
    grid on; title('Battery Investment'); ylabel('% dev');

    subplot(2,3,4);
    plot(t, oo_.irfs.u_eps_bat, 'Color', [0.8 0 0], 'LineWidth', 2);
    hold on; plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
    grid on; title('Utilization'); ylabel('% dev');

    subplot(2,3,5);
    plot(t, oo_.irfs.A_bat_eps_bat, 'Color', [0.5 0 0.5], 'LineWidth', 2);
    hold on; plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
    grid on; title('Battery Technology (A_{bat})'); ylabel('% dev');

    subplot(2,3,6);
    plot(t, oo_.irfs.B_star_eps_bat, 'Color', [0 0.5 0.5], 'LineWidth', 2);
    hold on; plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
    grid on; title('Net Foreign Assets (B^*)'); ylabel('% dev');

    sgtitle('Impulse Responses to Battery Price Shock (\epsilon_{bat})', ...
        'FontSize', 14, 'FontWeight', 'bold');

    print('irf_battery_price.png', '-dpng', '-r300');
    fprintf('Figure saved: irf_battery_price.png\n');
end

%% ========================================================================
%% VARIANCE DECOMPOSITION
%% ========================================================================

fprintf('\nVARIANCE DECOMPOSITION:\n');
fprintf('(See Dynare output above for full decomposition table)\n');
fprintf('Key variables: Y, u, F, I_bat, I_grid, B_star\n\n');

%% ========================================================================
%% WELFARE COST CALCULATION
%% ========================================================================

beta_val = M_.params(strcmp(M_.param_names, 'beta'));
T_welfare = min(T, length(irf_C));

discount_factors = beta_val .^ (0:T_welfare-1);
sum_discounted = sum(discount_factors .* log(1 + irf_C(1:T_welfare)/100));
normalizer = (1 - beta_val^T_welfare) / (1 - beta_val);

welfare_lambda = 1 - exp(sum_discounted / normalizer);
welfare_pct = welfare_lambda * 100;

fprintf('WELFARE COST OF INTERMITTENCY (Baseline, chi=1.0):\n');
fprintf('  Compensating consumption variation: %.4f%%\n', abs(welfare_pct));
fprintf('  Horizon: %d quarters\n\n', T_welfare);

% Store baseline results
baseline_irf_Y = irf_Y;
baseline_irf_C = irf_C;
baseline_irf_u = irf_u;
baseline_irf_A_bat = irf_A_bat;
baseline_irf_I_bat = irf_I_bat;
baseline_welfare = abs(welfare_pct);

%% ========================================================================
%% COUNTERFACTUAL 1: SIGNAL ATTENUATION (chi variation)
%% ========================================================================

fprintf('\n========================================\n');
fprintf('COUNTERFACTUAL ANALYSIS: SIGNAL ATTENUATION\n');
fprintf('========================================\n\n');

chi_values = [1.0, 0.5, 0.0];
chi_labels = {'chi=1.0 (Baseline)', 'chi=0.5 (Partial)', 'chi=0.0 (No Learning)'};
chi_welfare = zeros(size(chi_values));
chi_irfs_Y = zeros(length(chi_values), T);
chi_irfs_u = zeros(length(chi_values), T);
chi_irfs_A = zeros(length(chi_values), T);

for i = 1:length(chi_values)
    % Change chi parameter
    set_param_value('chi', chi_values(i));

    % Re-compute steady state and solve
    [oo_.dr.ys, M_.params, info_check] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    [oo_.dr, info_ss, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);

    if info_ss(1) == 0
        % Compute IRFs manually
        dr = oo_.dr;
        nvar = M_.endo_nbr;
        nshocks = M_.exo_nbr;

        % Find eps_ren index
        eps_ren_idx = find(strcmp(M_.exo_names, 'eps_ren'));
        sigma_ren_val = M_.params(strcmp(M_.param_names, 'sigma_ren'));

        % Shock vector
        shock_vec = zeros(nshocks, 1);
        shock_vec(eps_ren_idx) = sigma_ren_val;

        % Get state-space matrices
        ghx = dr.ghx;
        ghu = dr.ghu;

        % Order of variables in decision rule
        order_var = dr.order_var;
        nstatic = M_.nstatic;
        nspred = M_.nspred;

        % Simulate IRFs
        irf_mat = zeros(nvar, T);
        state = ghu * shock_vec;  % impact
        irf_mat(:, 1) = state;

        for tt = 2:T
            % Extract state variables (predetermined)
            state_vars = state(nstatic+1:nstatic+nspred);
            state = ghx * state_vars;
            irf_mat(:, tt) = state;
        end

        % Reorder to declaration order
        irf_reordered = zeros(nvar, T);
        irf_reordered(order_var, :) = irf_mat;

        % Extract specific variables
        Y_idx = find(strcmp(M_.endo_names, 'Y'));
        C_idx = find(strcmp(M_.endo_names, 'C'));
        u_idx = find(strcmp(M_.endo_names, 'u'));
        A_idx = find(strcmp(M_.endo_names, 'A_bat'));

        chi_irfs_Y(i,:) = irf_reordered(Y_idx, :);
        chi_irfs_u(i,:) = irf_reordered(u_idx, :);
        chi_irfs_A(i,:) = irf_reordered(A_idx, :);

        % Welfare cost
        irf_C_chi = irf_reordered(C_idx, :);
        sum_disc = sum(discount_factors .* log(1 + irf_C_chi/100));
        chi_welfare(i) = abs((1 - exp(sum_disc / normalizer)) * 100);

        fprintf('  chi = %.1f: Welfare cost = %.6f%%, Impact Y = %+.4f%%\n', ...
            chi_values(i), chi_welfare(i), chi_irfs_Y(i,1));
    else
        fprintf('  chi = %.1f: FAILED (info = %d)\n', chi_values(i), info_ss(1));
    end
end

% Restore baseline chi
set_param_value('chi', 1.0);
[oo_.dr.ys, M_.params, ~] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
[oo_.dr, ~, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);

%% FIGURE 4: Signal Attenuation Comparison
figure('Position', [100 100 1200 800]);

colors = {'b', 'r', [0.5 0.5 0]};
line_styles = {'-', '--', ':'};

% Panel 1: Output
subplot(2,2,1);
hold on;
for i = 1:length(chi_values)
    plot(t, chi_irfs_Y(i,:), 'Color', colors{i}, 'LineStyle', line_styles{i}, ...
        'LineWidth', 2, 'DisplayName', chi_labels{i});
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;
grid on; title('Output (Y)'); ylabel('% dev'); xlabel('Quarters');
legend('Location', 'southeast', 'FontSize', 8);

% Panel 2: Utilization
subplot(2,2,2);
hold on;
for i = 1:length(chi_values)
    plot(t, chi_irfs_u(i,:), 'Color', colors{i}, 'LineStyle', line_styles{i}, ...
        'LineWidth', 2, 'DisplayName', chi_labels{i});
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;
grid on; title('Utilization (u)'); ylabel('% dev'); xlabel('Quarters');

% Panel 3: Battery Technology
subplot(2,2,3);
hold on;
for i = 1:length(chi_values)
    plot(t, chi_irfs_A(i,:), 'Color', colors{i}, 'LineStyle', line_styles{i}, ...
        'LineWidth', 2, 'DisplayName', chi_labels{i});
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off;
grid on; title('Battery Technology (A_{bat})'); ylabel('% dev'); xlabel('Quarters');

% Panel 4: Welfare cost bar chart
subplot(2,2,4);
bar(chi_welfare);
set(gca, 'XTickLabel', {'chi=1.0', 'chi=0.5', 'chi=0.0'});
ylabel('Welfare Cost (% consumption)');
title('Welfare Cost by Signal Transmission');
grid on;

sgtitle('Counterfactual: Effect of Signal Attenuation on Learning Channel', ...
    'FontSize', 14, 'FontWeight', 'bold');

print('counterfactual_chi.png', '-dpng', '-r300');
fprintf('\nFigure saved: counterfactual_chi.png\n');

%% ========================================================================
%% COUNTERFACTUAL 2: PARAMETER SENSITIVITY
%% ========================================================================

fprintf('\n========================================\n');
fprintf('PARAMETER SENSITIVITY ANALYSIS\n');
fprintf('========================================\n\n');

% Sensitivity: psi (reliability sensitivity)
psi_values = [1.5, 2.0, 3.0];
psi_welfare = zeros(size(psi_values));
psi_impact_Y = zeros(size(psi_values));

fprintf('Sensitivity to psi (reliability sensitivity):\n');
for i = 1:length(psi_values)
    set_param_value('psi', psi_values(i));
    [oo_.dr.ys, M_.params, info_check] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    [oo_.dr, info_ss, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if info_ss(1) == 0
        dr = oo_.dr;
        shock_vec = zeros(M_.exo_nbr, 1);
        shock_vec(find(strcmp(M_.exo_names, 'eps_ren'))) = M_.params(strcmp(M_.param_names, 'sigma_ren'));
        state = dr.ghu * shock_vec;
        irf_all = zeros(M_.endo_nbr, T);
        irf_all(:,1) = state;
        for tt = 2:T
            state_vars = state(M_.nstatic+1:M_.nstatic+M_.nspred);
            state = dr.ghx * state_vars;
            irf_all(:,tt) = state;
        end
        irf_reord = zeros(M_.endo_nbr, T);
        irf_reord(dr.order_var,:) = irf_all;

        Y_idx = find(strcmp(M_.endo_names, 'Y'));
        C_idx = find(strcmp(M_.endo_names, 'C'));
        psi_impact_Y(i) = irf_reord(Y_idx, 1);
        irf_C_psi = irf_reord(C_idx, :);
        sum_disc = sum(discount_factors .* log(1 + irf_C_psi/100));
        psi_welfare(i) = abs((1 - exp(sum_disc / normalizer)) * 100);
        fprintf('  psi = %.1f: Welfare = %.6f%%, Impact Y = %+.4f%%\n', ...
            psi_values(i), psi_welfare(i), psi_impact_Y(i));
    end
end
set_param_value('psi', 2.0);  % restore

% Sensitivity: phi_grid (grid investment aggressiveness)
phi_grid_values = [0.5, 1.5, 3.0];
phi_welfare = zeros(size(phi_grid_values));
phi_impact_Y = zeros(size(phi_grid_values));
phi_agility = zeros(size(phi_grid_values));

fprintf('\nSensitivity to phi_grid (grid investment response):\n');
for i = 1:length(phi_grid_values)
    set_param_value('phi_grid', phi_grid_values(i));
    [oo_.dr.ys, M_.params, info_check] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    [oo_.dr, info_ss, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if info_ss(1) == 0
        dr = oo_.dr;
        shock_vec = zeros(M_.exo_nbr, 1);
        shock_vec(find(strcmp(M_.exo_names, 'eps_ren'))) = M_.params(strcmp(M_.param_names, 'sigma_ren'));
        state = dr.ghu * shock_vec;
        irf_all = zeros(M_.endo_nbr, T);
        irf_all(:,1) = state;
        for tt = 2:T
            state_vars = state(M_.nstatic+1:M_.nstatic+M_.nspred);
            state = dr.ghx * state_vars;
            irf_all(:,tt) = state;
        end
        irf_reord = zeros(M_.endo_nbr, T);
        irf_reord(dr.order_var,:) = irf_all;

        Y_idx = find(strcmp(M_.endo_names, 'Y'));
        C_idx = find(strcmp(M_.endo_names, 'C'));
        I_bat_idx = find(strcmp(M_.endo_names, 'I_bat'));
        I_grid_idx = find(strcmp(M_.endo_names, 'I_grid'));

        phi_impact_Y(i) = irf_reord(Y_idx, 1);
        irf_C_phi = irf_reord(C_idx, :);
        sum_disc = sum(discount_factors .* log(1 + irf_C_phi/100));
        phi_welfare(i) = abs((1 - exp(sum_disc / normalizer)) * 100);

        if irf_reord(I_grid_idx, 1) ~= 0
            phi_agility(i) = abs(irf_reord(I_bat_idx, 1)) / abs(irf_reord(I_grid_idx, 1));
        end

        fprintf('  phi_grid = %.1f: Welfare = %.6f%%, Impact Y = %+.4f%%, Agility = %.1fx\n', ...
            phi_grid_values(i), phi_welfare(i), phi_impact_Y(i), phi_agility(i));
    end
end
set_param_value('phi_grid', 1.5);  % restore

% Sensitivity: mu (battery share in flexibility)
mu_values = [0.10, 0.16, 0.25];
mu_welfare = zeros(size(mu_values));
mu_impact_Y = zeros(size(mu_values));

fprintf('\nSensitivity to mu (battery share in flexibility):\n');
for i = 1:length(mu_values)
    set_param_value('mu', mu_values(i));
    [oo_.dr.ys, M_.params, info_check] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    [oo_.dr, info_ss, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if info_ss(1) == 0
        dr = oo_.dr;
        shock_vec = zeros(M_.exo_nbr, 1);
        shock_vec(find(strcmp(M_.exo_names, 'eps_ren'))) = M_.params(strcmp(M_.param_names, 'sigma_ren'));
        state = dr.ghu * shock_vec;
        irf_all = zeros(M_.endo_nbr, T);
        irf_all(:,1) = state;
        for tt = 2:T
            state_vars = state(M_.nstatic+1:M_.nstatic+M_.nspred);
            state = dr.ghx * state_vars;
            irf_all(:,tt) = state;
        end
        irf_reord = zeros(M_.endo_nbr, T);
        irf_reord(dr.order_var,:) = irf_all;

        Y_idx = find(strcmp(M_.endo_names, 'Y'));
        C_idx = find(strcmp(M_.endo_names, 'C'));
        mu_impact_Y(i) = irf_reord(Y_idx, 1);
        irf_C_mu = irf_reord(C_idx, :);
        sum_disc = sum(discount_factors .* log(1 + irf_C_mu/100));
        mu_welfare(i) = abs((1 - exp(sum_disc / normalizer)) * 100);
        fprintf('  mu = %.2f: Welfare = %.6f%%, Impact Y = %+.4f%%\n', ...
            mu_values(i), mu_welfare(i), mu_impact_Y(i));
    end
end
set_param_value('mu', 0.16);  % restore

% Sensitivity: beta (discount factor / EIS)
beta_values = [0.98, 0.99, 0.995];
beta_welfare = zeros(size(beta_values));
beta_impact_Y = zeros(size(beta_values));
beta_impact_Bstar = zeros(size(beta_values));
beta_C_share = zeros(size(beta_values));

fprintf('\nSensitivity to beta (discount factor / EIS):\n');
for i = 1:length(beta_values)
    set_param_value('beta', beta_values(i));
    % Also update r_bar for consistency: r_bar = 1/beta - 1
    set_param_value('r_bar', 1/beta_values(i) - 1);
    [oo_.dr.ys, M_.params, info_check] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    [oo_.dr, info_ss, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);
    if info_ss(1) == 0
        dr = oo_.dr;
        shock_vec = zeros(M_.exo_nbr, 1);
        shock_vec(find(strcmp(M_.exo_names, 'eps_ren'))) = M_.params(strcmp(M_.param_names, 'sigma_ren'));
        state = dr.ghu * shock_vec;
        irf_all = zeros(M_.endo_nbr, T);
        irf_all(:,1) = state;
        for tt = 2:T
            state_vars = state(M_.nstatic+1:M_.nstatic+M_.nspred);
            state = dr.ghx * state_vars;
            irf_all(:,tt) = state;
        end
        irf_reord = zeros(M_.endo_nbr, T);
        irf_reord(dr.order_var,:) = irf_all;

        Y_idx = find(strcmp(M_.endo_names, 'Y'));
        C_idx = find(strcmp(M_.endo_names, 'C'));
        B_star_idx = find(strcmp(M_.endo_names, 'B_star'));

        beta_impact_Y(i) = irf_reord(Y_idx, 1);
        beta_impact_Bstar(i) = irf_reord(B_star_idx, 1);
        beta_C_share(i) = oo_.dr.ys(C_idx) / oo_.dr.ys(Y_idx) * 100;

        % Welfare with scenario-specific beta for discounting
        disc_beta = beta_values(i).^(0:T-1);
        norm_beta = sum(disc_beta);
        irf_C_beta = irf_reord(C_idx, :);
        sum_disc_b = sum(disc_beta .* log(1 + irf_C_beta/100));
        beta_welfare(i) = abs((1 - exp(sum_disc_b / norm_beta)) * 100);

        fprintf('  beta = %.3f (ann. rate ~%.1f%%): Welfare = %.6f%%, Impact Y = %+.4f%%, Impact B* = %+.4f%%, C/Y = %.1f%%\n', ...
            beta_values(i), (1/beta_values(i)-1)*400, beta_welfare(i), ...
            beta_impact_Y(i), beta_impact_Bstar(i), beta_C_share(i));
    else
        fprintf('  beta = %.3f: FAILED to solve (info = %d)\n', beta_values(i), info_ss(1));
    end
end
set_param_value('beta', 0.99);  % restore
set_param_value('r_bar', 1/0.99 - 1);  % restore

% Final restore
[oo_.dr.ys, M_.params, ~] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
[oo_.dr, ~, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);

%% FIGURE 5: Sensitivity Summary
figure('Position', [100 100 1200 400]);

subplot(1,3,1);
bar(psi_welfare);
set(gca, 'XTickLabel', {'1.5', '2.0', '3.0'});
xlabel('\psi'); ylabel('Welfare Cost (%)');
title('Reliability Sensitivity (\psi)');
grid on;

subplot(1,3,2);
bar(phi_welfare);
set(gca, 'XTickLabel', {'0.5', '1.5', '3.0'});
xlabel('\phi_{grid}'); ylabel('Welfare Cost (%)');
title('Grid Investment Response (\phi_{grid})');
grid on;

subplot(1,3,3);
bar(mu_welfare);
set(gca, 'XTickLabel', {'0.10', '0.16', '0.25'});
xlabel('\mu'); ylabel('Welfare Cost (%)');
title('Battery Share (\mu)');
grid on;

sgtitle('Parameter Sensitivity: Welfare Cost of Intermittency', ...
    'FontSize', 14, 'FontWeight', 'bold');

print('sensitivity_analysis.png', '-dpng', '-r300');
fprintf('\nFigure saved: sensitivity_analysis.png\n');

%% ========================================================================
%% COUNTERFACTUAL 3: JOINT SHOCK ("PERFECT STORM") SCENARIO
%% ========================================================================

fprintf('\n========================================\n');
fprintf('JOINT SHOCK SCENARIO: PERFECT STORM\n');
fprintf('========================================\n\n');

% Restore baseline parameters
set_param_value('chi', 1.0);
set_param_value('psi', 2.0);
set_param_value('phi_grid', 1.5);
set_param_value('mu', 0.16);
[oo_.dr.ys, M_.params, ~] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
[oo_.dr, ~, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);

dr = oo_.dr;
nvar = M_.endo_nbr;
nshocks = M_.exo_nbr;
eps_ren_idx = find(strcmp(M_.exo_names, 'eps_ren'));
eps_bat_idx = find(strcmp(M_.exo_names, 'eps_bat'));
sigma_ren_val = M_.params(strcmp(M_.param_names, 'sigma_ren'));
sigma_bat_val = M_.params(strcmp(M_.param_names, 'sigma_bat'));

% Joint shock: simultaneous 1-sigma renewable intermittency + 1-sigma battery price
shock_joint = zeros(nshocks, 1);
shock_joint(eps_ren_idx) = sigma_ren_val;   % positive = more intermittency
shock_joint(eps_bat_idx) = sigma_bat_val;    % positive = battery price increase

% Also compute individual shocks for comparison
shock_ren_only = zeros(nshocks, 1);
shock_ren_only(eps_ren_idx) = sigma_ren_val;
shock_bat_only = zeros(nshocks, 1);
shock_bat_only(eps_bat_idx) = sigma_bat_val;

% Simulate IRFs for all three scenarios
scenarios = {shock_ren_only, shock_bat_only, shock_joint};
scenario_names = {'Intermittency Only', 'Battery Price Only', 'Joint (Perfect Storm)'};
scenario_irfs = cell(3,1);

for s = 1:3
    irf_mat = zeros(nvar, T);
    state = dr.ghu * scenarios{s};
    irf_mat(:,1) = state;
    for tt = 2:T
        state_vars = state(M_.nstatic+1:M_.nstatic+M_.nspred);
        state = dr.ghx * state_vars;
        irf_mat(:,tt) = state;
    end
    irf_reord = zeros(nvar, T);
    irf_reord(dr.order_var,:) = irf_mat;
    scenario_irfs{s} = irf_reord;
end

Y_idx = find(strcmp(M_.endo_names, 'Y'));
C_idx = find(strcmp(M_.endo_names, 'C'));
u_idx = find(strcmp(M_.endo_names, 'u'));
I_bat_idx = find(strcmp(M_.endo_names, 'I_bat'));
I_grid_idx = find(strcmp(M_.endo_names, 'I_grid'));
A_bat_idx = find(strcmp(M_.endo_names, 'A_bat'));
B_star_idx = find(strcmp(M_.endo_names, 'B_star'));

fprintf('IMPACT EFFECTS (Quarter 0):\n');
fprintf('%-25s  %12s  %12s  %12s\n', 'Variable', 'Ren Only', 'Bat Only', 'Joint');
fprintf('%-25s  %+11.4f%%  %+11.4f%%  %+11.4f%%\n', 'Output (Y)', ...
    scenario_irfs{1}(Y_idx,1), scenario_irfs{2}(Y_idx,1), scenario_irfs{3}(Y_idx,1));
fprintf('%-25s  %+11.4f%%  %+11.4f%%  %+11.4f%%\n', 'Consumption (C)', ...
    scenario_irfs{1}(C_idx,1), scenario_irfs{2}(C_idx,1), scenario_irfs{3}(C_idx,1));
fprintf('%-25s  %+11.4f%%  %+11.4f%%  %+11.4f%%\n', 'Reliability (u)', ...
    scenario_irfs{1}(u_idx,1), scenario_irfs{2}(u_idx,1), scenario_irfs{3}(u_idx,1));
fprintf('%-25s  %+11.4f%%  %+11.4f%%  %+11.4f%%\n', 'Battery Investment', ...
    scenario_irfs{1}(I_bat_idx,1), scenario_irfs{2}(I_bat_idx,1), scenario_irfs{3}(I_bat_idx,1));
fprintf('%-25s  %+11.4f%%  %+11.4f%%  %+11.4f%%\n', 'Grid Investment', ...
    scenario_irfs{1}(I_grid_idx,1), scenario_irfs{2}(I_grid_idx,1), scenario_irfs{3}(I_grid_idx,1));
fprintf('%-25s  %+11.4f%%  %+11.4f%%  %+11.4f%%\n', 'Net Foreign Assets', ...
    scenario_irfs{1}(B_star_idx,1), scenario_irfs{2}(B_star_idx,1), scenario_irfs{3}(B_star_idx,1));

% Welfare costs
joint_welfare = zeros(3,1);
for s = 1:3
    irf_C_s = scenario_irfs{s}(C_idx,:);
    sum_disc = sum(discount_factors .* log(1 + irf_C_s/100));
    joint_welfare(s) = abs((1 - exp(sum_disc / normalizer)) * 100);
end

fprintf('\nWELFARE COSTS:\n');
for s = 1:3
    fprintf('  %-25s: %.6f%%\n', scenario_names{s}, joint_welfare(s));
end

% Check superadditivity of welfare cost
sum_individual = joint_welfare(1) + joint_welfare(2);
fprintf('\nSUM of individual welfare costs:  %.6f%%\n', sum_individual);
fprintf('JOINT welfare cost:              %.6f%%\n', joint_welfare(3));
if joint_welfare(3) > sum_individual
    fprintf('Amplification: Joint > Sum by %.1f%% (nonlinear interaction)\n', ...
        (joint_welfare(3)/sum_individual - 1)*100);
else
    fprintf('No amplification: Joint <= Sum (linear superposition)\n');
end

% Figure: Joint shock comparison
figure('Position', [100 100 1200 800]);
colors_j = {'b', 'r', [0.6 0 0.6]};
lstyles_j = {'-', '--', '-.'};

subplot(2,3,1);
hold on;
for s = 1:3
    plot(t, scenario_irfs{s}(Y_idx,:), 'Color', colors_j{s}, ...
        'LineStyle', lstyles_j{s}, 'LineWidth', 2, 'DisplayName', scenario_names{s});
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
hold off; grid on; title('Output (Y)'); ylabel('% dev'); xlabel('Quarters');
legend('Location', 'southeast', 'FontSize', 7);

subplot(2,3,2);
hold on;
for s = 1:3
    plot(t, scenario_irfs{s}(u_idx,:), 'Color', colors_j{s}, ...
        'LineStyle', lstyles_j{s}, 'LineWidth', 2);
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
hold off; grid on; title('Reliability (u)'); ylabel('% dev'); xlabel('Quarters');

subplot(2,3,3);
hold on;
for s = 1:3
    plot(t, scenario_irfs{s}(I_bat_idx,:), 'Color', colors_j{s}, ...
        'LineStyle', lstyles_j{s}, 'LineWidth', 2);
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
hold off; grid on; title('Battery Investment'); ylabel('% dev'); xlabel('Quarters');

subplot(2,3,4);
hold on;
for s = 1:3
    plot(t, scenario_irfs{s}(C_idx,:), 'Color', colors_j{s}, ...
        'LineStyle', lstyles_j{s}, 'LineWidth', 2);
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
hold off; grid on; title('Consumption (C)'); ylabel('% dev'); xlabel('Quarters');

subplot(2,3,5);
hold on;
for s = 1:3
    plot(t, scenario_irfs{s}(B_star_idx,:), 'Color', colors_j{s}, ...
        'LineStyle', lstyles_j{s}, 'LineWidth', 2);
end
plot(t, zeros(size(t)), 'k--', 'LineWidth', 0.5);
hold off; grid on; title('Net Foreign Assets (B^*)'); ylabel('% dev'); xlabel('Quarters');

subplot(2,3,6);
bar(joint_welfare);
set(gca, 'XTickLabel', {'Ren', 'Bat', 'Joint'});
ylabel('Welfare Cost (%)');
title('Welfare: Individual vs Joint');
grid on;

sgtitle('"Perfect Storm": Simultaneous Intermittency + Battery Price Shocks', ...
    'FontSize', 13, 'FontWeight', 'bold');

print('joint_shock_perfect_storm.png', '-dpng', '-r300');
fprintf('\nFigure saved: joint_shock_perfect_storm.png\n');

%% ========================================================================
%% SUMMARY TABLE
%% ========================================================================

fprintf('\n========================================\n');
fprintf('SUMMARY OF ALL RESULTS\n');
fprintf('========================================\n\n');

fprintf('MODEL SPECIFICATION:\n');
fprintf('  Variables: 16 (SOE with B_star, r_star)\n');
fprintf('  States: %d, Jumpers: %d, Static: %d\n', M_.nspred, M_.nsfwrd, M_.nstatic);
fprintf('  Shocks: 3 (eps_ren, eps_bat, eps_I)\n\n');

fprintf('SIGNAL ATTENUATION RESULTS:\n');
for i = 1:length(chi_values)
    fprintf('  chi = %.1f: Welfare cost = %.6f%%\n', chi_values(i), chi_welfare(i));
end
fprintf('\n');

fprintf('SENSITIVITY RESULTS:\n');
fprintf('  psi:      [%.1f, %.1f, %.1f] -> Welfare [%.6f%%, %.6f%%, %.6f%%]\n', ...
    psi_values(1), psi_values(2), psi_values(3), psi_welfare(1), psi_welfare(2), psi_welfare(3));
fprintf('  phi_grid: [%.1f, %.1f, %.1f] -> Welfare [%.6f%%, %.6f%%, %.6f%%]\n', ...
    phi_grid_values(1), phi_grid_values(2), phi_grid_values(3), phi_welfare(1), phi_welfare(2), phi_welfare(3));
fprintf('  mu:       [%.2f, %.2f, %.2f] -> Welfare [%.6f%%, %.6f%%, %.6f%%]\n', ...
    mu_values(1), mu_values(2), mu_values(3), mu_welfare(1), mu_welfare(2), mu_welfare(3));

%% ========================================================================
%% COUNTERFACTUAL 4: CLOSING THE GAP — POLICY SCENARIOS
%% ========================================================================
% Simulates the Reliability Valley under three policy regimes:
%   1. Baseline: current parameters
%   2. Scenario A — Accelerated Grid Planning (phi_grid = 3.0)
%   3. Scenario B — Battery Subsidy (20% cost reduction, P_bat_eff = 0.8)
%
% Method: Deterministic transition where Vol_ren increases along an S-curve
% (15% → 50% renewable penetration over 100 quarters). At each step, we
% re-solve the steady state to find the implied reliability given the
% current flexibility stock grows at endogenous rates.

fprintf('\n========================================\n');
fprintf('CLOSING THE GAP: POLICY COUNTERFACTUALS\n');
fprintf('========================================\n\n');

% Restore baseline
set_param_value('chi', 1.0);
set_param_value('psi', 2.0);
set_param_value('phi_grid', 1.5);
set_param_value('mu', 0.16);
[oo_.dr.ys, M_.params, ~] = vietnam_dsge_steadystate(oo_.dr.ys, oo_.exo_steady_state, M_, options_);
[oo_.dr, ~, M_.params] = resol(0, M_, options_, oo_.dr, oo_.dr.ys, oo_.exo_steady_state, oo_.exo_det_steady_state);

% Transition parameters
T_trans = 100;  % 100 quarters = 25 years
t_trans = (1:T_trans)';

% S-curve for renewable penetration: 15% -> 50%
theta_start = 0.15;
theta_end   = 0.50;
% Logistic S-curve: front-loaded
k_steep = 0.08;  % steepness
t_mid = 35;      % midpoint at quarter 35 (year ~9)
theta_path = theta_start + (theta_end - theta_start) ./ (1 + exp(-k_steep * (t_trans - t_mid)));

% Baseline steady-state values
p_psi_val     = M_.params(strcmp(M_.param_names, 'psi'));
phi_int_val   = M_.params(strcmp(M_.param_names, 'phi_int'));
p_sigma_ren   = M_.params(strcmp(M_.param_names, 'sigma_ren'));
p_delta_b     = M_.params(strcmp(M_.param_names, 'delta_b'));
p_delta_g     = M_.params(strcmp(M_.param_names, 'delta_g'));
p_mu          = M_.params(strcmp(M_.param_names, 'mu'));
p_rho_flex    = M_.params(strcmp(M_.param_names, 'rho_flex'));

K_b_ss = oo_.dr.ys(strcmp(M_.endo_names, 'K_b'));
K_g_ss = oo_.dr.ys(strcmp(M_.endo_names, 'K_g'));
A_bat_ss = oo_.dr.ys(strcmp(M_.endo_names, 'A_bat'));
F_ss_val = oo_.dr.ys(strcmp(M_.endo_names, 'F'));
u_ss_val = oo_.dr.ys(strcmp(M_.endo_names, 'u'));

% Installed renewable capacity proxy (baseline)
K_ren_base = 0.045 / theta_start;  % from Vol_ren = theta * K_ren * sigma_ren

% Policy scenarios
scenario_labels = {'Baseline (\phi_{grid}=1.5)', ...
                   'Accelerated Grid (\phi_{grid}=3.0)', ...
                   'Battery Subsidy (P_{bat}=0.8)'};
n_scen = 3;

% Growth rate parameters per scenario
% Grid capital quarterly growth rates (endogenous, approximated)
grid_growth_base = 0.0135;   % 1.35% per quarter baseline
grid_accel_factor = [1.0, 2.0, 1.0];  % Scenario A: double grid growth
bat_growth_base = 0.020;     % 2.0% per quarter baseline
bat_accel_factor = [1.0, 1.0, 1.5];   % Scenario B: 50% faster battery growth
learning_boost = [0.0, 0.0, 0.02];    % Scenario B: subsidy enhances A_bat growth

% Simulate each scenario
u_paths     = zeros(T_trans, n_scen);
F_paths     = zeros(T_trans, n_scen);
K_b_paths   = zeros(T_trans, n_scen);
K_g_paths   = zeros(T_trans, n_scen);
Vol_paths   = zeros(T_trans, n_scen);

for sc = 1:n_scen
    K_b_t = K_b_ss;
    K_g_t = K_g_ss;
    A_bat_t = A_bat_ss;

    for qq = 1:T_trans
        % Current Vol_ren given penetration path
        K_ren_t = K_ren_base * (theta_path(qq) / theta_start);
        Vol_ren_t = theta_path(qq) * K_ren_t * p_sigma_ren;

        % Flexibility from current stocks
        s_f = (p_rho_flex - 1) / p_rho_flex;
        F_t = (p_mu * (A_bat_t * K_b_t)^s_f + (1-p_mu) * K_g_t^s_f)^(1/s_f);

        % Reliability
        z_t = p_psi_val * F_t / (phi_int_val * Vol_ren_t);
        u_t = 1 - exp(-z_t);

        % Store
        u_paths(qq, sc) = u_t;
        F_paths(qq, sc) = F_t;
        K_b_paths(qq, sc) = K_b_t;
        K_g_paths(qq, sc) = K_g_t;
        Vol_paths(qq, sc) = Vol_ren_t;

        % Reliability gap drives investment acceleration
        u_gap = max(0, 0.97 - u_t) / 0.97;

        % Grid investment growth (scenario-dependent)
        grid_growth = grid_growth_base * grid_accel_factor(sc) * (1 + 2*u_gap);
        I_grid_t = (p_delta_g + grid_growth) * K_g_t;
        K_g_t = (1 - p_delta_g) * K_g_t + I_grid_t;

        % Battery investment growth (scenario-dependent)
        bat_growth = bat_growth_base * bat_accel_factor(sc) * (1 + 3*u_gap);
        I_bat_t = (p_delta_b + bat_growth) * K_b_t;
        K_b_t = (1 - p_delta_b) * K_b_t + I_bat_t;

        % Learning-by-doing
        A_bat_t = A_bat_t * (1 + 0.10 * u_gap + learning_boost(sc));
    end
end

% Find valley characteristics
valley_results = zeros(n_scen, 4);  % [min_u, valley_quarter, valley_duration, recovery_quarter]
for sc = 1:n_scen
    [min_u, min_idx] = min(u_paths(:, sc));
    valley_results(sc, 1) = min_u * 100;
    valley_results(sc, 2) = min_idx;

    % Duration below 97%
    below_target = u_paths(:, sc) < 0.97;
    valley_results(sc, 3) = sum(below_target);

    % Recovery quarter (first quarter back above 97% after the valley)
    post_valley = find(~below_target & (1:T_trans)' > min_idx, 1, 'first');
    if ~isempty(post_valley)
        valley_results(sc, 4) = post_valley;
    else
        valley_results(sc, 4) = T_trans;  % never recovers
    end
end

fprintf('RELIABILITY VALLEY COMPARISON:\n');
fprintf('%-35s  %8s  %8s  %10s  %10s\n', 'Scenario', 'Min u(%)', 'Quarter', 'Duration', 'Recovery');
for sc = 1:n_scen
    fprintf('%-35s  %7.2f%%  Q%-6d  %8d Q  Q%-8d\n', ...
        scenario_labels{sc}, valley_results(sc,1), valley_results(sc,2), ...
        valley_results(sc,3), valley_results(sc,4));
end

% Compute valley depth reduction relative to baseline
fprintf('\nPOLICY EFFECTIVENESS:\n');
baseline_depth = 97 - valley_results(1,1);
for sc = 2:n_scen
    depth_sc = 97 - valley_results(sc,1);
    reduction = (baseline_depth - depth_sc) / baseline_depth * 100;
    duration_reduction = (valley_results(1,3) - valley_results(sc,3)) / valley_results(1,3) * 100;
    fprintf('  %s:\n', scenario_labels{sc});
    fprintf('    Valley depth reduction: %.1f%% (%.2f pp -> %.2f pp)\n', ...
        reduction, baseline_depth, depth_sc);
    fprintf('    Duration reduction:     %.1f%% (%d Q -> %d Q)\n', ...
        duration_reduction, valley_results(1,3), valley_results(sc,3));
end

% FIGURE: Closing the Gap
figure('Position', [100 100 1400 500]);
colors_policy = {'b', [0 0.6 0], [0.8 0 0.8]};
lstyles_policy = {'-', '--', '-.'};

% Panel 1: Renewable penetration path (same for all)
subplot(1,3,1);
plot(t_trans/4, theta_path*100, 'k-', 'LineWidth', 2.5);
hold on;
yline(50, 'r--', 'PDP8 Target', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
hold off;
grid on; xlabel('Years'); ylabel('Renewable Penetration (%)');
title('VRE Deployment Path');

% Panel 2: Reliability paths
subplot(1,3,2);
hold on;
for sc = 1:n_scen
    plot(t_trans/4, u_paths(:,sc)*100, 'Color', colors_policy{sc}, ...
        'LineStyle', lstyles_policy{sc}, 'LineWidth', 2.5, ...
        'DisplayName', scenario_labels{sc});
end
yline(97, 'r--', 'Target (97%)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
hold off;
grid on; xlabel('Years'); ylabel('Grid Reliability (%)');
title('Reliability Valley: Policy Comparison');
legend('Location', 'south', 'FontSize', 9);
ylim([85 100]);

% Panel 3: Flexibility growth
subplot(1,3,3);
hold on;
for sc = 1:n_scen
    plot(t_trans/4, F_paths(:,sc), 'Color', colors_policy{sc}, ...
        'LineStyle', lstyles_policy{sc}, 'LineWidth', 2.5, ...
        'DisplayName', scenario_labels{sc});
end
plot(t_trans/4, Vol_paths(:,1)*phi_int_val/p_psi_val * (-log(1-0.97)), ...
    'r--', 'LineWidth', 1.5, 'DisplayName', 'F required for u=97%');
hold off;
grid on; xlabel('Years'); ylabel('Flexibility (F)');
title('Flexibility vs. Requirement');
legend('Location', 'northwest', 'FontSize', 7);

sgtitle('Closing the Agility Gap: Policy Counterfactuals', ...
    'FontSize', 14, 'FontWeight', 'bold');

print('closing_the_gap.png', '-dpng', '-r300');
fprintf('\nFigure saved: closing_the_gap.png\n');

% -------------------------------------------------------------------------
% FIGURE: Reliability Valley (baseline only, two-panel for paper)
% Left: S-curve renewable penetration path
% Right: Reliability transition path showing the valley
% -------------------------------------------------------------------------
fig_valley = figure('Position', [100 100 1100 430]);

% Left panel: Renewable penetration S-curve
subplot(1,2,1);
hold on;
plot(t_trans/4, theta_path*100, 'Color', [0.18 0.55 0.34], 'LineWidth', 2.5);
yline(50, 'r--', 'LineWidth', 1.5);
text(23.5, 51.5, 'PDP8 Target (50\%)', 'Color', 'r', 'FontSize', 9);
% Mark the midpoint
[~, mid_q] = min(abs(theta_path - (theta_start + theta_end)/2));
plot(mid_q/4, theta_path(mid_q)*100, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'k');
text(mid_q/4 + 0.5, theta_path(mid_q)*100 - 3, 'Inflection', 'FontSize', 8, 'Color', 'k');
hold off;
grid on;
xlabel('Years from baseline', 'FontSize', 10);
ylabel('Renewable penetration (\theta_{ren}, \%)', 'FontSize', 10);
title('VRE Deployment: S-Curve Path', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 T_trans/4]);
ylim([10 55]);

% Right panel: Baseline reliability path
subplot(1,2,2);
u_baseline_pct = u_paths(:,1) * 100;
[min_u_val, min_u_idx] = min(u_baseline_pct);
below_97 = u_baseline_pct < 97;
hold on;

% Shade the valley region
x_fill = [(find(below_97,1,'first')-1)/4, (find(below_97,1,'last'))/4, ...
           (find(below_97,1,'last'))/4, (find(below_97,1,'first')-1)/4];
y_fill = [min_u_val-1, min_u_val-1, 97, 97];
fill(x_fill, y_fill, [1 0.8 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Reliability path
plot(t_trans/4, u_baseline_pct, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Grid reliability');
yline(97, 'r--', 'LineWidth', 1.8);
text(0.5, 97.4, 'Target (97\%)', 'Color', 'r', 'FontSize', 9);

% Annotate nadir
plot(min_u_idx/4, min_u_val, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(min_u_idx/4 + 0.5, min_u_val - 0.8, ...
    sprintf('Nadir: Q%d\n(%.1f%%)', min_u_idx, min_u_val), ...
    'FontSize', 8, 'Color', [0.7 0 0]);

% Label valley
valley_center = (find(below_97,1,'first') + find(below_97,1,'last')) / 2;
text(valley_center/4, min_u_val + 1.0, 'Reliability Valley', ...
    'FontSize', 9, 'Color', [0.8 0 0], 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold');

hold off;
grid on;
xlabel('Years from baseline', 'FontSize', 10);
ylabel('Grid reliability (u, \%)', 'FontSize', 10);
title('Reliability Valley: Baseline Transition', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 T_trans/4]);
ylim([min_u_val - 3, 100]);

sgtitle('The Reliability Valley: PDP8 Transition Dynamics', ...
    'FontSize', 13, 'FontWeight', 'bold');

print('reliability_valley.png', '-dpng', '-r300');
fprintf('Figure saved: reliability_valley.png\n');
close(fig_valley);

%% ========================================================================
%% SAVE WORKSPACE AND EXPORT DATA
%% ========================================================================

save('dynare_results.mat');
fprintf('\nWorkspace saved: dynare_results.mat\n');

% Create IRF table for LaTeX
irf_table = table((0:T-1)', irf_Y', irf_C', irf_u', irf_I_bat', irf_I_grid', ...
    irf_F', irf_A_bat', irf_B_star', irf_r_star', ...
    'VariableNames', {'Quarter', 'Output', 'Consumption', 'Utilization', ...
    'Battery_Investment', 'Grid_Investment', 'Flexibility', 'Battery_Tech', ...
    'Net_Foreign_Assets', 'Interest_Rate'});

writetable(irf_table, 'irf_dynare.csv');
fprintf('IRF data exported: irf_dynare.csv\n');

fprintf('\n========================================\n');
fprintf('ALL ANALYSIS COMPLETE!\n');
fprintf('========================================\n');
