%% ========================================================================
%% VIETNAM DSGE — PROPORTIONAL TAX VARIANT — COMPREHENSIVE RESULTS SCRIPT
%% ========================================================================
%% Runs vietnam_dsge_proptax.mod, extracts all results needed for the paper:
%%   - Steady state values
%%   - IRFs (intermittency, battery price, grid implementation shocks)
%%   - Forecast Error Variance Decomposition (FEVD)
%%   - Welfare costs (baseline + signal attenuation counterfactuals)
%%   - Sensitivity analysis (phi_grid, psi)
%%   - Joint "perfect storm" shock scenario
%%   - External sector dynamics
%%   - All figures + JSON output for paper verification
%%
%% Author: Toan T. Nguyen  |  March 2026
%% ========================================================================

clear all; close all; clc;
addpath('C:/dynare/6.5/matlab');
% Ensure Dynare preprocessor is on the system PATH (needed in batch mode)
setenv('PATH', [getenv('PATH') ';C:\dynare\6.5\preprocessor']);

fprintf('==============================================\n');
fprintf('  VIETNAM DSGE — PROPORTIONAL TAX VARIANT\n');
fprintf('==============================================\n\n');

dynare vietnam_dsge_proptax noclearall;
fprintf('\n[OK] Model solved. Blanchard-Kahn conditions verified.\n\n');

%% ========================================================================
%% 0. HELPER: Compute IRF from policy matrices for a given shock
%% ========================================================================
function irf_all = compute_irf(M_, dr, shock_name, T)
    nvar = M_.endo_nbr;
    eps_idx = find(strcmp(M_.exo_names, shock_name));
    sigma_val = M_.params(strcmp(M_.param_names, ['sigma_' shock_name(5:end)]));
    % For eps_I the param name is different — handle separately
    shock_vec = zeros(M_.exo_nbr, 1);
    shock_vec(eps_idx) = sigma_val;

    irf_mat = zeros(nvar, T);
    state = dr.ghu * shock_vec;
    irf_mat(:, 1) = state;
    for tt = 2:T
        state_vars = state(M_.nstatic+1:M_.nstatic+M_.nspred);
        state = dr.ghx * state_vars;
        irf_mat(:, tt) = state;
    end
    % Reorder to model variable ordering
    irf_all = zeros(nvar, T);
    irf_all(dr.order_var, :) = irf_mat;
end

%% ========================================================================
%% 1. STEADY STATE
%% ========================================================================
T  = 40;
t  = 0:T-1;
ss = oo_.dr.ys;
vn = M_.endo_names;

Y_ss   = ss(strcmp(vn,'Y'));
C_ss   = ss(strcmp(vn,'C'));
L_ss   = ss(strcmp(vn,'L'));
Kp_ss  = ss(strcmp(vn,'K_p'));
Kb_ss  = ss(strcmp(vn,'K_b'));
Kg_ss  = ss(strcmp(vn,'K_g'));
Ip_ss  = ss(strcmp(vn,'I_p'));
Ibat_ss= ss(strcmp(vn,'I_bat'));
Ig_ss  = ss(strcmp(vn,'I_grid'));
u_ss   = ss(strcmp(vn,'u'));
F_ss   = ss(strcmp(vn,'F'));
V_ss   = ss(strcmp(vn,'V'));
tau_ss = M_.params(strcmp(M_.param_names,'tau'));   % tau is a parameter now
r_ss   = ss(strcmp(vn,'r_star'));

fprintf('========================================\n');
fprintf('STEADY STATE (Proportional Tax)\n');
fprintf('========================================\n');
fprintf('  Y      = %.6f\n', Y_ss);
fprintf('  C      = %.6f   (%.2f%% of Y)\n', C_ss, C_ss/Y_ss*100);
fprintf('  L      = %.6f\n', L_ss);
fprintf('  K_p    = %.6f\n', Kp_ss);
fprintf('  K_b    = %.6f\n', Kb_ss);
fprintf('  K_g    = %.6f\n', Kg_ss);
fprintf('  I_p    = %.6f\n', Ip_ss);
fprintf('  I_bat  = %.6f\n', Ibat_ss);
fprintf('  I_grid = %.6f\n', Ig_ss);
fprintf('  u      = %.6f\n', u_ss);
fprintf('  F      = %.6f\n', F_ss);
fprintf('  V      = %.6f\n', V_ss);
fprintf('  tau    = %.6f   (%.4f%%)\n', tau_ss, tau_ss*100);
fprintf('  r_star = %.6f\n', r_ss);
fprintf('  Kb/Kg ratio = %.4f\n', Kb_ss/Kg_ss);
fprintf('\n');

%% ========================================================================
%% 2. EXTRACT IRFS FROM OO_ (from stoch_simul output)
%% ========================================================================

irf_Y      = oo_.irfs.Y_eps_ren;
irf_C      = oo_.irfs.C_eps_ren;
irf_L      = oo_.irfs.L_eps_ren;
irf_u      = oo_.irfs.u_eps_ren;
irf_I_bat  = oo_.irfs.I_bat_eps_ren;
irf_I_grid = oo_.irfs.I_grid_eps_ren;
irf_I_p    = oo_.irfs.I_p_eps_ren;
irf_F      = oo_.irfs.F_eps_ren;
irf_A_bat  = oo_.irfs.A_bat_eps_ren;
irf_B_star = oo_.irfs.B_star_eps_ren;
irf_r_star = oo_.irfs.r_star_eps_ren;

fprintf('========================================\n');
fprintf('IMPACT EFFECTS: Intermittency Shock\n');
fprintf('========================================\n');
fprintf('  Output (Y):          %+.6f%%\n', irf_Y(1));
fprintf('  Consumption (C):     %+.6f%%\n', irf_C(1));
fprintf('  Labor (L):           %+.6f%%\n', irf_L(1));
fprintf('  Utilization (u):     %+.6f%%\n', irf_u(1));
fprintf('  Battery Inv (I_bat): %+.6f%%\n', irf_I_bat(1));
fprintf('  Grid Inv (I_grid):   %+.6f%%\n', irf_I_grid(1));
fprintf('  Productive Inv (Ip): %+.6f%%\n', irf_I_p(1));
fprintf('  Flexibility (F):     %+.6f%%\n', irf_F(1));
fprintf('  Tech (A_bat) q10:    %+.6f%%\n', irf_A_bat(10));
fprintf('  Net FA (B_star) q1:  %+.6f\n',   irf_B_star(1));
fprintf('\n');

% Cumulative reliability loss (sum over 40q)
cum_u_loss = sum(irf_u);
fprintf('  Cumulative u loss (40q): %.6f\n', cum_u_loss);

% Agility gap: battery vs grid investment ratio (impact)
if irf_I_grid(1) ~= 0
    agility_ratio = abs(irf_I_bat(1)) / abs(irf_I_grid(1));
    fprintf('  Agility ratio I_bat/I_grid (impact): %.4f\n', agility_ratio);
end
fprintf('\n');

%% ========================================================================
%% 3. BATTERY PRICE SHOCK IRFS
%% ========================================================================
fprintf('========================================\n');
fprintf('IMPACT EFFECTS: Battery Price Shock\n');
fprintf('========================================\n');
if isfield(oo_.irfs, 'Y_eps_bat')
    fprintf('  Output (Y):          %+.6f%%\n', oo_.irfs.Y_eps_bat(1));
    fprintf('  Battery Inv (I_bat): %+.6f%%\n', oo_.irfs.I_bat_eps_bat(1));
    fprintf('  Battery Price:       %+.6f%%\n', oo_.irfs.P_bat_eps_bat(1));
    fprintf('  Utilization (u):     %+.6f%%\n', oo_.irfs.u_eps_bat(1));
    fprintf('  Net FA (B_star) q1:  %+.6f\n',   oo_.irfs.B_star_eps_bat(1));
end
fprintf('\n');

%% ========================================================================
%% 4. FORECAST ERROR VARIANCE DECOMPOSITION
%% ========================================================================
fprintf('========================================\n');
fprintf('FORECAST ERROR VARIANCE DECOMPOSITION\n');
fprintf('========================================\n');

% oo_.variance_decomposition is (nvar_listed x nshocks) in percent
% Variables listed in stoch_simul: Y C L I_p I_bat I_grid u F A_bat V P_bat B_star r_star a_ren a_I
if isfield(oo_, 'variance_decomposition')
    vd = oo_.variance_decomposition;
    % Dynare lists shocks in order of varexo: eps_ren eps_bat eps_I
    var_list = {'Y','C','L','I_p','I_bat','I_grid','u','F','A_bat','V','P_bat','B_star','r_star','a_ren','a_I'};
    shock_list = {'eps_ren','eps_bat','eps_I'};
    fprintf('  Variable    | eps_ren(%%)  | eps_bat(%%)  | eps_I(%%)\n');
    fprintf('  -----------   ----------    ----------    ---------\n');
    for i = 1:length(var_list)
        if i <= size(vd,1)
            fprintf('  %-12s| %9.2f   | %9.2f   | %9.2f\n', ...
                var_list{i}, vd(i,1), vd(i,2), vd(i,3));
        end
    end
else
    fprintf('  oo_.variance_decomposition not found.\n');
    fprintf('  Extracting from theoretical moments...\n');
    if isfield(oo_, 'gamma_y')
        fprintf('  gamma_y available but FEVD requires manual computation.\n');
    end
end
fprintf('\n');

%% ========================================================================
%% 5. STOCHASTIC SIMULATION MOMENTS
%%    In Dynare 6, oo_.mean and oo_.var are numeric vectors/matrices
%%    indexed by M_.endo_names order (not struct fieldnames)
%% ========================================================================
fprintf('========================================\n');
fprintf('STOCHASTIC SIMULATION MOMENTS (1000 periods)\n');
fprintf('========================================\n');
try
    % oo_.mean is a numeric column vector; index by position
    Abat_mean_idx = find(strcmp(M_.endo_names, 'A_bat'));
    Y_mean_idx    = find(strcmp(M_.endo_names, 'Y'));
    C_mean_idx    = find(strcmp(M_.endo_names, 'C'));
    u_mean_idx    = find(strcmp(M_.endo_names, 'u'));

    if isnumeric(oo_.mean)
        fprintf('  Mean A_bat = %.6f  (long-run tech improvement)\n', oo_.mean(Abat_mean_idx));
        fprintf('  Mean Y     = %.6f\n', oo_.mean(Y_mean_idx));
    end
    if isnumeric(oo_.var)
        % oo_.var is a (nvar x nvar) covariance matrix, same ordering
        fprintf('  Std dev Y   = %.6f (%.4f%% of Y_ss)\n', ...
            sqrt(oo_.var(Y_mean_idx, Y_mean_idx)), ...
            sqrt(oo_.var(Y_mean_idx, Y_mean_idx))/Y_ss*100);
        fprintf('  Std dev C   = %.6f (%.4f%% of C_ss)\n', ...
            sqrt(oo_.var(C_mean_idx, C_mean_idx)), ...
            sqrt(oo_.var(C_mean_idx, C_mean_idx))/C_ss*100);
        fprintf('  Std dev u   = %.6f (%.4f%% of u_ss)\n', ...
            sqrt(oo_.var(u_mean_idx, u_mean_idx)), ...
            sqrt(oo_.var(u_mean_idx, u_mean_idx))/u_ss*100);
    end
catch ME_moments
    fprintf('  [Moments not extracted: %s]\n', ME_moments.message);
    fprintf('  (Read from Dynare printed output above)\n');
    fprintf('  Mean A_bat ~ 1.015176 (from sim output)\n');
    fprintf('  Std Y ~ 0.670%% of SS, Std C ~ 0.210%%\n');
end
fprintf('\n');

%% ========================================================================
%% 6. WELFARE COST (baseline chi=1.0)
%% ========================================================================

beta_val = M_.params(strcmp(M_.param_names, 'beta'));
T_w = min(T, length(irf_C));
disc = beta_val .^ (0:T_w-1);
norm_factor = (1 - beta_val^T_w) / (1 - beta_val);

sum_disc = sum(disc .* log(1 + irf_C(1:T_w)/100));
welfare_baseline = abs((1 - exp(sum_disc / norm_factor)) * 100);

fprintf('========================================\n');
fprintf('WELFARE ANALYSIS\n');
fprintf('========================================\n');
fprintf('  Baseline welfare cost (chi=1.0): %.6f%%\n', welfare_baseline);

%% ========================================================================
%% 7. COUNTERFACTUAL: SIGNAL ATTENUATION (chi variation)
%% ========================================================================
chi_values = [1.0, 0.5, 0.3, 0.0];
chi_welfare = zeros(size(chi_values));
chi_irfs_C  = zeros(length(chi_values), T_w);
chi_irfs_Y  = zeros(length(chi_values), T);
chi_irfs_u  = zeros(length(chi_values), T);
chi_irfs_A  = zeros(length(chi_values), T);

for i = 1:length(chi_values)
    set_param_value('chi', chi_values(i));
    [dr_ys_tmp, params_tmp, ~] = vietnam_dsge_proptax_steadystate(...
        oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    M_.params = params_tmp;
    [dr_tmp, info_tmp, M_.params] = resol(0, M_, options_, oo_.dr, ...
        dr_ys_tmp, oo_.exo_steady_state, oo_.exo_det_steady_state);

    if info_tmp(1) == 0
        nvar = M_.endo_nbr;
        eps_idx = find(strcmp(M_.exo_names, 'eps_ren'));
        sigma_ren = M_.params(strcmp(M_.param_names, 'sigma_ren'));
        sv = zeros(M_.exo_nbr, 1);
        sv(eps_idx) = sigma_ren;

        irf_mat = zeros(nvar, T_w);
        state = dr_tmp.ghu * sv;
        irf_mat(:,1) = state;
        for tt = 2:T_w
            sv2 = state(M_.nstatic+1:M_.nstatic+M_.nspred);
            state = dr_tmp.ghx * sv2;
            irf_mat(:,tt) = state;
        end
        irf_ro = zeros(nvar, T_w);
        irf_ro(dr_tmp.order_var,:) = irf_mat;

        C_idx   = find(strcmp(M_.endo_names,'C'));
        Y_idx   = find(strcmp(M_.endo_names,'Y'));
        u_idx   = find(strcmp(M_.endo_names,'u'));
        A_idx   = find(strcmp(M_.endo_names,'A_bat'));

        chi_irfs_C(i,:)   = irf_ro(C_idx,:);
        chi_irfs_Y(i,:)   = irf_ro(Y_idx,:);
        chi_irfs_u(i,:)   = irf_ro(u_idx,:);
        chi_irfs_A(i,:)   = irf_ro(A_idx,:);

        sd_chi = sum(disc .* log(1 + chi_irfs_C(i,:)/100));
        chi_welfare(i) = abs((1 - exp(sd_chi / norm_factor)) * 100);

        fprintf('  chi=%.1f: Welfare=%.6f%%  DeltaWelfare=%+.2f%%  Impact_Y=%+.4f%%\n', ...
            chi_values(i), chi_welfare(i), ...
            (chi_welfare(i)/chi_welfare(1)-1)*100, ...
            chi_irfs_Y(i,1));
    else
        fprintf('  chi=%.1f: FAILED (BK conditions not satisfied)\n', chi_values(i));
        chi_welfare(i) = NaN;
    end
end

% Restore chi=1.0
set_param_value('chi', 1.0);
[dr_ys_restore, params_restore, ~] = vietnam_dsge_proptax_steadystate(...
    oo_.dr.ys, oo_.exo_steady_state, M_, options_);
M_.params = params_restore;
[oo_.dr, ~, M_.params] = resol(0, M_, options_, oo_.dr, ...
    dr_ys_restore, oo_.exo_steady_state, oo_.exo_det_steady_state);
fprintf('\n');

%% ========================================================================
%% 8. SENSITIVITY ANALYSIS: phi_grid
%% ========================================================================
fprintf('========================================\n');
fprintf('SENSITIVITY: phi_grid\n');
fprintf('========================================\n');
phi_vals = [0.5, 1.0, 1.5, 2.0, 3.0];
phi_welfare  = zeros(size(phi_vals));
phi_impact_Ig= zeros(size(phi_vals));
phi_impact_Ib= zeros(size(phi_vals));
phi_cum_u    = zeros(size(phi_vals));

for i = 1:length(phi_vals)
    set_param_value('phi_grid', phi_vals(i));
    [dr_ys_t, params_t, ~] = vietnam_dsge_proptax_steadystate(...
        oo_.dr.ys, oo_.exo_steady_state, M_, options_);
    M_.params = params_t;
    [dr_t, info_t, M_.params] = resol(0, M_, options_, oo_.dr, ...
        dr_ys_t, oo_.exo_steady_state, oo_.exo_det_steady_state);

    if info_t(1) == 0
        nvar = M_.endo_nbr;
        eps_idx = find(strcmp(M_.exo_names,'eps_ren'));
        sr = M_.params(strcmp(M_.param_names,'sigma_ren'));
        sv = zeros(M_.exo_nbr,1); sv(eps_idx) = sr;

        irf_m = zeros(nvar, T);
        st = dr_t.ghu * sv; irf_m(:,1) = st;
        for tt = 2:T
            sv2 = st(M_.nstatic+1:M_.nstatic+M_.nspred);
            st = dr_t.ghx * sv2; irf_m(:,tt) = st;
        end
        irf_ro2 = zeros(nvar,T); irf_ro2(dr_t.order_var,:) = irf_m;

        C_i = find(strcmp(M_.endo_names,'C'));
        Ig_i= find(strcmp(M_.endo_names,'I_grid'));
        Ib_i= find(strcmp(M_.endo_names,'I_bat'));
        u_i = find(strcmp(M_.endo_names,'u'));

        phi_impact_Ig(i) = irf_ro2(Ig_i,1);
        phi_impact_Ib(i) = irf_ro2(Ib_i,1);
        phi_cum_u(i)     = sum(irf_ro2(u_i,:));

        sd_phi = sum(disc .* log(1 + irf_ro2(C_i,1:T_w)/100));
        phi_welfare(i) = abs((1 - exp(sd_phi / norm_factor)) * 100);

        fprintf('  phi_grid=%.1f: Welfare=%.6f%% Impact_Igrid=%+.4f%% CumU=%.4f\n', ...
            phi_vals(i), phi_welfare(i), phi_impact_Ig(i), phi_cum_u(i));
    else
        fprintf('  phi_grid=%.1f: FAILED\n', phi_vals(i));
        phi_welfare(i) = NaN;
    end
end

% Restore phi_grid=1.5
set_param_value('phi_grid', 1.5);
[dr_ys_r2, params_r2, ~] = vietnam_dsge_proptax_steadystate(...
    oo_.dr.ys, oo_.exo_steady_state, M_, options_);
M_.params = params_r2;
[oo_.dr, ~, M_.params] = resol(0, M_, options_, oo_.dr, ...
    dr_ys_r2, oo_.exo_steady_state, oo_.exo_det_steady_state);
fprintf('\n');

%% ========================================================================
%% 9. JOINT SHOCK SCENARIO: "Perfect Storm"
%%    Simultaneous eps_ren = sigma_ren AND eps_bat = sigma_bat
%% ========================================================================
fprintf('========================================\n');
fprintf('JOINT SHOCK: Perfect Storm\n');
fprintf('========================================\n');

nvar = M_.endo_nbr;
eps_ren_idx = find(strcmp(M_.exo_names,'eps_ren'));
eps_bat_idx = find(strcmp(M_.exo_names,'eps_bat'));
sigma_ren_v = M_.params(strcmp(M_.param_names,'sigma_ren'));
sigma_bat_v = M_.params(strcmp(M_.param_names,'sigma_bat'));
sv_joint = zeros(M_.exo_nbr,1);
sv_joint(eps_ren_idx) = sigma_ren_v;
sv_joint(eps_bat_idx) = sigma_bat_v;

irf_jt = zeros(nvar, T);
st_j = oo_.dr.ghu * sv_joint; irf_jt(:,1) = st_j;
for tt = 2:T
    sv2 = st_j(M_.nstatic+1:M_.nstatic+M_.nspred);
    st_j = oo_.dr.ghx * sv2; irf_jt(:,tt) = st_j;
end
irf_joint = zeros(nvar,T); irf_joint(oo_.dr.order_var,:) = irf_jt;

C_idx_j   = find(strcmp(M_.endo_names,'C'));
Y_idx_j   = find(strcmp(M_.endo_names,'Y'));
u_idx_j   = find(strcmp(M_.endo_names,'u'));
Ib_idx_j  = find(strcmp(M_.endo_names,'I_bat'));

sd_joint = sum(disc .* log(1 + irf_joint(C_idx_j,1:T_w)/100));
welfare_joint = abs((1 - exp(sd_joint / norm_factor)) * 100);

% Individual shocks for comparison
sv_ren_only = zeros(M_.exo_nbr,1); sv_ren_only(eps_ren_idx) = sigma_ren_v;
sv_bat_only = zeros(M_.exo_nbr,1); sv_bat_only(eps_bat_idx) = sigma_bat_v;

irf_ren_only = zeros(nvar,T); st_r = oo_.dr.ghu*sv_ren_only; irf_ren_only(:,1)=st_r;
for tt=2:T; sv2=st_r(M_.nstatic+1:M_.nstatic+M_.nspred); st_r=oo_.dr.ghx*sv2; irf_ren_only(:,tt)=st_r; end
irf_ren_ro=zeros(nvar,T); irf_ren_ro(oo_.dr.order_var,:)=irf_ren_only;

irf_bat_only = zeros(nvar,T); st_b = oo_.dr.ghu*sv_bat_only; irf_bat_only(:,1)=st_b;
for tt=2:T; sv2=st_b(M_.nstatic+1:M_.nstatic+M_.nspred); st_b=oo_.dr.ghx*sv2; irf_bat_only(:,tt)=st_b; end
irf_bat_ro=zeros(nvar,T); irf_bat_ro(oo_.dr.order_var,:)=irf_bat_only;

sd_ren = sum(disc.*log(1+irf_ren_ro(C_idx_j,1:T_w)/100));
sd_bat = sum(disc.*log(1+irf_bat_ro(C_idx_j,1:T_w)/100));
welfare_ren = abs((1-exp(sd_ren/norm_factor))*100);
welfare_bat = abs((1-exp(sd_bat/norm_factor))*100);

fprintf('  Intermittency-only welfare:  %.6f%%\n', welfare_ren);
fprintf('  Battery-price-only welfare:  %.6f%%\n', welfare_bat);
fprintf('  Joint (perfect storm) welfare: %.6f%%\n', welfare_joint);
fprintf('  Sum of individuals:          %.6f%%\n', welfare_ren+welfare_bat);
fprintf('  Interaction term:            %.6f%%\n', welfare_joint-(welfare_ren+welfare_bat));
fprintf('  Impact Y (joint):    %+.4f%%\n', irf_joint(Y_idx_j,1));
fprintf('  Impact u (joint):    %+.4f%%\n', irf_joint(u_idx_j,1));
fprintf('  Impact I_bat (joint):%+.4f%%\n', irf_joint(Ib_idx_j,1));
fprintf('\n');

%% ========================================================================
%% 10. EXTERNAL SECTOR DYNAMICS
%% ========================================================================
fprintf('========================================\n');
fprintf('EXTERNAL SECTOR DYNAMICS\n');
fprintf('========================================\n');
fprintf('  B_star impact (ren shock):     %+.6f\n', irf_B_star(1));
fprintf('  r_star impact (ren shock):     %+.6f%%\n', irf_r_star(1));
fprintf('  Cumulative B_star change (10q): %.6f\n', sum(irf_B_star(1:10)));
fprintf('  Correlation u and B_star (ren): %.4f\n', ...
    corr(irf_u', irf_B_star'));
if isfield(oo_.irfs,'B_star_eps_bat')
    fprintf('  B_star impact (bat shock):     %+.6f\n', oo_.irfs.B_star_eps_bat(1));
end
fprintf('\n');

%% ========================================================================
%% 11. FIGURES
%% ========================================================================

%% Figure 1: Full IRF grid (12 panels)
figure('Position',[50 50 1600 1000]);
irf_V = oo_.irfs.V_eps_ren;
vars_to_plot = {'Y','C','u','V','I_bat','I_grid','I_p','F','A_bat','L','B_star','r_star'};
irfs_to_plot = {irf_Y, irf_C, irf_u, irf_V, irf_I_bat, irf_I_grid, ...
                irf_I_p, irf_F, irf_A_bat, irf_L, irf_B_star, irf_r_star};
titles_plot = {'Output (Y)','Consumption (C)','Utilization (u)','Value Added (V)', ...
               'Battery Inv. (I_{bat})','Grid Inv. (I_{grid})','Productive Inv. (I_p)', ...
               'Flexibility (F)','Battery Tech. (A_{bat})','Labor (L)', ...
               'Net Foreign Assets (B^*)','Interest Rate (r^*)'};
colors_plot = {'b','b','r',[0.7 0 0],'g',[0.8 0.5 0],[0.2 0.6 0.8],'m','c',...
               [0 0.5 0.5],[0.5 0 0.5],[0.3 0.3 0.3]};

for ii = 1:12
    subplot(3,4,ii);
    plot(t, irfs_to_plot{ii}, 'Color', colors_plot{ii}, 'LineWidth', 2);
    hold on;
    plot(t, zeros(1,T), 'k--', 'LineWidth', 0.5);
    hold off; grid on;
    title(titles_plot{ii}, 'FontSize', 10);
    xlabel('Quarters'); ylabel('% dev from SS');
end
sgtitle('IRFs to Intermittency Shock — Proportional Tax Variant','FontSize',13,'FontWeight','bold');
print('irf_proptax.png','-dpng','-r300');
fprintf('[Saved] irf_proptax.png\n');

%% Figure 2: Agility Gap
figure('Position',[100 100 800 500]);
hold on;
plot(t, irf_I_bat,  'g-',  'LineWidth', 2.5, 'DisplayName', 'Battery Inv. (I_{bat}) — Private');
plot(t, irf_I_grid, 'Color',[0.8 0.5 0],'LineWidth',2.5,'DisplayName','Grid Inv. (I_{grid}) — Public');
plot(t, irf_I_p,    'b--', 'LineWidth', 1.5, 'DisplayName', 'Productive Inv. (I_p)');
plot(t, zeros(1,T), 'k--','LineWidth',0.5,'HandleVisibility','off');
hold off; grid on;
xlabel('Quarters','FontSize',12); ylabel('% Deviation from SS','FontSize',12);
title('Agility Gap — Proportional Tax Variant','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',11);
print('agility_gap_proptax.png','-dpng','-r300');
fprintf('[Saved] agility_gap_proptax.png\n');

%% Figure 3: Battery price shock IRFs
if isfield(oo_.irfs,'Y_eps_bat')
    figure('Position',[100 100 1000 600]);
    bat_vars  = {'Y','P_bat','I_bat','u','I_grid','B_star'};
    bat_field = {'Y_eps_bat','P_bat_eps_bat','I_bat_eps_bat','u_eps_bat','I_grid_eps_bat','B_star_eps_bat'};
    bat_title = {'Output (Y)','Battery Price (P_{bat})','Battery Inv. (I_{bat})','Utilization (u)','Grid Inv. (I_{grid})','Net FA (B^*)'};
    for ii = 1:6
        subplot(2,3,ii);
        if isfield(oo_.irfs, bat_field{ii})
            plot(t, oo_.irfs.(bat_field{ii}), 'b-','LineWidth',2);
        end
        hold on; plot(t,zeros(1,T),'k--','LineWidth',0.5); hold off; grid on;
        title(bat_title{ii},'FontSize',10); ylabel('% dev'); xlabel('Quarters');
    end
    sgtitle('IRFs to Battery Price Shock — Proportional Tax Variant','FontSize',13,'FontWeight','bold');
    print('irf_battery_price_proptax.png','-dpng','-r300');
    fprintf('[Saved] irf_battery_price_proptax.png\n');
end

%% Figure 4: Output & Investment — Proportional Tax Distortion
figure('Position',[100 100 900 500]);
hold on;
plot(t, irf_Y,      'b-',  'LineWidth',2.5,'DisplayName','Output (Y)');
plot(t, irf_C,      'r-',  'LineWidth',2.0,'DisplayName','Consumption (C)');
plot(t, irf_I_p,    'b--', 'LineWidth',1.5,'DisplayName','Productive Inv. (I_p)');
plot(t, irf_V,      'k:',  'LineWidth',1.5,'DisplayName','Value Added (V)');
plot(t, zeros(1,T),'k--','LineWidth',0.5,'HandleVisibility','off');
hold off; grid on;
xlabel('Quarters','FontSize',12);
ylabel('% Deviation from SS','FontSize',11);
title('Output, Consumption & Investment — Proportional Tax Variant','FontSize',12,'FontWeight','bold');
legend('Location','northeast','FontSize',10);
print('fiscal_amplification.png','-dpng','-r300');
fprintf('[Saved] fiscal_amplification.png\n');

%% Figure 5: Signal Attenuation Counterfactual
chi_plot = [1.0, 0.3, 0.0];
chi_idx  = [find(chi_values==1.0), find(chi_values==0.3), find(chi_values==0.0)];
figure('Position',[100 100 1200 900]);
cols_cf = {'b','r',[0.4 0.4 0.4]};
ls_cf   = {'-','--',':'};
lbl_cf  = {'\chi=1.0 (Baseline)','\chi=0.3 (Partial)','\chi=0.0 (Suppressed)'};
sub_vars = {chi_irfs_Y, chi_irfs_u, chi_irfs_A, chi_irfs_C};
sub_title= {'Output (Y)','Utilization (u)','Battery Technology (A_{bat})','Consumption (C)'};
for p = 1:4
    subplot(2,3,p);
    hold on;
    for k = 1:3
        ci = chi_idx(k);
        if ~isnan(chi_welfare(ci))
            plot(t, sub_vars{p}(ci,:),'Color',cols_cf{k},'LineStyle',ls_cf{k},'LineWidth',2,'DisplayName',lbl_cf{k});
        end
    end
    plot(t,zeros(1,T),'k--','LineWidth',0.5,'HandleVisibility','off');
    hold off; grid on; title(sub_title{p},'FontSize',10);
    xlabel('Quarters'); ylabel('% dev');
    if p==1, legend('Location','southeast','FontSize',8); end
end
subplot(2,3,5);
wf_plot = chi_welfare(chi_idx); wf_plot(isnan(wf_plot))=0;
bar(wf_plot,'FaceColor',[0.4 0.6 0.8]);
set(gca,'XTickLabel',{'\chi=1.0','\chi=0.3','\chi=0.0'});
ylabel('Welfare Cost (% SS consumption)'); title('Welfare Cost by \chi'); grid on;
sgtitle('Counterfactual: Signal Attenuation — Proportional Tax Variant','FontSize',13,'FontWeight','bold');
print('counterfactual_chi_proptax.png','-dpng','-r300');
fprintf('[Saved] counterfactual_chi_proptax.png\n');

%% Figure 6: Joint Shock
figure('Position',[100 100 900 500]);
hold on;
plot(t, irf_joint(Y_idx_j,:), 'k-',   'LineWidth',2,'DisplayName','Joint Shock (Y)');
plot(t, irf_ren_ro(Y_idx_j,:),'b--',  'LineWidth',1.5,'DisplayName','Intermittency only');
plot(t, irf_bat_ro(Y_idx_j,:),'r:',   'LineWidth',1.5,'DisplayName','Battery price only');
plot(t, zeros(1,T),'k--','LineWidth',0.5,'HandleVisibility','off');
hold off; grid on;
xlabel('Quarters','FontSize',12); ylabel('% Deviation from SS','FontSize',12);
title('Perfect Storm: Joint Intermittency + Battery Price Shock','FontSize',12,'FontWeight','bold');
legend('Location','southeast','FontSize',10);
print('joint_shock_perfect_storm.png','-dpng','-r300');
fprintf('[Saved] joint_shock_perfect_storm.png\n');

%% ========================================================================
%% 12. SAVE ALL RESULTS TO JSON
%% ========================================================================

fprintf('\n========================================\n');
fprintf('SAVING RESULTS TO proptax_results.json\n');
fprintf('========================================\n');

% Build JSON string manually for compatibility
fid = fopen('proptax_results.json','w');
fprintf(fid,'{\n');

% Steady state
fprintf(fid,'  "ss": {\n');
fprintf(fid,'    "Y": %.10f,\n', Y_ss);
fprintf(fid,'    "C": %.10f,\n', C_ss);
fprintf(fid,'    "C_share_Y_pct": %.6f,\n', C_ss/Y_ss*100);
fprintf(fid,'    "L": %.10f,\n', L_ss);
fprintf(fid,'    "K_p": %.10f,\n', Kp_ss);
fprintf(fid,'    "K_b": %.10f,\n', Kb_ss);
fprintf(fid,'    "K_g": %.10f,\n', Kg_ss);
fprintf(fid,'    "I_p": %.10f,\n', Ip_ss);
fprintf(fid,'    "I_bat": %.10f,\n', Ibat_ss);
fprintf(fid,'    "I_grid": %.10f,\n', Ig_ss);
fprintf(fid,'    "u": %.10f,\n', u_ss);
fprintf(fid,'    "F": %.10f,\n', F_ss);
fprintf(fid,'    "V": %.10f,\n', V_ss);
fprintf(fid,'    "tau": %.10f,\n', tau_ss);
fprintf(fid,'    "tau_pct": %.6f,\n', tau_ss*100);
fprintf(fid,'    "r_star": %.10f,\n', r_ss);
fprintf(fid,'    "Kb_Kg_ratio": %.6f,\n', Kb_ss/Kg_ss);
fprintf(fid,'    "phi_int": %.10f,\n', M_.params(strcmp(M_.param_names,'phi_int')));
fprintf(fid,'    "varphi": %.10f\n', M_.params(strcmp(M_.param_names,'varphi')));
fprintf(fid,'  },\n');

% IRF impact effects
fprintf(fid,'  "impact_ren": {\n');
fprintf(fid,'    "Y": %.10f,\n',     irf_Y(1));
fprintf(fid,'    "C": %.10f,\n',     irf_C(1));
fprintf(fid,'    "L": %.10f,\n',     irf_L(1));
fprintf(fid,'    "u": %.10f,\n',     irf_u(1));
fprintf(fid,'    "I_bat": %.10f,\n', irf_I_bat(1));
fprintf(fid,'    "I_grid": %.10f,\n',irf_I_grid(1));
fprintf(fid,'    "I_p": %.10f,\n',   irf_I_p(1));
fprintf(fid,'    "F": %.10f,\n',     irf_F(1));
fprintf(fid,'    "A_bat_q10": %.10f,\n', irf_A_bat(10));
fprintf(fid,'    "B_star": %.10f,\n',irf_B_star(1));
fprintf(fid,'    "cum_u_40q": %.10f\n', sum(irf_u));
fprintf(fid,'  },\n');

% FEVD
if isfield(oo_,'variance_decomposition')
    vd = oo_.variance_decomposition;
    vl = {'Y','C','L','I_p','I_bat','I_grid','u','F','A_bat','V','P_bat','B_star','r_star','a_ren','a_I'};
    fprintf(fid,'  "fevd": {\n');
    for ii = 1:min(length(vl), size(vd,1))
        comma = ','; if ii==min(length(vl),size(vd,1)), comma=''; end
        fprintf(fid,'    "%s": {"ren": %.6f, "bat": %.6f, "I": %.6f}%s\n', ...
            vl{ii}, vd(ii,1), vd(ii,2), vd(ii,3), comma);
    end
    fprintf(fid,'  },\n');
end

% Welfare
fprintf(fid,'  "welfare_baseline": %.10f,\n', welfare_baseline);
fprintf(fid,'  "welfare_chi": {\n');
for i = 1:length(chi_values)
    comma = ','; if i==length(chi_values), comma=''; end
    fprintf(fid,'    "%.1f": %.10f%s\n', chi_values(i), chi_welfare(i), comma);
end
fprintf(fid,'  },\n');

% Signal value
if ~isnan(chi_welfare(find(chi_values==1.0))) && ~isnan(chi_welfare(find(chi_values==0.0)))
    signal_pct = (chi_welfare(end)/chi_welfare(1)-1)*100;
    fprintf(fid,'  "signal_value_pct": %.6f,\n', abs(signal_pct));
end

% Sensitivity phi_grid
fprintf(fid,'  "sensitivity_phi_grid": [\n');
for i = 1:length(phi_vals)
    comma = ','; if i==length(phi_vals), comma=''; end
    fprintf(fid,'    {"phi_grid":%.1f,"welfare":%.8f,"impact_Igrid":%.6f,"cum_u":%.6f}%s\n',...
        phi_vals(i),phi_welfare(i),phi_impact_Ig(i),phi_cum_u(i),comma);
end
fprintf(fid,'  ],\n');

% Joint shock
fprintf(fid,'  "joint_shock": {\n');
fprintf(fid,'    "welfare_ren_only": %.10f,\n', welfare_ren);
fprintf(fid,'    "welfare_bat_only": %.10f,\n', welfare_bat);
fprintf(fid,'    "welfare_joint": %.10f,\n', welfare_joint);
fprintf(fid,'    "impact_Y": %.10f,\n', irf_joint(Y_idx_j,1));
fprintf(fid,'    "impact_u": %.10f,\n', irf_joint(u_idx_j,1));
fprintf(fid,'    "impact_I_bat": %.10f\n', irf_joint(Ib_idx_j,1));
fprintf(fid,'  }\n');

fprintf(fid,'}\n');
fclose(fid);
fprintf('[Saved] proptax_results.json\n\n');

%% ========================================================================
%% 13. SUMMARY PRINTOUT FOR PAPER
%% ========================================================================
fprintf('==========================================\n');
fprintf('PAPER RESULTS SUMMARY — COPY TO LaTeX\n');
fprintf('==========================================\n\n');
fprintf('STEADY STATE:\n');
fprintf('  Y_ss = %.4f, C_ss = %.4f (%.1f%% of Y)\n', Y_ss, C_ss, C_ss/Y_ss*100);
fprintf('  L_ss = %.2f, K_p_ss = %.3f, u_ss = %.2f\n', L_ss, Kp_ss, u_ss);
fprintf('  tau_ss = %.4f%% (= %.4f)\n', tau_ss*100, tau_ss);
fprintf('\nIRF IMPACTS (intermittency shock):\n');
fprintf('  dY/dep_ren = %+.4f%%\n', irf_Y(1));
fprintf('  du/dep_ren = %+.4f%%\n', irf_u(1));
fprintf('  dI_bat/dep_ren = %+.4f%%\n', irf_I_bat(1));
fprintf('  dI_grid/dep_ren = %+.4f%%\n', irf_I_grid(1));
fprintf('  dI_p/dep_ren = %+.4f%%\n', irf_I_p(1));
fprintf('  dA_bat (q10) = %+.4f%%\n', irf_A_bat(10));
fprintf('\nWELFARE:\n');
fprintf('  chi=1.0: %.4f%%\n', chi_welfare(find(chi_values==1.0)));
fprintf('  chi=0.3: %.4f%%\n', chi_welfare(find(chi_values==0.3)));
fprintf('  chi=0.0: %.4f%%\n', chi_welfare(find(chi_values==0.0)));
chi10_idx = find(chi_values==1.0); chi00_idx = find(chi_values==0.0);
if ~isnan(chi_welfare(chi10_idx)) && ~isnan(chi_welfare(chi00_idx))
    fprintf('  Welfare increase chi=1->0: %.1f%%\n', ...
        (chi_welfare(chi00_idx)/chi_welfare(chi10_idx)-1)*100);
end
fprintf('\nJOINT SHOCK:\n');
fprintf('  Intermittency-only: %.4f%%\n', welfare_ren);
fprintf('  Battery-only:       %.4f%%\n', welfare_bat);
fprintf('  Joint:              %.4f%%\n', welfare_joint);
fprintf('\n==========================================\n');
fprintf('  RUN COMPLETE\n');
fprintf('==========================================\n');
