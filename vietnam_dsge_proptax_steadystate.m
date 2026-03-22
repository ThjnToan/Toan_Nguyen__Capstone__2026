function [ys,params,check] = vietnam_dsge_proptax_steadystate(ys,exo,M_,options_)
% Steady state file for vietnam_dsge_proptax.mod (17-variable SOE version)
% Proportional tax variant: tau_t = I_grid_t / Y_t
% Calibrates varphi and phi_int to hit L=0.33 and u=0.97
% Key difference from baseline: varphi calibration uses (1-tau_ss) wedge
% in labor FOC: (1-tau_ss)*(1-alpha)*V = varphi*C*L^(1+sigma_L)

check = 0;
params = M_.params;

% ---- Extract parameters ----
p_beta     = params(strcmp(M_.param_names, 'beta'));
p_alpha    = params(strcmp(M_.param_names, 'alpha'));
p_delta_p  = params(strcmp(M_.param_names, 'delta_p'));
p_sigma_L  = params(strcmp(M_.param_names, 'sigma_L'));
p_sigma_E  = params(strcmp(M_.param_names, 'sigma_E'));
p_omega_E  = params(strcmp(M_.param_names, 'omega_E'));
p_psi      = params(strcmp(M_.param_names, 'psi'));
p_theta_ren = params(strcmp(M_.param_names, 'theta_ren'));
p_sigma_ren = params(strcmp(M_.param_names, 'sigma_ren'));
p_mu       = params(strcmp(M_.param_names, 'mu'));
p_rho_flex = params(strcmp(M_.param_names, 'rho_flex'));
p_delta_g  = params(strcmp(M_.param_names, 'delta_g'));
p_delta_b  = params(strcmp(M_.param_names, 'delta_b'));
p_u_target = params(strcmp(M_.param_names, 'u_target'));
p_r_bar    = params(strcmp(M_.param_names, 'r_bar'));
p_phi_b    = params(strcmp(M_.param_names, 'phi_b'));
p_B_star_ss = params(strcmp(M_.param_names, 'B_star_ss'));

% ---- Fixed values ----
A_bat = 1.0;
P_bat = 1.0;
u = p_u_target;  % 0.97
L = 0.33;

E_bar = 0.15;
Vol_ren_bar = p_theta_ren * (0.045/p_theta_ren) * p_sigma_ren;  % 0.0054

% External sector steady state
B_star = p_B_star_ss;  % 0
r_star = 1/p_beta - 1; % exact consistency with bond Euler
p_r_bar = r_star;

% Flexibility assets (same as baseline)
K_b = 0.19;
K_g = 1.14;
rf = p_rho_flex;
F = (p_mu * (A_bat * K_b)^((rf-1)/rf) + (1-p_mu) * K_g^((rf-1)/rf))^(rf/(rf-1));

% Calibrate phi_int from reliability equation
phi_int_val = -p_psi * F / (Vol_ren_bar * log(1 - u));

% Capital from no-arbitrage: (1-tau_ss)*alpha*Y/K_p + 1 - delta_p = 1 + r_star
% => (1-tau_ss)*alpha*Y/K_p = r_star + delta_p
% We iterate jointly since tau_ss = I_grid_ss/Y_ss = delta_g*K_g/Y
R_target = r_star + p_delta_p;
I_grid_ss = p_delta_g * K_g;  % = 0.01425

% Initial guess: K_p ignoring tau
K_p = (p_alpha * u^(p_alpha) * L^(1-p_alpha) / R_target)^(1/(1-p_alpha));

for iter = 1:200
    V = (u * K_p)^p_alpha * L^(1-p_alpha);
    Y = ((1-p_omega_E) * V^((p_sigma_E-1)/p_sigma_E) + p_omega_E * E_bar^((p_sigma_E-1)/p_sigma_E))^(p_sigma_E/(p_sigma_E-1));
    tau_ss = I_grid_ss / Y;
    % Capital Euler at SS: (1-tau_ss)*alpha*Y/K_p = r_star + delta_p
    K_p_new = (1 - tau_ss) * p_alpha * Y / R_target;
    if abs(K_p_new - K_p) < 1e-12
        break;
    end
    K_p = K_p_new;
end
K_p = K_p_new;
V = (u * K_p)^p_alpha * L^(1-p_alpha);
Y = ((1-p_omega_E) * V^((p_sigma_E-1)/p_sigma_E) + p_omega_E * E_bar^((p_sigma_E-1)/p_sigma_E))^(p_sigma_E/(p_sigma_E-1));
tau_ss = I_grid_ss / Y;

% Investment flows
I_p = p_delta_p * K_p;
I_bat = p_delta_b * K_b;
I_grid = I_grid_ss;

% Consumption from resource constraint
% At SS with B_star=0: Y = C + I_p + P_bat*I_bat + I_grid
C = Y - I_p - P_bat * I_bat - I_grid;

if C <= 0
    check = 1;
    return;
end

% Calibrate varphi from MODIFIED labor FOC:
% (1-tau_ss)*(1-alpha)*V = varphi*C*L^(1+sigma_L)
p_varphi = (1 - tau_ss) * (1-p_alpha) * V / (C * L^(1+p_sigma_L));

% ---- Write back calibrated parameters ----
params(strcmp(M_.param_names, 'phi_int')) = phi_int_val;
params(strcmp(M_.param_names, 'varphi')) = p_varphi;
params(strcmp(M_.param_names, 'E_bar')) = E_bar;
params(strcmp(M_.param_names, 'Vol_ren_bar')) = Vol_ren_bar;
params(strcmp(M_.param_names, 'r_bar')) = p_r_bar;

% ---- Assign steady state vector (19 variables: 17 original + a_ren, a_I) ----
ys(strcmp(M_.endo_names, 'Y')) = Y;
ys(strcmp(M_.endo_names, 'C')) = C;
ys(strcmp(M_.endo_names, 'L')) = L;
ys(strcmp(M_.endo_names, 'K_p')) = K_p;
ys(strcmp(M_.endo_names, 'K_b')) = K_b;
ys(strcmp(M_.endo_names, 'K_g')) = K_g;
ys(strcmp(M_.endo_names, 'I_p')) = I_p;
ys(strcmp(M_.endo_names, 'I_bat')) = I_bat;
ys(strcmp(M_.endo_names, 'I_grid')) = I_grid;
ys(strcmp(M_.endo_names, 'u')) = u;
ys(strcmp(M_.endo_names, 'F')) = F;
ys(strcmp(M_.endo_names, 'A_bat')) = A_bat;
ys(strcmp(M_.endo_names, 'V')) = V;
ys(strcmp(M_.endo_names, 'P_bat')) = P_bat;
ys(strcmp(M_.endo_names, 'B_star')) = B_star;
ys(strcmp(M_.endo_names, 'r_star')) = r_star;
ys(strcmp(M_.endo_names, 'tau')) = tau_ss;
ys(strcmp(M_.endo_names, 'a_ren')) = 0;   % intermittency AR(1) = 0 in SS
ys(strcmp(M_.endo_names, 'a_I'))   = 0;   % implementation AR(1) = 0 in SS

end
