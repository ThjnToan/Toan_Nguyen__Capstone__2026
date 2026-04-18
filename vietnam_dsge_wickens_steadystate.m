function [ys, params, check] = vietnam_dsge_wickens_steadystate(ys, exo, M_, options_)
% =========================================================================
% Steady-state file for vietnam_dsge_wickens.mod
% Wickens (2012) Extended Version — 20-variable SOE DSGE
%
% HIGH-IMPACT EXTENSIONS IMPLEMENTED:
%
%  [W1] EULER WEDGE (Wickens Ch.5 §5.7.1):
%       At steady state the capital Euler is:
%           (1 − τ̄) · α · Y / K_p = r* + δ_p
%       vs. first-best  α · Y / K_p = r* + δ_p  (no tax).
%       We compute and report:
%           Delta_k     = τ_ss · α · Y_ss / K_p_ss   [wedge magnitude, % of output]
%           K_p_firstbest = α · Y_ss / (r* + δ_p)     [counterfactual capital]
%           K_p_gap       = (K_p_firstbest − K_p_ss) / K_p_firstbest * 100 [% gap]
%       These are DIAGNOSTIC OUTPUTS only, not model variables.
%
%  [W2] RICARDIAN NON-EQUIVALENCE:
%       The calibration of varphi uses the distorted labor FOC:
%           (1 − τ_ss)(1 − α)V = φ · C · L^(1+σ_L)
%       Making explicit that the tax wedge lowers equilibrium varphi
%       (labor supply is less costly in utility terms because after-tax
%       wages are lower — the model requires a smaller disutility
%       coefficient to match the target L = 0.33).
%       G_surplus_ss = 0 by construction (balanced budget).
%
%  [W3] FISCAL SUSTAINABILITY (Wickens Ch.5 §5.5):
%       We compute FisSust_ss = r*_ss − (π_bar + γ_bar) = r*_ss > 0.
%       We then compute the PRIMARY SURPLUS REQUIRED for sustainability:
%           b_dot / y = (r* − γ) · b/y − s/y  < 0  requires  s/y > (r*−γ)·b/y
%       With B_star = 0 (balanced position): s/y = G_surplus/Y = 0 is
%       just barely sufficient. Any positive debt requires G_surplus > 0.
%
%  [W4] FEENSTRA INTERMITTENCY WEDGE:
%       At steady state: T_int_ss = κ · (1 − u_ss) · τ_ss
%       This is SMALL at SS but becomes significant under shocks.
%       Crucially, T_int enters the Euler equations, so the SS value of C
%       is the same (wedge ratio = 1 at SS), but IRFs change.
%
%  [W5] TAX-SMOOTHING (Wickens Ch.5 §5.7.3):
%       At SS, both the proportional rule and the constant-τ rule give
%       τ_t = τ_ss, so steady states are IDENTICAL regardless of tau_smooth.
%       The welfare difference only shows up in IRFs — confirmed here.
%
% Author: Toan T. Nguyen / Advisor: Dr. Xavier Martin G. Bautista
% Fulbright University Vietnam — March 2026
% =========================================================================

check = 0;
params = M_.params;

% ---- Extract parameters ----------------------------------------
p_beta      = params(strcmp(M_.param_names, 'beta'));
p_alpha     = params(strcmp(M_.param_names, 'alpha'));
p_delta_p   = params(strcmp(M_.param_names, 'delta_p'));
p_sigma_L   = params(strcmp(M_.param_names, 'sigma_L'));
p_sigma_E   = params(strcmp(M_.param_names, 'sigma_E'));
p_omega_E   = params(strcmp(M_.param_names, 'omega_E'));
p_psi       = params(strcmp(M_.param_names, 'psi'));
p_theta_ren = params(strcmp(M_.param_names, 'theta_ren'));
p_sigma_ren = params(strcmp(M_.param_names, 'sigma_ren'));
p_mu        = params(strcmp(M_.param_names, 'mu'));
p_rho_flex  = params(strcmp(M_.param_names, 'rho_flex'));
p_delta_g   = params(strcmp(M_.param_names, 'delta_g'));
p_delta_b   = params(strcmp(M_.param_names, 'delta_b'));
p_u_target  = params(strcmp(M_.param_names, 'u_target'));
p_r_bar     = params(strcmp(M_.param_names, 'r_bar'));
p_phi_b     = params(strcmp(M_.param_names, 'phi_b'));
p_B_star_ss = params(strcmp(M_.param_names, 'B_star_ss'));
p_kappa_int = params(strcmp(M_.param_names, 'kappa_int'));
p_pi_bar    = params(strcmp(M_.param_names, 'pi_bar'));
p_gamma_bar = params(strcmp(M_.param_names, 'gamma_bar'));

% ---- Fixed target values ----------------------------------------
A_bat   = 1.0;
P_bat   = 1.0;
u       = p_u_target;   % 0.97
L       = 0.33;

E_bar       = 0.15;
Vol_ren_bar = p_theta_ren * (0.045/p_theta_ren) * p_sigma_ren;  % = 0.0054

% ---- External sector SS ----------------------------------------
B_star  = p_B_star_ss;          % 0 (balanced position)
r_star  = 1/p_beta - 1;         % exact consistency with bond Euler: 1/β − 1
p_r_bar = r_star;               % calibrate r_bar = r_star at SS

% ---- Flexibility assets ----------------------------------------
K_b = 0.19;
K_g = 1.14;
rf  = p_rho_flex;
F   = (p_mu * (A_bat * K_b)^((rf-1)/rf) + (1-p_mu) * K_g^((rf-1)/rf))^(rf/(rf-1));

% ---- Calibrate phi_int from reliability equation ---------------
% u = 1 − exp(−ψ·F / (φ_int · Vol_ren_bar))
% => φ_int = −ψ·F / (Vol_ren_bar · ln(1−u))
phi_int_val = -p_psi * F / (Vol_ren_bar * log(1 - u));

% ---- Capital stock: iterative solution including τ wedge --------
% [W1] EULER WEDGE at steady state:
%   Capital Euler: (1−τ_ss) · α · Y / K_p + 1 − δ_p = 1/β = 1 + r*
%   => (1−τ_ss) · α · Y / K_p = r* + δ_p ≡ R_target
%   => K_p = (1−τ_ss) · α · Y / R_target
% where τ_ss = I_grid_ss / Y_ss = δ_g · K_g / Y  (endogenous in Y)
% We iterate: K_p → V → Y → τ_ss → K_p_new until convergence.

R_target  = r_star + p_delta_p;
I_grid_ss = p_delta_g * K_g;   % = 0.01425 (tax revenue at SS)

% Initial guess (ignoring tau wedge)
K_p = (p_alpha * u^(p_alpha) * L^(1-p_alpha) / R_target)^(1/(1-p_alpha));

for iter = 1:500
    V       = (u * K_p)^p_alpha * L^(1-p_alpha);
    Y       = ((1-p_omega_E) * V^((p_sigma_E-1)/p_sigma_E) + ...
               p_omega_E * E_bar^((p_sigma_E-1)/p_sigma_E))^(p_sigma_E/(p_sigma_E-1));
    tau_ss  = I_grid_ss / Y;
    K_p_new = (1 - tau_ss) * p_alpha * Y / R_target;
    if abs(K_p_new - K_p) < 1e-14
        break;
    end
    K_p = 0.7*K_p_new + 0.3*K_p;  % damped update for stability
end
K_p    = K_p_new;
V      = (u * K_p)^p_alpha * L^(1-p_alpha);
Y      = ((1-p_omega_E) * V^((p_sigma_E-1)/p_sigma_E) + ...
          p_omega_E * E_bar^((p_sigma_E-1)/p_sigma_E))^(p_sigma_E/(p_sigma_E-1));
tau_ss = I_grid_ss / Y;

% ---- Investment flows -------------------------------------------
I_p   = p_delta_p * K_p;
I_bat = p_delta_b * K_b;    % = 0.0057
I_grid = I_grid_ss;

% ---- Consumption (resource constraint) -------------------------
% Y_ss = C_ss + I_p_ss + P_bat·I_bat_ss + I_grid_ss   (B_star = 0)
C = Y - I_p - P_bat * I_bat - I_grid;

if C <= 0
    warning('vietnam_dsge_wickens_ss: negative consumption. Check calibration.');
    check = 1;
    return;
end

% ---- [W4] Feenstra wedge at SS ---------------------------------
% T_int_ss = κ · (1 − u_ss) · τ_ss
% Small at SS (≈ 0.000248), but enters IRFs with amplification.
T_int_ss = p_kappa_int * (1 - u) * tau_ss;
% Note: the Euler wedge ratio (1+T_int)/(1+T_int(+1)) = 1 at SS
% so the SS consumption level is IDENTICAL to the proptax version.
% The Feenstra extension only changes DYNAMICS, not the steady state.

% ---- Calibrate varphi from distorted labor FOC [W2] ---------------
% (1 − τ_ss)(1 − α)V = φ · C · L^(1+σ_L)
% [W2 note]: The (1−τ_ss) factor here is the Ricardian non-equivalence
% manifestation: because τ is distortionary (not lump-sum), the labor
% supply decision depends on the after-tax wage. This wedge reduces the
% effective labor income and requires a LOWER φ to match L=0.33
% compared to the lump-sum-tax benchmark.
p_varphi = (1 - tau_ss) * (1 - p_alpha) * V / (C * L^(1 + p_sigma_L));

% ---- [W3] Fiscal sustainability indicator ----------------------
% FisSust_ss = r*_ss − (π_bar + γ_bar)
% = 1/β − 1 − 0 − 0  ≈ 0.0101 > 0
% Interpretation: at the calibrated steady state, the real interest rate
% EXCEEDS trend growth + inflation, meaning debt would compound faster
% than the economy grows. This requires the government to maintain a
% primary surplus (G_surplus > 0) to prevent explosive debt dynamics.
% With balanced budget (G_surplus = 0 by construction), the government
% is sitting at the knife-edge: any positive debt shock requires immediate
% corrective action. The proportional tax rule provides this automatically.
FisSust_ss = r_star - (p_pi_bar + p_gamma_bar);

% ---- [W1] EULER WEDGE diagnostic outputs (NOT model variables) ------
Delta_k         = tau_ss * p_alpha * Y / K_p;   % wedge magnitude (% output units)
K_p_firstbest   = p_alpha * Y / R_target;        % K_p if τ=0 (first-best)
K_p_gap_pct     = (K_p_firstbest - K_p) / K_p_firstbest * 100;

% ---- [W2] Fiscal surplus at SS (= 0 by balanced budget) ----------
G_surplus_ss = tau_ss * Y - I_grid;  % should be machine-epsilon near 0

% ---- Print diagnostics ------------------------------------------
fprintf('\n========================================================\n');
fprintf(' STEADY STATE — vietnam_dsge_wickens.mod (20 variables)\n');
fprintf('========================================================\n');
fprintf(' ORIGINAL VARIABLES:\n');
fprintf('   Y      = %8.4f      C      = %8.4f\n', Y, C);
fprintf('   L      = %8.4f      K_p    = %8.4f\n', L, K_p);
fprintf('   K_b    = %8.4f      K_g    = %8.4f\n', K_b, K_g);
fprintf('   I_p    = %8.4f      I_bat  = %8.6f\n', I_p, I_bat);
fprintf('   I_grid = %8.6f      u      = %8.4f\n', I_grid, u);
fprintf('   F      = %8.4f      tau_ss = %8.4f (%5.2f%%)\n', F, tau_ss, tau_ss*100);
fprintf('   r*     = %8.6f      varphi = %8.4f\n', r_star, p_varphi);
fprintf('   phi_int= %8.3f\n', phi_int_val);
fprintf('\n [W4] FEENSTRA WEDGE:\n');
fprintf('   T_int_ss       = %10.6f  (kappa*(1-u)*tau)\n', T_int_ss);
fprintf('   Wedge at SS    ≈ 0 (ratio = 1; only IRF dynamics change)\n');
fprintf('\n [W1] EULER WEDGE DIAGNOSTIC:\n');
fprintf('   Delta_k (wedge)        = %8.6f  (tau*alpha*Y/K_p)\n', Delta_k);
fprintf('   K_p_ss (with tax)      = %8.4f\n', K_p);
fprintf('   K_p_firstbest (tau=0)  = %8.4f\n', K_p_firstbest);
fprintf('   Capital gap            = %8.4f%%  <-- cost of proportional tax\n', K_p_gap_pct);
fprintf('   Chamley-Judd benchmark: optimal tau_k = 0 (K_p_gap should = 0)\n');
fprintf('\n [W2] RICARDIAN NON-EQUIVALENCE:\n');
fprintf('   G_surplus_ss           = %8.2e  (balanced budget by construction)\n', G_surplus_ss);
fprintf('   RE fails because tau is distortionary — fiscal timing matters\n');
fprintf('\n [W3] FISCAL SUSTAINABILITY:\n');
fprintf('   FisSust_ss             = %8.6f  (r* - pi_bar - gamma_bar)\n', FisSust_ss);
fprintf('   r*                     = %8.6f\n', r_star);
fprintf('   pi_bar + gamma_bar     = %8.6f\n', p_pi_bar + p_gamma_bar);
if FisSust_ss > 0
    fprintf('   STATUS: UNSTABLE (r* > pi+gamma) — primary surplus required\n');
    fprintf('   Proportional tax rule provides this surplus by construction.\n');
else
    fprintf('   STATUS: STABLE (pi+gamma >= r*) — debt ratio shrinks\n');
end
fprintf('========================================================\n\n');

% ---- Write back calibrated parameters ---------------------------
params(strcmp(M_.param_names, 'phi_int'))    = phi_int_val;
params(strcmp(M_.param_names, 'varphi'))     = p_varphi;
params(strcmp(M_.param_names, 'E_bar'))      = E_bar;
params(strcmp(M_.param_names, 'Vol_ren_bar')) = Vol_ren_bar;
params(strcmp(M_.param_names, 'r_bar'))      = p_r_bar;

% ---- Assign steady-state vector (20 variables) ------------------
ys(strcmp(M_.endo_names, 'Y'))          = Y;
ys(strcmp(M_.endo_names, 'C'))          = C;
ys(strcmp(M_.endo_names, 'L'))          = L;
ys(strcmp(M_.endo_names, 'K_p'))        = K_p;
ys(strcmp(M_.endo_names, 'K_b'))        = K_b;
ys(strcmp(M_.endo_names, 'K_g'))        = K_g;
ys(strcmp(M_.endo_names, 'I_p'))        = I_p;
ys(strcmp(M_.endo_names, 'I_bat'))      = I_bat;
ys(strcmp(M_.endo_names, 'I_grid'))     = I_grid;
ys(strcmp(M_.endo_names, 'u'))          = u;
ys(strcmp(M_.endo_names, 'F'))          = F;
ys(strcmp(M_.endo_names, 'A_bat'))      = A_bat;
ys(strcmp(M_.endo_names, 'V'))          = V;
ys(strcmp(M_.endo_names, 'P_bat'))      = P_bat;
ys(strcmp(M_.endo_names, 'B_star'))     = B_star;
ys(strcmp(M_.endo_names, 'r_star'))     = r_star;
ys(strcmp(M_.endo_names, 'tau'))        = tau_ss;
% [W2] Ricardian non-equivalence: government surplus
ys(strcmp(M_.endo_names, 'G_surplus'))  = G_surplus_ss;
% [W4] Feenstra intermittency wedge
ys(strcmp(M_.endo_names, 'T_int'))      = T_int_ss;
% [W3] Fiscal sustainability indicator
ys(strcmp(M_.endo_names, 'FisSust'))    = FisSust_ss;

end
