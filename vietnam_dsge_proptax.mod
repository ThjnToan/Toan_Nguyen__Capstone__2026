// =========================================================================
// VIETNAM RENEWABLE ENERGY TRANSITION DSGE MODEL
// Small Open Economy with Exponential Reliability + Scarcity-Induced Technological Change
// ** PROPORTIONAL TAX VARIANT **
// Fiscal financing via proportional output tax tau_t = I_grid_t / Y_t
// rather than lump-sum taxation. This introduces a fiscal amplification
// channel: intermittency shocks raise grid investment, which raises tau_t,
// which distorts both the capital Euler and labor supply conditions.
//
// Author:  Toan T. Nguyen
// Advisor: Dr. Xavier Martin G. Bautista
// Fulbright University Vietnam, February 2026
//
// =========================================================================
// DYNARE SYNTAX PRIMER (for better readings for non-Dynarers)
// -------------------------------------------------------------------------
// - X(-1) refers to X_{t-1}  (last period's value)
// - X(+1) refers to X_{t+1}  (next period's expected value)
// - X     refers to X_t      (current period)
// - steady_state(X) is the deterministic steady-state value of X
// - exp(eps_ren) converts a log-normal shock into a level multiplier
// - Parameters phi_int and varphi are set to 1.0 here as placeholders;
//   they are recalibrated endogenously by vietnam_dsge_steadystate.m
//   to hit the targets L_ss = 0.33 and u_ss = 0.97.
// =========================================================================

// =========================================================================
// ENDOGENOUS VARIABLES (19 total — proportional tax variant adds tau; AR(1) states add a_ren, a_I)
// =========================================================================
// See Table "Equilibrium Definition" in the paper for the full mapping.
var
    Y           // Y_t       : Aggregate output (CES of value-added and energy)
    C           // C_t       : Household consumption
    L           // L_t       : Labor supply
    K_p         // K_{p,t}   : Productive (physical) capital stock
    K_b         // K_{b,t}   : Battery storage capital stock
    K_g         // K_{g,t}   : Public grid infrastructure capital stock
    I_p         // I_{p,t}   : Investment in productive capital
    I_bat       // I_{bat,t} : Investment in battery storage
    I_grid      // I_{grid,t}: Public grid investment
    u           // u_t       : Grid utilization / reliability (0 to 1)
    F           // F_t       : Aggregate flexibility (CES of battery + grid)
    A_bat       // A_{bat,t} : Battery technology level (scarcity-induced technological change)
    V           // V_t       : Value added (effective capital x labor)
    P_bat       // P_{bat,t} : World battery price (exogenous AR(1))
    B_star      // B^*_t     : Net foreign assets (>0 = creditor)
    r_star      // r^*_t     : Country-specific real interest rate
    tau         // tau_t     : Proportional output tax rate (= I_grid / Y)
    a_ren       // a_{ren,t} : Renewable intermittency process (AR(1) state variable)
    a_I         // a_{I,t}   : Grid implementation shock process (AR(1) state variable)
;

// =========================================================================
// EXOGENOUS SHOCKS (3 total)
// =========================================================================
varexo
    eps_ren     // varepsilon_{ren,t} : Renewable intermittency shock
    eps_bat     // varepsilon_{bat,t} : Global battery price shock
    eps_I       // varepsilon_{I,t}   : Grid investment implementation shock
;

// =========================================================================
// PARAMETERS
// =========================================================================
parameters
    // --- Household preferences ---
    beta            // Discount factor (quarterly)
    sigma_L         // Inverse Frisch elasticity of labor supply
    varphi          // Labor disutility weight [calibrated endogenously]

    // --- Production technology ---
    alpha           // Capital share in Cobb-Douglas value-added
    delta_p         // Depreciation rate of productive capital (quarterly)
    sigma_E         // Elasticity of substitution: value-added vs. energy
    omega_E         // Energy share in CES output aggregator

    // --- Reliability and flexibility ---
    psi             // Sensitivity of reliability to the flexibility ratio
    phi_int         // Grid integration complexity [calibrated endogenously]
    mu              // Battery share in CES flexibility aggregator
    rho_flex        // Elasticity of substitution: battery vs. grid capital
    delta_b         // Depreciation rate of battery capital (quarterly)
    delta_g         // Depreciation rate of grid capital (quarterly)
    u_target        // Target reliability level (government objective)

    // --- Renewable volatility ---
    sigma_ren       // Std. dev. of intermittency shock
    theta_ren       // Renewable capacity share (used to compute Vol_ren)
    rho_ren         // Persistence of intermittency shock AR(1)
    rho_I           // Persistence of grid implementation shock AR(1)

    // --- Scarcity-Induced Technological Change ---
    eta_bat         // Innovation elasticity (speed of technology improvement)
    chi             // Regulatory signal transmission (1 = full, 0 = none)

    // --- Battery price dynamics ---
    rho_bat         // AR(1) persistence of world battery price
    sigma_bat       // Std. dev. of battery price shock

    // --- Government ---
    phi_grid        // Aggressiveness of grid investment response to reliability gap

    // --- External sector (Schmitt-Grohe & Uribe 2003) ---
    r_bar           // World real interest rate (quarterly)
    phi_b           // Elasticity of interest rate premium to debt position
    B_star_ss       // Steady-state net foreign assets (normalization)

    // --- Calibrated aggregates ---
    E_bar           // Energy endowment (exogenous, normalized)
    Vol_ren_bar     // Steady-state renewable output volatility
;

// =========================================================================
// PARAMETER VALUES
// =========================================================================
// See Section 4 (Calibration) in the paper for sources and justification.

// Household preferences
beta     = 0.99;        // => annual discount rate ~4%
alpha    = 0.35;        // capital share (Vietnam national accounts)
delta_p  = 0.025;       // => 10% annual depreciation
sigma_L  = 1.0;         // unit Frisch elasticity
sigma_E  = 0.6;         // energy-value added substitution < 1 (complements)
omega_E  = 0.045;       // energy share in output (~4.5%)
varphi   = 1.0;         // PLACEHOLDER: recalibrated in steadystate file

// Reliability and flexibility
psi      = 2.0;         // reliability sensitivity to flexibility ratio
phi_int  = 1.0;         // PLACEHOLDER: recalibrated in steadystate file
sigma_ren = 0.12;       // intermittency shock std. dev. (Vietnamese solar CV)
theta_ren = 0.30;       // renewable share of generation capacity
mu       = 0.16;        // battery weight in flexibility aggregator
rho_flex = 0.40;        // battery-grid substitution elasticity < 1 (complements)
delta_g  = 0.0125;      // grid depreciation (~5% annual, long-lived)
delta_b  = 0.030;       // battery depreciation (~12% annual, shorter-lived)

// Shock persistence and innovation
rho_ren  = 0.85;        // intermittency shock AR(1) persistence
rho_I    = 0.70;        // grid implementation shock AR(1) persistence
eta_bat  = 0.10;        // innovation elasticity (induced innovation; BloombergNEF experience curves)
chi      = 1.0;         // full signal transmission (baseline)
rho_bat  = 0.90;        // battery price shock persistence
sigma_bat = 0.08;       // battery price shock std. dev.

// Government
phi_grid = 1.5;         // grid investment response elasticity
u_target = 0.97;        // reliability target (97% utilization, EVN target)

// Exogenous calibrated aggregates
E_bar       = 0.15;     // normalized energy endowment
Vol_ren_bar = 0.0054;   // steady-state renewable volatility

// External sector (Schmitt-Grohe & Uribe 2003)
r_bar    = 1/beta - 1;  // world quarterly interest rate (~4% annual)
phi_b    = 0.001;       // debt-elastic premium (small, for stationarity)
B_star_ss = 0;          // balanced trade in steady state

// =========================================================================
// MODEL EQUATIONS (17 equations for 17 endogenous variables)
// =========================================================================
// Each equation is labeled with its counterpart in the LaTeX manuscript.
// The system is exactly identified: 17 equations, 17 unknowns.
// ** PROPORTIONAL TAX VARIANT: tau added as 17th variable/equation **
//
// Classification:
//   Static variables (7):  Y, L, u, F, V, I_bat, I_grid, tau
//                          (tau = I_grid/Y is purely static; no new
//                           eigenvalue — tau(+1) in Euler uses Y(+1)
//                           and I_grid(+1) which are already in the system)
//                          (A_bat is also static: responds to current u,
//                           but F uses A_bat(-1), so no simultaneity loop)
//   Predetermined (7):     K_p(-1), K_b(-1), K_g(-1), A_bat(-1),
//                          P_bat(-1), B_star(-1), r_star(-1 via B_star)
//   Forward-looking (3):   C, Y, r_star
//                          (C and r_star appear with (+1) leads in Euler eqs;
//                           Y appears as Y(+1) in the capital Euler eq [1],
//                           making it a jumper in Dynare's eigenvalue count.
//                           K_p is predetermined, not forward-looking.)
//
// Timing convention for A_bat:
//   F_t uses A_bat(-1): technology induced in t-1 becomes productive in t.
//   A_bat responds to current u_t (eq [10]), updating A_bat,t which is then
//   carried into t+1 as A_bat(-1). This one-period-delay convention
//   eliminates the within-period simultaneity loop A_bat -> F -> u -> A_bat
//   that was present when F used A_bat (current).
//
// Blanchard-Kahn: 3 eigenvalues > 1 for 3 forward-looking variables.
// =========================================================================

model;

// ------- HOUSEHOLD BLOCK -------

// [1] Capital Euler equation (eq:euler in paper) -- PROPORTIONAL TAX VARIANT
// After-tax return to capital: (1-tau(+1)) * alpha*Y(+1)/K_p + (1-delta_p).
// tau(+1) = I_grid(+1)/Y(+1) is forward-looking but determined by
// the grid investment rule and output, both already in the system.
// This creates fiscal amplification: intermittency -> higher I_grid(+1)
// -> higher tau(+1) -> lower after-tax MPK -> reduced investment today.
1/C = beta * ((1 - tau(+1)) * alpha * Y(+1) / K_p + 1 - delta_p) / C(+1);

// [2] Bond Euler equation (eq:bond_euler in paper)
// Household FOC for international bonds: equates domestic and foreign
// returns. Note r_star(+1): the return on bonds purchased today is
// realized next period, making r_star a jumper variable.
1/C = beta * (1 + r_star) / C(+1);

// [3] Labor supply / labor market clearing (eq:labor_supply in paper) -- PROPORTIONAL TAX VARIANT
// Proportional tax drives a wedge: household equates disutility of labor
// to the AFTER-TAX marginal product. The (1-tau) term reduces effective
// labor income, lowering equilibrium labor supply and amplifying
// the contractionary effect of intermittency shocks.
(1 - tau) * (1-alpha) * V = varphi * C * L^(1+sigma_L);

// ------- PRODUCTION BLOCK -------

// [4] CES output aggregator (eq:final_output in paper)
// Combines value-added V and exogenous energy endowment E_bar with
// elasticity of substitution sigma_E. Energy and value-added are
// complements (sigma_E = 0.6 < 1).
Y = ((1-omega_E) * V^((sigma_E-1)/sigma_E) + omega_E * E_bar^((sigma_E-1)/sigma_E))^(sigma_E/(sigma_E-1));

// [5] Cobb-Douglas value added (eq:value_added in paper)
// Effective capital = u * K_p(-1): reliability enters multiplicatively,
// so intermittency acts like a negative TFP shock.
// Note: K_p(-1) because capital installed last period is used today.
V = (u * K_p(-1))^alpha * L^(1-alpha);

// ------- RELIABILITY BLOCK -------

// [6] Exponential reliability function (eq:utilization in paper)
// Maps the ratio of flexibility assets to renewable volatility into
// grid utilization. The exponential form ensures u in (0,1) and
// captures the convexity: flexibility is most valuable when scarce.
// exp(a_ren) scales volatility up/down with the intermittency AR(1) process.
// NOTE: a_ren (not eps_ren) is used here so the AR(1) persistence rho_ren=0.85
// is correctly inherited. eps_ren is the i.i.d. innovation; a_ren carries memory.
u = 1 - exp(-psi * F / (phi_int * Vol_ren_bar * exp(a_ren)));

// [7] CES flexibility aggregator (eq:flexibility_bundle in paper)
// Combines battery storage (private, agile) and grid capital (public,
// inertial). Battery services = A_bat(-1) * K_b(-1) (technology-augmented).
// rho_flex < 1 implies batteries and grid are complements
// (Proposition 1 / superadditivity).
// NOTE: A_bat(-1) is used here (lagged technology), not A_bat (current).
// This eliminates the within-period simultaneity loop A_bat -> F -> u -> A_bat.
// Under this convention, scarcity-induced innovation in period t is reflected
// in A_bat,t (computed in eq [10] from current u_t), but only becomes
// productive in period t+1 when it enters F as A_bat(-1). This is a
// one-period-delay interpretation: innovation is induced today, operational tomorrow.
F = (mu * (A_bat(-1) * K_b(-1))^((rho_flex-1)/rho_flex) + (1-mu) * K_g(-1)^((rho_flex-1)/rho_flex))^(rho_flex/(rho_flex-1));

// ------- INVESTMENT BLOCK -------

// [8] Battery investment rule (eq:battery_investment in paper)
// Private battery investment responds to CURRENT-PERIOD reliability u
// and current battery price P_bat. This is internally consistent:
// I_bat only affects K_b next period (via eq [9]), so there is no
// within-period feedback from I_bat back into u. u is already
// determined before this equation runs, since F (and hence u) depends
// only on A_bat(-1) and K_b(-1) --- both predetermined.
I_bat / steady_state(I_bat) = (u_target / u)^phi_grid / P_bat;

// [9] Battery capital accumulation (eq:battery_capital in paper)
// Standard perpetual inventory method.
K_b = (1 - delta_b) * K_b(-1) + I_bat;

// [10] Scarcity-induced technological change in battery storage (eq:learning_by_doing in paper)
// Technology improves when reliability falls below target (u < u_target),
// creating a gap that drives innovation. chi controls how much of the
// scarcity signal passes through regulation (chi=0: no innovation).
// Timing: A_bat,t responds to current u_t, then becomes A_bat(-1) next period.
// Since F uses A_bat(-1) (not A_bat), the simultaneity loop is broken:
//   u_t is determined first from F_t = CES(A_bat(-1)*K_b(-1), K_g(-1)),
//   then A_bat,t is updated from the reliability gap (u_target - u_t).
//   A_bat,t only enters production in t+1 via F_{t+1}.
A_bat = A_bat(-1) * (1 + eta_bat * chi * (u_target - u) / u_target);

// [11] Grid investment rule (eq:grid_investment_rule in paper)
// Government follows a mechanical reliability-targeting rule (not Ramsey
// optimal). Investment accelerates when reliability falls below target.
// phi_grid governs aggressiveness; a_I captures implementation shocks
// (e.g., multi-year project delays, budget constraints) with AR(1) persistence rho_I.
I_grid / steady_state(I_grid) = (u_target / u)^phi_grid * exp(a_I);

// [12] Grid capital accumulation (eq:grid_capital in paper)
K_g = (1 - delta_g) * K_g(-1) + I_grid;

// ------- CAPITAL ACCUMULATION -------

// [13] Productive capital accumulation (eq:productive_capital in paper)
K_p = (1 - delta_p) * K_p(-1) + I_p;

// ------- EXTERNAL SECTOR (Schmitt-Grohe & Uribe 2003) -------

// [14] Current account / resource constraint (eq:current_account in paper)
// Open-economy budget constraint. B_star > 0 means net creditor.
// Trade balance = Y - C - I_p - P_bat*I_bat - I_grid.
// Net interest income = r_star(-1) * B_star(-1).
B_star = (1 + r_star(-1)) * B_star(-1) + Y - C - I_p - P_bat * I_bat - I_grid;

// [15] Debt-elastic interest rate premium (eq:interest_premium in paper)
// Ensures stationarity of net foreign assets. When the economy borrows
// (B_star < 0), r_star rises above r_bar, stabilizing the position.
// Without this, B_star would have a unit root (well-known SGU 2003 result).
r_star = r_bar + phi_b * (exp(B_star_ss - B_star) - 1);

// ------- EXOGENOUS PROCESSES -------

// [16] Battery price AR(1) process (eq:battery_price_shock in paper)
// World battery price follows a persistent stochastic process.
// P_bat = 1 in steady state (normalization via log specification).
log(P_bat) = rho_bat * log(P_bat(-1)) + eps_bat;

// [17] Proportional tax rate (eq:gov_budget in paper) -- NEW IN PROPTAX VARIANT
// Government balanced budget: tau_t * Y_t = I_grid_t.
// Tax rate adjusts endogenously each period to finance grid investment.
// In steady state: tau_ss = I_grid_ss / Y_ss = delta_g * K_g_ss / Y_ss (~1.5%).
// Dynamic movements in tau create the fiscal amplification channel.
tau = I_grid / Y;

// [18] Renewable intermittency AR(1) process
// a_ren is the persistent component of the intermittency shock.
// eps_ren is the i.i.d. innovation (std dev = sigma_ren = 0.12).
// In steady state: a_ren_ss = 0 (shock mean is zero).
// This is the standard way to implement AR(1) shocks in Dynare.
a_ren = rho_ren * a_ren(-1) + eps_ren;

// [19] Grid implementation shock AR(1) process
// a_I captures persistent delays in grid infrastructure deployment
// (multi-year permitting, construction schedules).
// eps_I is the i.i.d. innovation (std dev = 0.05).
// rho_I = 0.70: implementation shocks are moderately persistent.
a_I = rho_I * a_I(-1) + eps_I;

end;

// =========================================================================
// INITIAL VALUES AND STEADY STATE COMPUTATION
// =========================================================================
// These are initial guesses for the nonlinear solver. The actual steady
// state is computed by vietnam_dsge_steadystate.m, which:
//   1. Sets u = 0.97 and L = 0.33 as targets
//   2. Calibrates phi_int from the reliability equation (~55.6)
//   3. Calibrates varphi from the labor FOC (~9.7)
//   4. Iterates K_p <-> Y to convergence via the Euler equation
//   5. Derives C residually from the resource constraint
// =========================================================================

initval;
    Y = 0.95; C = 0.64; L = 0.33; K_p = 11.7;
    K_b = 0.19; K_g = 1.14; I_p = 0.29; I_bat = 0.0057;
    I_grid = 0.01425; u = 0.97; F = 0.526; A_bat = 1.0;
    V = 1.14; P_bat = 1.0;
    B_star = 0; r_star = 0.01;
    tau = 0.01425 / 0.95;  // tau_ss = I_grid_ss / Y_ss ~ 0.015
    a_ren = 0;              // intermittency process = 0 in steady state
    a_I = 0;                // implementation shock = 0 in steady state
end;

// Compute the deterministic steady state (calls vietnam_dsge_steadystate.m)
steady;

// Verify Blanchard-Kahn conditions:
// Expects 3 eigenvalues > 1 for 3 forward-looking variables (C, Y, r_star).
// C and r_star carry (+1) leads in the Euler equations; Y appears as Y(+1)
// in the capital Euler [1] and is thus classified as a jumper by Dynare.
check;

// =========================================================================
// SHOCK SPECIFICATION
// =========================================================================
// Standard deviations of the three exogenous innovations.
// Dynare convention: "stderr X" means the std. dev. of the shock is X.

shocks;
    var eps_ren; stderr sigma_ren;      // 0.12 (Vietnamese solar CV)
    var eps_bat; stderr sigma_bat;      // 0.08 (global lithium price vol)
    var eps_I;   stderr 0.05;           // implementation shock
end;

// =========================================================================
// SOLUTION AND SIMULATION
// =========================================================================
// order=1:    first-order perturbation (log-linearization around SS)
// irf=40:     impulse responses for 40 quarters (10 years)
// periods=1000: simulate 1000 periods for moments
// nograph:    suppress Dynare's default plots (we make custom ones in
//             run_dynare.m)

stoch_simul(order=1, irf=40, periods=1000, nograph)
    Y C L I_p I_bat I_grid u F A_bat V P_bat B_star r_star tau a_ren a_I;
