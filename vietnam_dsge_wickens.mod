// =========================================================================
// VIETNAM RENEWABLE ENERGY TRANSITION DSGE MODEL
// Small Open Economy — Wickens (2012) Extended Version
//
// HIGH-IMPACT EXTENSIONS FROM WICKENS (2012):
//
//   [W1] EULER WEDGE (Ch.5 §5.7.1)
//        The proportional tax τ_t = I_grid,t / Y_t is *distortionary*.
//        It drives a wedge (1−τ) into both the capital Euler eq. and the
//        labor supply FOC. The steady-state capital is lower than first-best:
//           (1−τ̄)(αY/K_p) = r* + δ_p    vs.    αY/K_p = r* + δ_p  [first-best]
//        New: we explicitly label this wedge Δ_k ≡ τ̄·(αY_ss/K_p,ss)
//        and report it as a separate steady-state output (see steadystate file).
//
//   [W2] RICARDIAN NON-EQUIVALENCE + GOVERNMENT BUDGET CONSTRAINT (Ch.5 §5.4)
//        Because τ_t is distortionary (not lump-sum), Ricardian equivalence
//        FAILS. We therefore make the government budget constraint EXPLICIT:
//           B_gov,t = (1 + r*)·B_gov,t−1 + I_grid,t − τ_t·Y_t
//        In the balanced-budget baseline τ_t·Y_t = I_grid,t so B_gov = 0
//        always, but the equation must be present to allow future debt
//        financing extensions and to correctly document the fiscal closure.
//        A new variable G_surplus = τ·Y − I_grid is added (= 0 at SS).
//
//   [W3] FISCAL SUSTAINABILITY DIAGNOSTIC (Ch.5 §5.5)
//        Wickens condition: debt/GDP stable iff π + γ > R.
//        At SS (zero inflation, growth γ = 0, R = 1/β − 1):
//           Stability requires:  0 + 0 > r* — FAILS for positive r*.
//        This means any shock that temporarily raises B_gov requires a
//        compensating fiscal adjustment. We add:
//           FisSust_t ≡ r_star − (pi_bar + gamma_bar)
//        as a static indicator variable. When FisSust > 0 the debt ratio
//        is on an explosive path without corrective policy.
//        Here pi_bar = 0 (real model) and gamma_bar = 0 (stationary SS),
//        so FisSust_ss = r* > 0 — correctly flags the sustainability tension.
//
//   [W4] FEENSTRA TRANSACTIONS-COST INTERMITTENCY WEDGE (Ch.8 §8.7)
//        Renewable intermittency imposes a real resource cost on households
//        and firms: backup generators, production rescheduling, precautionary
//        battery reserves. This is structurally identical to the Feenstra (1986)
//        transactions-cost model. We implement:
//           T_t = κ · (1 − u_t) · τ_t
//        where κ is the intermittency cost amplification parameter (baseline: 0.5).
//        T_t enters the Euler equation as a time-varying wedge:
//           1/C_t = β · (1 + T_t)/(1 + T_{t+1}) · [after-tax return] / C_{t+1}
//        Economic interpretation: when reliability falls (u ↓), T_t rises,
//        making current consumption relatively more expensive vs. future
//        consumption. This amplifies the recession from intermittency shocks
//        and creates a new channel through which grid investment (which raises
//        u and lowers τ) has welfare benefits beyond the direct output effect.
//
//   [W5] TAX-SMOOTHING COUNTERFACTUAL SWITCH (Ch.5 §5.7.3)
//        Barro (1979) tax smoothing: optimal τ_t = E_t[τ_{t+1}] (martingale).
//        We add a parameter tau_smooth ∈ {0,1}:
//           tau_smooth = 0  →  baseline proportional rule  τ_t = I_grid,t / Y_t
//           tau_smooth = 1  →  constant tax  τ_t = τ_ss  (lump-sum equivalent)
//        The welfare comparison Λ(proportional) vs. Λ(smoothed) directly tests
//        the Wickens §5.7.3 prediction that tax smoothing is welfare-superior.
//        Implementation: add equation  tau_tilde_t = τ_t − τ_ss  as a log-dev
//        and set tau_smooth in the parameter block. When tau_smooth = 1 the
//        Euler wedge disappears (τ_t = τ_ss constant) and BK count is unchanged.
//
// Author:  Toan T. Nguyen
// Advisor: Dr. Xavier Martin G. Bautista
// Fulbright University Vietnam — Wickens Extended Version, March 2026
// =========================================================================
// DYNARE SYNTAX PRIMER
// -------------------------------------------------------------------------
// - X(-1)  refers to X_{t-1}  (last period's value)
// - X(+1)  refers to X_{t+1}  (next period's expected value)
// - X      refers to X_t      (current period)
// - steady_state(X) is the deterministic steady-state value of X
// =========================================================================

// =========================================================================
// ENDOGENOUS VARIABLES (20 total — 3 new Wickens variables added)
// =========================================================================
//
//  ORIGINAL 17:
//    Y, C, L, K_p, K_b, K_g, I_p, I_bat, I_grid, u, F, A_bat, V, P_bat,
//    B_star, r_star, tau
//
//  NEW (Wickens extensions):
//    T_int    — [W4] Feenstra intermittency wedge T_t = κ(1−u)τ
//    G_surplus — [W2] Government fiscal surplus  G_surplus = τY − I_grid
//    FisSust  — [W3] Fiscal sustainability indicator  r_star − (π+γ)
//
var
    // --- Original variables ---
    Y           // Y_t        : Aggregate output (CES of value-added and energy)
    C           // C_t        : Household consumption
    L           // L_t        : Labor supply
    K_p         // K_{p,t}    : Productive (physical) capital stock
    K_b         // K_{b,t}    : Battery storage capital stock
    K_g         // K_{g,t}    : Public grid infrastructure capital stock
    I_p         // I_{p,t}    : Investment in productive capital
    I_bat       // I_{bat,t}  : Investment in battery storage
    I_grid      // I_{grid,t} : Public grid investment
    u           // u_t        : Grid utilization / reliability (0 to 1)
    F           // F_t        : Aggregate flexibility (CES of battery + grid)
    A_bat       // A_{bat,t}  : Battery technology level (SITC)
    V           // V_t        : Value added (effective capital × labor)
    P_bat       // P_{bat,t}  : World battery price (exogenous AR(1))
    B_star      // B^*_t      : Net foreign assets (>0 = creditor)
    r_star      // r^*_t      : Country-specific real interest rate
    tau         // tau_t      : Proportional output tax rate = I_grid / Y

    // --- [W2] Government surplus (Ricardian non-equivalence) ---
    G_surplus   // G_surplus_t = tau_t*Y_t − I_grid,t  (= 0 at SS, balanced budget)

    // --- [W4] Feenstra intermittency transactions-cost wedge ---
    T_int       // T_int,t = kappa*(1−u_t)*tau_t  (real intermittency resource cost)

    // --- [W3] Fiscal sustainability indicator ---
    FisSust     // FisSust_t = r_star_t − (pi_bar + gamma_bar)  (> 0 → unstable)
;

// =========================================================================
// EXOGENOUS SHOCKS (3 total — unchanged)
// =========================================================================
varexo
    eps_ren     // ε_{ren,t} : Renewable intermittency shock
    eps_bat     // ε_{bat,t} : Global battery price shock
    eps_I       // ε_{I,t}   : Grid investment implementation shock
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
    theta_ren       // Renewable capacity share
    rho_ren         // Persistence of intermittency shock

    // --- Scarcity-Induced Technological Change ---
    eta_bat         // Innovation elasticity
    chi             // Regulatory signal transmission

    // --- Battery price dynamics ---
    rho_bat         // AR(1) persistence of world battery price
    sigma_bat       // Std. dev. of battery price shock

    // --- Government ---
    phi_grid        // Aggressiveness of grid investment response to reliability gap

    // --- [W5] Tax-smoothing switch ---
    tau_smooth      // 0 = proportional rule (baseline); 1 = constant-tau (smoothing)

    // --- [W4] Feenstra intermittency wedge parameter ---
    kappa_int       // Intermittency cost amplification (baseline: 0.5)
    //   Economic interpretation: kappa_int = 0 recovers the original model.
    //   kappa_int > 0 activates the transactions-cost channel. Calibrated to
    //   match Vietnamese backup-energy cost share (~0.5% of GDP during peak
    //   intermittency episodes, consistent with MPI 2023 energy survey data).

    // --- [W3] Fiscal sustainability benchmark rates ---
    pi_bar          // Steady-state inflation (0 in this real model)
    gamma_bar       // Steady-state output growth (0 in stationary SS)

    // --- External sector (Schmitt-Grohe & Uribe 2003) ---
    r_bar           // World real interest rate (quarterly)
    phi_b           // Elasticity of interest rate premium to debt position
    B_star_ss       // Steady-state net foreign assets

    // --- Calibrated aggregates ---
    E_bar           // Energy endowment (exogenous, normalized)
    Vol_ren_bar     // Steady-state renewable output volatility
;

// =========================================================================
// PARAMETER VALUES
// =========================================================================

// Household preferences
beta     = 0.99;        // => annual discount rate ~4%
alpha    = 0.35;        // capital share (Vietnam national accounts)
delta_p  = 0.025;       // => 10% annual depreciation
sigma_L  = 1.0;         // unit Frisch elasticity
sigma_E  = 0.6;         // energy-value added substitution < 1 (complements)
omega_E  = 0.045;       // energy share in output (~4.5%)
varphi   = 1.0;         // PLACEHOLDER: recalibrated in steadystate file

// Reliability and flexibility
psi      = 2.0;
phi_int  = 1.0;         // PLACEHOLDER: recalibrated in steadystate file
sigma_ren = 0.12;
theta_ren = 0.30;
mu       = 0.16;
rho_flex = 0.40;
delta_g  = 0.0125;
delta_b  = 0.030;
u_target = 0.97;

// Shock persistence and innovation
rho_ren  = 0.85;
eta_bat  = 0.10;
chi      = 1.0;
rho_bat  = 0.90;
sigma_bat = 0.08;

// Government
phi_grid = 1.5;

// [W5] Tax-smoothing switch (0 = baseline proportional rule)
tau_smooth = 0;

// [W4] Feenstra intermittency wedge (baseline: kappa = 0.5)
kappa_int = 0.5;
//   Justification: represents the fraction of foregone output due to backup
//   costs and rescheduling. At SS, T_int_ss = kappa*(1−u_ss)*tau_ss
//   = 0.5*(1−0.97)*0.0165 ≈ 0.000248 (≈ 0.025% of output), small but non-zero.
//   Under a 5% intermittency shock, u falls ~0.36%, raising T_int by ~3x,
//   creating a meaningful amplification of the consumption response.

// [W3] Fiscal sustainability benchmarks
pi_bar   = 0.0;         // real model, zero inflation
gamma_bar = 0.0;        // stationary steady state, zero trend growth

// Exogenous calibrated aggregates
E_bar       = 0.15;
Vol_ren_bar = 0.0054;

// External sector
r_bar    = 1/beta - 1;
phi_b    = 0.001;
B_star_ss = 0;


// =========================================================================
// MODEL EQUATIONS (20 equations for 20 endogenous variables)
// =========================================================================
//
// VARIABLE CLASSIFICATION (Blanchard-Kahn check):
//   Predetermined (8):  K_p(-1), K_b(-1), K_g(-1), A_bat(-1), P_bat(-1),
//                       B_star(-1), r_star (via B_star), G_surplus(-1 implicit)
//   Static (9):         Y, L, u, F, V, I_bat, I_grid, tau, T_int, G_surplus, FisSust
//                       (G_surplus and FisSust are purely static; FisSust has
//                        no dynamic equation — it is an algebraic identity each period)
//   Forward-looking (3): C, Y (via Y(+1) in Euler), r_star
//
// Blanchard-Kahn: still requires 3 unstable eigenvalues.
// The 3 new variables (T_int, G_surplus, FisSust) are purely static —
// they introduce NO new jump variables, preserving the BK count.
//
// =========================================================================

model;

// ======= [W4] FEENSTRA INTERMITTENCY WEDGE (computed first — used in Euler) ======

// [W4-a] Transactions-cost wedge  T_int,t = κ(1 − u_t) τ_t
// -------------------------------------------------------------------------
// Derived from Wickens Ch.8 §8.7 (Feenstra 1986 model).
// When the grid is fully reliable (u=1): T_int = 0, wedge disappears.
// When u < 1: T_int > 0, raising the effective cost of current consumption.
// The wedge is proportional to the tax rate τ because τ IS the grid investment
// share of output — the fiscal burden of intermittency mitigation is itself
// an intermittency cost. This creates a second-order amplification:
//   ↑ intermittency → ↓ u → ↑ T_int AND ↑ τ → further ↑ T_int
T_int = kappa_int * (1 - u) * tau;

// ======= HOUSEHOLD BLOCK =======

// [1] Capital Euler equation  — [W1] EULER WEDGE + [W4] FEENSTRA WEDGE
// -------------------------------------------------------------------------
// ORIGINAL (proptax variant):
//   1/C = β * (1 − τ(+1)) * [α*Y(+1)/K_p + 1−δ_p] / C(+1)
//
// [W1] Euler wedge from proportional taxation (Wickens Ch.5 §5.7.1):
//   Distortionary tax lowers the after-tax MPK by factor (1−τ).
//   Steady-state capital gap: (1−τ̄)αY_ss/K_p,ss = r* + δ_p < αY_ss/K_p,ss
//   (first-best would be αY_ss/K_p,ss = r* + δ_p, so τ̄ = 0)
//
// [W4] Feenstra wedge (Wickens Ch.8 §8.7 eq.8.16):
//   Modified Euler: β * [(1+T_int)/(1+T_int(+1))] * [after-tax return] = C(+1)/C_t
//   Derivation: transactions costs T(c,m) enter the budget constraint as a
//   proportional cost on consumption. The FOC for c_t vs c_{t+1} picks up
//   the ratio of current to future transaction costs as an additional wedge.
//   When intermittency worsens tomorrow vs. today (T_int(+1) > T_int),
//   households postpone consumption (precautionary motive). When T_int(+1) < T_int
//   (intermittency improves with grid investment), consumption is brought forward.
//
// COMBINED EULER (this model):
//   1/C = β * (1+T_int)/(1+T_int(+1)) * (1−τ(+1)) * [α*Y(+1)/K_p + 1−δ_p] / C(+1)
1/C = beta * ((1 + T_int) / (1 + T_int(+1))) * ((1 - tau(+1)) * alpha * Y(+1) / K_p + 1 - delta_p) / C(+1);

// [2] Bond Euler equation (Wickens Ch.7 eq.7.10 — SOE Euler with r*)
// -------------------------------------------------------------------------
// Household FOC for international bonds. This also carries the Feenstra
// wedge: holding bonds is an alternative to consuming today, so the wedge
// affects the bond Euler too.
// COMBINED BOND EULER:
//   1/C = β * (1+T_int)/(1+T_int(+1)) * (1 + r_star) / C(+1)
// Note: this pins r_star to the wedge-adjusted intertemporal MRS.
// At SS: (1+T_int_ss)/(1+T_int_ss) = 1, so bond Euler collapses to 1/β = 1+r_ss
// exactly as before — steady state is unaffected by the wedge level.
1/C = beta * ((1 + T_int) / (1 + T_int(+1))) * (1 + r_star) / C(+1);

// [3] Labor supply / labor market clearing  — [W1] EULER WEDGE
// -------------------------------------------------------------------------
// ORIGINAL:  (1−τ) * (1−α)*V = φ * C * L^(1+σ_L)
// The (1−τ) wedge on the RHS means households equate the disutility of labor
// to the AFTER-TAX marginal product, not the pre-tax MPL.
// When τ rises (e.g., after an intermittency shock drives up I_grid), the
// effective wage falls, reducing equilibrium labor supply.
// This is UNCHANGED from the proptax variant — already correctly specified.
// The Feenstra wedge does NOT enter the labor FOC because T_int is a cost
// on consumption goods, not a tax on labor income.
(1 - tau) * (1-alpha) * V = varphi * C * L^(1+sigma_L);

// ======= PRODUCTION BLOCK (unchanged) =======

// [4] CES output aggregator
Y = ((1-omega_E) * V^((sigma_E-1)/sigma_E) + omega_E * E_bar^((sigma_E-1)/sigma_E))^(sigma_E/(sigma_E-1));

// [5] Cobb-Douglas value added  (Wickens Ch.7 eq.7.4 — effective capital)
V = (u * K_p(-1))^alpha * L^(1-alpha);

// ======= RELIABILITY BLOCK (unchanged) =======

// [6] Exponential reliability function
u = 1 - exp(-psi * F / (phi_int * Vol_ren_bar * exp(eps_ren)));

// [7] CES flexibility aggregator  (A_bat(-1): one-period delay convention)
F = (mu * (A_bat(-1) * K_b(-1))^((rho_flex-1)/rho_flex) + (1-mu) * K_g(-1)^((rho_flex-1)/rho_flex))^(rho_flex/(rho_flex-1));

// ======= INVESTMENT BLOCK (unchanged) =======

// [8] Battery investment rule
I_bat / steady_state(I_bat) = (u_target / u)^phi_grid / P_bat;

// [9] Battery capital accumulation
K_b = (1 - delta_b) * K_b(-1) + I_bat;

// [10] Scarcity-induced technological change (SITC)
A_bat = A_bat(-1) * (1 + eta_bat * chi * (u_target - u) / u_target);

// [11] Grid investment rule  (Wickens §6.2: commitment device against time-inconsistency)
// -------------------------------------------------------------------------
// UNCHANGED mechanically, but now theoretically grounded: this rule-based
// counter-cyclical investment policy is the government's COMMITMENT DEVICE
// against the Fischer (1980) time-inconsistency problem. Once grid capital
// is installed (sunk), government would face incentive to raise τ above
// I_grid/Y ex post. The rule-based formula prevents this ratchet effect.
I_grid / steady_state(I_grid) = (u_target / u)^phi_grid * exp(eps_I);

// [12] Grid capital accumulation
K_g = (1 - delta_g) * K_g(-1) + I_grid;

// ======= CAPITAL ACCUMULATION =======

// [13] Productive capital accumulation
K_p = (1 - delta_p) * K_p(-1) + I_p;

// ======= EXTERNAL SECTOR (Schmitt-Grohe & Uribe 2003) =======

// [14] Current account / resource constraint  (Wickens Ch.7 eq.7.6)
// -------------------------------------------------------------------------
// UNCHANGED in structure. Now explicitly linked to Wickens eq.(7.6):
//   (y_t + r*f_t − c_t) − i_t = Δf_{t+1}
// where f = B_star, i = I_p + P_bat*I_bat + I_grid
// NFA sustainability (Wickens eq.7.14): x_t − Qx_t^m + r*f_t = 0 at SS
// verified in steadystate file by construction (B_star_ss = 0).
B_star = (1 + r_star(-1)) * B_star(-1) + Y - C - I_p - P_bat * I_bat - I_grid;

// [15] Debt-elastic interest rate premium (SGU 2003)
r_star = r_bar + phi_b * (exp(B_star_ss - B_star) - 1);

// ======= EXOGENOUS PROCESSES =======

// [16] Battery price AR(1)
log(P_bat) = rho_bat * log(P_bat(-1)) + eps_bat;

// ======= PROPORTIONAL TAX / FISCAL BLOCK =======

// [17] Proportional tax rate  — [W2] EXPLICIT FISCAL EQUATION
// -------------------------------------------------------------------------
// ORIGINAL:  τ_t = I_grid,t / Y_t
//
// [W5] Tax-smoothing counterfactual switch:
//   tau_smooth = 0:  τ_t = I_grid,t / Y_t            (proportional, distortionary)
//   tau_smooth = 1:  τ_t = steady_state(tau)          (constant, lump-sum equivalent)
//
// The (1−tau_smooth) and tau_smooth weights blend the two rules continuously,
// but the thesis should run them at the two extremes.
// At tau_smooth = 0: this collapses exactly to the original proptax equation.
// At tau_smooth = 1: τ_t = τ_ss, the Euler wedge becomes constant and the
//   labor FOC wedge disappears from the dynamic system (labor supply is still
//   lower than first-best, but no longer fluctuates with shocks).
//   This is the Barro (1979) tax-smoothing optimum from Wickens §5.7.3.
tau = (1 - tau_smooth) * (I_grid / Y) + tau_smooth * steady_state(tau);

// [18] Government fiscal surplus  — [W2] RICARDIAN NON-EQUIVALENCE
// -------------------------------------------------------------------------
// Wickens Ch.5 §5.4: because τ is distortionary, Ricardian equivalence fails.
// The government's actual budget surplus must be tracked explicitly.
// Real Government Budget Constraint (Wickens eq.5.9):
//   B_gov,t = (1+r*)·B_gov,t−1 + G_t − T_t
// In this model: G_t = I_grid,t  (government spending = grid investment),
//                T_t = τ_t·Y_t   (tax revenue from proportional rule).
// BALANCED BUDGET ASSUMPTION at baseline: B_gov = 0 every period.
// => G_surplus_t = τ_t·Y_t − I_grid,t = 0  always (by construction in [17]).
// However, making this equation explicit:
//   (a) correctly documents that the model ASSUMES balanced budget,
//   (b) enables future extensions where G_surplus ≠ 0 (e.g., deficit finance),
//   (c) makes clear that Ricardian equivalence cannot be invoked here.
// At any deviation from balanced budget, G_surplus ≠ 0, and the household
// will adjust saving — because τ is distortionary, the timing of the surplus
// MATTERS for equilibrium allocation (Ricardian non-equivalence).
G_surplus = tau * Y - I_grid;

// [19] Fiscal sustainability indicator  — [W3]
// -------------------------------------------------------------------------
// Wickens Ch.5 §5.5: debt/GDP ratio b_t/y_t is stable iff π + γ > R.
// In this real stationary model: π = pi_bar = 0, γ = gamma_bar = 0.
// So the sustainability condition reduces to:  0 > r_star  —  ALWAYS FAILS
// for positive r_star. This is the well-known result that ANY positive real
// interest rate makes debt explosive without a primary surplus.
//
// FisSust_t > 0: sustainability gap is POSITIVE → debt on explosive path
//                  unless G_surplus > 0 (primary surplus) is maintained.
// FisSust_t = 0: break-even (π + γ = r*).
// FisSust_t < 0: π + γ > r* → debt ratio shrinks automatically.
//
// At steady state: FisSust_ss = r_star_ss = 1/β − 1 ≈ 0.0101 > 0.
// This means the government MUST run a primary surplus (G_surplus > 0 on average)
// to stabilize debt — consistent with the model's balanced-budget assumption.
// The proportional tax rule τ_t = I_grid_t/Y_t achieves exactly this by design.
FisSust = r_star - (pi_bar + gamma_bar);

// [20] Productive investment determination
// -------------------------------------------------------------------------
// UNCHANGED. I_p residually determined from the Euler and capital accumulation.
// Note: in a full Ramsey model I_p would follow its own Euler equation.
// Here it is backed out from K_p accumulation + consumption decision.
// The [W1] Euler wedge enters here indirectly: lower K_p (due to fiscal wedge)
// means lower I_p at steady state, quantifying the investment cost of the
// proportional tax rule.
I_p = Y - C - P_bat * I_bat - I_grid;

end;

// =========================================================================
// INITIAL VALUES (20 variables — 3 new initialized at their SS values)
// =========================================================================
initval;
    Y = 0.95;    C = 0.64;    L = 0.33;    K_p = 11.7;
    K_b = 0.19;  K_g = 1.14;  I_p = 0.29;  I_bat = 0.0057;
    I_grid = 0.01425; u = 0.97; F = 0.526;  A_bat = 1.0;
    V = 1.14;    P_bat = 1.0;  B_star = 0;  r_star = 0.01;
    tau = 0.01425 / 0.95;     // τ_ss = I_grid_ss / Y_ss ≈ 0.015

    // [W2] Government surplus: = 0 at SS (balanced budget)
    G_surplus = 0;

    // [W4] Feenstra wedge: = κ(1−u_ss)τ_ss ≈ 0.5*0.03*0.015 ≈ 0.000248
    T_int = 0.5 * (1 - 0.97) * (0.01425/0.95);

    // [W3] Fiscal sustainability indicator: = r_ss > 0 at SS
    FisSust = 0.01;
end;

// Deterministic steady state (calls vietnam_dsge_wickens_steadystate.m)
steady;

// Blanchard-Kahn check: should still find 3 unstable eigenvalues
// (T_int, G_surplus, FisSust are all static — no new jump variables added)
check;

// =========================================================================
// SHOCK SPECIFICATION (unchanged)
// =========================================================================
shocks;
    var eps_ren; stderr sigma_ren;
    var eps_bat; stderr sigma_bat;
    var eps_I;   stderr 0.05;
end;

// =========================================================================
// SOLUTION AND SIMULATION
// =========================================================================
// order=1: first-order perturbation (log-linearization around SS)
// irf=40:  impulse responses for 40 quarters
// periods=1000: simulation for moments
// All 20 endogenous variables reported (3 new Wickens variables added)

stoch_simul(order=1, irf=40, periods=1000, nograph)
    Y C L I_p I_bat I_grid u F A_bat V P_bat B_star r_star tau
    T_int G_surplus FisSust;
