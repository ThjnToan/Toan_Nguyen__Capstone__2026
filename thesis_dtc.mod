// ==================================================================
// FILE: thesis_dtc.mod (FINAL OPEN ECONOMY VERSION)
// SCENARIO: Small Open Economy with Imported Battery Tech
// DESCRIPTION: Core DTC-DSGE model capturing the "Agility Gap"
// ==================================================================

// 1. VARIABLES
var 
    Y C I K L w r_k u F 
    A_bat MPK_bat               
    K_g K_b I_g I_b             
    MPK_b_physical              
    tau_res G_cost          
    Welfare
    // --- OPEN ECONOMY VARIABLES ---
    Q           // Real Exchange Rate
    X           // Exports
    B_star      // Net Foreign Assets
    r_star_eff  // Effective Foreign Interest Rate
    P_tech      // Imported Battery Price (Global Shock)
    TradeBal;   // Trade Balance (X - Value of Imports)

varexo 
    e_ren       // Renewable Supply Shock (Domestic Intermittency)
    e_price;    // Global Battery Price Shock (Imported Greenflation)

// 2. PARAMETERS
parameters 
    beta delta alpha psi phi_int mu rho sigma_L A 
    eta_dtc delta_g delta_b phi_g
    phi_sub cost_scale chi_price
    // --- OPEN ECONOMY PARAMETERS ---
    theta_x     // Export elasticity
    phi_b       // Debt adjustment cost (risk premium)
    r_star      // World interest rate
    phi_adj_bat; // Battery Installation Adjustment Cost

// 3. CALIBRATION 
// --- Macro Fundamentals ---
beta    = 0.99;    
delta   = 0.025;
alpha   = 0.35;    
sigma_L = 1.0;     
A       = 1.0;

// --- Energy & Reliability ---
psi     = 2.0;     // Sensitivity of reliability to flexibility gap
phi_int = 0.3;     // Renewable Intermittency Coefficient
mu      = 0.16;    // Battery Share (Calibrated to PDP8)
rho     = 0.4;     // Complementarity (Low substitution)

// --- Innovation & Infrastructure ---
eta_dtc    = 0.5;     // Elasticity of innovation to price signals
delta_g    = 0.015;   // Depreciation: Grid (Slower)
delta_b    = 0.025;   // Depreciation: Batteries (Faster)
phi_g      = 2.0;     // Investment Inertia (Public Grid)
phi_sub    = 0.0;     // R&D Subsidy (Baseline = 0)
cost_scale = 0.05; 
chi_price  = 1.0;     // Market Friction (1.0 = Perfect Price Signal)

// --- Open Economy & Frictions ---
theta_x = 1.5;            // Trade elasticity
phi_b   = 0.001;          // Portfolio adjustment cost (stationarity)
r_star  = (1/beta) - 1;   // World interest rate matches domestic SS

// [CRITICAL FRICTION] 
// Absorptive Capacity Constraint: Slows down private battery installation
// to prevent unrealistic instantaneous stabilization via trade.
phi_adj_bat = 100.0;

// 4. THE MODEL
model;
    // --- HOUSEHOLDS & FIRMS ---
    // 1. Domestic Euler Equation
    1/C = beta * (1/C(+1)) * (r_k(+1) + 1 - delta);
    
    // 2. Foreign Euler (UIP Condition for Exchange Rate Q)
    1 = beta * (C / C(+1)) * (1 + r_star_eff) * (Q(+1) / Q);
    
    // 3. Effective Foreign Rate (Debt Elastic) 
    r_star_eff = r_star - phi_b * B_star;

    // 4. Labor Supply
    w = C * L^sigma_L;

    // 5. Production Function (With Reliability Wedge u)
    Y = A * ( (u * K(-1))^alpha ) * (L^(1-alpha));
    w = (1-alpha) * Y / L;
    r_k = alpha * Y / K(-1);

    // --- RELIABILITY BLOCK ---
    // Utilization depends on the gap between Flexibility (F) and Intermittency
    u = 1 - exp( -psi * (F / (phi_int * (1+e_ren))) );
    
    // Flexibility Aggregation (CES: Grid vs. Batteries)
    F = ( mu * (A_bat * K_b(-1))^((rho-1)/rho) + (1-mu) * K_g(-1)^((rho-1)/rho) )^(rho/(rho-1));

    // --- INVESTMENT & ACCUMULATION ---
    // Public Grid: Subject to Investment Inertia (Time-to-Build proxy)
    I_g = delta_g * 15 * (0.99 / u)^phi_g;
    K_g = (1 - delta_g) * K_g(-1) + I_g(-4); 

    // Battery Marginal Product (Physical)
    MPK_b_physical = (Y / F) * (1 - u) * (mu * (F / (A_bat*K_b(-1)))^(1/rho)) * A_bat;
    
    // 6. Battery Arbitrage (Open Economy with Adjustment Costs)
    // Private sector imports batteries (I_b) at global price (P_tech * Q).
    // The adjustment cost term [phi_adj_bat * (...)] creates the "Agility Gap".
    P_tech * Q * (1 + phi_adj_bat * (I_b/K_b(-1) - delta_b)) 
        = beta * (C / C(+1)) * (MPK_b_physical(+1) + P_tech(+1)*Q(+1)*(1 - delta_b));
    
    K_b = (1 - delta_b) * K_b(-1) + I_b;

    // --- INNOVATION (Directed Technical Change) ---
    tau_res = phi_sub * (0.99 - u);
    G_cost  = cost_scale * tau_res * Y;
    MPK_bat = MPK_b_physical;
    
    // Law of Motion for Battery Technology
    log(A_bat) = 0.95 * log(A_bat(-1)) + eta_dtc * (1 + tau_res) * chi_price * (log(MPK_bat) - log(MPK_bat(-1)));

    // --- EXTERNAL SECTOR ---
    // 7. Global Price Shock Process
    log(P_tech) = 0.9 * log(P_tech(-1)) + e_price; 

    // 8. Export Demand
    X = 0.2 * Y * Q^theta_x; 

    // 9. Balance of Payments
    // Imports = Consumption Imports (implied) + Battery Tech Imports (P_tech * I_b)
    Q * B_star = Q * B_star(-1) * (1 + r_star_eff(-1)) + X - Q * (P_tech * I_b);
    
    // 10. Trade Balance Definition
    TradeBal = X - Q * (P_tech * I_b);

    // --- CLOSING ---
    // 11. Market Clearing
    // Domestic Output Y covers C, Domestic I, Public I_g, R&D costs, and Exports.
    // Note: I_b is excluded here because it is imported.
    Y = C + I + I_g + G_cost + X;
    
    K = (1-delta)*K(-1) + I;
    Welfare = log(C) - (L^(1+sigma_L))/(1+sigma_L) + beta * Welfare(+1);
end;

// 5. STEADY STATE
initval;
    e_ren = 0;
    e_price = 0;
    P_tech = 1;
    Q = 1;
    B_star = 0;
    r_star_eff = (1/beta) - 1;
    
    A_bat = 1;
    L     = 0.93;
    u     = 0.99;    
    
    MPK_b_physical = (1/beta) - (1-delta_b);
    MPK_bat = MPK_b_physical;
    r_k = (1/beta) - (1-delta);
    
    K_g = 15;
    I_g = delta_g * K_g; 
    K_b = 2.0;
    I_b = delta_b * K_b;
    F   = 6.0; 
    K   = 30;
    I   = delta * K;
    tau_res = 0;
    G_cost = 0;
    Y = 3.2;
    
    // SOE Closing: Exports cover Battery Imports in SS
    X = I_b; 
    
    C = Y - I - I_g - G_cost - X;
    w = (1-alpha)*Y/L;
    Welfare = (log(C) - (L^(1+sigma_L))/(1+sigma_L)) / (1-beta);
end;
steady; 

// 6. SIMULATION
shocks;
    // STRESS TEST: 50% Negative Shock to Renewable Potential
    var e_ren; stderr 0.50; 
    // Global Price Shock (Inactive for baseline, activated in SOE script)
    var e_price; stderr 0.10; 
end;

// Run Impulse Response Functions (IRFs)
stoch_simul(order=1, irf=40, nograph);

// 7. PLOTTING
verbatim;
    T = 1:40;
    
    % --- FIGURE 1: THE RELIABILITY PENALTY ---
    figure('Name', 'Fig1_Reliability');
    subplot(1,3,1); 
    plot(T, oo_.irfs.Y_e_ren, 'k-', 'LineWidth', 2); 
    title('Output (Y)'); grid on; axis tight;
    
    subplot(1,3,2); 
    plot(T, oo_.irfs.u_e_ren, 'r-', 'LineWidth', 2); 
    title('Reliability (u)'); grid on; axis tight;
    
    subplot(1,3,3); 
    plot(T, oo_.irfs.F_e_ren, 'b-', 'LineWidth', 2); 
    title('Flexibility Stock (F)'); grid on; axis tight;

    % --- FIGURE 2: MACRO TRANSMISSION ---
    figure('Name', 'Fig2_Macro');
    subplot(2,2,1); 
    plot(T, oo_.irfs.I_e_ren, 'b-', 'LineWidth', 2); 
    title('Investment (I)'); grid on; axis tight;
    
    subplot(2,2,2); 
    plot(T, oo_.irfs.L_e_ren, 'k-', 'LineWidth', 2); 
    title('Labor (L)'); grid on; axis tight;
    
    subplot(2,2,3); 
    plot(T, oo_.irfs.C_e_ren, 'k--', 'LineWidth', 2); 
    title('Consumption (C)'); grid on; axis tight;
    
    subplot(2,2,4); 
    plot(T, oo_.irfs.w_e_ren, 'r-', 'LineWidth', 2); 
    title('Real Wages (w)'); grid on; axis tight;

    % --- FIGURE 3: INNOVATION DYNAMICS ---
    figure('Name', 'Fig3_Innovation');
    subplot(1,3,1); 
    plot(T, oo_.irfs.MPK_bat_e_ren, 'r-', 'LineWidth', 2); 
    title('Shadow Price (MPK)'); grid on; axis tight;
    
    subplot(1,3,2); 
    plot(T, oo_.irfs.A_bat_e_ren, 'g-', 'LineWidth', 2); 
    title('Battery Tech (A_{bat})'); grid on; axis tight;
    
    subplot(1,3,3); 
    plot(T, oo_.irfs.I_b_e_ren, 'b-', 'LineWidth', 2); 
    title('Battery Inv (I_b)'); grid on; axis tight;
end;
