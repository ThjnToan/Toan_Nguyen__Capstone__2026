// ==================================================================
// FILE: thesis_pdp8.mod (OPEN ECONOMY + RAMP + ADJUSTMENT COSTS)
// SCENARIO: Long-Run Transition in a Small Open Economy
// ==================================================================

// 1. VARIABLES
var 
    Y C I K L w r_k u F 
    A_bat MPK_bat               
    K_g K_b I_g I_b             
    MPK_b_physical              
    tau_res G_cost          
    Welfare
    // --- OPEN ECONOMY VARS ---
    Q           // Real Exchange Rate
    X           // Exports
    B_star      // Net Foreign Assets
    r_star_eff  // Effective Foreign Interest Rate
    P_tech      // Imported Battery Price
    TradeBal;   // Trade Balance

varexo 
    E_ren;      // The Policy Target (Exogenous path)

// 2. PARAMETERS
parameters 
    beta delta alpha psi phi_int mu rho sigma_L A 
    eta_dtc delta_g delta_b phi_g
    phi_sub cost_scale chi_price
    // --- OPEN ECONOMY PARAMS ---
    theta_x     // Export elasticity
    phi_b       // Debt adjustment cost
    r_star      // World interest rate
    phi_adj_bat; // NEW: Installation Cost

// 3. CALIBRATION 
beta    = 0.99;    
delta   = 0.025;   
alpha   = 0.35;
sigma_L = 1.0;     
A       = 1.0;

// Energy & Reliability
psi     = 2.0;     
phi_int = 0.3;     
mu      = 0.16;
rho     = 0.4;     

// Innovation
eta_dtc = 0.5;     
delta_g = 0.015;   
delta_b = 0.025;
phi_g   = 2.0;     
phi_sub    = 0.0;  
cost_scale = 0.05; 
chi_price  = 1.0;

// Open Economy
theta_x = 1.5;     
phi_b   = 0.001;   
r_star  = (1/beta) - 1; 

// --- CRITICAL UPDATE: Adjustment Cost ---
// Matches your thesis_dtc.mod to prevent "instant" transition
phi_adj_bat = 50.0; 

// 4. THE MODEL
model;
    // --- HOUSEHOLDS ---
    1/C = beta * (1/C(+1)) * (r_k(+1) + 1 - delta);
    
    // Foreign Euler (UIP)
    1 = beta * (C / C(+1)) * (1 + r_star_eff) * (Q(+1) / Q);
    r_star_eff = r_star - phi_b * B_star;

    w = C * L^sigma_L;

    // --- FIRMS ---
    Y = A * ( (u * K(-1))^alpha ) * (L^(1-alpha));
    w = (1-alpha) * Y / L;
    r_k = alpha * Y / K(-1);

    // --- RELIABILITY ---
    // E_ren is exogenous here (The PDP8 Target)
    u = 1 - exp( -psi * (F / (phi_int * E_ren)) );
    F = ( mu * (A_bat * K_b(-1))^((rho-1)/rho) + (1-mu) * K_g(-1)^((rho-1)/rho) )^(rho/(rho-1));

    // --- INVESTMENT ---
    I_g = delta_g * 15 * (0.99 / u)^phi_g;
    K_g = (1 - delta_g) * K_g(-1) + I_g(-4); 

    MPK_b_physical = (Y / F) * (1 - u) * (mu * (F / (A_bat*K_b(-1)))^(1/rho)) * A_bat;
    
    // --- UPDATED BATTERY ARBITRAGE (With Adjustment Cost) ---
    // This slows down the private sector response, creating the "Transition Trap"
    P_tech * Q * (1 + phi_adj_bat * (I_b/K_b(-1) - delta_b)) = beta * (C / C(+1)) * (MPK_b_physical(+1) + P_tech(+1)*Q(+1)*(1 - delta_b));
    
    K_b = (1 - delta_b) * K_b(-1) + I_b;

    // --- INNOVATION ---
    tau_res = phi_sub * (0.99 - u);
    G_cost  = cost_scale * tau_res * Y;
    MPK_bat = MPK_b_physical;
    log(A_bat) = 0.95 * log(A_bat(-1)) + eta_dtc * (1 + tau_res) * chi_price * (log(MPK_bat) - log(MPK_bat(-1)));
    
    // --- EXTERNAL SECTOR ---
    P_tech = 1; // Constant price for this scenario
    
    // Exports
    X = 0.2 * Y * Q^theta_x;
    
    // Balance of Payments
    Q * B_star = Q * B_star(-1) * (1 + r_star_eff(-1)) + X - Q * (P_tech * I_b);
    TradeBal = X - Q * (P_tech * I_b);

    // --- CLOSING ---
    Y = C + I + I_g + G_cost + X;
    K = (1-delta)*K(-1) + I;
    Welfare = log(C) - (L^(1+sigma_L))/(1+sigma_L) + beta * Welfare(+1);
end;

// 5. TRANSITION SETUP 

// A. Start (2025)
initval;
    E_ren = 10;
    P_tech = 1;
    Q = 1;
    B_star = 0;
    r_star_eff = (1/beta) - 1;

    A_bat = 1; L = 0.93; u = 0.99;
    MPK_b_physical = (1/beta) - (1-delta_b);
    MPK_bat = MPK_b_physical;
    r_k = (1/beta) - (1-delta);
    
    K_g = 15; I_g = delta_g * K_g; 
    K_b = 2.0; I_b = delta_b * K_b; F = 6.0; 
    K = 30; I = delta * K;
    tau_res = 0; G_cost = 0;
    Y = 3.2;
    
    X = I_b; 
    C = Y - I - I_g - G_cost - X;
    w = (1-alpha)*Y/L;
    Welfare = (log(C) - (L^(1+sigma_L))/(1+sigma_L)) / (1-beta);
end;
steady; 

// B. End (2050)
endval;
    E_ren = 15; 
end;
steady;

// 6. CREATE THE RAMP & SOLVE
perfect_foresight_setup(periods=100);

// --- OVERWRITE PATH WITH LINEAR RAMP ---
verbatim;
    ramp_path = linspace(10, 15, size(oo_.exo_simul, 1))';
    oo_.exo_simul(:, strmatch('E_ren', M_.exo_names, 'exact')) = ramp_path;
end;
// ---------------------------------------

perfect_foresight_solver(robust_lin_solve);
