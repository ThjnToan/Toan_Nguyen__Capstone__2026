// ==================================================================
// FILE: thesis_dtc.mod
// SCENARIO: Small Open Economy with Imported Battery Tech
// ==================================================================

// 1. VARIABLES
var 
    Y C I K L w r_k u F 
    A_bat MPK_bat               
    K_g K_b I_g I_b             
    MPK_b_physical              
    tau_res G_cost              
    Welfare
    // --- NEW OPEN ECONOMY VARIABLES ---
    Q           // Real Exchange Rate
    X           // Exports
    B_star      // Net Foreign Assets
    r_star_eff  // Effective Foreign Interest Rate
    P_tech      // Imported Battery Price (Shock Process)
    TradeBal;   // Trade Balance (X - Q*Imports)

varexo 
    e_ren;      // Renewable Supply Shock
    e_price;    // GLOBAL BATTERY PRICE SHOCK (New)

// 2. PARAMETERS
parameters 
    beta delta alpha psi phi_int mu rho sigma_L A 
    eta_dtc delta_g delta_b phi_g
    phi_sub cost_scale chi_price
    // --- NEW OPEN ECONOMY PARAMETERS ---
    theta_x     // Export elasticity
    phi_b       // Debt adjustment cost (risk premium)
    r_star;     // World interest rate

// 3. CALIBRATION (MUST MATCH PDP8 EXACTLY)
beta    = 0.99;    
delta   = 0.025;   
alpha   = 0.35;    
sigma_L = 1.0;     
A       = 1.0;     

// Energy & Reliability Parameters (The Critical Block)
psi     = 2.0;     // Matches PDP8 ramp
phi_int = 0.3;     
mu      = 0.16;    // Vietnam Calibration (Low Storage)
rho     = 0.4;     // Complementarity

// Innovation & Rigor
eta_dtc = 0.5;     
delta_g = 0.015;   
delta_b = 0.025;   
phi_g   = 2.0;     
phi_sub    = 0.0;  // Baseline: No Subsidy
cost_scale = 0.05; 
chi_price  = 1.0;  // Baseline: Perfect Market (We shock this in Ch 7)

// Open Economy Calibration
theta_x = 1.5;     // Standard trade elasticity
phi_b   = 0.001;   // Small friction to close the model
r_star  = (1/beta) - 1; // 0.0101

// 4. THE MODEL
model;
    // --- HOUSEHOLDS & FIRMS ---
    // 1. Domestic Euler
    1/C = beta * (1/C(+1)) * (r_k(+1) + 1 - delta);
    
    // 2. Foreign Euler (UIP Condition)
    // Determines Real Exchange Rate (Q)
    1 = beta * (C / C(+1)) * (1 + r_star_eff) * (Q(+1) / Q);
    
    // 3. Effective Foreign Rate (Debt Elastic) 
    r_star_eff = r_star - phi_b * B_star;

    // 4. Labor Supply
    w = C * L^sigma_L;

    // 5. Production
    Y = A * ( (u * K(-1))^alpha ) * (L^(1-alpha));
    w = (1-alpha) * Y / L;
    r_k = alpha * Y / K(-1);

    // --- RELIABILITY BLOCK ---
    u = 1 - exp( -psi * (F / (phi_int * (1+e_ren))) );
    F = ( mu * (A_bat * K_b(-1))^((rho-1)/rho) + (1-mu) * K_g(-1)^((rho-1)/rho) )^(rho/(rho-1));

    // --- INVESTMENT & ACCUMULATION ---
    I_g = delta_g * 15 * (0.99 / u)^phi_g;
    K_g = (1 - delta_g) * K_g(-1) + I_g(-4); 

    MPK_b_physical = (Y / F) * (1 - u) * (mu * (F / (A_bat*K_b(-1)))^(1/rho)) * A_bat;
    
    // 6. Battery Arbitrage (Modified for Imported Tech)
    // The cost of battery investment depends on Global Price (P_tech) and Exchange Rate (Q)
    P_tech * Q = beta * (C / C(+1)) * (MPK_b_physical(+1) + P_tech(+1)*Q(+1)*(1 - delta_b));
    K_b = (1 - delta_b) * K_b(-1) + I_b;

    // --- INNOVATION ---
    tau_res = phi_sub * (0.99 - u);
    G_cost  = cost_scale * tau_res * Y;
    MPK_bat = MPK_b_physical;
    log(A_bat) = 0.95 * log(A_bat(-1)) + eta_dtc * (1 + tau_res) * chi_price * (log(MPK_bat) - log(MPK_bat(-1)));

    // --- EXTERNAL SECTOR ---
    // 7. Global Price Shock Process
    log(P_tech) = 0.9 * log(P_tech(-1)) + e_price; 

    // 8. Export Demand [cite: 500]
    X = 0.2 * Y * Q^theta_x; // Calibrated so Exports are ~20% of GDP in SS

    // 9. Balance of Payments [cite: 506]
    // Change in Net Foreign Assets = Net Exports (Exports - Value of Imports)
    // We import Battery Investment (I_b).
    Q * B_star = Q * B_star(-1) * (1 + r_star_eff(-1)) + X - Q * (P_tech * I_b);
    
    // 10. Trade Balance Definition
    TradeBal = X - Q * (P_tech * I_b);

    // --- CLOSING ---
    // 11. Market Clearing (Open Economy)
    // Y is domestic output. It is consumed, invested in K, Grid, or Exported.
    // Note: I_b is NOT here because it is imported (foreign good).
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
    
    MPK_b_physical = (1/beta) - (1-delta_b); // Approx for initial guess
    MPK_bat = MPK_b_physical;
    r_k = (1/beta) - (1-delta);
    
    K_g = 15;
    I_g = delta_g * K_g; 
    K_b = 2.0;
    I_b = delta_b * K_b;
    F   = 6.0; 
    K = 30;
    I = delta * K;
    tau_res = 0;
    G_cost = 0;
    Y = 3.2;
    
    // SOE Closing
    X = I_b; // In SS with B*=0, Exports pay for Battery Imports
    // Adjust constant in export equation implicitly or let solver find Q.
    // Here we let Dynare solve it.
    
    C = Y - I - I_g - G_cost - X;
    w = (1-alpha)*Y/L;
    Welfare = (log(C) - (L^(1+sigma_L))/(1+sigma_L)) / (1-beta);
end;
steady; 

// 6. SIMULATION (Impulse Response)
shocks;
    var e_ren; stderr 0.10; // 10% Renewable Supply Shock
    var e_price; stderr 0.10; // 10% Global Price Shock capability
end;

// Only run IRFs for the baseline shock, suppress others until called
stoch_simul(order=1, irf=40, nograph);

// 7. PLOTTING (CORRECTED)
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
    
    % CHANGED: Plot Flexibility (F) instead of the Shock (E_ren)
    subplot(1,3,3); 
    plot(T, oo_.irfs.F_e_ren, 'b-', 'LineWidth', 2); 
    title('Flexibility Stock (F)'); grid on; axis tight;
    % Save as image_fig1.png when generated

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
    % Save as image_fig2.png when generated

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
    % Save as image_fig3.png when generated

end;
