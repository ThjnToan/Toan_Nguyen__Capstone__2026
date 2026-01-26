// ==================================================================
// FILE: thesis_dtc.mod (FINAL FIXED VERSION)
// ==================================================================

// 1. VARIABLES
var 
    Y C I K L w r_k u F 
    A_bat MPK_bat               
    K_g K_b I_g I_b             
    MPK_b_physical              
    tau_res G_cost          
    Welfare
    Q X B_star r_star_eff P_tech TradeBal
    // Persistent Shock Variable
    E_ren_process; 

varexo 
    e_ren_shock // Domestic "Stress Test" Shock
    e_price;    // Global Price Shock

// 2. PARAMETERS
parameters 
    beta delta alpha psi phi_int mu rho sigma_L A 
    eta_dtc delta_g delta_b phi_g
    phi_sub cost_scale chi_price
    theta_x phi_b r_star phi_adj_bat E_scale
    rho_ren_shock; // Persistence Parameter

// 3. CALIBRATION 
beta    = 0.99;    
delta   = 0.025;
alpha   = 0.35;    
sigma_L = 1.0;     
A       = 1.0;

psi     = 2.0;     
phi_int = 0.3;     
mu      = 0.16;    
rho     = 0.4;
E_scale = 10.0; // Scaling Fix

eta_dtc    = 0.5;     
delta_g    = 0.015;   
delta_b    = 0.025;   
phi_g      = 2.0;     
phi_sub    = 0.0;     
cost_scale = 0.05; 
chi_price  = 1.0;     

theta_x = 1.5;            
phi_b   = 0.001;          
r_star  = (1/beta) - 1;   
phi_adj_bat = 50.0;

// Persistence (0.9 = shock lasts ~1.5 years)
rho_ren_shock = 0.9; 

// 4. THE MODEL
model;
    // --- HOUSEHOLDS & FIRMS ---
    1/C = beta * (1/C(+1)) * (r_k(+1) + 1 - delta);
    1 = beta * (C / C(+1)) * (1 + r_star_eff) * (Q(+1) / Q);
    r_star_eff = r_star - phi_b * B_star;
    w = C * L^sigma_L;

    Y = A * ( (u * K(-1))^alpha ) * (L^(1-alpha));
    w = (1-alpha) * Y / L;
    r_k = alpha * Y / K(-1);

    // --- RELIABILITY BLOCK (FIXED) ---
    // 1. AR(1) Process for the Grid Stress
    E_ren_process = rho_ren_shock * E_ren_process(-1) + e_ren_shock;

    // 2. Reliability Function
    // FIX: Changed denominator to (1 + E_ren_process) to avoid division by zero.
    // Logic: A POSITIVE shock to E_ren_process increases the "Intermittency Load",
    // causing u (Reliability) to DROP.
    u = 1 - exp( -psi * (F / (phi_int * E_scale * (1 + E_ren_process))) );
    
    F = ( mu * (A_bat * K_b(-1))^((rho-1)/rho) + (1-mu) * K_g(-1)^((rho-1)/rho) )^(rho/(rho-1));

    // --- INVESTMENT ---
    I_g = delta_g * 15 * (0.99 / u)^phi_g;
    K_g = (1 - delta_g) * K_g(-1) + I_g(-4); 

    MPK_b_physical = (Y / F) * (1 - u) * (mu * (F / (A_bat*K_b(-1)))^(1/rho)) * A_bat;
    
    P_tech * Q * (1 + phi_adj_bat * (I_b/K_b(-1) - delta_b)) 
        = beta * (C / C(+1)) * (MPK_b_physical(+1) + P_tech(+1)*Q(+1)*(1 - delta_b));
    
    K_b = (1 - delta_b) * K_b(-1) + I_b;

    // --- INNOVATION ---
    tau_res = phi_sub * (0.99 - u);
    G_cost  = cost_scale * tau_res * Y;
    MPK_bat = MPK_b_physical;
    log(A_bat) = 0.95 * log(A_bat(-1)) + eta_dtc * (1 + tau_res) * chi_price * (log(MPK_bat) - log(MPK_bat(-1)));

    // --- EXTERNAL SECTOR ---
    log(P_tech) = 0.9 * log(P_tech(-1)) + e_price; 
    X = 0.2 * Y * Q^theta_x; 
    Q * B_star = Q * B_star(-1) * (1 + r_star_eff(-1)) + X - Q * (P_tech * I_b);
    TradeBal = X - Q * (P_tech * I_b);

    // --- CLOSING ---
    Y = C + I + I_g + G_cost + X;
    K = (1-delta)*K(-1) + I;
    Welfare = log(C) - (L^(1+sigma_L))/(1+sigma_L) + beta * Welfare(+1);
end;

// 5. STEADY STATE
initval;
    E_ren_process = 0; 
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
    
    X = I_b; 
    C = Y - I - I_g - G_cost - X;
    w = (1-alpha)*Y/L;
    Welfare = (log(C) - (L^(1+sigma_L))/(1+sigma_L)) / (1-beta);
end;
steady; 

// 6. SIMULATION
shocks;
    // STRESS TEST: Positive shock increases Intermittency Load
    // This causes Reliability (u) to fall.
    var e_ren_shock; stderr 0.10; 
    var e_price; stderr 0.10; 
end;

stoch_simul(order=1, irf=40, nograph);

// 7. PLOTTING
verbatim;
    T = 1:40;
    figure('Name', 'Fig1_Reliability');
    subplot(1,3,1); plot(T, oo_.irfs.Y_e_ren_shock, 'k-', 'LineWidth', 2); title('Output (Y)'); grid on; axis tight;
    subplot(1,3,2); plot(T, oo_.irfs.u_e_ren_shock, 'r-', 'LineWidth', 2); title('Reliability (u)'); grid on; axis tight;
    subplot(1,3,3); plot(T, oo_.irfs.F_e_ren_shock, 'b-', 'LineWidth', 2); title('Flexibility Stock (F)'); grid on; axis tight;
    
    figure('Name', 'Fig2_Macro');

    subplot(2,2,1); plot(T, oo_.irfs.I_e_ren_shock, 'b-', 'LineWidth', 2); title('Investment (I)'); grid on; axis tight;
    subplot(2,2,2); plot(T, oo_.irfs.L_e_ren_shock, 'k-', 'LineWidth', 2); title('Labor (L)'); grid on; axis tight;
    subplot(2,2,3); plot(T, oo_.irfs.C_e_ren_shock, 'k--', 'LineWidth', 2); title('Consumption (C)'); grid on; axis tight;
    subplot(2,2,4); plot(T, oo_.irfs.w_e_ren_shock, 'r-', 'LineWidth', 2); title('Real Wages (w)'); grid on; axis tight;
    
    figure('Name', 'Fig3_Innovation');
    subplot(1,3,1); plot(T, oo_.irfs.MPK_bat_e_ren_shock, 'r-', 'LineWidth', 2); title('Shadow Price (MPK)'); grid on; axis tight;
    subplot(1,3,2); plot(T, oo_.irfs.A_bat_e_ren_shock, 'g-', 'LineWidth', 2); title('Battery Tech (A_{bat})'); grid on; axis tight;
    subplot(1,3,3); plot(T, oo_.irfs.I_b_e_ren_shock, 'b-', 'LineWidth', 2); title('Battery Inv (I_b)'); grid on; axis tight;
end;
