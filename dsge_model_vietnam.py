"""
DSGE Model: Renewable Intermittency and Grid Reliability (Vietnam)
CORRECTED VERSION - Fixed Reliability Valley Dynamics

Author: Toan T. Nguyen
Advisor: Dr. Xavier Martin G. Bautista
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import pandas as pd

# =========================================================================
# 1. PARAMETERS
# =========================================================================

class Parameters:
    """Store all model parameters"""
    
    def __init__(self):
        # Preferences
        self.beta = 0.99
        self.sigma_L = 1.0
        
        # Production
        self.alpha = 0.35
        self.delta_p = 0.025
        self.delta_b = 0.03
        self.delta_g = 0.0125
        
        # Energy
        self.omega_E = 0.05
        self.sigma_E = 0.5
        
        # Reliability (micro-founded)
        self.psi = 2.0
        self.phi_int = 0.30
        self.E_scale = 1.0
        
        # Flexibility (from PDP8)
        self.mu = 0.16
        self.rho = 0.40
        
        # Innovation
        self.eta_bat = 0.50
        self.chi = 0.70
        
        # Grid investment rule
        self.phi_grid = 1.5
        self.u_target = 0.97
        
        # Shock persistence
        self.rho_ren = 0.80
        
        # Normalizations
        self.A_ss = 1.0

# =========================================================================
# 2. STEADY STATE (Simplified Calibration)
# =========================================================================

def compute_steady_state(params):
    """
    Compute steady state with direct calibration
    """
    
    ss = {}
    
    # Capital stocks (from PDP8)
    ss['K_g'] = 3.0
    ss['K_b'] = 0.5  # Ratio matches PDP8: 0.5/3.0 ≈ 0.167
    ss['K_p'] = 10.0
    ss['A_bat'] = 1.0
    
    # Flexibility
    ss['F'] = (params.mu * (ss['A_bat'] * ss['K_b'])**((params.rho-1)/params.rho) + 
               (1-params.mu) * ss['K_g']**((params.rho-1)/params.rho))**(params.rho/(params.rho-1))
    
    # Utilization (target)
    ss['u'] = params.u_target
    
    # Verify utilization is consistent with flexibility
    # u = 1 - exp(-psi * F / (phi_int * E_scale))
    # This gives us the required F for target u
    # -psi * F / phi_int = ln(1 - u)
    # F_required = -phi_int * ln(1 - u) / psi
    F_required = -params.phi_int * params.E_scale * np.log(1 - params.u_target) / params.psi
    
    # Scale factor to match
    scale = F_required / ss['F']
    print(f"Flexibility scale factor needed: {scale:.3f}")
    
    # Adjust capital stocks proportionally
    ss['K_g'] = ss['K_g'] * scale
    ss['K_b'] = ss['K_b'] * scale
    ss['F'] = F_required
    
    # Basic macro variables
    ss['L'] = 0.33
    ss['Y'] = 1.0
    ss['E'] = 0.15 * ss['Y']
    ss['C'] = 0.65 * ss['Y']
    
    # Investment flows
    ss['I_p'] = params.delta_p * ss['K_p']
    ss['I_bat'] = params.delta_b * ss['K_b']
    ss['I_grid'] = params.delta_g * ss['K_g']
    
    # Value added
    ss['V'] = params.A_ss * (ss['u'] * ss['K_p'])**params.alpha * ss['L']**(1-params.alpha)
    
    # Prices
    ss['r'] = 1/params.beta - 1
    ss['W'] = ss['L']**params.sigma_L * ss['C']
    ss['R_k'] = params.alpha * params.A_ss * (ss['u'] * ss['K_p'])**(params.alpha-1) * ss['L']**(1-params.alpha)
    
    # Shadow value
    ss['lambda_F'] = (params.alpha * ss['Y'] / ss['K_p']) * params.psi * \
                     np.exp(-params.psi * ss['F'] / (params.phi_int * params.E_scale)) / \
                     (params.phi_int * params.E_scale)
    
    print(f"\nSteady State Computed:")
    print(f"  K_b = {ss['K_b']:.3f}, K_g = {ss['K_g']:.3f}")
    print(f"  F = {ss['F']:.3f}, u = {ss['u']:.4f}")
    print(f"  Shadow value λ_F = {ss['lambda_F']:.6f}")
    
    return ss

# =========================================================================
# 3. IMPULSE RESPONSE FUNCTIONS (Corrected)
# =========================================================================

def compute_irf(params, ss, shock_size=0.05, T=40):
    """
    Compute IRF to intermittency shock with corrected dynamics
    """
    
    # Initialize storage
    irf = {var: np.zeros(T) for var in ['Y', 'C', 'u', 'F', 'I_bat', 'I_grid', 
                                          'A_bat', 'lambda_F', 'K_b', 'K_g']}
    
    # Intermittency shock (AR1)
    eps_ren = np.zeros(T)
    eps_ren[0] = shock_size
    for t in range(1, T):
        eps_ren[t] = params.rho_ren * eps_ren[t-1]
    
    # State variables (levels)
    K_b = np.zeros(T)
    K_g = np.zeros(T)
    A_bat = np.zeros(T)
    
    K_b[0] = ss['K_b']
    K_g[0] = ss['K_g']
    A_bat[0] = ss['A_bat']
    
    for t in range(T):
        # [1] Flexibility
        F_t = (params.mu * (A_bat[t] * K_b[t])**((params.rho-1)/params.rho) + 
               (1-params.mu) * K_g[t]**((params.rho-1)/params.rho))**(params.rho/(params.rho-1))
        irf['F'][t] = (F_t / ss['F'] - 1) * 100
        
        # [2] Utilization (key: shock INCREASES denominator, LOWERS utilization)
        denominator = params.phi_int * params.E_scale * (1 + eps_ren[t])
        u_t = 1 - np.exp(-params.psi * F_t / denominator)
        irf['u'][t] = (u_t / ss['u'] - 1) * 100
        
        # [3] Output (reduced by lower utilization)
        # Output elasticity: approx alpha * 0.5 (due to CES structure)
        Y_t = ss['Y'] * (u_t / ss['u'])**(params.alpha * 0.5)
        irf['Y'][t] = (Y_t / ss['Y'] - 1) * 100
        
        # [4] Shadow value of flexibility
        lambda_F_t = (params.alpha * Y_t / ss['K_p']) * \
                     (params.psi / denominator) * \
                     np.exp(-params.psi * F_t / denominator)
        irf['lambda_F'][t] = (lambda_F_t / ss['lambda_F'] - 1) * 100
        
        # [5] Grid investment (counter-cyclical)
        I_grid_t = ss['I_grid'] * (params.u_target / u_t)**params.phi_grid
        I_grid_t = np.clip(I_grid_t, ss['I_grid'] * 0.5, ss['I_grid'] * 2.5)
        irf['I_grid'][t] = (I_grid_t / ss['I_grid'] - 1) * 100
        
        # [6] Battery investment (elastic response to shadow value)
        battery_elasticity = 0.8
        I_bat_t = ss['I_bat'] * (lambda_F_t / ss['lambda_F'])**battery_elasticity
        I_bat_t = np.clip(I_bat_t, ss['I_bat'] * 0.3, ss['I_bat'] * 2.5)
        irf['I_bat'][t] = (I_bat_t / ss['I_bat'] - 1) * 100
        
        # [7] Consumption (resource constraint)
        C_t = Y_t - ss['I_p'] - I_bat_t - I_grid_t
        C_t = max(C_t, 0.3 * ss['C'])
        irf['C'][t] = (C_t / ss['C'] - 1) * 100
        
        # [8] Update states
        if t < T - 1:
            K_b[t+1] = (1 - params.delta_b) * K_b[t] + I_bat_t
            K_g[t+1] = (1 - params.delta_g) * K_g[t] + I_grid_t
            
            # Learning (bounded)
            learning = params.eta_bat * params.chi * (lambda_F_t / ss['lambda_F'] - 1)
            learning = np.clip(learning, -0.01, 0.03)
            A_bat[t+1] = A_bat[t] * (1 + learning)
            
            irf['K_b'][t+1] = (K_b[t+1] / ss['K_b'] - 1) * 100
            irf['K_g'][t+1] = (K_g[t+1] / ss['K_g'] - 1) * 100
            irf['A_bat'][t+1] = (A_bat[t+1] / ss['A_bat'] - 1) * 100
    
    return irf, eps_ren

# =========================================================================
# 4. RELIABILITY VALLEY (Corrected Transition Path)
# =========================================================================

def compute_reliability_valley(params, ss, T=100):
    """
    Simulate transition where renewable capacity grows faster than flexibility
    KEY: Front-loaded renewable deployment vs gradual flexibility buildout
    """
    
    # Renewable penetration: Rapid growth in years 2-8, then slower
    # Models aggressive PDP8 targets with front-loaded solar/wind deployment
    ren_initial = 0.15
    ren_final = 0.50
    
    # S-curve with early acceleration (models subsidy-driven solar boom)
    x = np.linspace(0, 1, T)
    ren_path = ren_initial + (ren_final - ren_initial) * (1 / (1 + np.exp(-10*(x - 0.35))))
    
    # Flexibility assets: slower initial growth (planning/construction lags)
    # Then accelerates as reliability concerns mount
    F_path = np.zeros(T)
    K_b_path = np.zeros(T)
    K_g_path = np.zeros(T)
    
    F_path[0] = ss['F']
    K_b_path[0] = ss['K_b']
    K_g_path[0] = ss['K_g']
    
    for t in range(1, T):
        # Grid capital: steady growth with moderate acceleration
        grid_growth = 0.013 if t < 40 else 0.017
        K_g_path[t] = K_g_path[t-1] * (1 + grid_growth)
        
        # Battery capital: strong growth
        battery_growth = 0.019
        K_b_path[t] = K_b_path[t-1] * (1 + battery_growth)
        
        # Aggregate flexibility (CES)
        F_path[t] = (params.mu * K_b_path[t]**((params.rho-1)/params.rho) + 
                     (1-params.mu) * K_g_path[t]**((params.rho-1)/params.rho))**(params.rho/(params.rho-1))
    
    # Utilization path
    u_path = np.zeros(T)
    
    for t in range(T):
        # Intermittency scales with renewable penetration
        # Higher VRE share → more volatility → harder to maintain reliability
        # Use power of 1.1 (not 1.2) for more moderate scaling
        intermittency_factor = params.phi_int * params.E_scale * (ren_path[t] / ren_initial)**1.1
        
        u_path[t] = 1 - np.exp(-params.psi * F_path[t] / intermittency_factor)
    
    # Find valley (minimum utilization)
    valley_idx = np.argmin(u_path)
    
    return ren_path, u_path, F_path, K_b_path, K_g_path, valley_idx

# =========================================================================
# 5. WELFARE CALCULATION
# =========================================================================

def compute_welfare_cost(irf, ss, params, T=40):
    """
    Consumption-equivalent welfare loss using proper utility metric
    """
    
    # Baseline utility (steady state)
    U_ss = np.log(ss['C']) - (ss['L']**(1 + params.sigma_L)) / (1 + params.sigma_L)
    
    # Lifetime utility with shocks
    U_shock = 0
    for t in range(T):
        C_t = ss['C'] * (1 + irf['C'][t] / 100)
        L_t = ss['L']  # Assume labor roughly constant
        U_t = np.log(C_t) - (L_t**(1 + params.sigma_L)) / (1 + params.sigma_L)
        U_shock += params.beta**t * U_t
    
    # Lifetime utility without shocks
    U_baseline = U_ss / (1 - params.beta)
    
    # Find consumption compensation: (1+lambda)*C gives same utility
    # Solve: log((1+lambda)*C) = U_shock implies lambda = exp(U_shock - U_baseline) - 1
    # Approximate for small changes: lambda ≈ (U_shock - U_baseline)
    
    welfare_loss_pct = (U_baseline - U_shock) / abs(U_baseline) * 100
    
    return abs(welfare_loss_pct)

# =========================================================================
# 6. MAIN EXECUTION
# =========================================================================

if __name__ == "__main__":
    
    print("="*70)
    print("Vietnam Renewable Energy DSGE Model - CORRECTED")
    print("="*70)
    
    params = Parameters()
    
    # Steady state
    print("\n[1] Computing Steady State...")
    ss = compute_steady_state(params)
    
    # IRF
    print("\n[2] Computing Impulse Responses...")
    irf, shock = compute_irf(params, ss, shock_size=0.05, T=40)
    
    print("\nPeak IRF Responses (% from SS):")
    for var in ['Y', 'C', 'u', 'F', 'I_grid', 'I_bat', 'lambda_F']:
        peak_neg = np.min(irf[var])
        peak_pos = np.max(irf[var])
        peak = peak_neg if abs(peak_neg) > abs(peak_pos) else peak_pos
        print(f"  {var:10s}: {peak:8.3f}%")
    
    # Welfare
    welfare_cost = compute_welfare_cost(irf, ss, params)
    print(f"\n[3] Welfare Cost: {abs(welfare_cost):.4f}% of lifetime consumption")
    
    # Reliability valley
    print("\n[4] Computing Reliability Valley...")
    ren_path, u_path, F_path, K_b_path, K_g_path, valley_idx = compute_reliability_valley(params, ss, T=100)
    
    valley_year = valley_idx / 4
    valley_u = u_path[valley_idx]
    valley_gap = (valley_u - params.u_target) / params.u_target * 100
    
    print(f"  Valley occurs at quarter {valley_idx} ({valley_year:.1f} years)")
    print(f"  Minimum utilization: {valley_u:.4f} ({valley_u*100:.2f}%)")
    print(f"  Gap from target: {valley_gap:.2f}%")
    
    # =====================================================================
    # PLOTTING
    # =====================================================================
    
    print("\n[5] Generating Figures...")
    
    # IRF Figure
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle('Impulse Response to 5% Renewable Intermittency Shock\n(Corrected Dynamics)', 
                 fontsize=14, fontweight='bold')
    
    variables = [
        ('Y', 'Output'),
        ('C', 'Consumption'),
        ('u', 'Utilization (Reliability)'),
        ('F', 'Flexibility Assets'),
        ('I_grid', 'Grid Investment'),
        ('I_bat', 'Battery Investment'),
        ('lambda_F', 'Shadow Value of Flexibility'),
        ('A_bat', 'Battery Technology'),
        ('K_b', 'Battery Capital Stock')
    ]
    
    for idx, (var, label) in enumerate(variables):
        ax = axes[idx // 3, idx % 3]
        ax.plot(irf[var], linewidth=2.5, color='darkblue')
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.set_title(label, fontweight='bold', fontsize=11)
        ax.set_xlabel('Quarters', fontsize=9)
        ax.set_ylabel('% Deviation from SS', fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, 40])
        
        # Highlight key features
        if var == 'u':
            ax.axhline(y=-3, color='red', linestyle=':', linewidth=1, alpha=0.5)
            ax.text(30, -3.5, 'Reliability drops', fontsize=8, color='red')
        elif var == 'I_bat':
            peak_t = np.argmax(irf[var])
            ax.scatter([peak_t], [irf[var][peak_t]], color='red', s=50, zorder=5)
            ax.text(peak_t+2, irf[var][peak_t], f'{irf[var][peak_t]:.1f}%', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/irf_corrected.png', dpi=300, bbox_inches='tight')
    print("  Saved: irf_corrected.png")
    
    # Reliability Valley Figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    quarters = np.arange(len(ren_path))
    years = quarters / 4
    
    ax1.plot(years, ren_path * 100, linewidth=2.5, label='Renewable Penetration', color='green')
    ax1.set_xlabel('Years', fontsize=11)
    ax1.set_ylabel('Renewable Share of Energy (%)', fontsize=11)
    ax1.set_title('PDP8 Transition Path:\nRenewable Deployment', fontweight='bold', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_ylim([10, 55])
    
    ax2.plot(years, u_path * 100, linewidth=2.5, label='Grid Reliability', color='darkblue')
    ax2.axhline(y=params.u_target * 100, color='red', linestyle='--', linewidth=2, 
                label=f'Target ({params.u_target*100:.0f}%)', alpha=0.7)
    ax2.axvline(x=valley_year, color='orange', linestyle=':', linewidth=2.5, alpha=0.7,
                label=f'Reliability Valley (Year {valley_year:.1f})')
    ax2.fill_between(years, params.u_target * 100, u_path * 100, 
                     where=(u_path < params.u_target), alpha=0.3, color='red',
                     label='Below Target')
    ax2.set_xlabel('Years', fontsize=11)
    ax2.set_ylabel('Grid Reliability (%)', fontsize=11)
    ax2.set_title('The "Reliability Valley"\nDuring Energy Transition', fontweight='bold', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    ax2.set_ylim([90, 99])
    
    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/valley_corrected.png', dpi=300, bbox_inches='tight')
    print("  Saved: valley_corrected.png")
    
    # Export data
    results_df = pd.DataFrame({
        'Quarter': range(40),
        'Output': irf['Y'],
        'Consumption': irf['C'],
        'Utilization': irf['u'],
        'Flexibility': irf['F'],
        'Grid_Investment': irf['I_grid'],
        'Battery_Investment': irf['I_bat'],
        'Shadow_Value': irf['lambda_F']
    })
    results_df.to_csv('/mnt/user-data/outputs/irf_corrected.csv', index=False)
    
    valley_df = pd.DataFrame({
        'Quarter': quarters,
        'Year': years,
        'Renewable_Penetration': ren_path,
        'Reliability': u_path,
        'Flexibility': F_path,
        'Battery_Capital': K_b_path,
        'Grid_Capital': K_g_path
    })
    valley_df.to_csv('/mnt/user-data/outputs/valley_corrected.csv', index=False)
    
    ss_df = pd.DataFrame([ss])
    ss_df.to_csv('/mnt/user-data/outputs/steady_state_corrected.csv', index=False)
    
    print("\n" + "="*70)
    print("Simulation Complete (Corrected Version)")
    print("="*70)
