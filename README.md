# Renewable Intermittency and Grid Reliability: A Directed Technical Change Model of Vietnam

**Author:** Toan T. Nguyen (Justin)  
**Advisor:** Dr. Xavier Martin G. Bautista  
**Institution:** Fulbright University Vietnam  
**Date:** Spring 2026

---

## 1. Project Overview

This repository contains the MATLAB/Dynare codebase and LaTeX manuscript for the capstone thesis investigating the macroeconomic challenges of Vietnam's clean energy transition.

The project develops a **Small Open Economy DSGE model** calibrated to Vietnam's **Power Development Plan 8 (PDP8)**. It examines the structural mismatch between rapid Variable Renewable Energy (VRE) deployment and the slower accumulation of grid flexibility assets.

### Key Economic Mechanisms

1. **Reliability Penalty (Endogenous TFP):** Renewable intermittency is modeled as a negative TFP shock governed by an exponential reliability constraint. When grid flexibility is insufficient, installed physical capital cannot be fully utilized.

2. **The "Agility Gap":** A structural friction where private battery adoption is agile but exposed to global price shocks, while public grid investment responds sluggishly due to bureaucratic "Time-to-Build" lags.

3. **Endogenous Learning-by-Doing:** Battery technology improves endogenously, driven by the "Shadow Price of Reliability"—the scarcity signal generated when the grid underperforms.

4. **Current Account Vulnerability:** Global battery price spikes act as a macroeconomic trilemma, simultaneously draining the trade balance and increasing sovereign borrowing premiums.

---

## 2. File Structure

```
/
├── main.tex                      # LaTeX manuscript
├── main.pdf                      # Compiled thesis
├── vietnam_dsge.mod             # Dynare structural model
├── run_dynare.m                 # Master execution script
├── vietnam_dsge_steadystate.m  # Analytical steady-state solver
└── save_to_github.bat          # Auto-commit helper
```

---

## 3. Core Model (`vietnam_dsge.mod`)

The Dynare model defines 16 linearized structural equations for the small open economy.

**Key equations:**
- **Steady state (lines 143-144):** `r_bar = 1/beta - 1` — Foreign interest rate tied to discount factor
- **Reliability (line 175):** `u = xi * (F - Vol_ren)` — Utilization falls when flexibility fails to keep pace with renewable volatility
- **Flexibility CES (line 176-177):** `F = s_b * (A_bat + K_b) + (1 - s_b) * K_g` — Aggregate flexibility from batteries and grid
- **Battery investment (line 183):** `I_bat = -phi_grid * u - P_bat` — Investment rises when reliability falls
- **Learning-by-Doing (line 185):** `A_bat - A_bat(-1) = -eta_bat * chi * u` — Innovation accelerates with reliability scarcity

---

## 4. Running the Code

### Prerequisites
- MATLAB (R2023a or newer)
- [Dynare](https://www.dynare.org/) (version 5.x or 6.x)
- Add Dynare to MATLAB path: `addpath('c:\dynare\6.1\matlab')`

### Execution
1. Open `run_dynare.m` in MATLAB
2. Run the script (takes ~5-10 seconds)
3. Outputs include:
   - Generated IRF figures
   - `agility_gap_dynare.png` — 25-year reliability valley simulation
   - `counterfactual_chi.png` — Welfare analysis across market liberalization scenarios
   - `joint_shock_perfect_storm.png` — Combined shock scenario

### Compile Manuscript
Run `xelatex main.tex` or use your preferred LaTeX editor.

---

## 5. Version Control

Double-click `save_to_github.bat` to auto-commit and push changes to the `main` branch.