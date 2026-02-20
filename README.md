# Renewable Intermittency and Grid Reliability: A Directed Technical Change Model of Vietnam

**Author:** Toan T. Nguyen (Justin)  
**Advisor:** Dr. Xavier Martin G. Bautista  
**Institution:** Fulbright University Vietnam  
**Date:** January 2026

## 1. Project Overview
This repository contains the Dynare/MATLAB codebase for the capstone thesis **"Renewable Intermittency and Grid Reliability: A Directed Technical Change Model of Vietnam."**

The project develops a **Small Open Economy DSGE model** calibrated to Vietnam's **Power Development Plan 8 (PDP8)**. It investigates the macroeconomic costs of the green transition ("Greenflation") arising from the structural mismatch between variable renewable supply and grid flexibility.

### Key Economic Mechanisms
1.  **Endogenous Reliability Penalty:** Renewable intermittency acts as a negative Total Factor Productivity (TFP) shock when grid flexibility is insufficient.
2.  **The "Agility Gap":** A structural friction where private battery adoption is fast but expensive (imported), while public grid investment is cheap (domestic) but slow due to "Time-to-Build" lags.
3.  **Directed Technical Change (DTC):** Innovation in battery technology is endogenously driven by the "Shadow Price of Reliability".
4.  **Small Open Economy:** Vietnam imports battery technology, linking energy security to the trade balance and global price shocks.

---

## 2. Core Model Files (`.mod`)

These two files are the "engines" of the thesis. They share the same underlying economic structure but are solved under different environments (Stochastic vs. Deterministic).

### A. `thesis_dtc.mod` (The Stochastic Model)
* **Purpose:** Simulates the economy's reaction to **"Surprise" Shocks** (e.g., a sudden renewable drought or a global battery price spike). Used for Impulse Response Functions (IRFs) and Welfare Analysis.
* **Key Features:**
    * **Unit Correction (`E_scale = 10.0`):** Normalizes energy units to match the scale of economic capital ($K_g \approx 15$) to ensure the reliability function is numerically stable.
    * **AR(1) Shock Process:** Renewable deficits are modeled as persistent shocks ($\rho = 0.9$) to capture the 1.5-year duration of **ENSO climate cycles** (El Niño) in Vietnam, rather than white noise.
    * **The Reliability Function:**
        $$u_t = 1 - \exp\left(-\psi \frac{F_t}{\phi_{int} \cdot E_{scale} \cdot (1 + \hat{E}_{t})}\right)$$
        This equation serves as the "TFP Link" between energy physics and the macroeconomy.

### B. `thesis_pdp8.mod` (The Transition Model)
* **Purpose:** Simulates the **Deterministic Transition Path (2025–2050)**. It solves for the perfect foresight path as Vietnam implements the PDP8 renewable targets.
* **Key Features:**
    * **The Linear Ramp:** Forces the renewable target `E_ren` to rise linearly from 10 to 15 over 100 quarters (25 years).
    * **The "Transition Trap":** Demonstrates that reliability ($u$) dips during the mid-transition phase because grid investment ($K_g$) lags behind the renewable ramp due to the 4-quarter construction delay.

---

## 3. Analysis Scripts (`.m`)

Run these MATLAB scripts to reproduce the specific figures and tables in the thesis.  
**Note:** These scripts automatically load `thesis_dtc.mod`.

| Script Name | Output | Description & Thesis Link |
| :--- | :--- | :--- |
| **`run_welfare_analysis.m`** | **CEV Calculation** | Calculates the **Consumption Equivalent Variation**. It compares lifetime utility during a crisis vs. steady state. Outputs the % of consumption households would pay to avoid the "Reliability Valley". |
| **`run_sensitivity_rho.m`** | **Figure A.2** | **"Agility Gap" Robustness Check.** Compares the investment response when Grid/Batteries are Complements ($\rho=0.4$) vs. Substitutes ($\rho=0.9$). Proves that complementarity drives the volatility. |
| **`run_optimal_subsidy.m`** | **Figure 4.7** | **Optimal R&D Policy.** Loops through subsidy levels ($\phi_{sub}$) to find the peak welfare point. Shows the trade-off between tech gain and tax distortion. |
| **`run_second_best_policy.m`** | **Figure 8** | **Fiscal Impotence Analysis.** Compares optimal subsidies under Perfect Markets ($\chi=1$) vs. Regulated Markets ($\chi=0.3$). Shows that subsidies fail if price signals are jammed. |
| **`run_soe_shock.m`** | **Figure 6** | **Imported Greenflation.** Simulates a 10% global battery price shock. Shows how external tech dependence hurts domestic GDP. |
| **`run_transition_pdp8.m`** | **Chapter 6 Figures** | Plots the full 25-year transition path, highlighting the "Reliability Valley" and the Trade Balance deficit. |
| **`run_optimal_mix.m`** | **Figure 5** | **Optimal Flexibility Mix.** Varying the battery share $\mu$ to find the optimal portfolio. Shows Vietnam's current mix ($\mu=0.16$) is suboptimal. |
| **`run_counterfactual.m`** | **Figure 4** | **Value of Innovation.** Compares the recession depth with and without endogenous DTC. Highlights the "Innovation Rescue Effect." |

---

## 4. "Cheatsheet" of Critical Parameters

If examining the calibration, these are the most impactful parameters defining the Vietnam scenario:

| Parameter | Value | Description |
| :--- | :--- | :--- |
| **`mu` (Battery Weight)** | **0.16** | Derived from PDP8 investment data (\$2.4B Battery vs. \$14.9B Grid). Reflects the grid's heavy reliance on transmission over storage. |
| **`rho` (Substitution Elasticity)** | **0.4** | Sets Grid and Batteries as **Complements**, not substitutes. This is the mathematical root of the "Agility Gap". |
| **`rho_ren_shock` (Persistence)** | **0.9** | Matches the half-life of **El Niño** events (~1.5 years), ensuring the stress test is physically realistic for Vietnam. |
| **`phi_adj_bat` (Friction)** | **50.0** | Represents installation bottlenecks. Prevents the private sector from instantly fixing the grid, creating the "Agility Gap" window. |

---

## 5. How to Run the Code

1.  **Prerequisites:** Ensure MATLAB is installed with **Dynare 4.6** or higher. Add the Dynare `matlab` folder to your path.
2.  **To Run a Specific Analysis:**
    * Open the desired `.m` script (e.g., `run_welfare_analysis.m`).
    * Click **Run**. The script will automatically clean up old files, call Dynare, and generate the plots.
3.  **To Run the Raw Models:**
    * **Stochastic:** Type `dynare thesis_dtc` in the Command Window.
    * **Transition:** Type `dynare thesis_pdp8` in the Command Window.
