# Renewable Intermittency and Grid Reliability: A Directed Technical Change Model of Vietnam

**Author:** Toan T. Nguyen (Justin)  
**Advisor:** Dr. Xavier Martin G. Bautista  
**Institution:** Fulbright University Vietnam  
**Date:** Spring 2026

---

## 1. Project Abstract & Theoretical Contribution

This repository contains the complete MATLAB/Dynare codebase and LaTeX manuscript for the capstone thesis: **"Renewable Intermittency and Grid Reliability: A Directed Technical Change Model of Vietnam."**

The project develops a novel **Small Open Economy DSGE model** calibrated to Vietnam's **Power Development Plan 8 (PDP8)**. It investigates the macroeconomic challenges of the clean energy transition, specifically the structural mismatch between the rapid deployment of Variable Renewable Energy (VRE) and the sluggish accumulation of grid flexibility assets. 

### Key Economic Mechanisms Modeled
1. **The Reliability Penalty (Endogenous TFP):** Renewable intermittency is not modeled simply as a cost shock, but as a negative TFP shock governed by an exponential reliability constraint. When grid flexibility is insufficient, installed physical capital cannot be fully utilized.
2. **The "Agility Gap":** The model captures a structural friction in investment dynamics. Private battery adoption is agile but exposed to global price shocks, while public grid transmission investment is insulated from international prices but responds sluggishly due to bureaucratic "Time-to-Build" lags.
3. **Endogenous Learning-by-Doing (Directed Technical Change):** Battery technology improves endogenously. Crucially, the speed of innovation is driven by the "Shadow Price of Reliability"—the scarcity signal generated when the grid underperforms.
4. **Current Account Vulnerability:** Vietnam is modeled as a technology importer. Global battery price spikes (e.g., lithium shortages) don't just delay the green transition; they act as a macroeconomic trilemma by simultaneously draining the trade balance and increasing the sovereign borrowing premium.

---

## 2. Codebase Architecture & File Structure

The project has been aggressively consolidated for elegance and reproducibility. The entire analysis runs through a core Dynare model and a single MATLAB execution script.

### Directory Overview
```text
/
├── main.tex                  # The comprehensive LaTeX manuscript (compiles to main.pdf)
├── main.pdf                  # The final compiled thesis document
├── vietnam_dsge.mod          # The core Dynare structural model
├── run_dynare.m              # The master MATLAB execution script
├── vietnam_dsge_steadystate.m# Analytical steady-state solver for Dynare
└── save_to_github.bat        # Auto-commit helper script
```

---

## 3. Deep Dive: Core Files Explained for Advisors

### A. `vietnam_dsge.mod` (The Structural DSGE Model)
This is the mathematical heart of the project. It defines the endogenous variables, parameters, and the 16 linearized structural equations that govern the small open economy. 

**Key Code Snippets & Modeling Choices:**
* **Lines 143-144 (Steady State Constraints):** 
  ```matlab
  r_bar = 1/beta - 1; 
  ```
  The steady-state foreign interest rate is tied analytically to the discount factor to close the small open economy and prevent explosive net foreign asset paths.
* **Line 175 (The Reliability Wedge / Agility Gap):** 
  ```matlab
  u = xi * (F - Vol_ren);
  ```
  This is the linearized version of the exponential reliability function. It dictates that utilization ($u$) falls precisely when the required flexibility ($F$) fails to keep pace with renewable volatility ($Vol_ren$).
* **Line 176-177 (Aggregate Flexibility CES):** 
  ```matlab
  F = s_b * (A_bat + K_b) + (1 - s_b) * K_g;
  ```
  Flexibility is a share-weighted composite of private batteries ($K_b$) and public grid transmission ($K_g$). Note that technology ($A_{bat}$) enters additively in logs, augmenting physical battery capital perfectly.
* **Line 183 (Reduced-Form Battery Investment Rule):** 
  ```matlab
  I_bat = -phi_grid * u - P_bat;
  ```
  *Advisor Note:* A fully endogenous, forward-looking optimization of battery capital subject to the non-linear exponential reliability constraint would violate Blanchard-Kahn conditions for local approximation. This rule mathematically proxies the gradient of the shadow price: investment rises when reliability falls (negative $u$ deviation), but is dampened by world battery prices ($P_{bat}$).
* **Line 185 (Endogenous Learning-by-Doing):** 
  ```matlab
  A_bat - A_bat(-1) = -eta_bat * chi * u;
  ```
  Innovation ($A_{bat}$) accelerates when reliability drops. The parameter `chi` ($\chi$) represents the transmission of the scarcity signal. For Vietnam (with regulated tariffs), $\chi$ is close to 0, choking off private learning.

### B. `run_dynare.m` (The Master Execution Script)
This MATLAB script automates the entire analytical pipeline. You do not need to interact with the Dynare command-line; running this script does everything.

**Execution Flow:**
1. **Environment Cleanup:** Clears out old Dynare build folders (`+vietnam_dsge`) to prevent messy conflicts.
2. **Dynare Execution:** Calls `dynare vietnam_dsge.mod` programmatically to solve the rational expectations equilibrium and generate Impulse Response Functions (IRFs).
3. **Data Extraction:** Pulls the IRF data from `oo_.irfs` for the three main shocks:
   * `eps_ren` (Domestic renewable intermittency shock)
   * `eps_bat` (Imported battery price shock)
   * `eps_I` (Government grid construction delay)
4. **Counterfactual & Qualitative Analysis:** The script generates customized, publication-ready MATLAB plots that Dynare cannot produce natively. 
   * **`agility_gap_dynare.png`**: Simulates the "Reliability Valley" over a 25-year transition, proving the timing mismatch between fast solar growth and slow grid expansion.
   * **`counterfactual_chi.png`**: Solves the model iteratively for different values of $\chi$, proving that liberalizing Vietnam's power market (moving to capacity pricing) significantly reduces the welfare costs of intermittency.
   * **`joint_shock_perfect_storm.png`**: Evaluates the "Perfect Storm" scenario where an ENSO drought coincides with a global battery price spike.
5. **Console Reporting:** Prints nicely formatted text to the MATLAB console detailing the precise Welfare Cost (in Consumption Equivalent Variation) and the Variance Decomposition.

---

## 4. How to Reproduce the Thesis Results

To rebuild the data, regenerate all figures, and view the analysis:

1. **Prerequisites:**
   * Install MATLAB (R2023a or newer recommended).
   * Install [Dynare](https://www.dynare.org/) (version 5.x or 6.x).
   * Ensure Dynare is added to your MATLAB path. (e.g., `addpath('c:\dynare\6.1\matlab')`).
2. **Run Execution Script:**
   * Open `run_dynare.m` in MATLAB.
   * Click **Run** (or type `run_dynare` in the command window).
   * The script takes ~5-10 seconds. It will populate the repository with the updated `.png` figures and print the core quantitative findings to the console.
3. **Compile the Manuscript:**
   * If you wish to rebuild the PDF, run `xelatex main.tex` (or use your preferred LaTeX editor like Overleaf, TeXworks, or VS Code LaTeX Workshop).

---

## 5. Version Control ("Save to GitHub")

For the author: If you make changes to the `.mod` or `.tex` files, simply double-click the `save_to_github.bat` script in the Windows directory. It will automatically stage all changes, commit them to the `main` branch, and push them to this exact repository seamlessly.
