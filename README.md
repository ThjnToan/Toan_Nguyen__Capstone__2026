# Renewable Intermittency and Grid Investment Dynamics

**Scarcity-Induced Battery Innovation in Vietnam's Energy Transition**

**Author:** Toan T. Nguyen (Justin) · **Advisor:** Dr. Xavier Martin G. Bautista
**Institution:** Fulbright University Vietnam · **Date:** Spring 2025

---

## Overview

This repository contains the MATLAB/Dynare codebase and LaTeX manuscript for a capstone thesis investigating the macroeconomic challenges of Vietnam's clean energy transition. The project develops a **Small Open Economy DSGE model** with endogenous scarcity-induced technological change in battery storage and a grid reliability constraint, calibrated to Vietnam's **Power Development Plan 8 (PDP8)**.

### Key Mechanisms

| Mechanism | Description |
|---|---|
| **Reliability Penalty** | Renewable intermittency reduces capital utilization through grid instability — an endogenous TFP wedge |
| **Agility Gap** | Private battery investment responds quickly but is dominated by global price shocks; public grid investment accumulates slowly through fiscal rules |
| **Scarcity-Induced Innovation** | Battery technology improves endogenously when the grid reliability gap signals scarcity |
| **Dual External Exposure** | Simultaneous vulnerability to global battery supply chains and domestic fiscal constraints |

### Main Findings

- Intermittency shocks account for **95.4% of output volatility**
- Global battery-price shocks drive **100% of battery investment volatility**
- Suppressing scarcity signals raises the welfare cost of intermittency by **20%**
- Battery and grid capital are **imperfect complements** — neither alone solves the reliability problem

---

## Repository Structure

```
├── fulbright_thesis.tex              # Main thesis (compiles thesis.pdf)
├── fulbrightthesis.cls               # Thesis document class
├── fuvmacro.sty                      # Macro/style package
├── main_body.tex                     # Thesis body (introduction through appendices)
├── main.tex                          # Standalone article variant
│
├── vietnam_dsge.mod                  # Baseline Dynare model
├── vietnam_dsge_steadystate.m        # Steady-state solver (baseline)
├── run_dynare.m                      # Script to run baseline model
├── vietnam_dsge_proptax.mod          # Proportional-tax variant
├── vietnam_dsge_proptax_steadystate.m# Steady-state solver (proptax)
├── run_dynare_proptax.m              # Script to run proptax model
│
├── vietnam_dsge_wickens.mod          # Wickens (2008) extension variant
├── vietnam_dsge_wickens_steadystate.m
│
├── test_ss.m                         # Steady-state validation script
├── references.bib                    # Bibliography
├── figures/                          # Generated figures
├── data/                             # Calibration data
├── scripts/                          # Additional MATLAB scripts
├── docs/                             # Documentation
│
├── vietnam_dsge/                     # Dynare output (baseline)
├── vietnam_dsge_proptax/             # Dynare output (proptax)
│
├── *.png                             # Generated IRF/valley figures
└── *.docx                            # Advisory memos
```

---

## Compiling the Thesis

```bash
pdflatex fulbright_thesis.tex
bibtex fulbright_thesis
pdflatex fulbright_thesis.tex
pdflatex fulbright_thesis.tex
```

## Running the Model

1. Open MATLAB with Dynare 6.x on the path
2. Run `run_dynare.m` (baseline) or `run_dynare_proptax.m` (proportional tax variant)
3. Generated figures are saved as PNGs in the root directory

---

## Citation

```bibtex
@thesis{nguyen2025renewable,
  author  = {Toan T. Nguyen},
  title   = {Renewable Intermittency and Grid Investment Dynamics:
             The Role of Scarcity-Induced Battery Innovation in Vietnam},
  school  = {Fulbright University Vietnam},
  year    = {2025}
}
```
