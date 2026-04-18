# Deep Calibration Review: Vietnam DSGE Model
**Date:** April 4, 2026
**Objective:** Identify and implement opportunities for better parameter calibration
**Methodology:** Data-driven discipline similar to φ_b calibration approach

---

## Executive Summary

This review examines all 30+ model parameters against available Vietnamese macroeconomic and sectoral data. Following the successful φ_b calibration (which used real data to discipline an order-of-magnitude estimate), we identify **five parameters that can be improved** with data-driven approaches:

1. **β (discount factor)** — Currently narrative (4% from VGB yields); can use actual historical real-rate series
2. **α (capital share)** — Currently narrative (inverse of labor share); can be anchored to national accounts labor compensation
3. **σ_E (energy-value added elasticity)** — Currently literature-based; can be estimated from Vietnam energy-GDP covariance
4. **ρ_ren (intermittency persistence)** — Currently calibrated from IEA; can be validated with EVN monthly generation data
5. **φ_grid (grid investment response)** — Currently narrative (EVN 2023 shock response); can be strengthened with multi-year quasi-experimental data

Additionally, we flag **three parameters for sensitivity analysis** that are currently fixed:

6. **μ (battery weight in flexibility)** — Critical lever on steady-state K_b/K_g ratio; sensitivity around 0.12–0.20 warranted
7. **σ_E (substitution elasticity)** — Controls energy intensity of output; ±0.2 swing materially affects reliability wedge
8. **ρ_flex (battery-grid complementarity)** — Currently 0.40; range [0.30, 0.60] includes both "weak" and "strong" complementarity regimes

---

## Part I: Parameters That CAN Be Better Calibrated

### 1. Discount Factor β → Real Interest Rate Series

**Current calibration:**
- β = 0.99 (quarterly)
- Narrative: "4% annualized real return on Vietnamese Government Bonds"
- Source: GSO 2024 (stated without specific year or bond tenor)

**Data available:**
- Vietnam real interest rate series (1993–2023) in workspace: `vietnam_real_interest_rate_1993_2023.csv`
- 28 annual observations (5 missing values)
- Recent window (2013–2023): mean 5.24% annual, implies β = 0.9502
- Window (2016–2023): mean 4.72% annual, implies β = 0.9549

**Why this matters:**
The discount factor directly governs:
- Capital accumulation dynamics: K_p,ss = (1-τ)αβ/(1-β(1-δ_p))
- Steady-state consumption share: higher β → lower impatience → higher C/Y
- Welfare costs: β weights future losses in compensating variation
- External borrowing premium: r* = 1/β - 1, so β change directly shifts r_bar

**Recommended improvement:**

Choose β based on a defensible window rather than narrative. Options:

- **2000–2023 (full modern era):** 3.14% annual → β = 0.9691. Includes post-WTO accession.
- **2013–2023 (recent decade):** 5.24% annual → β = 0.9502. Reflects current monetary regime.
- **2016–2023 (post-inflation):** 4.72% annual → β = 0.9549. Excludes 2010–2015 volatility.
- **Narrative (VGB target):** 4.00% annual → β = 0.9604. Status quo, defensible for fiscal planning.

**Implementation:**
- Select one window and report sensitivity ±50 bps
- Calculate steady state for three β values to show impact on K_p, C, and welfare cost
- Results are robust if sensitivity shows <5% changes in key IRFs

**Estimated effort:** 30 minutes

---

### 2. Capital Share α → Labor Compensation Data

**Current calibration:**
- α = 0.35 (capital share)
- Narrative: "Inverse of 0.65 labor share; GSO national accounts"
- Source: Cited but not specific (which year? which series?)

**Why this matters:**
Capital share α directly controls:
- Steady-state capital-output ratio: K_p,ss/Y = (1-τ)αβ/(1-β(1-δ_p))
- Output elasticity to capital shocks
- Calibration consistency with long-run growth accounting

**Recommended improvement:**

Retrieve GSO labor compensation share from 2015–2024 national accounts. Expected range: 45–50% labor share → α ∈ [0.50, 0.55].

**Steps:**
1. GSO direct download of compensation of employees / GDP (2015–2024)
2. Penn World Table check: if available, verify with `labsh` for Vietnam
3. WIOD verification: Input-Output tables for labor-capital shares

**If α = 0.35 is correct:**
- Document with specific GSO table reference and year(s)
- State: "Labor compensation share in [sector/economy-wide] was [X]% in [year], yielding α = 0.35"

**If α ≠ 0.35:**
- Recalibrate steady state
- Sensitivity: K_p,ss proportional to α, so 10% change in α → ~10% change in K_p,ss

**Estimated effort:** 1–2 hours

---

### 3. Energy-Value Added Elasticity σ_E → Data-Driven Estimate

**Current calibration:**
- σ_E = 0.6 (elasticity of substitution)
- Source: Heutel (2012) US range 0.4–0.7; "consistent with IEA (2023) Southeast Asia estimates"
- **Currently literature-based, not data-grounded for Vietnam**

**Why this matters:**
- Controls energy-capital complementarity
- Directly affects reliability function
- Higher σ_E → easier substitution → weaker reliability constraint

**Recommended improvement:**

Estimate σ_E from covariance of log(energy services) and log(value added) using 2015–2023 data.

**Specification:**
```
log(E_t) = α + β log(Y_t) + ε_t
σ_E = 1 / (1 - β)
```

If β ≈ 0.40–0.50, then σ_E ≈ [1.67, 2.0], suggesting energy is a complement (not σ_E < 1).

**Caveat:** High β may reflect rapid electrification rather than true substitution elasticity. If so, use literature benchmark and flag as policy-dependent.

**Estimated effort:** 2–3 hours

---

### 4. Intermittency Persistence ρ_ren → EVN Generation Data

**Current calibration:**
- ρ_ren = 0.85 (quarterly)
- Source: "IEA quarterly AR(1) matches seasonal autocorrelation"
- Validated against monsoon cycle

**Why this matters:**
- Directly controls impulse-response duration
- Higher ρ_ren → longer-lived shocks → larger welfare costs

**Recommended improvement:**

Retrieve EVN monthly solar + wind generation data (2015–2023, ~100 observations) and estimate AR(1) on deseasonalized series.

**Steps:**
1. Collect monthly renewable generation from EVN public reports
2. Deseasonalize using seasonal dummies
3. Regress log(residual)_t on log(residual)_{t-3} to estimate quarterly ρ_ren
4. Compare to current calibration

**If ρ_ren ≈ 0.85:**
- Validates current choice
- Report with source: "EVN 2015–2023 data, point estimate ρ_ren = X with 90% CI [X, X]"

**If ρ_ren ≠ 0.85:**
- Recalibrate and run sensitivity on IRFs
- Flag if discrepancy due to increasing solar penetration

**Estimated effort:** 3–4 hours

---

### 5. Grid Investment Response φ_grid → Quasi-Experimental Data

**Current calibration:**
- φ_grid = 1.5 (semi-elasticity)
- Source: "June 2023 EVN power shortage: reliability dropped 1.5pp, government increased 2024 investment by 12–15%"
- Interpretation: φ_grid = 12% / 1.5% ≈ 8; conservatively set to 1.5 for implementation lag

**Why this matters:**
- Controls responsiveness of public investment to reliability shocks
- Directly affects Reliability Valley depth
- Determines welfare cost of implementation delays

**Recommended improvement:**

Compile multi-year data on EVN reliability vs. actual investment, 2015–2023:

- Retrieve EVN annual reports (2015–2023) for reliability metrics (SAIFI, SAIDI)
- Retrieve Ministry of Finance energy investment budgets (2015–2023)
- Estimate semi-elasticity: φ_grid = (ΔI_grid / I_baseline) / (Δu)
- Compare to current value 1.5 and check variance

**Expected result:** Multi-year average φ_grid likely 2–4 (more responsive), but high variance due to budget execution lags.

**Estimated effort:** 4–5 hours

---

## Part II: Parameters for Sensitivity Analysis

These parameters are currently **not varied** in sensitivity tables but materially affect steady state or dynamics.

### 6. Battery Weight in Flexibility μ = 0.16

**Why variable:**
- K_b/K_g ratio directly proportional to μ / (1-μ)
- Currently μ = 0.16 → K_b/K_g ≈ 19%
- Small changes materially affect reliability wedge

**Sensitivity recommendation:**
- Add 3-point sensitivity: μ ∈ {0.12, 0.16, 0.20}
- Re-run steady state for each to show K_b, K_g, F, u_ss
- Higher μ → larger K_b → tighter reliability wedge → lower welfare cost

**Estimated effort:** 1 hour

---

### 7. Energy-Value Added Substitution σ_E = 0.6

**Why variable:**
- Currently fixed; described as "gross complements"
- Literature range: 0.4–0.7
- σ_E < 1 → energy is a bottleneck → larger reliability wedge

**Sensitivity recommendation:**
- Add 3-point sensitivity: σ_E ∈ {0.4, 0.6, 0.8}
- Recalculate steady state and reliability elasticity for each
- Show output response to energy shortage under three complementarity regimes

**Estimated impact:**
- σ_E = 0.4: stronger energy-output co-movement → larger reliability shock
- σ_E = 0.8: weaker complementarity → smaller welfare cost
- ±0.2 range implies ±15–20% variation in reliability wedge

**Estimated effort:** 1 hour

---

### 8. Battery-Grid Substitution Elasticity ρ_flex = 0.40

**Why variable:**
- Controls complementarity strength
- Hart (2019): 0.30 (Leontief); others: 0.50–0.60 (loose)
- Currently fixed at 0.40 with no sensitivity

**Sensitivity recommendation:**
- Add 3-point sensitivity: ρ_flex ∈ {0.30, 0.40, 0.60}
- Run IRF to intermittency shock for each regime
- Show I_bat and I_grid responses: weak vs. strong complementarity changes agility gap

**Estimated impact:**
- ρ_flex = 0.30: I_bat and I_grid move proportionally → agility gap shrinks
- ρ_flex = 0.60: I_bat can increase without I_grid → larger agility gap
- Likely ±10–15% variance in welfare cost

**Estimated effort:** 1 hour

---

## Part III: Well-Calibrated Parameters (Do Not Change)

These parameters have strong data foundations and should remain fixed:

| Parameter | Value | Data anchor | Confidence |
|-----------|-------|-------------|-----------|
| β | 0.99 | Vietnam real rate (3.14% mean) | High if window documented |
| α | 0.35 | GSO labor share ≈ 65% | High if GSO cited |
| δ_p | 0.025 | Standard 10% annual | Very high |
| δ_g | 0.0125 | Vietnam transmission 5% annual | High |
| δ_b | 0.030 | Li-ion cells 12% annual | Very high |
| σ_L | 1.0 | Standard NK | Moderate |
| ρ_ren | 0.85 | IEA/monsoon data | Moderate-High |
| σ_ren | 0.12 | GSO solar CV | High |
| θ_ren | 0.30 | IEA/PDP8 | Moderate-High |
| φ_b | 0.039 | Vietnam real-rate regression | Moderate (R²=0.16) |

---

## Action Plan

**Immediate (Before Final Submission):** 1 hour
- Confirm β window (2000–2023 or 2016–2023)
- Cite specific GSO table for α

**Medium-term:** 12–15 hours
- σ_E data estimation (2–3 hours)
- ρ_ren validation (3–4 hours)
- φ_grid multi-year analysis (4–5 hours)
- Sensitivity tables (μ, σ_E, ρ_flex): 3 hours

**Total estimated effort:** 13–16 hours to implement all recommendations.

**Prioritization:** β, α documentation (immediate); σ_E, ρ_ren, φ_grid (next version).

