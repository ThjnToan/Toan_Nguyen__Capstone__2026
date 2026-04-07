# Future Validation and Robustness Analysis Roadmap

**Date:** April 4, 2026
**Status:** Intermediate (immediate priorities completed; medium-term roadmap established)

---

## Overview

This document outlines the strategic roadmap for strengthening the Vietnam DSGE model's calibration and robustness analysis beyond the immediate priorities already completed. The roadmap is organized into three phases:

1. **Phase 1 (Completed):** Data-anchored immediate priorities (β, α documentation) — ~1 hour
2. **Phase 2 (Recommended Next):** Medium-term parameter validation (3-4 months to full implementation)
3. **Phase 3 (Optional):** Long-term extensions and advanced sensitivity analysis

---

## Phase 1: ✅ Completed (April 4, 2026)

### 1.1 Discount Factor β = 0.99

**Task:** Anchor to Vietnamese historical real interest rates
**Status:** ✅ DONE

- Data source: World Bank real interest rate series (1993-2023)
- Window selected: Post-WTO modern era (2000-2023)
- Sample statistics: Mean 3.14% annual; implied β = 0.9691
- Documentation: Main.tex lines 843, 866, 928-930
- Robustness note: Sensitivity to β ∈ [0.95, 0.99] mentioned but not yet computed

### 1.2 Capital Share α = 0.35

**Task:** Document source of labor share data from GSO
**Status:** ✅ DONE

- Data source: GSO national accounts (growth accounting)
- Measurement: Labor compensation share ≈ 0.65 across Vietnamese economy
- Window: 2015-2024
- Documentation: Main.tex lines 867, 932
- Method: Standard inverse (capital share = 1 - labor share)

---

## Phase 2: Medium-Term Parameter Validation (Estimated 12-15 hours total)

### 2.1 Renewable Intermittency Persistence ρ_ren = 0.85

**Priority:** HIGH (foundational to reliability mechanism)
**Effort:** 3-4 hours
**Current status:** Data-referenced but not computed from raw data

#### Immediate Actions (Next 1-2 weeks):

1. **Acquire EVN generation data**
   - Source: Vietnam Electricity (EVN) National Load Dispatch Centre (NLDC) historical generation records
   - Frequency: Monthly utility-scale solar and wind generation (2015-2024 recommended)
   - Alternative: IEA Vietnam Energy Outlook Annex tables (publicly available)
   - Data cleaning: Deseasonalize using 12-month moving average or seasonal decomposition

2. **Compute sample autocorrelation**
   ```python
   # Pseudocode workflow
   monthly_capacity_factor = (monthly_generation / installed_capacity)
   deseasonalized = monthly_capacity_factor - seasonal_component
   quarterly_ar1 = scipy.stats.pearsonr(deseasonalized[:-1], deseasonalized[1:])
   # Should yield ρ_ren ≈ 0.85 at quarterly frequency
   ```

3. **Document conversion from monthly to quarterly**
   - If IEA data is monthly, compute quarterly averages before AR(1) estimation
   - Report exact sample autocorrelation with 95% confidence interval
   - Compare to model value ρ_ren = 0.85

4. **Add to main.tex**
   - Section "Renewable Intermittency Validation" (new subsection after calibration table)
   - Report: Sample autocorrelation coefficient, data source, time period, confidence bounds
   - Sensitivity table: Impact effects for ρ_ren ∈ {0.70, 0.85, 0.95} on output and investment volatility

#### Key Questions:
- Does actual deseasonalized autocorrelation support ρ_ren = 0.85, or should the parameter be adjusted?
- Are there time-varying patterns (2015-2020 vs 2020-2024) suggesting changing intermittency profiles?
- How sensitive are key welfare results to ±0.15 swing in ρ_ren?

---

### 2.2 Energy Elasticity of Substitution σ_E = 0.6

**Priority:** MEDIUM (affects energy-capital complementarity; already has literature support)
**Effort:** 2-3 hours
**Current status:** Cited from Heutel (2012) and IEA; not estimated from Vietnamese data

#### Immediate Actions:

1. **Assemble Vietnam energy-GDP data**
   - GSO: Annual energy consumption (MWh) from 1990-2024
   - GSO: Annual value-added (constant 2015 VND) from 1990-2024
   - WB: Energy consumption (kg oil equivalent) as robustness check
   - Deflate nominal series; compute real energy-value-added ratio

2. **Estimate energy elasticity**
   ```python
   # Log-linear regression: ln(E/VA) = a + σ_E * ln(P_E/P_VA) + ε
   # Two approaches:

   # Approach 1 (Short-run): Use year-over-year changes
   # Δln(E/VA) on Δln(P_E/P_VA)

   # Approach 2 (Long-run): Detrend and compute trend elasticity
   # VAR/cointegration if using 1990-2024 non-stationary series
   ```

3. **Compare to Heutel (2012) range**
   - Heutel US estimate: σ_E ∈ [0.4, 0.7]
   - IEA Southeast Asia range: σ_E ≈ 0.5-0.7
   - Model choice: σ_E = 0.6 (midpoint, gross complements σ_E < 1)

4. **Add to main.tex**
   - Subsection "Energy-Value Added Elasticity Validation"
   - Report: Estimated elasticity with 95% CI, data source, time period
   - Sensitivity table: Impulse responses for σ_E ∈ {0.4, 0.6, 0.8}

#### Key Questions:
- Does Vietnamese energy-GDP data support σ_E = 0.6, or is Vietnam more complementary (σ_E < 0.6)?
- Are there structural breaks (e.g., 2008 crisis, industrialization phases)?
- How sensitive are battery/grid investment ratios to σ_E variations?

---

### 2.3 Grid Investment Response φ_grid = 1.5

**Priority:** MEDIUM (fiscal policy channel; already calibrated to 2023 incident)
**Effort:** 4-5 hours
**Current status:** Calibrated to single 2023 EVN shortage event; lacks multi-year evidence

#### Immediate Actions:

1. **Assemble EVN fiscal response data**
   - Ministry of Finance supplementary budgets (2015-2024)
   - EVN capital investment disbursements (MoF budget execution reports)
   - Grid reliability metrics (EVN operational reports, quarterly)
   - Build event history: Grid reliability shortfalls and lagged investment responses

2. **Conduct quasi-experimental analysis**
   ```python
   # Event study approach: Identify grid reliability dips (Δu < -0.005)
   # Measure investment response 1-4 quarters ahead
   # Regression: I_{grid,t+k} = α + β(Δu_t) + ε
   # Estimate impulse response across multiple events, not just 2023
   ```

3. **Validate reactive vs anticipatory behavior**
   - Test for anticipation: Does investment rise *before* reliability drops?
   - If yes, revise φ_grid to capture forward-looking behavior (may be higher)
   - If no, current reactive calibration is correct

4. **Account for implementation delays**
   - Distinguish between budget *allocation* (immediate) and *disbursement* (2-3 year lag)
   - Model φ_grid captures allocation response; accumulation follows investment rule

5. **Add to main.tex**
   - Subsection "Government Investment Response Validation"
   - Report: Multi-year event study results, coefficient estimates, event windows
   - Sensitivity table: Grid-investment elasticity for φ_grid ∈ {1.0, 1.5, 2.0}

#### Key Questions:
- Is φ_grid = 1.5 representative, or does it cherry-pick the unusual 2023 response?
- Have there been other reliability shortfalls (2010-2012) with measurable investment responses?
- How much of the 2023 response was anticipatory (pre-shock) vs reactive (post-shock)?

---

## Phase 3: Robustness and Sensitivity Analysis (Estimated 3-5 hours)

### 3.1 Multi-Parameter Sensitivity Tables

**Priority:** MEDIUM (documentation of parameter uncertainty)
**Effort:** 1-2 hours

Add four new tables to the sensitivity analysis appendix:

1. **Battery Weight μ ∈ {0.12, 0.16, 0.20}**
   - Impact on I_bat vs I_grid investment ratios
   - Impact on welfare costs across χ values
   - Rationale: PDP8 trajectory shows battery share rising from ~14% (current) to ~64% (2027-2030)

2. **Energy Elasticity σ_E ∈ {0.4, 0.6, 0.8}**
   - Gross complements (σ_E = 0.4), baseline (0.6), substitutes (0.8)
   - Impact on energy-GDP response to reliability shocks
   - Impact on K_p accumulation (energy-complementary vs substitutable regimes)

3. **Intermittency Persistence ρ_ren ∈ {0.70, 0.85, 0.95}**
   - Impact on output and utilization volatility
   - Impact on battery adoption (learning-by-doing channel)
   - Interpretation: seasonal vs trend-driven intermittency

4. **Battery Depreciation δ_b ∈ {0.02, 0.03, 0.04}**
   - Currently 0.03 (12% annual)
   - Robustness: accounts for uncertainty in battery cycle life (8-10 vs 12+ years)
   - Impact on K_b/K_g ratio and steady-state capital structure

### 3.2 Welfare Robustness Across Multiple Dimensions

**Priority:** LOW (nice-to-have for completeness)
**Effort:** 2-3 hours

Create a welfare sensitivity matrix:
- Rows: χ ∈ {0.0, 0.3, 1.0} (regulatory transmission)
- Columns: ρ_ren ∈ {0.70, 0.85, 0.95} (intermittency persistence)
- Values: Welfare cost Λ (%)
- Interpretation: All qualitative conclusions robust if welfare rankings unchanged

---

## Implementation Timeline

### Immediate (Before final submission) — ~1 week
- [x] Complete Phase 1 (β, α calibration anchoring)
- [ ] Add ρ_ren validation note to main.tex (30 min)
- [ ] Recompile PDF with validation roadmap embedded

### Medium-term (Next 4-8 weeks) — if pursuing journal submission
- [ ] Execute Phase 2.1: Acquire EVN data, compute ρ_ren empirically (3-4 hours)
- [ ] Execute Phase 2.2: Estimate σ_E from Vietnam energy-GDP data (2-3 hours)
- [ ] Execute Phase 2.3: Conduct multi-year φ_grid event study (4-5 hours)
- [ ] Create Phase 3.1 sensitivity tables (1-2 hours)
- [ ] Update main.tex with validation subsections and new tables

### Long-term (Post-submission feedback cycle)
- [ ] Phase 3.2 welfare robustness matrix (if referees request)
- [ ] Deep dive on specific parameter if sensitivity analysis flags high elasticity
- [ ] Consider model extensions (e.g., endogenous reliability targets, non-linear depreciation)

---

## Data Sources and Access Notes

### Required for Phase 2 Implementation

| Parameter | Data Source | Access Level | Estimated Time |
|-----------|------------|--------------|-----------------|
| ρ_ren     | EVN NLDC or IEA Vietnam Energy Outlook | Public | 2-3 hours to clean |
| σ_E       | GSO annual energy & VA; WB energy data | Public | 1-2 hours |
| φ_grid    | MoF budget execution reports; EVN reports | Public (with registration) | 3-4 hours |

### Existing Calibration Data (Already Assembled)

- vietnam_wb_indicators.csv (World Bank 1985-2024)
- vietnam_real_interest_rate_1993_2023.csv (real rates + β calibration)
- PDP8 investment tables (battery & grid planned investment)

---

## Quality Assurance Checkpoints

After each Phase 2 parameter is validated, verify:

1. **Does data-implied value support model choice?**
   - If point estimate differs by >20% from model value, explain discrepancy
   - If confidence interval excludes model value, consider re-calibration

2. **Do sensitivity results change qualitative conclusions?**
   - FEVD rankings should remain robust (intermittency > battery price > implementation)
   - Welfare rankings under different χ should remain invariant

3. **Are results policy-relevant?**
   - Does parameter uncertainty material change policy recommendations?
   - If yes, highlight in updated Conclusions

---

## Success Criteria

Phase 2 implementation will be considered complete when:

- ✅ All five parameters (β, α, ρ_ren, σ_E, φ_grid) are documented with data sources
- ✅ Sample statistics (means, CIs, autocorrelations) are reported in main text
- ✅ Sensitivity tables show qualitative robustness across plausible parameter ranges
- ✅ Discrepancies between data-implied and model values are explained
- ✅ Appendix contains complete parameter validation subsections with references to raw data

This roadmap provides a structured path to full data-anchored calibration while keeping the core paper submission-ready.
