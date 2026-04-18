# Calibration Data Inventory for Vietnam DSGE

This file separates calibration inputs into three groups:

1. Directly observable public data
2. Observable with additional collection work
3. Structural parameters that are not directly observable and must be calibrated or benchmarked

## A. Directly Observable Public Data Collected Now

| Model object | Variable / target | Vietnam value found | Source | Notes |
|---|---|---:|---|---|
| Current account balance to GDP | CA/GDP | 2024: 5.9% | World Bank indicator `BN.CAB.XOKA.GD.ZS` | Public WB page states Viet Nam 2024 = 5.9 |
| Household consumption share | C/Y proxy | 2024: 54.0% of GDP | World Bank indicator `NE.CON.PRVT.ZS` | Household + NPISH final consumption expenditure |
| Government consumption share | G/Y proxy | 2024: 8.9% of GDP | World Bank indicator `NE.CON.GOVT.ZS` | Useful to reconcile model consumption residual |
| GDP level | Y (nominal anchor) | 2024: USD 476.388 bn | World Bank indicator `NY.GDP.MKTP.CD` | Page reports 476,388,230.31 (thousand US$) |
| CPI inflation | Inflation / ex-post real-rate adjustment | 2024: 3.6% | World Bank indicator `FP.CPI.TOTL.ZG` | Can help convert nominal bond yield to real |
| Renewable electricity share | Renewable output share of total electricity | 2021: 41.17% | World Bank indicator `EG.ELC.RNEW.ZS` (IEA source) | Useful external validation for energy block |
| Renewable electricity excluding hydro | Non-hydro renewable output | 2021: 25,778,100 thousand kWh | World Bank indicator `EG.ELC.RNWX.KH` (IEA source) | Equals about 25.78 TWh |
| Electricity access | Power system coverage proxy | 2023: 99.8% | World Bank indicator `EG.ELC.ACCS.ZS` | Not reliability, but useful sector context |

## B. Data Already Entered in Workspace

| Model object | Variable / target | Workspace source | Notes |
|---|---|---|---|
| Annual CA/GDP series | CA/GDP, 2016-2025 | `vietnam_ca_2016_2025.csv` | User-provided series used for phi_b calibration |
| Annual real interest rate series | Real rate (%), 1993-2023 | `vietnam_real_interest_rate_1993_2023.csv` | User-provided historical series with missing years left blank |
| Real-rate calibration summary | Implied annual/quarterly real rate and beta by window | `vietnam_real_rate_summary.json` | Computed from user-provided real-rate series |
| CA target moments | mean/std/ac1 of CA/GDP | `vietnam_ca_targets.json` | Built from user-provided annual series |
| SGU-style pure CA estimate | phi_b (conditional) | `smm_phi_b_sgu_ca_only.json` | Annual CA-only mode |
| Exact 1D phi solve | phi_b (conditional, high precision) | `phi_b_exact_vietnam.json` | Uses fixed auxiliary parameters |
| Joint profile estimate | phi_b band | `phi_b_profile_joint_vietnam.json` | Best phi about 0.0454, 10% band about 0.0356-0.0603 |
| Joint bootstrap check | phi_b uncertainty (fast run) | `smm_phi_b_vietnam_ca_override_bootstrap_mini.json` | Quick CI run, not final high-accuracy run |

## C. Observable with Additional Collection Work

These are measurable, but not yet collected into the workspace.

| Parameter / target | What to collect | Best source candidates | Why needed |
|---|---|---|---|
| `beta`, `r_bar` | 10Y Vietnam government bond yield series, monthly or quarterly | AsianBondsOnline, HNX, MOF, IMF IFS | To replace narrative 4% real-rate anchor with actual series |
| `alpha` | Labor compensation share or labor income share | Penn World Table (`labsh`), GSO national accounts, WIOD | Capital share can be set as `1 - labor share` |
| Investment share | Gross capital formation / GDP series | World Bank `NE.GDI.TOTL.ZS`, GSO national accounts | Current paper cites GSO; WB fetch page timed out during retrieval |
| `theta_ren` | Renewable share of electricity generation/capacity | IEA, Ember, EVN, PDP8 | Needed for renewable penetration anchor |
| `sigma_ren`, `rho_ren` | Monthly or quarterly solar/wind output volatility and persistence | IEA monthly electricity, Ember, EVN dispatch/generation data | Needed for intermittency shock process |
| `mu` | Battery share in total battery+grid investment | PDP8 planned investment tables, MOIT, EVN | Needed for flexibility CES weighting |
| `phi_grid`, `rho_I`, `sigma_I` | Annual or quarterly grid-investment disbursement and project delay data | Ministry of Finance execution reports, World Bank project documents, EVN annual reports | Needed for fiscal implementation shock block |
| `u_target`, reliability validation | SAIDI, SAIFI, ENS, outage rate, transmission reliability | EVN annual reports, ERAV, MOIT | Needed to discipline reliability wedge |
| `phi_b` market anchor | Sovereign spread / CDS / EMBI / external borrowing spread | Bloomberg, JP Morgan EMBI, IMF, World Bank IDS, central bank external debt bulletins | Needed to give market-based external-finance discipline |

## D. Parameters Not Directly Observable

These cannot be "found" as raw data series. They require calibration, estimation, or external literature priors.

| Parameter | Why it is not directly observable | Recommended discipline |
|---|---|---|
| `phi_b` | Elasticity of borrowing premium to debt position is a model object, not a published national statistic | Estimate from moments + market spread evidence + sensitivity band |
| `psi` | Reliability elasticity of output/utilization to flexibility ratio is a structural mapping | Calibrate to outage/reliability targets and engineering ranges |
| `rho_flex` | Substitution/complementarity between batteries and grid is not directly reported by statistical agencies | Literature benchmark + sensitivity analysis |
| `eta_bat` | Innovation elasticity from reliability gaps to battery productivity is structural | Literature on induced innovation / learning curves |
| `chi` | Regulatory transmission of scarcity signals is institutional and latent | Counterfactual scenario parameter, not data series |
| `delta_g`, `delta_b`, `delta_p` | Economic depreciation is model-based, not a single observed official number | Use engineering lifetimes and accounting evidence |
| `sigma_E` | CES elasticity between energy and value added is a structural production parameter | External literature benchmark |

## E. Practical Recommendation

Use the following hierarchy in the thesis:

1. Replace all possible macro calibration targets with official public series from World Bank / GSO / EVN / PDP8.
2. For sectoral and engineering objects, build a source appendix listing the exact report/table used.
3. For structural parameters, explicitly label them as "calibrated" or "externally benchmarked," not "measured."
4. For `phi_b`, report:
   - SGU benchmark: 0.000742 or 0.001
   - Vietnam data-driven conditional estimate: around 0.009 in pure annual CA mode
   - Vietnam joint external-block estimate/profile region: about 0.045 with band roughly 0.036-0.060

5. For `beta`/`r_bar`, use the user-provided Vietnam real-rate series as a data anchor and report window sensitivity:
   - 2013-2023 mean real rate: 5.24% annual => 1.285% quarterly => implied beta about 0.9502
   - 2016-2023 mean real rate: 4.72% annual => 1.161% quarterly => implied beta about 0.9549

## F. Immediate Missing Items to Collect Next

1. 10Y Vietnam government bond yield history
2. Labor share series for Vietnam
3. Gross capital formation / GDP official annual series
4. PDP8 battery and grid investment tables
5. EVN reliability and outage indicators
6. Monthly renewable generation data for volatility and persistence
7. External borrowing spread / sovereign CDS or EMBI proxy
