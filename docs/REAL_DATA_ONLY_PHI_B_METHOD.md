# Real-Data-Only Calibration of Debt-Elastic Premium (phi_b)

This note documents a strict real-data-only procedure for calibrating `phi_b`.

## Update: observed debt-stock method added

In addition to the CA-recursion debt proxy method, a second strict real-data-only method is now implemented using observed external debt stock data from the World Bank.

## Data used (only observed series)

1. `vietnam_ca_2016_2025.csv` (annual CA/GDP)
2. `vietnam_real_interest_rate_1993_2023.csv` (annual real interest rate, %)

Overlap window used by the script: 2016-2023 (8 annual observations).

## Estimation equation

1. Compute baseline real rate over overlap:
   `r_bar = mean(r_t)`
2. Compute observed spread:
   `spread_t = r_t - r_bar`
3. Construct debt proxy from CA identity (debt-positive convention):
   `d_t = (1 + r_{t-1}) d_{t-1} - ca_t`, with normalized `d_0 = 0`
4. Estimate:
   `spread_t = phi_b * (exp(d_t - d_bar) - 1)`

`phi_b` is estimated by constrained OLS:

`phi_b = max(0, sum(x_t*y_t)/sum(x_t^2))`, where `x_t = exp(d_t-d_bar)-1` and `y_t=spread_t`.

## Output

Results are written to:

- `phi_b_realdata_only.json`
- `phi_b_realdata_debtstock.json`
- `vietnam_wb_indicators.csv`
- `vietnam_wb_indicators_meta.json`

Current run result (from your data):

- unconstrained estimate: `-0.1924`
- constrained estimate: `0.0`

Observed debt-stock method result (preferred strict real-data-only anchor):

- overlap sample: 1993-2023 with non-missing real-rate years (`n=28`)
- debt proxy replaced by observed external debt stock ratio (`DT.DOD.DECT.GN.ZS`)
- estimated `phi_b`: `0.03866`
- fit RMSE (spread equation): `0.05664`

## Interpretation

A constrained estimate of zero means this specific observed sample does not support a positive debt-elastic premium under this reduced-form mapping.

This does NOT imply the model should set `phi_b = 0` mechanically. It implies:

1. The 2016-2023 overlap is short.
2. Debt is proxied (not directly observed NFA market valuation series).
3. Real-rate movements reflect many factors beyond debt premium.

Use this as a strict data anchor/check, then combine with structural model discipline and robustness analysis.

With observed debt stock included, the strict real-data-only estimate is positive and economically close to the profile-based baseline range used in the DSGE exercises.
