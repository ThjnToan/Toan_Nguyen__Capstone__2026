# SMM Setup for Debt-Elastic Premium (phi_b)

This workflow estimates the debt-elastic risk-premium parameter phi_b using Simulated Method of Moments (SMM).

## 1) Files added

- estimate_phi_b_smm.py: End-to-end SMM estimator.
- smm_phi_b_results.json: Output file written by the estimator.

## 2) Data format expected

Provide a CSV with quarterly observations and at least these columns:

- nfa_to_gdp
- ca_to_gdp
- spread

Column aliases are supported in the script (for example: b_star, current_account_to_gdp, interest_rate).

If ca_to_gdp is missing but NFA and spread are present, the script reconstructs CA approximately.

## 3) Baseline run command

python estimate_phi_b_smm.py --data your_data.csv --out smm_phi_b_results.json

Optional useful flags:

- --bootstrap 100
- --maxiter 300
- --n-sim 6000

## 4) Built-in self-test (recommended first run)

Use this to validate the estimator before using real data.

python estimate_phi_b_smm.py --data dummy.csv --self-test --self-test-out smm_self_test_data.csv --out smm_self_test_results.json

The self-test generates synthetic data from known parameters and checks whether estimation recovers them.

## 5) Recommended empirical moments for thesis reporting

Report at minimum:

- mean, std, and ac1 of NFA-to-GDP
- std and ac1 of CA-to-GDP
- mean and std of spread

Then report:

- estimated phi_b
- objective value
- bootstrap confidence interval for phi_b
- sensitivity to alternative moment sets and sample windows

## 6) Interpretation notes

- A larger phi_b implies faster risk-premium increases when external debt rises.
- Check units carefully: if B is in GDP-ratio units, phi_b should be interpreted in that scale.
- If objective value is very high, verify stationarity and data preprocessing (trend breaks, structural regime changes, or inconsistent spread units).

## 7) Using Vietnam CA/GDP moment overrides

If you only have a trusted CA/GDP series (for example 2016-2025 annual data), you can override CA moments in SMM while keeping other moments from the base dataset.

1. Build moment targets from annual CA data:

python prepare_sgu_ca_targets.py --in vietnam_ca_2016_2025.csv --out vietnam_ca_targets.json

2. Re-estimate with CA moment overrides:

python estimate_phi_b_smm.py --data your_data.csv --moments-override vietnam_ca_targets.json --out smm_phi_b_vietnam_ca_override.json

The override file injects `std_ca` and `ac1_ca` directly into the moment vector. This is useful for SGU-style anchoring when CA moments are more credible than NFA/spread moments from the same source.

## 8) Pure SGU-CA mode (annual-frequency consistent)

If you want full alignment with annual CA data, use pure SGU-CA mode. This mode targets only annual CA moments (`std` and `ac1`) computed from non-overlapping 4-quarter blocks in simulation.

Example:

python estimate_phi_b_smm.py --data your_data.csv --moments-override vietnam_ca_targets.json --sgu-ca-only --out smm_phi_b_sgu_ca_only.json

Optional fixed parameters in this mode:

- --rho-ca-fixed (default 0.57)
- --gamma-r-fixed (default 0.80)
- --sigma-ca-fixed (default 0.003)
- --annual-block-size (default 4)

Interpretation: this mode estimates `phi_b` conditional on the fixed non-`phi_b` parameters. It is designed for transparent SGU-style CA-moment anchoring, not for fully joint structural estimation.

## 9) High-precision one-dimensional `phi_b` solve

For a reproducible scalar solve (tight tolerance) targeting annual CA moments:

python calibrate_phi_b_exact.py --targets vietnam_ca_targets.json --base-data your_data.csv --anchor smm_phi_b_results.json --out phi_b_exact_vietnam.json

This computes a high-precision `phi_b` conditional on fixed (`rho_ca`, `gamma_r`, `sigma_ca`) and exports objective profiling around the optimum.

Important caveat: if the objective is very flat in `phi_b`, your moments do not strongly identify `phi_b` under the fixed-parameter configuration. In that case, report an identification range and perform joint calibration over (`phi_b`, `rho_ca`, `gamma_r`, `sigma_ca`) or add spread/NFA moments.

## 10) Joint-identification next steps (recommended)

1. Joint fit with bootstrap CI (fast configuration):

python estimate_phi_b_smm.py --data your_data.csv --moments-override vietnam_ca_targets.json --out smm_phi_b_vietnam_ca_override_bootstrap_mini.json --bootstrap 10 --maxiter 40 --n-sim 1200 --burn-in 200

2. Profile `phi_b` while re-optimizing (`rho_ca`, `gamma_r`, `sigma_ca`):

python profile_phi_b_joint.py --data your_data.csv --moments-override vietnam_ca_targets.json --phi-min 0.001 --phi-max 0.08 --n-grid 17 --n-sim 2500 --burn-in 400 --maxiter 70 --out phi_b_profile_joint_vietnam.json

Interpretation workflow:

- Use `smm_phi_b_vietnam_ca_override_bootstrap_mini.json` for a quick uncertainty interval around the joint estimate.
- Use `phi_b_profile_joint_vietnam.json` to report a profile-based identification band (the file includes a 10% objective-inflation band).
