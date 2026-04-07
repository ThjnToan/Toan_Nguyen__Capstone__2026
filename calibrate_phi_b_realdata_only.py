"""
Calibrate debt-elastic premium phi_b using only real observed data.

Inputs (workspace files):
- vietnam_ca_2016_2025.csv: annual CA/GDP series (decimal)
- vietnam_real_interest_rate_1993_2023.csv: annual real interest rate (%)

Method:
1) Overlap years where both CA and real rate exist.
2) Build observed spread as r_t - r_bar, where r_bar is sample mean over overlap.
3) Construct debt proxy recursively from CA identity (debt-positive convention):
   d_t = (1 + r_{t-1}) d_{t-1} - ca_t
4) Estimate phi_b in spread_t = phi_b * (exp(d_t - d_bar) - 1) by constrained OLS:
   phi_b = max(0, sum(x_t y_t) / sum(x_t^2))

This is a real-data-only reduced-form calibration.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


def load_ca(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            y = int(r["year"])
            v = r["ca_to_gdp"].strip()
            if v:
                out[y] = float(v)
    return out


def load_real_rate(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            y = int(r["year"])
            v = r["real_interest_rate_pct"].strip()
            if v:
                out[y] = float(v) / 100.0
    return out


def main() -> None:
    ca = load_ca(Path("vietnam_ca_2016_2025.csv"))
    rr = load_real_rate(Path("vietnam_real_interest_rate_1993_2023.csv"))

    years = sorted(set(ca.keys()).intersection(rr.keys()))
    if len(years) < 5:
        raise ValueError("Need at least 5 overlap years between CA and real-rate data.")

    # Observed real rates and baseline (sample mean over overlap window)
    r_obs = [rr[y] for y in years]
    r_bar = sum(r_obs) / len(r_obs)
    spread = [rr[y] - r_bar for y in years]

    # Debt proxy recursion with debt-positive convention
    debt = []
    d_prev = 0.0
    for i, y in enumerate(years):
        if i == 0:
            d_t = d_prev - ca[y]
        else:
            r_lag = rr[years[i - 1]]
            d_t = (1.0 + r_lag) * d_prev - ca[y]
        debt.append(d_t)
        d_prev = d_t

    d_bar = sum(debt) / len(debt)
    x = [math.exp(d - d_bar) - 1.0 for d in debt]
    y = spread

    den = sum(v * v for v in x)
    if den <= 1e-14:
        raise ValueError("Insufficient variation in debt proxy transform to identify phi_b.")

    phi_hat_unconstrained = sum(a * b for a, b in zip(x, y)) / den
    phi_hat = max(0.0, phi_hat_unconstrained)

    y_hat = [phi_hat * v for v in x]
    resid = [a - b for a, b in zip(y, y_hat)]
    sse = sum(e * e for e in resid)
    mse = sse / len(resid)
    rmse = math.sqrt(mse)

    out = {
        "mode": "real-data-only",
        "overlap_years": years,
        "n_obs": len(years),
        "inputs": {
            "ca_file": "vietnam_ca_2016_2025.csv",
            "real_rate_file": "vietnam_real_interest_rate_1993_2023.csv",
        },
        "constructed_series": {
            "r_bar_decimal": r_bar,
            "spread_decimal_by_year": {str(y): s for y, s in zip(years, spread)},
            "debt_proxy_by_year": {str(y): d for y, d in zip(years, debt)},
            "debt_proxy_mean": d_bar,
        },
        "estimate": {
            "phi_b_hat": phi_hat,
            "phi_b_unconstrained": phi_hat_unconstrained,
            "rmse_spread_fit": rmse,
        },
        "equation": "spread_t = phi_b * (exp(d_t - d_bar) - 1)",
        "notes": [
            "Strictly uses observed CA and observed real rates.",
            "Debt is a proxy reconstructed from CA identity with normalized initial debt d_0=0.",
            "Small sample; interpret as reduced-form calibration anchor, not structural truth."
        ],
    }

    Path("phi_b_realdata_only.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
