"""
Real-data-only calibration of phi_b using observed external debt stock.

Uses:
- vietnam_wb_indicators.csv (contains ext_debt_to_gni)
- vietnam_real_interest_rate_1993_2023.csv (real rate, %)

Equation estimated:
spread_t = phi_b * (exp(d_t - d_bar) - 1)
where spread_t = r_t - r_bar and d_t is observed external debt ratio proxy.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent


def load_real_rate(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            v = r["real_interest_rate_pct"].strip()
            if v:
                out[int(r["year"])] = float(v) / 100.0
    return out


def load_debt_ratio(path: Path) -> dict[int, float]:
    out: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            v = str(r.get("ext_debt_to_gni", "")).strip()
            if v:
                # convert percent to decimal ratio
                out[int(r["year"])] = float(v) / 100.0
    return out


def main() -> None:
    rr = load_real_rate(BASE_DIR / "vietnam_real_interest_rate_1993_2023.csv")
    debt = load_debt_ratio(BASE_DIR / "vietnam_wb_indicators.csv")

    years = sorted(set(rr.keys()).intersection(debt.keys()))
    if len(years) < 8:
        raise ValueError("Need at least 8 overlap years for debt-stock calibration.")

    r_vals = [rr[y] for y in years]
    r_bar = sum(r_vals) / len(r_vals)
    spread = [rr[y] - r_bar for y in years]

    d_vals = [debt[y] for y in years]
    d_bar = sum(d_vals) / len(d_vals)

    x = [math.exp(d - d_bar) - 1.0 for d in d_vals]
    y = spread

    den = sum(v * v for v in x)
    if den <= 1e-14:
        raise ValueError("Insufficient variation in debt stock to identify phi_b")

    phi_uncon = sum(a * b for a, b in zip(x, y)) / den
    phi_hat = max(0.0, phi_uncon)

    y_hat = [phi_hat * v for v in x]
    resid = [a - b for a, b in zip(y, y_hat)]
    rmse = math.sqrt(sum(e * e for e in resid) / len(resid))

    out = {
        "mode": "real-data-only-observed-debt-stock",
        "overlap_years": years,
        "n_obs": len(years),
        "inputs": {
            "real_rate_file": "vietnam_real_interest_rate_1993_2023.csv",
            "debt_stock_file": "vietnam_wb_indicators.csv",
            "debt_variable": "ext_debt_to_gni (World Bank DT.DOD.DECT.GN.ZS)",
        },
        "constructed": {
            "r_bar_decimal": r_bar,
            "d_bar_decimal": d_bar,
            "spread_by_year": {str(y): s for y, s in zip(years, spread)},
            "debt_ratio_by_year": {str(y): d for y, d in zip(years, d_vals)},
        },
        "estimate": {
            "phi_b_hat": phi_hat,
            "phi_b_unconstrained": phi_uncon,
            "rmse_spread_fit": rmse,
        },
        "equation": "spread_t = phi_b * (exp(d_t - d_bar) - 1)",
        "notes": [
            "Strictly uses observed debt stock ratio and observed real rates.",
            "Debt variable is gross external debt to GNI proxy, not NIIP.",
            "Interpret as reduced-form data anchor."
        ],
    }

    (BASE_DIR / "phi_b_realdata_debtstock.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
