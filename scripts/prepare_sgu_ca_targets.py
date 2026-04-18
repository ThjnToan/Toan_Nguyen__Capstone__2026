"""
Prepare SGU-style CA/GDP moment targets from annual Vietnam data.

Usage:
    python prepare_sgu_ca_targets.py --in vietnam_ca_2016_2025.csv --out vietnam_ca_targets.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def ac1(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    if x.size < 3:
        return 0.0
    x0 = x[:-1] - np.mean(x[:-1])
    x1 = x[1:] - np.mean(x[1:])
    denom = np.sqrt(np.dot(x0, x0) * np.dot(x1, x1))
    if denom <= 1e-14:
        return 0.0
    return float(np.dot(x0, x1) / denom)


def load_ca_series(path: Path) -> np.ndarray:
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}

    if "ca_to_gdp" in cols:
        return df[cols["ca_to_gdp"]].astype(float).to_numpy()
    if "ca_to_gdp_pct" in cols:
        return (df[cols["ca_to_gdp_pct"]].astype(float) / 100.0).to_numpy()

    raise ValueError("Input must contain ca_to_gdp (decimal) or ca_to_gdp_pct (percent).")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build CA/GDP moment targets for SGU-style calibration.")
    parser.add_argument("--in", dest="infile", required=True, type=str, help="Input CSV with annual CA/GDP.")
    parser.add_argument("--out", dest="outfile", default="vietnam_ca_targets.json", type=str, help="Output JSON path.")
    args = parser.parse_args()

    ca = load_ca_series(Path(args.infile))
    ca = ca[np.isfinite(ca)]
    if ca.size < 5:
        raise ValueError("Need at least 5 observations.")

    mean_ann = float(np.mean(ca))
    std_ann = float(np.std(ca, ddof=1))
    ac1_ann = float(ac1(ca))

    # Frequency conversion proxy under i.i.d. within-year aggregation.
    std_qtr_proxy = float(std_ann / np.sqrt(4.0))

    result = {
        "series_info": {
            "frequency": "annual",
            "n_obs": int(ca.size),
            "years": "2016-2025",
        },
        "targets_annual": {
            "mean_ca_to_gdp": mean_ann,
            "std_ca_to_gdp": std_ann,
            "ac1_ca_to_gdp": ac1_ann,
        },
        "targets_quarterly_proxy": {
            "std_ca_to_gdp": std_qtr_proxy,
            "notes": "Proxy using annual std divided by sqrt(4). Use with caution.",
        },
        "smm_moment_overrides": {
            "std_ca": std_ann,
            "ac1_ca": ac1_ann,
        },
    }

    out_path = Path(args.outfile)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
