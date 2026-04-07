"""
Fetch selected Vietnam macro indicators from World Bank API and save CSV.
"""

from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

INDICATORS = {
    "ca_to_gdp": "BN.CAB.XOKA.GD.ZS",      # current account balance (% GDP)
    "ext_debt_to_gni": "DT.DOD.DECT.GN.ZS",# external debt stocks (% GNI)
    "gdp_current_usd": "NY.GDP.MKTP.CD",   # GDP (current US$)
    "inv_to_gdp": "NE.GDI.TOTL.ZS",        # gross capital formation (% GDP)
}


def fetch_indicator(indicator: str) -> dict[int, float]:
    url = f"https://api.worldbank.org/v2/country/VNM/indicator/{indicator}?format=json&per_page=200"
    with urllib.request.urlopen(url, timeout=30) as r:
        payload = json.loads(r.read().decode("utf-8"))
    out: dict[int, float] = {}
    for row in payload[1]:
        y = int(row["date"])
        v = row.get("value")
        if v is not None:
            out[y] = float(v)
    return out


def main() -> None:
    series = {name: fetch_indicator(code) for name, code in INDICATORS.items()}
    all_years = sorted(set().union(*[set(s.keys()) for s in series.values()]))

    out_path = BASE_DIR / "vietnam_wb_indicators.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["year", *series.keys()])
        for y in all_years:
            row = [y]
            for name in series.keys():
                v = series[name].get(y)
                row.append("" if v is None else v)
            w.writerow(row)

    meta = {
        "indicators": INDICATORS,
        "n_years": len(all_years),
        "latest": {
            name: {
                "year": max(s.keys()) if s else None,
                "value": s[max(s.keys())] if s else None,
            }
            for name, s in series.items()
        },
    }
    (BASE_DIR / "vietnam_wb_indicators_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(json.dumps(meta, indent=2))
    print(f"saved: {out_path}")


if __name__ == "__main__":
    main()
