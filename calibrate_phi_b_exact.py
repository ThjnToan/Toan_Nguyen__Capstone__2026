"""
High-precision one-dimensional calibration of phi_b to annual CA moments.

This script is designed to make phi_b calibration reproducible and transparent:
- Fix non-phi parameters (rho_ca, gamma_r, sigma_ca)
- Simulate the external block with long samples
- Match annual CA moments (std and optional ac1) with scalar optimization

Example:
  python calibrate_phi_b_exact.py \
    --targets vietnam_ca_targets.json \
    --base-data your_data.csv \
    --anchor smm_phi_b_results.json \
    --out phi_b_exact_vietnam.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from estimate_phi_b_smm import (  # type: ignore
    SmmConfig,
    load_data,
    simulate_external_paths,
    compute_ca_annual_moments,
)


def read_targets(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    ann = payload.get("targets_annual", payload)
    if "std_ca_to_gdp" not in ann:
        raise ValueError("targets file must contain targets_annual.std_ca_to_gdp")
    out = {"std_ca_ann": float(ann["std_ca_to_gdp"])}
    if "ac1_ca_to_gdp" in ann:
        out["ac1_ca_ann"] = float(ann["ac1_ca_to_gdp"])
    return out


def read_anchor(path: Path, rho: float | None, gamma: float | None, sigma: float | None) -> tuple[float, float, float]:
    if rho is not None and gamma is not None and sigma is not None:
        return float(rho), float(gamma), float(sigma)

    payload = json.loads(path.read_text(encoding="utf-8"))
    th = payload.get("theta_hat", {})
    if not all(k in th for k in ("rho_ca", "gamma_r", "sigma_ca")):
        raise ValueError("anchor file must contain theta_hat.rho_ca, gamma_r, sigma_ca")
    return float(th["rho_ca"]), float(th["gamma_r"]), float(th["sigma_ca"])


def main() -> None:
    p = argparse.ArgumentParser(description="Calibrate phi_b with high-precision scalar optimization.")
    p.add_argument("--targets", required=True, type=str, help="JSON with annual CA targets.")
    p.add_argument("--base-data", required=True, type=str, help="CSV for b_bar anchor (mean_b).")
    p.add_argument("--anchor", default="smm_phi_b_results.json", type=str, help="JSON with anchor params rho_ca, gamma_r, sigma_ca.")
    p.add_argument("--rho-ca", default=None, type=float, help="Optional manual rho_ca override.")
    p.add_argument("--gamma-r", default=None, type=float, help="Optional manual gamma_r override.")
    p.add_argument("--sigma-ca", default=None, type=float, help="Optional manual sigma_ca override.")
    p.add_argument("--w-std", default=1.0, type=float, help="Weight on std target.")
    p.add_argument("--w-ac1", default=1.25, type=float, help="Weight on ac1 target.")
    p.add_argument("--beta", default=0.99, type=float, help="Discount factor.")
    p.add_argument("--seed", default=42, type=int, help="Simulation seed.")
    p.add_argument("--n-sim", default=60000, type=int, help="Simulation length for precise objective.")
    p.add_argument("--burn-in", default=2000, type=int, help="Burn-in periods.")
    p.add_argument("--annual-block-size", default=4, type=int, help="Quarterly-to-annual aggregation block.")
    p.add_argument("--phi-min", default=1e-5, type=float, help="Lower bound for phi_b.")
    p.add_argument("--phi-max", default=0.25, type=float, help="Upper bound for phi_b.")
    p.add_argument("--out", default="phi_b_exact_vietnam.json", type=str, help="Output JSON path.")
    args = p.parse_args()

    targets = read_targets(Path(args.targets))
    rho_ca, gamma_r, sigma_ca = read_anchor(Path(args.anchor), args.rho_ca, args.gamma_r, args.sigma_ca)

    data = load_data(Path(args.base_data))
    mean_b = float(np.mean(data["b"]))

    cfg = SmmConfig(
        beta=float(args.beta),
        seed=int(args.seed),
        n_sim=int(args.n_sim),
        burn_in=int(args.burn_in),
        maxiter=200,
    )

    std_t = float(targets["std_ca_ann"])
    ac1_t = float(targets.get("ac1_ca_ann", np.nan))

    def objective(phi: float) -> float:
        theta = np.array([phi, rho_ca, gamma_r, sigma_ca], dtype=float)
        sim = simulate_external_paths(theta, b_bar=mean_b, cfg=cfg)
        m = compute_ca_annual_moments(sim["ca"], block_size=int(args.annual_block_size))

        loss = 0.0
        d_std = (m["std_ca_ann"] - std_t) / max(abs(std_t), 1e-8)
        loss += float(args.w_std) * d_std * d_std

        if np.isfinite(ac1_t):
            d_ac1 = (m["ac1_ca_ann"] - ac1_t) / max(abs(ac1_t), 1e-8)
            loss += float(args.w_ac1) * d_ac1 * d_ac1

        return float(loss)

    res = minimize_scalar(
        objective,
        method="bounded",
        bounds=(float(args.phi_min), float(args.phi_max)),
        options={"xatol": 1e-10, "maxiter": 500},
    )

    phi_hat = float(res.x)
    theta_hat = np.array([phi_hat, rho_ca, gamma_r, sigma_ca], dtype=float)
    sim_hat = simulate_external_paths(theta_hat, b_bar=mean_b, cfg=cfg)
    model_ann = compute_ca_annual_moments(sim_hat["ca"], block_size=int(args.annual_block_size))

    # Small profile around optimum for transparency.
    lo = max(float(args.phi_min), phi_hat * 0.5)
    hi = min(float(args.phi_max), phi_hat * 1.5)
    grid = np.linspace(lo, hi, 21)
    profile = [{"phi_b": float(x), "objective": float(objective(float(x)))} for x in grid]

    out = {
        "phi_b_hat": phi_hat,
        "objective_value": float(res.fun),
        "success": bool(res.success),
        "message": str(res.message),
        "mode": "phi-exact-1d",
        "targets_annual": targets,
        "model_moments_annual": model_ann,
        "fixed_parameters": {
            "rho_ca": rho_ca,
            "gamma_r": gamma_r,
            "sigma_ca": sigma_ca,
            "mean_b": mean_b,
        },
        "weights": {
            "w_std": float(args.w_std),
            "w_ac1": float(args.w_ac1),
        },
        "config": {
            "beta": cfg.beta,
            "seed": cfg.seed,
            "n_sim": cfg.n_sim,
            "burn_in": cfg.burn_in,
            "annual_block_size": int(args.annual_block_size),
            "phi_min": float(args.phi_min),
            "phi_max": float(args.phi_max),
            "targets_path": str(args.targets),
            "base_data_path": str(args.base_data),
            "anchor_path": str(args.anchor),
        },
        "objective_profile": profile,
    }

    out_path = Path(args.out)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print(f"phi_b_hat: {phi_hat:.12f}")
    print(f"objective: {float(res.fun):.12f}")
    print(f"saved: {out_path}")


if __name__ == "__main__":
    main()
