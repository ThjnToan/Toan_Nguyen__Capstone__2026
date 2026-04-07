"""
Profile-identification curve for phi_b in the external block.

For each fixed phi_b on a grid, optimize (rho_ca, gamma_r, sigma_ca)
against the selected moment set and store the minimized objective.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict

import numpy as np

from scipy.optimize import differential_evolution

from estimate_phi_b_smm import (  # type: ignore
    SmmConfig,
    load_data,
    compute_moments,
    apply_moment_overrides,
    objective,
    simulate_external_block,
)


def optimize_conditional(phi_b: float, data_moments: Dict[str, float], cfg: SmmConfig, weights: Dict[str, float]) -> Dict[str, float]:
    bounds = [
        (0.00, 0.98),   # rho_ca
        (0.00, 2.00),   # gamma_r
        (1e-6, 0.05),   # sigma_ca
    ]

    def obj_rest(x: np.ndarray) -> float:
        theta = np.array([phi_b, x[0], x[1], x[2]], dtype=float)
        return objective(theta, data_moments, weights, cfg)

    res = differential_evolution(
        func=obj_rest,
        bounds=bounds,
        strategy="best1bin",
        maxiter=cfg.maxiter,
        popsize=14,
        tol=1e-7,
        polish=True,
        seed=cfg.seed,
        updating="deferred",
        workers=1,
    )

    x = res.x
    return {
        "objective": float(res.fun),
        "rho_ca": float(x[0]),
        "gamma_r": float(x[1]),
        "sigma_ca": float(x[2]),
        "success": bool(res.success),
        "message": str(res.message),
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Profile phi_b by re-optimizing non-phi parameters.")
    p.add_argument("--data", required=True, type=str)
    p.add_argument("--moments-override", default="", type=str)
    p.add_argument("--phi-min", default=1e-5, type=float)
    p.add_argument("--phi-max", default=0.20, type=float)
    p.add_argument("--n-grid", default=21, type=int)
    p.add_argument("--beta", default=0.99, type=float)
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--n-sim", default=3000, type=int)
    p.add_argument("--burn-in", default=500, type=int)
    p.add_argument("--maxiter", default=90, type=int)
    p.add_argument("--out", default="phi_b_profile_joint.json", type=str)
    args = p.parse_args()

    data = load_data(Path(args.data))
    data_moments = compute_moments(data["b"], data["ca"], data["spread"])

    if args.moments_override:
        over = Path(args.moments_override)
        if not over.exists():
            raise FileNotFoundError(f"Override file not found: {over}")
        data_moments = apply_moment_overrides(data_moments, over)

    cfg = SmmConfig(
        beta=float(args.beta),
        seed=int(args.seed),
        n_sim=int(args.n_sim),
        burn_in=int(args.burn_in),
        maxiter=int(args.maxiter),
    )

    # Same weighting scheme as estimator.
    weights = {
        "mean_b": 0.75,
        "std_b": 1.00,
        "ac1_b": 1.25,
        "std_ca": 1.00,
        "ac1_ca": 1.25,
        "mean_spread": 1.25,
        "std_spread": 1.25,
    }

    phi_grid = np.linspace(float(args.phi_min), float(args.phi_max), int(args.n_grid))
    rows = []

    for phi in phi_grid:
        cond = optimize_conditional(float(phi), data_moments, cfg, weights)
        theta = np.array([phi, cond["rho_ca"], cond["gamma_r"], cond["sigma_ca"]], dtype=float)
        mm = simulate_external_block(theta, data_moments, cfg)
        rows.append({
            "phi_b": float(phi),
            **cond,
            "model_std_ca": float(mm["std_ca"]),
            "model_ac1_ca": float(mm["ac1_ca"]),
        })

    best = min(rows, key=lambda r: r["objective"]) if rows else None

    # Identification band based on objective inflation <= 10% from minimum.
    ident_band = []
    if best is not None:
        thr = best["objective"] * 1.10
        ident_band = [r["phi_b"] for r in rows if r["objective"] <= thr]

    out = {
        "mode": "phi-profile-joint",
        "data_moments": data_moments,
        "config": {
            "beta": cfg.beta,
            "seed": cfg.seed,
            "n_sim": cfg.n_sim,
            "burn_in": cfg.burn_in,
            "maxiter": cfg.maxiter,
            "phi_min": float(args.phi_min),
            "phi_max": float(args.phi_max),
            "n_grid": int(args.n_grid),
            "moments_override": str(args.moments_override),
        },
        "weights": weights,
        "profile": rows,
        "best": best,
        "identification_band_10pct_obj": {
            "phi_min": min(ident_band) if ident_band else None,
            "phi_max": max(ident_band) if ident_band else None,
        },
    }

    Path(args.out).write_text(json.dumps(out, indent=2), encoding="utf-8")
    if best is not None:
        print(f"best phi_b: {best['phi_b']:.6f}")
        print(f"best objective: {best['objective']:.6f}")
    print(f"saved: {args.out}")


if __name__ == "__main__":
    main()
