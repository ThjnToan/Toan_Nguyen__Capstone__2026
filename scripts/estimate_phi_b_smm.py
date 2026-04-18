"""
Estimate debt-elastic premium phi_b with Simulated Method of Moments (SMM).

This script is intentionally lightweight so it can run directly in the thesis
workspace without Dynare dependencies.

Model block (quarterly):
    r_t = r_bar + phi_b * (exp(B_bar - B_t) - 1)
    CA_t = rho_ca * CA_{t-1} - gamma_r * (r_t - r_bar) + eps_t
    B_t = (1 + r_{t-1}) * B_{t-1} + CA_t

Moments targeted:
    - mean(B), std(B), ac1(B)
    - std(CA), ac1(CA)
    - mean(spread), std(spread)

Input data:
    Preferred CSV columns (case-insensitive aliases supported):
      nfa_to_gdp, ca_to_gdp, spread
    Optional:
      If ca_to_gdp is missing but nfa/spread exist, CA is reconstructed from
      CA_t = B_t - (1 + r_{t-1}) * B_{t-1}.

Usage examples:
    python estimate_phi_b_smm.py --data irf_dynare.csv --out smm_phi_b_results.json
    python estimate_phi_b_smm.py --data your_data.csv --bootstrap 100
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution


# ------------------------------
# Utilities
# ------------------------------


def ac1(x: np.ndarray) -> float:
    """Lag-1 autocorrelation with safe fallback."""
    x = np.asarray(x, dtype=float)
    if x.size < 3:
        return 0.0
    x0 = x[:-1] - np.mean(x[:-1])
    x1 = x[1:] - np.mean(x[1:])
    denom = np.sqrt(np.dot(x0, x0) * np.dot(x1, x1))
    if denom <= 1e-14:
        return 0.0
    return float(np.dot(x0, x1) / denom)


def _match_column(df: pd.DataFrame, aliases: Tuple[str, ...]) -> str | None:
    normalized = {c.lower().strip(): c for c in df.columns}
    for a in aliases:
        if a in normalized:
            return normalized[a]
    return None


def load_data(data_path: Path) -> Dict[str, np.ndarray]:
    """Load external-sector data and return cleaned arrays."""
    df = pd.read_csv(data_path)

    nfa_col = _match_column(
        df,
        (
            "nfa_to_gdp",
            "net_foreign_assets",
            "net_foreign_assets_to_gdp",
            "b_star",
            "b",
        ),
    )
    ca_col = _match_column(
        df,
        (
            "ca_to_gdp",
            "current_account_to_gdp",
            "current_account",
            "trade_balance_to_gdp",
            "tb_to_gdp",
            "ca",
            "tb",
        ),
    )
    spread_col = _match_column(
        df,
        (
            "spread",
            "interest_spread",
            "sovereign_spread",
            "interest_rate",
            "r_star",
            "r",
        ),
    )

    if nfa_col is None:
        raise ValueError("Could not find an NFA column. Expected aliases include nfa_to_gdp or b_star.")

    b = df[nfa_col].astype(float).to_numpy()

    if spread_col is None:
        raise ValueError("Could not find spread/interest column. Expected aliases include spread or interest_rate.")

    spread_raw = df[spread_col].astype(float).to_numpy()

    # If column looks like a gross interest rate level in decimals, keep as is.
    # If it looks like basis points (e.g., values > 1 on average), convert to decimal.
    if np.nanmean(np.abs(spread_raw)) > 1.0:
        spread = spread_raw / 10000.0
    else:
        spread = spread_raw

    if ca_col is not None:
        ca = df[ca_col].astype(float).to_numpy()
    else:
        # Reconstruct CA using a conservative assumption that spread proxies r_{t-1}-r_bar.
        # Here we assume r_bar is near zero in short samples and reconstruct approximately.
        ca = np.zeros_like(b)
        for t in range(1, len(b)):
            ca[t] = b[t] - (1.0 + spread[t - 1]) * b[t - 1]

    mask = np.isfinite(b) & np.isfinite(ca) & np.isfinite(spread)
    b = b[mask]
    ca = ca[mask]
    spread = spread[mask]

    if b.size < 20:
        raise ValueError("Need at least 20 usable observations after cleaning for SMM.")

    return {"b": b, "ca": ca, "spread": spread}


# ------------------------------
# Moment construction
# ------------------------------


def compute_moments(b: np.ndarray, ca: np.ndarray, spread: np.ndarray) -> Dict[str, float]:
    return {
        "mean_b": float(np.mean(b)),
        "std_b": float(np.std(b, ddof=1)),
        "ac1_b": float(ac1(b)),
        "std_ca": float(np.std(ca, ddof=1)),
        "ac1_ca": float(ac1(ca)),
        "mean_spread": float(np.mean(spread)),
        "std_spread": float(np.std(spread, ddof=1)),
    }


def apply_moment_overrides(data_moments: Dict[str, float], override_path: Path) -> Dict[str, float]:
    """Override selected target moments from a JSON file.

    Supported keys (top-level or nested under 'smm_moment_overrides'):
      mean_b, std_b, ac1_b, std_ca, ac1_ca, mean_spread, std_spread
    """
    payload = json.loads(override_path.read_text(encoding="utf-8"))
    overrides = payload.get("smm_moment_overrides", payload)

    out = dict(data_moments)
    for k in out.keys():
        if k in overrides:
            out[k] = float(overrides[k])
    return out


def to_annual_series(x: np.ndarray, block_size: int = 4) -> np.ndarray:
    """Convert a quarterly ratio series to annual averages by non-overlapping blocks."""
    x = np.asarray(x, dtype=float)
    n_blocks = x.size // block_size
    if n_blocks < 3:
        return np.asarray([], dtype=float)
    x_trim = x[: n_blocks * block_size]
    return x_trim.reshape(n_blocks, block_size).mean(axis=1)


def compute_ca_annual_moments(ca_qtr: np.ndarray, block_size: int = 4) -> Dict[str, float]:
    ca_ann = to_annual_series(ca_qtr, block_size=block_size)
    if ca_ann.size < 3:
        return {"std_ca_ann": 1e6, "ac1_ca_ann": 1e6}
    return {
        "std_ca_ann": float(np.std(ca_ann, ddof=1)),
        "ac1_ca_ann": float(ac1(ca_ann)),
    }


@dataclass
class SmmConfig:
    beta: float = 0.99
    seed: int = 42
    n_sim: int = 4000
    burn_in: int = 500
    maxiter: int = 200


# ------------------------------
# Simulation model
# ------------------------------


def simulate_external_block(theta: np.ndarray, data_targets: Dict[str, float], cfg: SmmConfig) -> Dict[str, float]:
    """Simulate model-implied moments for parameter vector theta.

    theta = [phi_b, rho_ca, gamma_r, sigma_ca]
    """
    phi_b, rho_ca, gamma_r, sigma_ca = theta

    # Parameter guards for numerical safety
    if phi_b <= 0.0 or sigma_ca <= 0.0 or abs(rho_ca) >= 0.999 or gamma_r < 0.0:
        return {
            "mean_b": 1e6,
            "std_b": 1e6,
            "ac1_b": 1e6,
            "std_ca": 1e6,
            "ac1_ca": 1e6,
            "mean_spread": 1e6,
            "std_spread": 1e6,
        }

    rng = np.random.default_rng(cfg.seed)

    b_bar = data_targets["mean_b"]
    r_bar = 1.0 / cfg.beta - 1.0

    T = cfg.n_sim + cfg.burn_in
    b = np.zeros(T)
    ca = np.zeros(T)
    spread = np.zeros(T)

    b[0] = b_bar

    for t in range(1, T):
        # Risk premium (spread over r_bar)
        spread[t - 1] = phi_b * (np.exp(b_bar - b[t - 1]) - 1.0)
        r_t = r_bar + spread[t - 1]

        eps_t = rng.normal(0.0, sigma_ca)
        ca[t] = rho_ca * ca[t - 1] - gamma_r * spread[t - 1] + eps_t

        b[t] = (1.0 + r_t) * b[t - 1] + ca[t]

        # Soft clipping keeps simulation stable in bad regions of parameter space.
        if abs(b[t]) > 5.0:
            b[t] = np.sign(b[t]) * 5.0

    spread[-1] = phi_b * (np.exp(b_bar - b[-1]) - 1.0)

    b_s = b[cfg.burn_in :]
    ca_s = ca[cfg.burn_in :]
    sp_s = spread[cfg.burn_in :]

    return compute_moments(b_s, ca_s, sp_s)


def simulate_external_paths(theta: np.ndarray, b_bar: float, cfg: SmmConfig) -> Dict[str, np.ndarray]:
    """Simulate quarterly paths for B, CA, and spread under theta."""
    phi_b, rho_ca, gamma_r, sigma_ca = theta

    if phi_b <= 0.0 or sigma_ca <= 0.0 or abs(rho_ca) >= 0.999 or gamma_r < 0.0:
        return {
            "b": np.full(10, 1e6, dtype=float),
            "ca": np.full(10, 1e6, dtype=float),
            "spread": np.full(10, 1e6, dtype=float),
        }

    rng = np.random.default_rng(cfg.seed)
    r_bar = 1.0 / cfg.beta - 1.0

    T = cfg.n_sim + cfg.burn_in
    b = np.zeros(T)
    ca = np.zeros(T)
    spread = np.zeros(T)

    b[0] = float(b_bar)

    for t in range(1, T):
        spread[t - 1] = phi_b * (np.exp(b_bar - b[t - 1]) - 1.0)
        r_t = r_bar + spread[t - 1]

        eps_t = rng.normal(0.0, sigma_ca)
        ca[t] = rho_ca * ca[t - 1] - gamma_r * spread[t - 1] + eps_t
        b[t] = (1.0 + r_t) * b[t - 1] + ca[t]

        if abs(b[t]) > 5.0:
            b[t] = np.sign(b[t]) * 5.0

    spread[-1] = phi_b * (np.exp(b_bar - b[-1]) - 1.0)

    return {
        "b": b[cfg.burn_in :],
        "ca": ca[cfg.burn_in :],
        "spread": spread[cfg.burn_in :],
    }


# ------------------------------
# SMM objective and estimation
# ------------------------------


def objective(theta: np.ndarray, data_moments: Dict[str, float], weights: Dict[str, float], cfg: SmmConfig) -> float:
    model_m = simulate_external_block(theta, data_moments, cfg)

    loss = 0.0
    for k in data_moments:
        scale = max(abs(data_moments[k]), 1e-6)
        diff = (model_m[k] - data_moments[k]) / scale
        loss += weights[k] * diff * diff
    return float(loss)


def estimate_phi_b(data_moments: Dict[str, float], cfg: SmmConfig) -> Dict[str, object]:
    # Weight autocorrelations and spread moments slightly more for identification.
    weights = {
        "mean_b": 0.75,
        "std_b": 1.00,
        "ac1_b": 1.25,
        "std_ca": 1.00,
        "ac1_ca": 1.25,
        "mean_spread": 1.25,
        "std_spread": 1.25,
    }

    bounds = [
        (1e-5, 0.25),   # phi_b
        (0.00, 0.98),   # rho_ca
        (0.00, 2.00),   # gamma_r
        (1e-6, 0.05),   # sigma_ca
    ]

    result = differential_evolution(
        func=lambda th: objective(th, data_moments, weights, cfg),
        bounds=bounds,
        strategy="best1bin",
        maxiter=cfg.maxiter,
        popsize=16,
        tol=1e-7,
        polish=True,
        seed=cfg.seed,
        updating="deferred",
        workers=1,
    )

    theta_hat = result.x
    model_mom = simulate_external_block(theta_hat, data_moments, cfg)

    return {
        "theta_hat": {
            "phi_b": float(theta_hat[0]),
            "rho_ca": float(theta_hat[1]),
            "gamma_r": float(theta_hat[2]),
            "sigma_ca": float(theta_hat[3]),
        },
        "objective_value": float(result.fun),
        "success": bool(result.success),
        "message": str(result.message),
        "data_moments": data_moments,
        "model_moments": model_mom,
        "weights": weights,
    }


def estimate_phi_b_sgu_ca_only(
    targets_ann: Dict[str, float],
    b_bar: float,
    cfg: SmmConfig,
    rho_ca_fixed: float,
    gamma_r_fixed: float,
    sigma_ca_fixed: float,
    block_size: int = 4,
) -> Dict[str, object]:
    """Estimate phi_b by matching only annual CA moments (SGU-style CA mode)."""

    if not (0.0 <= rho_ca_fixed < 0.999):
        raise ValueError("rho_ca_fixed must satisfy 0 <= rho_ca_fixed < 0.999")
    if gamma_r_fixed < 0.0:
        raise ValueError("gamma_r_fixed must be nonnegative")
    if sigma_ca_fixed <= 0.0:
        raise ValueError("sigma_ca_fixed must be strictly positive")

    w_std = 1.0
    w_ac1 = 1.25

    std_target = float(targets_ann["std_ca_ann"])
    ac1_target = float(targets_ann["ac1_ca_ann"])

    def obj_phi(phi_vec: np.ndarray) -> float:
        phi = float(phi_vec[0])
        theta = np.array([phi, rho_ca_fixed, gamma_r_fixed, sigma_ca_fixed], dtype=float)
        sim = simulate_external_paths(theta, b_bar=b_bar, cfg=cfg)
        m = compute_ca_annual_moments(sim["ca"], block_size=block_size)

        diff_std = (m["std_ca_ann"] - std_target) / max(abs(std_target), 1e-6)
        diff_ac1 = (m["ac1_ca_ann"] - ac1_target) / max(abs(ac1_target), 1e-6)
        return float(w_std * diff_std * diff_std + w_ac1 * diff_ac1 * diff_ac1)

    res = differential_evolution(
        func=obj_phi,
        bounds=[(1e-5, 0.25)],
        strategy="best1bin",
        maxiter=cfg.maxiter,
        popsize=16,
        tol=1e-7,
        polish=True,
        seed=cfg.seed,
        updating="deferred",
        workers=1,
    )

    phi_hat = float(res.x[0])
    theta_hat = np.array([phi_hat, rho_ca_fixed, gamma_r_fixed, sigma_ca_fixed], dtype=float)
    sim_hat = simulate_external_paths(theta_hat, b_bar=b_bar, cfg=cfg)
    model_ann = compute_ca_annual_moments(sim_hat["ca"], block_size=block_size)

    return {
        "theta_hat": {
            "phi_b": phi_hat,
            "rho_ca": float(rho_ca_fixed),
            "gamma_r": float(gamma_r_fixed),
            "sigma_ca": float(sigma_ca_fixed),
        },
        "objective_value": float(res.fun),
        "success": bool(res.success),
        "message": str(res.message),
        "mode": "sgu-ca-only",
        "data_moments_annual": {
            "std_ca_ann": std_target,
            "ac1_ca_ann": ac1_target,
        },
        "model_moments_annual": model_ann,
        "weights": {
            "std_ca_ann": w_std,
            "ac1_ca_ann": w_ac1,
        },
    }


def bootstrap_phi_b(data: Dict[str, np.ndarray], cfg: SmmConfig, n_bootstrap: int) -> Dict[str, float]:
    """Simple block bootstrap over rows to quantify phi_b uncertainty."""
    if n_bootstrap <= 0:
        return {}

    b = data["b"]
    ca = data["ca"]
    spread = data["spread"]
    n = len(b)

    # Quarterly data: 4-quarter blocks preserve short persistence.
    block = 4
    n_blocks = int(np.ceil(n / block))

    rng = np.random.default_rng(cfg.seed + 100)
    estimates = []

    for _ in range(n_bootstrap):
        starts = rng.integers(0, max(1, n - block + 1), size=n_blocks)
        idx = []
        for s in starts:
            idx.extend(range(s, min(s + block, n)))
        idx = np.array(idx[:n], dtype=int)

        b_bs = b[idx]
        ca_bs = ca[idx]
        sp_bs = spread[idx]

        dm = compute_moments(b_bs, ca_bs, sp_bs)
        est = estimate_phi_b(dm, cfg)
        estimates.append(est["theta_hat"]["phi_b"])

    arr = np.array(estimates, dtype=float)
    return {
        "phi_b_bootstrap_mean": float(np.mean(arr)),
        "phi_b_bootstrap_std": float(np.std(arr, ddof=1)),
        "phi_b_bootstrap_p05": float(np.quantile(arr, 0.05)),
        "phi_b_bootstrap_p50": float(np.quantile(arr, 0.50)),
        "phi_b_bootstrap_p95": float(np.quantile(arr, 0.95)),
    }


def generate_synthetic_data(theta_true: np.ndarray, cfg: SmmConfig, n_obs: int = 240) -> Dict[str, np.ndarray]:
    """Generate synthetic quarterly data for estimator validation."""
    dm_proxy = {
        "mean_b": 0.0,
        "std_b": 0.01,
        "ac1_b": 0.5,
        "std_ca": 0.005,
        "ac1_ca": 0.5,
        "mean_spread": 0.001,
        "std_spread": 0.001,
    }
    m = simulate_external_block(theta_true, dm_proxy, cfg)

    # Regenerate the full simulated paths with the same DGP as the objective.
    phi_b, rho_ca, gamma_r, sigma_ca = theta_true
    rng = np.random.default_rng(cfg.seed + 7)
    r_bar = 1.0 / cfg.beta - 1.0
    b_bar = 0.0

    T = n_obs + cfg.burn_in
    b = np.zeros(T)
    ca = np.zeros(T)
    spread = np.zeros(T)

    for t in range(1, T):
        spread[t - 1] = phi_b * (np.exp(b_bar - b[t - 1]) - 1.0)
        r_t = r_bar + spread[t - 1]
        eps_t = rng.normal(0.0, sigma_ca)
        ca[t] = rho_ca * ca[t - 1] - gamma_r * spread[t - 1] + eps_t
        b[t] = (1.0 + r_t) * b[t - 1] + ca[t]
        if abs(b[t]) > 5.0:
            b[t] = np.sign(b[t]) * 5.0

    spread[-1] = phi_b * (np.exp(b_bar - b[-1]) - 1.0)

    return {
        "b": b[cfg.burn_in :],
        "ca": ca[cfg.burn_in :],
        "spread": spread[cfg.burn_in :],
        "reference_moments": m,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Estimate debt-elastic premium phi_b using SMM.")
    parser.add_argument("--data", required=True, type=str, help="Path to CSV data file.")
    parser.add_argument("--out", default="smm_phi_b_results.json", type=str, help="Output JSON file.")
    parser.add_argument("--bootstrap", default=0, type=int, help="Number of bootstrap replications.")
    parser.add_argument("--beta", default=0.99, type=float, help="Discount factor for r_bar = 1/beta - 1.")
    parser.add_argument("--seed", default=42, type=int, help="Random seed.")
    parser.add_argument("--n-sim", default=4000, type=int, help="Simulation length per objective evaluation.")
    parser.add_argument("--burn-in", default=500, type=int, help="Burn-in periods for simulation.")
    parser.add_argument("--maxiter", default=200, type=int, help="Max optimizer iterations.")
    parser.add_argument("--self-test", action="store_true", help="Run estimation on synthetic data generated from known parameters.")
    parser.add_argument("--self-test-out", default="smm_self_test_data.csv", type=str, help="Synthetic CSV path used when --self-test is enabled.")
    parser.add_argument("--moments-override", default="", type=str, help="Optional JSON path with moment overrides (e.g., std_ca, ac1_ca).")
    parser.add_argument("--sgu-ca-only", action="store_true", help="Estimate phi_b by matching only annual CA moments.")
    parser.add_argument("--rho-ca-fixed", default=0.57, type=float, help="Fixed rho_ca used in --sgu-ca-only mode.")
    parser.add_argument("--gamma-r-fixed", default=0.80, type=float, help="Fixed gamma_r used in --sgu-ca-only mode.")
    parser.add_argument("--sigma-ca-fixed", default=0.003, type=float, help="Fixed sigma_ca used in --sgu-ca-only mode.")
    parser.add_argument("--annual-block-size", default=4, type=int, help="Quarterly-to-annual block size for --sgu-ca-only mode.")
    args = parser.parse_args()

    data_path = Path(args.data)
    out_path = Path(args.out)

    if (not args.self_test) and (not data_path.exists()):
        raise FileNotFoundError(f"Data file not found: {data_path}")

    cfg = SmmConfig(
        beta=float(args.beta),
        seed=int(args.seed),
        n_sim=int(args.n_sim),
        burn_in=int(args.burn_in),
        maxiter=int(args.maxiter),
    )

    if args.self_test:
        theta_true = np.array([0.05, 0.65, 0.40, 0.003], dtype=float)
        synth = generate_synthetic_data(theta_true, cfg)
        synth_df = pd.DataFrame(
            {
                "nfa_to_gdp": synth["b"],
                "ca_to_gdp": synth["ca"],
                "spread": synth["spread"],
            }
        )
        synth_path = Path(args.self_test_out)
        synth_df.to_csv(synth_path, index=False)
        data = {"b": synth["b"], "ca": synth["ca"], "spread": synth["spread"]}
    else:
        data = load_data(data_path)

    data_moments = compute_moments(data["b"], data["ca"], data["spread"])

    if args.sgu_ca_only:
        if not args.moments_override:
            raise ValueError("--sgu-ca-only requires --moments-override with annual CA targets.")
        override_path = Path(args.moments_override)
        if not override_path.exists():
            raise FileNotFoundError(f"Moments override file not found: {override_path}")

        payload = json.loads(override_path.read_text(encoding="utf-8"))
        ann = payload.get("targets_annual", {})
        if "std_ca_to_gdp" not in ann or "ac1_ca_to_gdp" not in ann:
            raise ValueError("Override file must contain targets_annual.std_ca_to_gdp and targets_annual.ac1_ca_to_gdp")

        targets_ann = {
            "std_ca_ann": float(ann["std_ca_to_gdp"]),
            "ac1_ca_ann": float(ann["ac1_ca_to_gdp"]),
        }

        result = estimate_phi_b_sgu_ca_only(
            targets_ann=targets_ann,
            b_bar=float(data_moments["mean_b"]),
            cfg=cfg,
            rho_ca_fixed=float(args.rho_ca_fixed),
            gamma_r_fixed=float(args.gamma_r_fixed),
            sigma_ca_fixed=float(args.sigma_ca_fixed),
            block_size=int(args.annual_block_size),
        )
    else:
        if args.moments_override:
            override_path = Path(args.moments_override)
            if not override_path.exists():
                raise FileNotFoundError(f"Moments override file not found: {override_path}")
            data_moments = apply_moment_overrides(data_moments, override_path)

        result = estimate_phi_b(data_moments, cfg)

    if args.bootstrap > 0:
        result["bootstrap"] = bootstrap_phi_b(data, cfg, int(args.bootstrap))

    result["config"] = {
        "beta": cfg.beta,
        "seed": cfg.seed,
        "n_sim": cfg.n_sim,
        "burn_in": cfg.burn_in,
        "maxiter": cfg.maxiter,
        "data_path": str(data_path),
        "self_test": bool(args.self_test),
        "moments_override": str(args.moments_override) if args.moments_override else "",
        "sgu_ca_only": bool(args.sgu_ca_only),
        "rho_ca_fixed": float(args.rho_ca_fixed),
        "gamma_r_fixed": float(args.gamma_r_fixed),
        "sigma_ca_fixed": float(args.sigma_ca_fixed),
        "annual_block_size": int(args.annual_block_size),
    }

    if args.self_test:
        result["self_test_truth"] = {
            "phi_b": 0.05,
            "rho_ca": 0.65,
            "gamma_r": 0.40,
            "sigma_ca": 0.003,
        }
        result["config"]["self_test_data_path"] = str(Path(args.self_test_out))

    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("SMM estimation completed.")
    print(f"Estimated phi_b: {result['theta_hat']['phi_b']:.6f}")
    print(f"Objective value: {result['objective_value']:.6f}")
    print(f"Saved results: {out_path}")


if __name__ == "__main__":
    main()
