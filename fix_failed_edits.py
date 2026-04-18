#!/usr/bin/env python3
"""
Apply 18 failed text replacements to main_body.tex.
Uses slice-based replacement: find start+end positions, replace region directly.
"""

import os

FILE = r"C:\Users\Laptop K1\OneDrive\Documents\MATLAB\main_body.tex"

with open(FILE, "rb") as f:
    raw = f.read()

content = raw.decode("utf-8")
content = content.replace("\r\n", "\n")

results = {}

def slice_replace(content, start_phrase, end_phrase, new_text, label):
    """Find text between start_phrase and end_phrase (inclusive) and replace with new_text."""
    s = content.find(start_phrase)
    if s < 0:
        return content, f"FAILED - start phrase not found: {repr(start_phrase[:60])}"
    e = content.find(end_phrase, s)
    if e < 0:
        return content, f"FAILED - end phrase not found: {repr(end_phrase[:60])}"
    e += len(end_phrase)
    new_content = content[:s] + new_text + content[e:]
    return new_content, "SUCCESS"


# ============================================================
# METH-1
# ============================================================
content, results["METH-1"] = slice_replace(
    content,
    "paragraph{Battery Investment Decision.}",
    "reflecting the scarcity of storage relative to grid capital. The reliability gradient ($\\partial u / \\partial F = 0.20$) is moderate because the system operates near the flat portion of the exponential curve at $u_{ss} = 0.97$.",
    r"""\paragraph{Battery Investment Decision.}
	A fully endogenous optimization of battery capital subject to the non-linear reliability constraint would create non-convexities that complicate the Blanchard-Kahn conditions. The investment rule below proxies the gradient of the structural shadow price while preserving local saddle-path stability; the full shadow-price derivation is provided in Appendix~\ref{app:deriv_euler}.

	Table~\ref{tab:mpk_decomp} evaluates the shadow rental rate at steady state, yielding $R_{b,ss} = 10.93\%$ quarterly. The \textit{flexibility bottleneck} ($\partial F / \partial(A_{bat} K_b) = 2.04$) is the dominant force, reflecting battery scarcity in the CES aggregator, while the reliability gradient is moderate because the system operates near $u_{ss} = 0.97$.""",
    "METH-1"
)

# ============================================================
# METH-3
# ============================================================
content, results["METH-3"] = slice_replace(
    content,
    "The spending item $I_{p,t}$ in the budget constraint",
    "The Euler equation \\eqref{eq:euler} below follows directly from the first-order condition with respect to $K_{p,t+1}$.",
    r"The household chooses $\{C_t, L_t, K_{p,t+1}, B^*_t\}$ to maximize utility subject to the budget constraint \eqref{eq:hh_budget} and the standard capital accumulation law \eqref{eq:productive_capital}; the Lagrangian derivation is in Appendix~\ref{app:deriv_euler}.",
    "METH-3"
)

# ============================================================
# METH-6
# ============================================================
content, results["METH-6"] = slice_replace(
    content,
    "This exponential formulation ensures that:",
    "capturing the convex integration costs documented in the empirical literature.",
    r"This exponential form is standard in reliability engineering \citep{zakeri2015, lund2015}: $u_t \in (0,1)$ for all $F_t > 0$, reliability is bounded between grid collapse and perfection, and marginal returns to flexibility are strictly diminishing.",
    "METH-6"
)

# ============================================================
# METH-8
# ============================================================
content, results["METH-8"] = slice_replace(
    content,
    "The model is embedded in a small open economy framework",
    "in equilibrium, the domestic marginal product of capital equals the external borrowing rate.",
    r"""The model follows the debt-elastic interest rate approach of \citet{schmittgrohe2003}, giving the household access to international bond markets at a country-specific rate.

	\subsubsection{Bond Euler Equation}

	The household's optimality condition for international bond holdings yields:
	\begin{equation}
		\frac{1}{C_t} = \beta \, \mathbb{E}_t \left[ \frac{1 + r^*_t}{C_{t+1}} \right]
		\label{eq:bond_euler}
	\end{equation}

	Combining with \eqref{eq:euler} delivers the no-arbitrage condition $(1-\tau)\alpha Y_{t+1}/K_{p,t} + 1 - \delta_p = 1 + r^*_t$, confirming that the domestic after-tax return equals the external borrowing rate in equilibrium.""",
    "METH-8"
)

# ============================================================
# METH-10
# ============================================================
content, results["METH-10"] = slice_replace(
    content,
    "\\subsubsection{Current Account and Budget Constraint}",
    "resolving the well-known unit root problem in small open economy models \\citep{schmittgrohe2003}.",
    r"""\subsubsection{Debt-Elastic Interest Rate Premium}

	The country-specific interest rate rises with external debt, following \citet{schmittgrohe2003}:
	\begin{equation}
		r^*_t = \bar{r} + \phi_b \left( \exp(\bar{B}^* - B^*_t) - 1 \right)
		\label{eq:interest_premium}
	\end{equation}

	With $\bar{r} = 1/\beta - 1$ and $\bar{B}^* = 0$, this premium stabilizes the net foreign asset process and resolves the unit-root problem standard in small open economy models.""",
    "METH-10"
)

# ============================================================
# CAL-6
# ============================================================
content, results["CAL-6"] = slice_replace(
    content,
    "In particular, under the interpretation in Equation (\\ref{eq:utilization}),",
    "the structural necessity of the transition.",
    r"	The model's 97\% steady-state utilization implies 3\% energy not served, within EVN's reported 3--5\% range. Qualitative results are robust for $\psi \in [1.5, 3.0]$ and $\mu \in [0.20, 0.35]$.",
    "CAL-6"
)

# ============================================================
# CAL-7
# ============================================================
content, results["CAL-7"] = slice_replace(
    content,
    "The technology adoption elasticity is set to $\\eta_{bat} = 0.10$ (quarterly). This parameter is anchored",
    "place $\\eta_{bat} = 0.10$ on solid empirical footing for both dimensions of the mechanism.",
    r"	The technology adoption elasticity is set to $\eta_{bat} = 0.10$ (quarterly), within the 0.05--0.35 induced-innovation elasticity range estimated by \citet{popp2002} and consistent with the 18--22\% battery experience curve of \citet{kittner2017}. A sustained 10\% reliability gap in the model yields approximately 1\% quarterly battery efficiency improvement, matching the secular pace of cost reduction documented in \citet{way2022}.",
    "CAL-7"
)

# ============================================================
# BASE-1
# ============================================================
content, results["BASE-1"] = slice_replace(
    content,
    "On impact, grid utilization (reliability) declines by 1.30\\% relative to steady state. While a single shock may appear transitory,",
    "the reliability penalty persists until the underlying flexibility deficit is addressed through investment.",
    r"On impact, grid utilization (reliability) declines by 1.30\% relative to steady state, which directly reduces the effective capital services available for production ($V_t = Z(u_t K_{p,t-1})^\alpha L_t^{1-\alpha}$) and generates an output loss of approximately 0.54\%---as shown in the top-left panel of Figure~\ref{fig:irf_baseline}. This \textit{reliability penalty} is not a standard demand shock but a binding physical constraint: when renewable generation falls unexpectedly and storage is insufficient, installed capital cannot be fully operated. The penalty operates through two channels---direct reduction in effective capital services and an implicit cost wedge from backup generation---and persists until the underlying flexibility deficit is addressed through investment.",
    "BASE-1"
)

# ============================================================
# BASE-2
# ============================================================
content, results["BASE-2"] = slice_replace(
    content,
    "The variance decomposition reveals that battery investment volatility is 100\\% driven by global battery price shocks",
    "improving welfare outcomes.",
    r"	Battery investment volatility is 100\% driven by global battery price shocks, while grid investment follows the public reliability-targeting rule (equation~\ref{eq:grid_investment_rule}) with parameter $\phi_{grid} = 1.5$. The proportional tax $\tau_{ss} = 1.43\%$ is fixed: its effect operates through the steady-state channel (reducing $K_{p,ss}$ to 9.83 vs.\ 9.97 under lump-sum taxation), not through dynamic amplification. The equal proportional response of battery and grid investment (1.95\% each) masks an important asymmetry in levels: grid investment rises 2.5$\times$ more in absolute units due to its larger steady-state base.",
    "BASE-2"
)

# ============================================================
# BASE-3
# ============================================================
content, results["BASE-3"] = slice_replace(
    content,
    "The reliability gap also activates the model's scarcity-induced technological change mechanism. By quarter 10,",
    "this innovation channel collapses entirely.",
    r"The reliability gap activates scarcity-induced technological change: by quarter 10, battery technology ($A_{bat}$) has improved by approximately 0.65\% above baseline, with the innovation governed by $A_{bat,t} = A_{bat,t-1}[1 + \eta_{bat}\chi(u^*-u_t)/u^*]$. The one-period delay means improvements only enter flexibility from $t+1$ onward. The effect is cumulative: the 1000-quarter ergodic mean of $A_{bat}$ reaches 1.015 (autocorrelation 0.997), a 1.5\% endogenous technology gain. When regulatory suppression drives $\chi \to 0$, this innovation channel collapses entirely.",
    "BASE-3"
)

# ============================================================
# BASE-5
# ============================================================
content, results["BASE-5"] = slice_replace(
    content,
    "\\paragraph{(1) Front-Loaded Renewable Deployment.}",
    "while the complementary flexibility infrastructure---both public grid and private storage---accumulates incrementally.",
    r"""		Three structural features produce the valley. First, Vietnam's PDP8 front-loads renewable deployment (the Feed-in Tariff boom nearly doubled solar capacity in 18 months), so intermittency rises faster than flexibility can respond. Second, public grid infrastructure follows a bureaucratic investment rule with 5--7 year commissioning lags: even aggressive fiscal responses cannot offset the accumulated deficit within the transition window. Third, battery deployment, though more agile, is bounded by global supply chains and financing constraints (calibrated at 2.0\% quarterly growth); and because storage and transmission are complements ($\rho = 0.4$), neither alone can substitute for the other. The valley thus represents a predictable consequence of unbalanced transition planning.""",
    "BASE-5"
)

# ============================================================
# Write back
# ============================================================
with open(FILE, "wb") as f:
    f.write(content.encode("utf-8"))

successes = sum(1 for v in results.values() if v == "SUCCESS")
failures = sum(1 for v in results.values() if v != "SUCCESS")

print(f"\n{'='*60}")
print(f"RESULTS: {successes} succeeded, {failures} failed")
print(f"{'='*60}")
for k, v in results.items():
    status = "OK" if v == "SUCCESS" else "FAIL"
    print(f"  [{status}] {k}: {v}")

file_size = os.path.getsize(FILE)
print(f"\nNew file size: {file_size:,} bytes ({file_size/1024:.1f} KB)")
