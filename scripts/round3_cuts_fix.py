#!/usr/bin/env python3
"""
round3_cuts_fix.py  — Fix remaining failed cuts from round3_cuts.py
"""

FILE = r"C:\Users\Laptop K1\OneDrive\Documents\MATLAB\main_body.tex"

def slice_replace(content, start_anchor, end_anchor, new_text, cut_name=""):
    try:
        start = content.index(start_anchor)
    except ValueError:
        print(f"FAIL [{cut_name}]: start anchor not found: {start_anchor[:80]!r}")
        return content, False
    try:
        end = content.index(end_anchor, start) + len(end_anchor)
    except ValueError:
        print(f"FAIL [{cut_name}]: end anchor not found: {end_anchor[:80]!r}")
        return content, False
    content = content[:start] + new_text + content[end:]
    print(f"OK   [{cut_name}]")
    return content, True

with open(FILE, "r", encoding="utf-8") as f:
    content = f.read()

original_len = len(content)
ok_count = 0
fail_count = 0

# ============================================================
# FIX 3-D1: The cut broke the equation - the begin{equation} + Y_t line
# is still there but missing the intro. Fix: add the missing \begin{equation}\n\tY_t =
# ============================================================
content, ok = slice_replace(
    content,
    "Output equals domestic absorption plus net capital outflows: + P_{bat,t} I_{bat,t} + I_{grid,t} + \\underbrace{(B^*_t - (1+r^*_{t-1})B^*_{t-1})}_{\\text{net capital outflow}}\n\t\t\\label{eq:goods_market_clearing}\n\t\\end{equation}",
    "\n\t\n\t\\noindent This is identical to equation~\\eqref{eq:current_account}",
    "Output equals domestic absorption plus net capital outflows:\n\t\\begin{equation}\n\t\tY_t = C_t + I_{p,t} + P_{bat,t} I_{bat,t} + I_{grid,t} + \\underbrace{(B^*_t - (1+r^*_{t-1})B^*_{t-1})}_{\\text{net capital outflow}}\n\t\t\\label{eq:goods_market_clearing}\n\t\\end{equation}",
    "FIX-D1: restore broken equation environment"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-A1: Remove duplicate FOC paragraph (current state: single \t indent)
# ============================================================
content, ok = slice_replace(
    content,
    "The household chooses $\\{C_t, L_t, K_{p,t+1}, B^*_t\\}$ to maximize utility subject to the budget constraint \\eqref{eq:hh_budget} and the standard capital accumulation law \\eqref{eq:productive_capital}; the Lagrangian derivation is in Appendix~\\ref{app:deriv_euler}.\n\t\n\t\n\\paragraph{Euler equation for capital.}",
    "\\paragraph{Euler equation for capital.}",
    "\\paragraph{Euler equation for capital.}",
    "3-A1: remove duplicate FOC paragraph"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-A2: Trim budget constraint footnote prose
# ============================================================
content, ok = slice_replace(
    content,
    "The household rents out the capital stock carried into the period, $K_{p,t-1}$, so that capital income is $R_t K_{p,t-1}$. From the firm's profit maximization, $w_t = (1-\\alpha)V_t/L_t$ and $R_t = \\alpha V_t / K_{p,t-1}$, so that total factor income equals $w_t L_t + R_t K_{p,t-1} = (1-\\alpha)V_t + \\alpha V_t = V_t$.\\footnote{\\label{fn:timing_convention}Capital installed by the end of period $t-1$ is productive in period $t$, earning rental rate $R_t$; investment $I_{p,t}$ does not contribute until $t+1$. Consequently, the value-added function uses the predetermined stock $K_{p,t-1}$: $V_t = Z(u_t K_{p,t-1})^\\alpha L_t^{1-\\alpha}$.}",
    "\n\t\n\tThe household chooses $\\{C_t, L_t, K_{p,t+1}, B^*_t\\}$ to maximize utility subject to the budget constraint \\eqref{eq:hh_budget} and the standard capital accumulation law \\eqref{eq:productive_capital}; the Lagrangian derivation is in Appendix~\\ref{app:deriv_euler}.",
    "The household earns labour income $(1-\\tau)w_t L_t$ and capital rental income $R_t K_{p,t-1}$, invests in productive capital $I_{p,t}$ and battery capital $P_{bat,t}I_{bat,t}$, and saves via foreign bonds $B^*_t$.\\footnote{\\label{fn:timing_convention}Capital installed by the end of period $t-1$ is productive in period $t$; investment $I_{p,t}$ does not contribute until $t+1$.}",
    "3-A2: trim budget constraint prose"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-B1: Cut the long CES explanation
# ============================================================
content, ok = slice_replace(
    content,
    "Output is a weighted aggregate of value-added $V_t$ and energy $E_t$, where neither input can fully substitute for the other when $\\sigma_E < 1$: shutting down either drives output toward zero, making energy and productive activity essential complements.\n\t\n\tThis formulation follows the environmental macroeconomics literature \\citep{heutel2012, airaudo2025green} and captures the limited short-run substitutability between energy and non-energy inputs. When the elasticity of substitution $\\sigma_E$ is low (calibrated at 0.6), disruptions to energy supply translate into disproportionately large output losses, a feature that is empirically relevant for energy-intensive emerging economies.\n\n\tValue-added production combines labor and effective capital services:",
    "\n\n\t\\begin{equation}\n\t\tV_t = Z\\,(u_t K_{p,t-1})^\\alpha L_t^{1-\\alpha}",
    "Output is a CES aggregate of value-added $V_t$ and fixed energy endowment $\\bar{E}$; when $\\sigma_E = 0.6 < 1$, energy and productive activity are gross complements \\citep{heutel2012}.\n\n\tValue-added production combines labor and effective capital services:",
    "3-B1: cut CES explanation paragraph"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-B2: Trim verbose value-added explanation
# ============================================================
content, ok = slice_replace(
    content,
    "Value-added is produced from ``effective capital'' $u_t K_{p,t-1}$---the fraction of \\textit{previously installed} capacity that is actually operational---combined with labor. A factory with machines that cannot run due to power outages produces less, even if the machines physically exist. The timing $K_{p,t-1}$ reflects the standard capital-in-advance convention: capital installed by the end of period $t-1$ is the stock available for production in period $t$ (see footnote~\\ref{fn:timing_convention}).\n\t\n\tHere, $K_{p,t-1}$ is the physical stock of productive capital carried into period $t$, and $u_t \\in (0,1)$ represents the fraction of that stock that can be operated given grid conditions. Unlike standard utilization models (e.g., \\citet{greenwood1988}), $u_t$ is not chosen optimally by firms and does not influence capital depreciation. Instead, it captures system-wide outages, curtailment, and power shortages that prevent installed capital from being used. In this sense, grid instability manifests as an endogenous, state-dependent productivity wedge: a 3\\% reliability shortfall ($u_t = 0.97$) is equivalent to a TFP loss of approximately $\\alpha \\times 3\\% \\approx 1\\%$.",
    "\n\t\n\tIn the implemented model, energy enters the CES aggregator",
    "Value-added uses ``effective capital'' $u_t K_{p,t-1}$: the share of installed capacity operational under current grid conditions. A 3\\% reliability shortfall acts as a $\\sim$1\\% TFP loss.",
    "3-B2: trim value-added explanation"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-B3: Remove "In the implemented model..." and "The household's optimality..." paragraphs
# ============================================================
content, ok = slice_replace(
    content,
    "\n\t\n\tIn the implemented model, energy enters the CES aggregator as a fixed endowment $\\bar{E}$, reflecting the assumption that total energy supply is determined by installed capacity and policy targets (PDP8) rather than optimized period-by-period. This simplification allows the model to focus on the key channel of interest---how renewable intermittency affects grid reliability and, through it, effective capital utilization---without requiring a full energy market clearing block. The intermittency shock enters through the reliability function (equation \\eqref{eq:utilization} below) rather than through energy supply directly.\n\n\tThe household's optimality conditions---the capital Euler equation \\eqref{eq:euler} and labor supply condition \\eqref{eq:labor_supply}, as derived in Section~\\ref{sec:households}---together with the production technology determine the equilibrium allocation of capital and labor.",
    "\n\n\t\\subsection{Grid Reliability and Flexibility Assets}",
    "\n\n\t\\subsection{Grid Reliability and Flexibility Assets}",
    "3-B3: remove energy endowment and FOC paragraphs"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-C2: Trim post-itemize explanation (after C1 already ran)
# ============================================================
content, ok = slice_replace(
    content,
    "Technology growth is zero when reliability is at target and accelerates proportionally to the shortfall. The compounding structure captures path-dependence and saturation documented in battery experience curves \\citep{kittner2017, way2022}.",
    "\n\n\t\\textbf{Timing and causal ordering.}",
    "Technology growth is zero when reliability is at target and accelerates proportionally to the shortfall \\citep{kittner2017, way2022}.",
    "3-C2: trim duplicate sentence"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 3-C3: Remove the long DTC comparison paragraph (still present)
# ============================================================
content, ok = slice_replace(
    content,
    "This formulation differs from Directed Technical Change in three key respects. First, there is no R\\&D sector or explicit allocation of scientists. Second, innovation is directed by the \\textit{reliability gap}---a physical grid signal---rather than by inter-sectoral price competition. Third, technological progress compounds from an inherited knowledge base, capturing the path-dependence and saturation properties of empirical battery learning curves \\citep{kittner2017, way2022}. The $\\chi$ parameter formalizes the empirical finding that policy-induced scarcity signals---such as carbon markets or capacity payments---are necessary to fully activate the innovation response \\citep{calel2016, aghion2016}. Underlying all three investment and innovation channels is a single structural property of the exponential utilization function: its marginal effect of flexibility on reliability,\n\t\\begin{equation}\n\t\t\\frac{\\partial u_t}{\\partial F_t} = (1 - u_t) \\cdot \\frac{\\psi}{\\phi_{int} \\cdot \\text{Vol}_{ren,t}},\n\t\\end{equation}\n\tis proportional to the current reliability shortfall $(1-u_t)$. When reliability is already high, additional flexibility yields diminishing returns; when the grid is stressed, flexibility becomes critically valuable. The battery investment rule \\eqref{eq:battery_investment} and the innovation equation \\eqref{eq:learning_by_doing} both respond to this scarcity through the ratio $u^*/u_t$, which amplifies each response precisely when the shortfall is largest.",
    "\n\n\n\t\\subsection{Government} \\label{sec:government}",
    "\n\n\t\\subsection{Government} \\label{sec:government}",
    "3-C3: remove DTC comparison and structural-property paragraphs"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 4-A: Replace verbose transmission channels with summary
# (note: "subsection" missing backslash in file)
# ============================================================
content, ok = slice_replace(
    content,
    "subsection{Transmission Channels} \\label{subsec:transmission}\n\t\n\tThe linearized system reveals four distinct transmission channels through which intermittency shocks affect the macroeconomy. The first three are present in the lump-sum baseline; the fourth is unique to the proportional tax variant. These channels interact dynamically, and their relative importance varies with the structural parameters of the model.",
    "\n\n\t\t\\paragraph{Channel 1: The Reliability Wedge.}",
    "\\subsection{Transmission Channels} \\label{subsec:transmission}\nSee Appendix~\\ref{app:log_lin} for the complete log-linearized system. The four channels are summarised in Section~\\ref{sec:methodology} (Mechanism Summary); their quantitative importance is documented in the FEVD (Table~\\ref{tab:fevd}).",
    "4-A: replace verbose transmission channels"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 4-A2: Remove the remaining Channel paragraphs
# ============================================================
content, ok = slice_replace(
    content,
    "\n\n\t\t\\paragraph{Channel 1: The Reliability Wedge.}\n\t\tAn adverse intermittency shock ($\\varepsilon^{ren}_t > 0$) directly reduces utilization through Equation~\\eqref{eq:reliability_lin}. Via Equation~\\eqref{eq:va_lin}, this lowers value-added with elasticity $\\alpha$. Because energy enters the CES aggregator as a complement to value-added ($\\sigma_E < 1$), the output contraction is amplified through the nested production structure. The immediate impact on output is:\n\t\t\\begin{equation}\n\t\t\t\\frac{\\partial \\hat{Y}_t}{\\partial \\varepsilon^{ren}_t} \\approx -\\alpha \\xi s_V^{(\\sigma_E-1)/\\sigma_E}(1-\\omega_E) < 0\n\t\t\\end{equation}\n\t\t\n\t\tThis is the ``greenflation'' channel: supply-side volatility from renewables generates simultaneous output losses and energy price increases. The mechanism is structurally analogous to the adverse supply shocks analyzed in \\citet{blanchard2007}, who show that the macroeconomic consequences of energy supply disruptions depend critically on the degree of real wage rigidity and the substitutability between energy and other inputs. In the present model, the low calibrated value of $\\sigma_E$ ensures that intermittency shocks behave as ``unfavorable supply shocks'' in the Blanchard-Gal\\'{i} sense, generating a policy-relevant trade-off between output stabilization and energy price stabilization.\n\n\t\t\\paragraph{Channel 2: The Innovation Response.}\n\t\tThe reliability decline widens the gap $(u^* - u_t)/u^*$, which accelerates scarcity-induced technological change in storage (Equation~\\ref{eq:abat_lin}). The resulting improvement in $A_{bat,t}$ increases the effective services per unit of battery capital in period $t+1$, gradually restoring flexibility and utilization. Note the one-period delay: innovation induced by the reliability shortfall in period $t$ raises $A_{bat,t}$, which enters $F_{t+1}$ as $A_{bat,t-1}$ and improves reliability only from period $t+1$ onward. The speed of this response depends critically on three parameters: the innovation elasticity $\\eta_{bat}$, the regulatory transmission coefficient $\\chi$, and the battery share in flexibility $s_b$. When any of these is small, the innovation channel is weak, and adjustment relies on costly capital accumulation rather than efficient technological improvement.\n\n\n\t\t\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}\n\t\tThe government investment rule (Equation~\\ref{eq:gridrule_lin}) responds to reliability shortfalls, but with two sources of delay. First, the investment rule operates with a one-period implementation lag inherent in the capital accumulation equation. Second, institutional inertia (captured by $\\phi_{grid}$) limits the magnitude of the fiscal response.\n\n\t\t\t\\paragraph{Channel 4: The Fiscal Wedge Channel.}\n\t\tThis channel identifies the level distortion introduced by proportional taxation, which is absent in the lump-sum baseline. While the tax rate $\\tau$ remains fixed during the business cycle, it imposes a permanent wedge on the incentives to work and invest (equations~\\ref{eq:euler_lin}--\\ref{eq:labor_lin}). The steady-state capital stock $K_{p,ss}$ is approximately 1.4\\% lower in the proportional tax variant compared to the lump-sum benchmark, reflecting the lower after-tax marginal product of capital. This means that for a given percentage intermittency shock, the economy operates from a more constrained initial condition, making each percentage decline in reliability costlier in terms of final output units. This channel implies that proportional taxation reduces the economy's structural resilience to intermittency shocks even without providing dynamic amplification.\n\t\t\n\t\t\\paragraph{Interaction: The ``Reliability Valley.''}\n\t\tThe interaction of Channels 1--4 generates the ``Reliability Valley'' identified in the abstract. During a transition to high renewable penetration, the deterministic path of increasing intermittency outpaces both the sluggish accumulation of public grid capital (Channel 3) and the gradual improvement in storage technology (Channel 2). The resulting reliability shortfall persists for several years, during which households bear welfare losses through reduced consumption and increased energy price volatility. The depth and duration of the valley depend on the relative speed of private innovation versus public infrastructure accumulation---the ``Agility Gap.''",
    "\n\n\n\t\\section{Calibration}",
    "",
    "4-A2: remove Channel 1-4 paragraphs"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 4-B: Replace verbose Equilibrium Definition
# ============================================================
content, ok = slice_replace(
    content,
    "\\subsection{Equilibrium Definition}\n\n\tThe log-linearized equilibrium conditions are the linearized counterparts of the 16-equation system defined in Section~\\ref{sec:methodology}; the full system is collected in Appendix~\\ref{app:log_lin}. The linearized system satisfies the Blanchard-Kahn conditions \\citep{blanchard1980} for a unique saddle-path stable solution, with 3 eigenvalues outside the unit circle matching the 3 forward-looking variables ($C_t$, $Y_t$, and $r^*_t$).\\footnote{$C_t$ and $r^*_t$ are forward-looking because they appear with a lead in the capital Euler and bond Euler equations respectively. $Y_t$ enters the capital Euler as $Y_{t+1}$, making it a jumper in Dynare\\'s eigenvalue classification. The six predetermined state variables are $K_{p,t-1}$, $K_{b,t-1}$, $K_{g,t-1}$, $A_{bat,t-1}$, $P_{bat,t-1}$, and $B^*_{t-1}$.} Seven state variables and 6 static variables complete the system.",
    "\n\n\tsubsection{Transmission Channels} \\label{subsec:transmission}",
    "\\subsection{Equilibrium Definition}\nThe recursive competitive equilibrium is formally stated in Appendix~\\ref{app:equilibrium}; the Blanchard-Kahn conditions are satisfied with 3 unstable eigenvalues matching the 3 forward-looking variables ($C_t$, $Y_t$, $r^*_t$).",
    "4-B: replace equilibrium definition with cross-reference"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 5-A1: Trim external sector regression narrative
# ============================================================
content, ok = slice_replace(
    content,
    "The debt-elastic premium follows the SGU stationarity device. Rather than adopting the SGU-standard value of $0.001$ (calibrated to a generic small open economy), I discipline the curvature using Vietnamese data as an order-of-magnitude anchor. Let $r_t$ denote Vietnam's observed annual real interest rate (World Bank, 1993--2023) and let $d_t$ denote gross external debt as a share of GNI (World Bank series DT.DOD.DECT.GN.ZS). I construct the demeaned spread and debt terms $y_t = r_t - \\bar{r}$ and $x_t = \\exp(d_t - \\bar{d}) - 1$, where bars denote sample means over the 28 overlapping annual observations (excluding years with missing data). The premium curvature is then estimated by the no-intercept least-squares slope\n\t\t\\[\n\t\t\\hat{\\phi}_b = \\max\\!\\left(0, \\frac{\\sum_t x_t y_t}{\\sum_t x_t^2}\\right),\n\t\t\\]\n\t\tyielding $\\hat{\\phi}_b = 0.03866$, which I round to $\\phi_b = 0.039$.\n\n\tThe regression provides an order-of-magnitude anchor ($R^2 = 0.16$) that rules out the frictionless extreme ($\\phi_b \\approx 0$) and implausibly large frictions ($\\phi_b \\gg 0.10$); the primary justification is economic: the SGU-standard $\\phi_b = 0.001$ generates an implausible $-21.8\\%$ quarterly investment collapse while $\\phi_b = 0.039$ yields the credible $-3.5\\%$ response consistent with Vietnamese national accounts data (see Table~\\ref{tab:phi_b_sensitivity}).",
    "\n\n\t\\begin{table}[H]\n\t\t\\centering\n\t\t\\caption[Structural Parameter Calibration",
    "Following \\citet{schmittgrohe2003}, I calibrate $\\phi_b$ from a no-intercept regression of demeaned Vietnamese real interest rates on demeaned external-debt/GNI (World Bank 1993--2023), yielding $\\hat{\\phi}_b = 0.039$; the SGU-standard $\\phi_b = 0.001$ implies an implausible $-21.8\\%$ investment collapse (Table~\\ref{tab:phi_b_sensitivity}).",
    "5-A1: trim external sector regression narrative"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 5-B1: Trim calibration validation intro
# ============================================================
content, ok = slice_replace(
    content,
    "To validate the calibration, I check that the steady state matches key Vietnamese macroeconomic moments. Table~\\ref{tab:validation} distinguishes between \\textit{calibration targets} (moments directly used to pin parameters) and \\textit{free checks} (moments not targeted during calibration but implied by the model's structure):",
    "\n\n\t\\begin{table}[H]\n\t\t\\centering\n\t\t\\caption[Calibration Validation",
    "Table~\\ref{tab:validation} validates the calibration against key Vietnamese moments, distinguishing calibration targets (T) from free checks (F):",
    "5-B1: trim calibration validation intro"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 5-C: Trim the calibration targets paragraph (has a literal \t prefix)
# ============================================================
content, ok = slice_replace(
    content,
    "\\\tThe calibration targets are closely matched: $\\beta = 0.99$ reflects Vietnam's 3.14\\% mean real rate (World Bank 2000--2023); $\\alpha = 0.35$ is the residual capital share given a labor income share of 0.65 (GSO); and $\\phi_b = 0.039$ is the data-anchored debt-elastic premium documented in Section~\\ref{subsec:external_sector}, with sensitivity reported in Table~\\ref{tab:phi_b_sensitivity} below.",
    "\n\n\t\t\\begin{table}[H]\n\t\t\t\\centering\n\t\t\t\\caption[Sensitivity to Debt-Elastic Premium",
    "",
    "5-C: trim calibration targets prose paragraph"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 6-A: Replace large steady-state table with compact version
# ============================================================
content, ok = slice_replace(
    content,
    "Table \\ref{tab:steady_state} reports the model's steady state. The calibration targets a baseline reliability rate of 97.0\\%, consistent with Vietnam's official grid performance objective under PDP8. The battery-to-grid capital ratio of 0.167 reflects the structural imbalance in Vietnam's current investment allocation: approximately \\$2.4 billion earmarked for Battery Energy Storage Systems (BESS) versus \\$14.9 billion for transmission and grid infrastructure over the 2021--2030 period.\n\n\t\t\\begin{table}[H]\n\t\t\t\\centering\n\t\t\t\\caption[Model Steady State]{\\textbf{Model Steady State}}\n\t\t\t\\label{tab:steady_state}\n\t\t\t\\renewcommand{\\arraystretch}{1.3}\n\t\t\t\\begin{tabular}{l c p{6.5cm}}\n\t\t\t\t\\toprule\n\t\t\t\t\\textbf{Variable} & \\textbf{Value} & \\textbf{Interpretation} \\\\\n\t\t\t\t\\midrule\n\t\t\t\tOutput ($Y$) & 1.000 & Single output normalization: $Y_{ss}=1$ enforced by TFP scale $Z=1.131$ \\\\\n\t\t\t\tConsumption ($C$) & 0.734 & 73.4\\% of output \\\\\n\t\t\t\tLabor ($L$) & 0.330 & Calibration target: $\\tfrac{1}{3}$ of unit time endowment \\\\\n\t\t\t\tProductive Capital ($K_p$) & 9.829 & Production capacity \\\\\n\t\t\t\tBattery Capital ($K_b$) & 0.190 & Private storage capacity \\\\\n\t\t\t\tGrid Capital ($K_g$) & 1.140 & Public transmission infrastructure \\\\\n\t\t\t\tValue Added ($V$) & 1.211 & $Z\\,(u \\cdot K_p)^\\alpha L^{1-\\alpha}$ \\\\\n\t\t\t\tUtilization ($u$) & 0.970 & Grid reliability (97.0\\%) \\\\\n\t\t\t\tFlexibility ($F$) & 0.526 & Aggregate smoothing capacity \\\\\n\t\t\t\tBattery Technology ($A_{bat}$) & 1.000 & Initial technology index; free scale choice (only growth rate $\\eta_{bat}$ has economic content)$^{\\ddagger}$ \\\\\n\t\t\t\tBattery Price ($P_{bat}$) & 1.000 & AR(1) unconditional mean; $\\mathbb{E}[\\log P_{bat}]=0$ structurally pins $P_{bat,ss}=1$ (not a free choice)$^{\\ddagger}$ \\\\\n\t\t\t\tNet Foreign Assets ($B^*$) & 0.000 & Zero NFA position; consistent with Vietnam's mean CA/GDP $\\approx -0.03\\%$ over 2000--2024 \\\\\n\t\t\t\tWorld Interest Rate ($r^*$) & 0.0101 & $= 1/\\beta - 1$ (quarterly) \\\\\n\t\t\t\t\\midrule\n\t\t\t\t\\multicolumn{3}{l}{\\textit{Steady-State Investment Flows}} \\\\\n\t\t\t\tProductive Investment ($I_p$) & 0.246 & $= \\delta_p \\cdot K_p$ \\\\\n\t\t\t\tBattery Investment ($I_{bat}$) & 0.0057 & $= \\delta_b \\cdot K_b$ \\\\\n\t\t\t\tGrid Investment ($I_{grid}$) & 0.0143 & $= \\delta_g \\cdot K_g$ \\\\\n\t\t\t\tTax Rate ($\\tau$) & 0.01425 & Grid investment / output ($=\\delta_g K_{g,ss}$ since $Y_{ss}=1$) \\\\\n\t\t\t\t\\bottomrule\n\t\t\t\\end{tabular}\n\t\t\t\\vspace{0.2cm}\n\t\t\t\\parbox{\\textwidth}{\\footnotesize \\textbf{Note:} Initial flexibility ratios are derived from the Power Development Plan 8 (PDP8) investment targets. The consumption share of 73.4\\% and the zero net foreign asset position ($B^*_{ss} = 0$) are standard for the small open economy steady state. The labor disutility weight $\\varphi$ is calibrated endogenously such that households supply $L_{ss} = 0.33$ units of labor. The steady-state tax rate $\\tau_{ss}=1.43\\%$ is pinned by the grid maintenance requirement $\\delta_g K_{g,ss}$ relative to output $Y_{ss}=1$. $^\\ddagger$The technology index $A_{bat} = 1$ is a normalization, while $P_{bat} = 1$ is the structural mean of the price process.}\n\t\t\\end{table}",
    "\n\n\t\t\\subsubsection{The Efficiency Wedge of Renewable Reliability}",
    "\\begin{table}[H]\n\t\t\\centering\n\t\t\\caption[Key Steady-State Values]{\\textbf{Key Steady-State Values}}\n\t\t\\label{tab:steady_state}\n\t\t\\renewcommand{\\arraystretch}{1.2}\n\t\t\\small\n\t\t\\begin{tabular}{lclc}\n\t\t\t\\toprule\n\t\t\t\\textbf{Variable} & \\textbf{Value} & \\textbf{Variable} & \\textbf{Value} \\\\\n\t\t\t\\midrule\n\t\t\tOutput ($Y_{ss}$) & 1.000 & Utilization ($u_{ss}$) & 0.970 \\\\\n\t\t\tConsumption share & 73.4\\% & Flexibility ($F_{ss}$) & 0.526 \\\\\n\t\t\tLabor ($L_{ss}$) & 0.330 & $K_b/K_g$ ratio & 0.167 \\\\\n\t\t\tProd.\\ capital ($K_p$) & 9.829 & World rate ($r^*_{ss}$) & 0.0101 \\\\\n\t\t\tTax rate ($\\tau_{ss}$) & 1.43\\% & Net foreign assets & 0.000 \\\\\n\t\t\t\\bottomrule\n\t\t\\end{tabular}\n\t\t\\vspace{0.1cm}\n\t\t\\parbox{\\textwidth}{\\footnotesize \\textbf{Note:} Full 16-variable steady state in Appendix~\\ref{app:steadystate}. All calibration targets are met: $u_{ss}=0.97$, $L_{ss}=0.33$, $I/Y=26.6\\%$.}\n\t\t\\end{table}\n\n\t\t\\subsubsection{The Efficiency Wedge of Renewable Reliability}",
    "6-A: replace large SS table with compact version"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 6-B: Trim FEVD intro paragraph
# ============================================================
content, ok = slice_replace(
    content,
    "The variance decomposition identifies the dominant sources of volatility in each endogenous variable. Table~\\ref{tab:fevd} reports the contribution of each shock at the 40-quarter horizon.",
    "\n\n\t\t\\begin{table}[H]\n\t\t\t\\centering\n\t\t\t\\caption[Forecast Error Variance Decomposition",
    "Table~\\ref{tab:fevd} reports each shock's contribution to forecast error variance at the 40-quarter horizon.",
    "6-B: trim FEVD intro paragraph"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 8-F: Remove "Second" bold finding from Section 7 summary (single \t indent)
# ============================================================
content, ok = slice_replace(
    content,
    "\n\n\t\\textbf{Second, battery investment is dominated by global price dynamics.} Variance decomposition reveals that 100\\% of battery investment volatility is driven by global battery price shocks, while intermittency accounts for only 1\\%. This highlights Vietnam's price-taking position in global battery markets and suggests that domestic reliability policies alone cannot determine the pace of storage deployment.",
    "\n\n\t\\textbf{Third, scarcity-induced technological change provides",
    "\n\n\t\\textbf{Third, scarcity-induced technological change provides",
    "8-F: remove duplicate Second finding from Section 7 summary"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# 8-F2: Remove "Fourth" bold finding from Section 7 summary
# ============================================================
content, ok = slice_replace(
    content,
    "\n\n\t\\textbf{Fourth, the model's transmission channels are quantitatively validated.} The Dynare simulation confirms expected signs: output and reliability are highly correlated (0.986), and battery price shocks affect investment negatively ($\\rho = -0.979$ between $P_{bat}$ and $I_{bat}$). The proportional tax is fixed, so its role appears entirely through steady-state comparisons: the constant wedge $\\tau_{ss}=1.43\\%$ shifts equilibrium quantities, not their dynamic comovement. These correlations validate the model's internal consistency.",
    "\n\n\t\n\n\t\\section{Conclusion}",
    "\n\n\t\\section{Conclusion}",
    "8-F2: remove Fourth finding from Section 7 summary"
)
ok_count += ok; fail_count += (not ok)

# ============================================================
# Write output
# ============================================================
with open(FILE, "w", encoding="utf-8", newline="\n") as f:
    f.write(content)

new_len = len(content)
print(f"\nDone. {ok_count} cuts OK, {fail_count} cuts FAILED.")
print(f"File size: {original_len} -> {new_len} chars (removed {original_len - new_len} chars, "
      f"{(original_len - new_len)/original_len*100:.1f}%)")
