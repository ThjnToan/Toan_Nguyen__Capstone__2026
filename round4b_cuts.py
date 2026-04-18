content = open('main_body.tex', encoding='utf-8').read()
original_len = len(content)
applied = 0; failed = []

def S(label, start, end, new):
    global content, applied
    try:
        s = content.index(start)
        e = content.index(end, s) + len(end)
        content = content[:s] + new + content[e:]
        applied += 1; print(f'  OK: {label}')
    except ValueError as ex:
        failed.append(label); print(f'  FAIL: {label} - {ex}')

# FIX: dangling colon in government section
S('gov-fix dangling colon',
    ', calibrated to balance the steady-state budget:\n\n\t\\subsection{Battery Price Dynamics}',
    '\n\n\t\\subsection{Battery Price Dynamics}',
    ' ($\\tau_{ss} = \\delta_g K_{g,ss}/Y_{ss} = 0.0143$), so tax revenue exactly covers replacement grid investment.\n\n\t\\subsection{Battery Price Dynamics}')

# CUT 3-3: Vol_ren - 3 paragraphs -> 2 sentences
S('3-3 Vol_ren paragraphs',
    'where $\\theta_{ren}$ is the renewable penetration share (dimensionless), $K_{ren}$ is installed renewable capacity (GW), and $\\sigma_{ren}$ is the coefficient of variation of the intermittency shock (dimensionless). As defined, $\\overline{\\text{Vol}}_{ren}$ carries units of GW. Dimensional consistency in the reliability function \\eqref{eq:utilization} is preserved by the scaling parameter $\\phi_{int}$, which has units of GW per unit of flexibility stock, so that the ratio $F_t / (\\phi_{int} \\cdot \\text{Vol}_{ren,t})$ is dimensionless as required by the exponential.',
    'is allowed to increase along the policy-determined VRE expansion path.',
    'where $\\theta_{ren}$, $K_{ren}$, and $\\sigma_{ren}$ are the penetration share, installed capacity, and CV of the intermittency shock. The shock $\\varepsilon_{ren,t}$ enters multiplicatively as $\\overline{\\text{Vol}}_{ren}\\cdot\\exp(\\varepsilon_{ren,t})$; in transition simulations $\\overline{\\text{Vol}}_{ren}$ rises with VRE expansion.')

# CUT 3-4b: Post-flexibility-equation 3 paragraphs -> 1 sentence
S('3-4b flexibility post-eq',
    'Flexibility is a weighted combination of technology-augmented battery capacity ($A_{bat,t-1} K_{b,t-1}$) and grid infrastructure ($K_{g,t-1}$), where neither can fully replace the other. A grid with abundant batteries but no transmission lines cannot deliver stored power; a grid with extensive transmission but no storage cannot smooth intermittent supply. The negative exponent $s_F < 0$ (since $\\rho = 0.4 < 1$) ensures that flexibility collapses if either input approaches zero---a mathematical expression of this engineering constraint.',
    'See Proposition~\\ref{prop:superadditivity} in Appendix~\\ref{app:proof_superadditivity} for the statement and quantitative characterization.',
    'The negative exponent $s_F<0$ ($\\rho=0.4<1$) enforces imperfect complementarity; all inputs are predetermined from period $t-1$. See Proposition~\\ref{prop:superadditivity} (Appendix~\\ref{app:proof_superadditivity}).')

# CUT 3-5b: Battery investment post-eq paragraph (the one that wasn't removed)
S('3-5b battery post-eq paragraph',
    'The exponent $\\phi_{grid}$ governs the elasticity of the investment response to the reliability gap---higher values imply more aggressive battery deployment when the grid is stressed. The price term $1/P_{bat,t}$ reflects the cost side: Vietnam is a price-taker in global battery markets, so higher world battery prices reduce real investment even when the return to storage is high. At steady state, the rule collapses to $I_{bat,ss}/\\bar{I}_{bat} = 1$, consistent with capital replacement at $I_{bat,ss} = \\delta_b K_{b,ss}$.',
    '\n\n\t\\subsection{Scarcity-Induced Technological Change in Battery Storage}',
    'Higher $\\phi_{grid}$ implies more aggressive deployment when the grid is stressed; the price term $1/P_{bat,t}$ reflects Vietnam\'s price-taking position in global markets.\n\n\t\\subsection{Scarcity-Induced Technological Change in Battery Storage}')

# CUT 3-6: SITC itemized + 3 paragraphs -> 1 sentence
S('3-6 SITC itemized + paragraphs',
    '\n\twhere:\n\t\\begin{itemize}\n\t\t\\item $\\eta_{bat} > 0$ is the technology adoption elasticity governing how strongly the reliability scarcity signal translates into innovation\n\t\t\\item $\\chi \\in (0,1]$ captures regulatory distortions in price signals\n\t\t\\item $u^*$ is the target reliability level\n\t\t\\item $u_t$ is the current reliability\n\t\\end{itemize}',
    'Unlike Directed Technical Change, there is no R\\&D sector: innovation compounds from an inherited knowledge base directed by the physical reliability gap, activating the scarcity signals documented in \\citet{calel2016, aghion2016}.',
    '\nwhere $\\eta_{bat}$ is the innovation elasticity, $\\chi\\in(0,1]$ captures regulatory signal attenuation, and $u^*$ is the reliability target. Innovation is zero at $u_t=u^*$, accelerates when reliability falls below target, and enters $F_{t+1}$ with a one-period delay eliminating within-period simultaneity.')

# CUT 3-8: Greenflation paragraph after Ch1 equation
S('3-8 greenflation paragraph',
    'This is the ``greenflation\'\' channel: supply-side volatility from renewables generates simultaneous output losses and energy price increases. The mechanism is structurally analogous to the adverse supply shocks analyzed in \\citet{blanchard2007}, who show that the macroeconomic consequences of energy supply disruptions depend critically on the degree of real wage rigidity and the substitutability between energy and other inputs. In the present model, the low calibrated value of $\\sigma_E$ ensures that intermittency shocks behave as ``unfavorable supply shocks\'\' in the Blanchard-Gal\\\'\\{\\i\\} sense, generating a policy-relevant trade-off between output stabilization and energy price stabilization.',
    '\n\n\t\\paragraph{Channel 2: The Innovation Response.}',
    '\n\n\t\\paragraph{Channel 2: The Innovation Response.}')

# CUT 3-9: Channel 2 long -> 2 sentences
S('3-9 channel2 long paragraph',
    'The reliability decline widens the gap $(u^* - u_t)/u^*$, which accelerates scarcity-induced technological change in storage (Equation~\\ref{eq:abat_lin}). The resulting improvement in $A_{bat,t}$ increases the effective services per unit of battery capital in period $t+1$, gradually restoring flexibility and utilization. Note the one-period delay: innovation induced by the reliability shortfall in period $t$ raises $A_{bat,t}$, which enters $F_{t+1}$ as $A_{bat,t-1}$ and improves reliability only from period $t+1$ onward. The speed of this response depends critically on three parameters: the innovation elasticity $\\eta_{bat}$, the regulatory transmission coefficient $\\chi$, and the battery share in flexibility $s_b$. When any of these is small, the innovation channel is weak, and adjustment relies on costly capital accumulation rather than efficient technological improvement.',
    '\n\n\n\t\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}',
    'The reliability decline widens the gap $(u^*-u_t)/u^*$, accelerating innovation in storage (Eq.~\\ref{eq:abat_lin}) and gradually restoring flexibility. Innovation in period $t$ enters $F_{t+1}$ with a one-period delay.\n\n\t\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}')

# CUT 3-10: Channel 4 long -> 2 sentences
S('3-10 channel4 paragraph',
    'This channel identifies the level distortion introduced by proportional taxation, which is absent in the lump-sum baseline. While the tax rate $\\tau$ remains fixed during the business cycle, it imposes a permanent wedge on the incentives to work and invest (equations~\\ref{eq:euler_lin}--\\ref{eq:labor_lin}). The steady-state capital stock $K_{p,ss}$ is approximately 1.4\\% lower in the proportional tax variant compared to the lump-sum benchmark, reflecting the lower after-tax marginal product of capital. This means that for a given percentage intermittency shock, the economy operates from a more constrained initial condition, making each percentage decline in reliability costlier in terms of final output units. This channel implies that proportional taxation reduces the economy\'s structural resilience to intermittency shocks even without providing dynamic amplification.',
    '\n\n\t\\paragraph{Interaction: The ``Reliability Valley.\'\'}',
    'The fixed tax $\\tau$ imposes a permanent wedge: $K_{p,ss}$ is $\\approx$1.4\\% below the lump-sum benchmark, reducing structural resilience without dynamic amplification.\n\n\t\\paragraph{Interaction: The ``Reliability Valley.\'\'}')

# CUT 5-1: "core contribution" sentence
S('5-1 core contribution',
    '\n\tA core contribution of this capstone is the rigorous calibration of the reliability constraint, which is micro-founded on the physical properties of the Vietnamese grid rather than arbitrary assumptions.\n\n\t\\subsubsection{Renewable Intermittency}',
    '\n\t\\subsubsection{Renewable Intermittency}',
    '\n\t\\subsubsection{Renewable Intermittency}')

# CUT 5-3: sigma_f 2 paragraphs -> 1 sentence
S('5-3 sigma_f paragraphs',
    '\\textbf{2. Elasticity of Substitution ($\\sigma_f = 0.4$):}\n\tStandard macro-energy models often assume high substitutability between energy capital types. However, \\citet{hart2019} and \\citet{hirth2013} demonstrate that transmission and storage are \\textit{imperfect complements} in the presence of intermittency. Transmission moves energy across space, while storage moves energy across time. When solar generation is zero (nighttime), infinite transmission capacity yields zero marginal utility.',
    'This is consistent with \\citet{hart2019}\'s empirical finding that the cross-price elasticity between storage and transmission investment is negative.',
    '\\textbf{2. Elasticity of Substitution ($\\sigma_f = 0.4$):} Battery and grid are gross complements ($\\sigma_f<1$): transmission moves power across space, storage across time, and neither alone delivers system reliability \\citep{hart2019,hirth2013}.')

# CUT 5-5: chi paragraph -> 2 sentences
S('5-5 chi paragraph',
    'The price signal transmission parameter is set to $\\chi = 1.0$ in the baseline (full price transmission). The regulated-market counterfactual uses $\\chi = 0.3$, anchored to the empirical retail electricity price pass-through in Vietnam. Under EVN\'s administered tariff system, retail electricity prices are set by the Ministry of Industry and Trade and updated infrequently; IEEFA (2022) documents that Vietnam\'s retail tariffs have been held below cost-recovery levels for extended periods, with effective pass-through of wholesale scarcity cost signals estimated at approximately 25--35\\% \\citep{ieefa2022}. The value $\\chi = 0.3$ is calibrated to the midpoint of this range, representing a regime in which only 30 cents of every dollar of reliability-gap signal reaches the innovating sector --- consistent with Vietnam\'s fixed-tariff electricity market structure prior to the 2023 tariff reform.',
    '\n\n\t\\subsection{Technology Adoption Elasticity and Innovation}',
    'The baseline sets $\\chi=1.0$ (full signal transmission); the regulated counterfactual uses $\\chi=0.3$, calibrated to Vietnam\'s documented retail tariff pass-through of 25--35\\% \\citep{ieefa2022}.\n\n\t\\subsection{Technology Adoption Elasticity and Innovation}')

open('main_body.tex', 'w', encoding='utf-8').write(content)
print(f'\nApplied {applied}/{applied+len(failed)} cuts')
if failed: print('Failed:', failed)
print(f'Size: {original_len} -> {len(content)} bytes')
