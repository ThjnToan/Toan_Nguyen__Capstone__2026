import re, sys

content = open('main_body.tex', encoding='utf-8').read()
original_len = len(content)
applied = 0
failed = []

def slice_replace(label, start_anchor, end_anchor, new_text):
    global content, applied
    try:
        s = content.index(start_anchor)
        e = content.index(end_anchor, s) + len(end_anchor)
        content = content[:s] + new_text + content[e:]
        applied += 1
        print(f'  OK: {label}')
    except ValueError as ex:
        failed.append(label)
        print(f'  FAIL: {label} — {ex}')

def simple_replace(label, old, new):
    global content, applied
    if old in content:
        content = content.replace(old, new, 1)
        applied += 1
        print(f'  OK: {label}')
    else:
        failed.append(label)
        print(f'  FAIL: {label}')

# ─── SECTION 3: Reliability Constraint intro (8 sentences → 1) ───
slice_replace('3-1 reliability intro',
    'Grid reliability is modeled as a reduced-form approximation',
    'The reliability constraint is formalized as:',
    'The reliability constraint maps the flexibility-to-volatility ratio into grid utilization:\n\n\tThe reliability constraint is formalized as:')

# ─── SECTION 3: Variable itemized list → 1 sentence ───
slice_replace('3-2 reliability variable list',
    '\n\tThe variables and parameters are:\n\t\\begin{itemize}',
    '\\end{itemize}',
    '\nwhere $F_t$ is total flexibility, $\\phi_{int}\\approx 55.6$ scales volatility to flexibility requirements, and $\\psi=2.0$ governs sensitivity.')

# ─── SECTION 3: Three paragraphs about Vol_ren (lines 292-298) → 2 sentences ───
slice_replace('3-3 Vol_ren paragraphs',
    'where $\\theta_{ren}$ is the renewable penetration share (dimensionless), $K_{ren}$ is installed renewable capacity (GW), and $\\sigma_{ren}$ is the coefficient of variation of the intermittency shock (dimensionless). As defined',
    'increasing $\\overline{\\text{Vol}}_{ren}$ is allowed to increase along the policy-determined VRE expansion path.',
    'where $\\theta_{ren}$, $K_{ren}$, and $\\sigma_{ren}$ are the penetration share, installed capacity, and CV of the intermittency shock. The intermittency shock $\\varepsilon_{ren,t}$ enters multiplicatively as $\\overline{\\text{Vol}}_{ren}\\cdot\\exp(\\varepsilon_{ren,t})$; in transition simulations $\\overline{\\text{Vol}}_{ren}$ rises with VRE expansion.')

# ─── SECTION 3: Flexibility Assets intro paragraph → 1 sentence ───
slice_replace('3-4a flexibility assets intro',
    'Total system flexibility $F_t$ is produced by combining two distinct assets: private battery storage and public grid infrastructure. These assets serve complementary roles. Battery storage shifts energy \\textit{across time} (charging when supply exceeds demand, discharging when it falls short), while grid infrastructure shifts energy \\textit{across space} (transmitting power from surplus to deficit regions). Using the shorthand $s_F \\equiv (\\rho - 1)/\\rho$ for the flexibility substitution parameter (where $\\rho$ is the elasticity of substitution between battery and grid capital), the flexibility bundle is:',
    '\n\t\\begin{equation}',
    'Total system flexibility $F_t$ combines private battery storage (time-shifting) and public grid (space-shifting) via a CES aggregator:\n\n\t\\begin{equation}')

# ─── SECTION 3: Post-flexibility-equation verbose paragraph → 1 sentence ───
slice_replace('3-4b flexibility post-eq paragraph',
    'Flexibility is a weighted combination of technology-augmented battery capacity ($A_{bat,t-1} K_{b,t-1}$) and grid infrastructure ($K_{g,t-1}$), where neither can fully replace the other. A grid with abundant batteries but no transmission lines cannot deliver stored power; a grid with extensive transmission but no storage cannot smooth intermittent supply. The negative exponent $s_F < 0$ (since $\\rho = 0.4 < 1$) ensures that flexibility collapses if either input approaches zero---a mathematical expression of this engineering constraint.\n\n\t\\textbf{Timing convention.} All three inputs to the flexibility aggregator are predetermined: capital stocks $K_{b,t-1}$ and $K_{g,t-1}$ installed by the end of period $t-1$, and battery technology $A_{bat,t-1}$ induced in period $t-1$. Scarcity-induced innovation in period $t$ updates $A_{bat,t}$ (equation~\\ref{eq:learning_by_doing}), but this improvement only enters flexibility in period $t+1$. This one-period-delay convention means that innovative effort is directed today but becomes operationally productive tomorrow, eliminating any within-period simultaneity between the innovation and reliability equations.\n\n\tThe elasticity of substitution $\\rho = 0.4$ is calibrated below unity to enforce this complementarity. The parameter $\\mu = 0.16$ governs the relative weight of battery storage in the aggregate, reflecting Vietnam\'s current early-stage storage deployment. Innovation affects flexibility through $A_{bat,t-1}$: technology improvements induced by last period\'s reliability gap raise the usable energy-shifting capacity per unit of installed storage in the current period.\n\n\tThe complementarity parameter $\\rho < 1$ has a precise formal consequence for the returns to coordinated flexibility investment. See Proposition~\\ref{prop:superadditivity} in Appendix~\\ref{app:proof_superadditivity} for the statement and quantitative characterization.',
    '\n\n\t\\subsubsection{Laws of Motion}',
    'The negative exponent $s_F<0$ ($\\rho=0.4<1$) enforces imperfect complementarity; all inputs are predetermined from period $t-1$, eliminating within-period simultaneity. See Proposition~\\ref{prop:superadditivity} (Appendix~\\ref{app:proof_superadditivity}) for the superadditivity result.\n\n\t\\subsubsection{Laws of Motion}')

# ─── SECTION 3: Battery investment — long paragraph before the rule equation → 1 sentence ───
slice_replace('3-5 battery investment pre-eq paragraph',
    'The rental rate $R_{b,t}$ depends on the reliability gradient $\\partial u / \\partial F = \\psi(1-u_t)/(\\phi_{int} \\cdot \\text{Vol}_{ren,t})$, which is proportional to the reliability gap $(1-u_t)$. When reliability falls below target, $(1-u_t)$ rises, steepening the gradient and raising $R_{b,t}$ above its steady-state value. The ratio $u^*/u_t$ is a monotone transformation of this gap: at steady state $u_t = u^* = 0.97$ and the ratio equals unity; when $u_t$ drops below $u^*$, the ratio rises, driving investment above replacement levels. This motivates the reduced-form investment rule:',
    '\n\t\\begin{equation}',
    'Rising reliability gap $(1-u_t)$ steepens the gradient and raises $R_{b,t}$, motivating:\n\n\t\\begin{equation}')

# ─── SECTION 3: Scarcity-induced tech — itemized list + 2 long paragraphs → 1 sentence ───
slice_replace('3-6 SITC itemized + paragraphs',
    '\n\twhere:\n\t\\begin{itemize}\n\t\t\\item $\\eta_{bat} > 0$ is the technology adoption elasticity governing how strongly the reliability scarcity signal translates into innovation\n\t\t\\item $\\chi \\in (0,1]$ captures regulatory distortions in price signals\n\t\t\\item $u^*$ is the target reliability level\n\t\t\\item $u_t$ is the current reliability\n\t\\end{itemize}\n\n\tThis specification ensures that technological progress is zero when reliability is at target ($u_t = u^*$) and accelerates proportionally when reliability falls below target. The mechanism is distinct from both pure learning-by-doing (where the driver is cumulative deployment or production volume) and pure Directed Technical Change (where scientists allocate R\\&D effort across sectors). Instead, it captures the empirically documented phenomenon that factor scarcity raises the private return to innovation in that factor, directing technological capabilities toward the binding bottleneck \\citep{popp2002, acemoglu2012environment, aghion2016, hassler2021}. The accumulation structure --- compounding from $A_{bat,t-1}$ with bounded growth --- ensures that improvements persist and exhibit the saturation properties documented in battery experience curves \\citep{kittner2017, way2022}.\n\n\t\\textbf{Timing and causal ordering.} The update $A_{bat,t}$ is determined within period $t$ from the current reliability gap, but it only enters the flexibility aggregator \\eqref{eq:flexibility_bundle} in period $t+1$ as $A_{bat,t-1}$. This one-period delay means the sequence within any period $t$ is strictly recursive and free of simultaneity: (i) $F_t$ is computed from predetermined $A_{bat,t-1}$ and $K_{b,t-1}$; (ii) $u_t$ follows from $F_t$ and the current intermittency shock; (iii) $A_{bat,t}$ is then updated from the reliability gap $(u^* - u_t)/u^*$. Economically, this captures the reality that scarcity signals are processed and innovation decisions are made within the current quarter, but the resulting improvements to battery systems become operationally available only in the following period.\n\n\tThe parameter $\\chi$ captures regulatory distortions that dampen the transmission of reliability scarcity into private incentives. In regulated markets with incomplete remuneration for ancillary services (e.g., fixed tariffs, lack of capacity payments, or missing real-time pricing), $\\chi < 1$, weakening innovation responses even when flexibility is scarce. In fully liberalized markets with scarcity pricing, $\\chi = 1$.\n\n\tUnlike Directed Technical Change, there is no R\\&D sector: innovation compounds from an inherited knowledge base directed by the physical reliability gap, activating the scarcity signals documented in \\citet{calel2016, aghion2016}.',
    '\n\n\\subsection{Government}',
    '\nwhere $\\eta_{bat}$ is the innovation elasticity, $\\chi\\in(0,1]$ captures regulatory signal attenuation ($\\chi=1$ baseline; $\\chi=0.3$ regulated counterfactual), and $u^*$ is the reliability target. Innovation is zero at $u_t=u^*$, accelerates when reliability falls below target, and applies with a one-period delay into $F_{t+1}$.\n\n\\subsection{Government}')

# ─── SECTION 3: Government — long proportional tax paragraph → 2 sentences ───
slice_replace('3-7 govt tax paragraph',
    'The government finances grid investment through a proportional tax $\\tau$ on household output income. The tax rate is a \\textit{fixed parameter}, calibrated in steady state to satisfy the government\'s budget constraint:',
    'Since $\\tau$ is fixed, it does not respond to cyclical fluctuations in grid investment --- it acts as a constant fiscal wedge that affects the steady-state level of capital and output but does not amplify business cycle dynamics.',
    'The government finances grid investment through a fixed proportional tax $\\tau_{ss}=I_{grid,ss}/Y_{ss}=1.43\\%$, calibrated to balance the steady-state budget:')

# ─── SECTION 3: Model Mechanics — remove the long "greenflation" paragraph after Ch 1 equation ───
slice_replace('3-8 greenflation paragraph',
    'This is the ``greenflation\'\' channel: supply-side volatility from renewables generates simultaneous output losses and energy price increases. The mechanism is structurally analogous to the adverse supply shocks analyzed in \\citet{blanchard2007}, who show that the macroeconomic consequences of energy supply disruptions depend critically on the degree of real wage rigidity and the substitutability between energy and other inputs. In the present model, the low calibrated value of $\\sigma_E$ ensures that intermittency shocks behave as ``unfavorable supply shocks\'\' in the Blanchard-Gal\\\'\\{\\i\\} sense, generating a policy-relevant trade-off between output stabilization and energy price stabilization.',
    '\n\n\t\\paragraph{Channel 2: The Innovation Response.}',
    '\n\n\t\\paragraph{Channel 2: The Innovation Response.}')

# ─── SECTION 3: Channel 2 long paragraph → 2 sentences ───
slice_replace('3-9 channel2 paragraph',
    'The reliability decline widens the gap $(u^* - u_t)/u^*$, which accelerates scarcity-induced technological change in storage (Equation~\\ref{eq:abat_lin}). The resulting improvement in $A_{bat,t}$ increases the effective services per unit of battery capital in period $t+1$, gradually restoring flexibility and utilization. Note the one-period delay: innovation induced by the reliability shortfall in period $t$ raises $A_{bat,t}$, which enters $F_{t+1}$ as $A_{bat,t-1}$ and improves reliability only from period $t+1$ onward. The speed of this response depends critically on three parameters: the innovation elasticity $\\eta_{bat}$, the regulatory transmission coefficient $\\chi$, and the battery share in flexibility $s_b$. When any of these is small, the innovation channel is weak, and adjustment relies on costly capital accumulation rather than efficient technological improvement.',
    '\n\n\n\t\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}',
    'The reliability decline widens the gap $(u^*-u_t)/u^*$, accelerating scarcity-induced innovation in storage (Eq.~\\ref{eq:abat_lin}) and gradually restoring flexibility. Innovation in period $t$ enters $F_{t+1}$ with a one-period delay.\n\n\t\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}')

# ─── SECTION 3: Channel 4 Fiscal Wedge — cut long paragraph → 2 sentences ───
slice_replace('3-10 channel4 paragraph',
    'This channel identifies the level distortion introduced by proportional taxation, which is absent in the lump-sum baseline. While the tax rate $\\tau$ remains fixed during the business cycle, it imposes a permanent wedge on the incentives to work and invest (equations~\\ref{eq:euler_lin}--\\ref{eq:labor_lin}). The steady-state capital stock $K_{p,ss}$ is approximately 1.4\\% lower in the proportional tax variant compared to the lump-sum benchmark, reflecting the lower after-tax marginal product of capital. This means that for a given percentage intermittency shock, the economy operates from a more constrained initial condition, making each percentage decline in reliability costlier in terms of final output units. This channel implies that proportional taxation reduces the economy\'s structural resilience to intermittency shocks even without providing dynamic amplification.',
    '\n\n\t\\paragraph{Interaction: The ``Reliability Valley.\'\'}',
    'The fixed tax $\\tau$ imposes a permanent wedge on capital returns: $K_{p,ss}$ is $\\approx1.4\\%$ lower than the lump-sum benchmark, reducing structural resilience to intermittency shocks without dynamic amplification.\n\n\t\\paragraph{Interaction: The ``Reliability Valley.\'\'}')

# ─── SECTION 5: "Core contribution" sentence → delete ───
simple_replace('5-1 core contribution sentence',
    '\n\tA core contribution of this capstone is the rigorous calibration of the reliability constraint, which is micro-founded on the physical properties of the Vietnamese grid rather than arbitrary assumptions.\n\n\t',
    '\n\t')

# ─── SECTION 5: mu calibration — 3 paragraphs with 2 equations → 1 sentence ───
slice_replace('5-2 mu calibration verbose',
    'The parameter $\\mu$ in the CES flexibility aggregator governs the relative importance of battery storage versus grid transmission in producing system flexibility. Rather than calibrating to a single point-in-time snapshot, I discipline $\\mu$ to reflect the forward-looking battery investment share across the PDP8 planning horizon (2021--2030), capturing Vietnam\'s policy trajectory toward accelerated energy storage deployment.\n\n\t\\textit{Near-term period (2021--2026):}',
    'Averaging near-term ($\\mu\\approx 0.14$) and late-period ($\\mu\\approx 0.64$) shares, I adopt $\\mu = 0.16$ as a conservative near-terminal value consistent with PDP8 commitments; results are robust to $\\mu \\in [0.12, 0.20]$.',
    '$\\mu=0.16$ is calibrated from PDP8 investment flows: the near-term battery share is $0.24/(0.24+1.49)\\approx 0.14$, rising to $0.64$ in 2027--2030; the baseline adopts a conservative near-terminal value, robust to $\\mu\\in[0.12,0.20]$.')

# ─── SECTION 5: sigma_f — 2 paragraphs → 1 sentence ───
slice_replace('5-3 sigma_f explanation',
    '\\textbf{2. Elasticity of Substitution ($\\sigma_f = 0.4$):}\n\tStandard macro-energy models often assume high substitutability between energy capital types. However, \\citet{hart2019} and \\citet{hirth2013} demonstrate that transmission and storage are \\textit{imperfect complements} in the presence of intermittency. Transmission moves energy across space, while storage moves energy across time. When solar generation is zero (nighttime), infinite transmission capacity yields zero marginal utility.\n\n\tIn the model, the CES flexibility aggregator is parameterized as $F_t = \\bigl(\\mu\\,(A_{bat}K_{b})^{(\\sigma_f-1)/\\sigma_f} + (1-\\mu)\\,K_{g}^{(\\sigma_f-1)/\\sigma_f}\\bigr)^{\\sigma_f/(\\sigma_f-1)}$, where $\\sigma_f$ is the \\textit{elasticity of substitution} between battery storage and grid capital. Setting $\\sigma_f = 0.4 < 1$ places battery and grid in the gross-complements region ($\\sigma_f < 1$), consistent with the physical argument that neither input alone can deliver system reliability without the other. The CES exponent $(\\sigma_f - 1)/\\sigma_f = -1.5$ and outer power $\\sigma_f/(\\sigma_f-1) = -0.67$ enforce this complementarity in the Dynare model. This is consistent with \\citet{hart2019}\'s empirical finding that the cross-price elasticity between storage and transmission investment is negative.',
    '\n\n\t\\textbf{3. Reliability Function Parameters:}',
    '\\textbf{2. Elasticity of Substitution ($\\sigma_f = 0.4$):} Battery and grid are gross complements ($\\sigma_f<1$): transmission moves energy across space while storage moves it across time, and neither alone delivers reliability \\citep{hart2019}.\n\n\t\\textbf{3. Reliability Function Parameters:}')

# ─── SECTION 5: Reliability function parameters — long itemized block → 2 sentences ───
slice_replace('5-4 reliability param block',
    '\\begin{itemize}\n\t\t\\item $\\psi = 2.0$: Sensitivity parameter governing how rapidly reliability degrades when the flexibility-to-volatility ratio declines. $\\psi$ and $\\phi_{int}$ are \\textit{jointly} identified from two independent empirical moments, avoiding the under-identification that would arise from the level condition alone.\n\n\t\t\\textit{First moment (level):} $u_{ss} = 0.97$',
    '$\\phi_{int}$ is not a free parameter once $\\psi$ is identified --- it is uniquely determined by the two moments and the calibrated flexibility stocks.\n\t\\end{itemize}',
    '$\\psi=2.0$ is identified from the slope moment $\\Delta u/(\\Delta\\ln F\\cdot(1-u_{ss}))\\approx 1.85$ (EVN data: reliability rose 1.5~pp following a 20\\% transmission expansion, 2016--2020 \\citep{worldbank2022}); $\\phi_{int}\\approx 55.6$ is then pinned by the level condition $u_{ss}=0.97$. Both are robust for $\\psi\\in[1.5,3.0]$.')

# ─── SECTION 5: chi paragraph — long paragraph → 2 sentences ───
slice_replace('5-5 chi paragraph',
    'The price signal transmission parameter is set to $\\chi = 1.0$ in the baseline (full price transmission). The regulated-market counterfactual uses $\\chi = 0.3$, anchored to the empirical retail electricity price pass-through in Vietnam. Under EVN\'s administered tariff system, retail electricity prices are set by the Ministry of Industry and Trade and updated infrequently; IEEFA (2022) documents that Vietnam\'s retail tariffs have been held below cost-recovery levels for extended periods, with effective pass-through of wholesale scarcity cost signals estimated at approximately 25--35\\% \\citep{ieefa2022}. The value $\\chi = 0.3$ is calibrated to the midpoint of this range, representing a regime in which only 30 cents of every dollar of reliability-gap signal reaches the innovating sector --- consistent with Vietnam\'s fixed-tariff electricity market structure prior to the 2023 tariff reform.',
    '\n\n\t\\subsection{Technology Adoption Elasticity and Innovation}',
    'The baseline sets $\\chi=1.0$ (full signal transmission); the regulated counterfactual uses $\\chi=0.3$, calibrated to Vietnam\'s documented retail tariff pass-through of 25--35\\% \\citep{ieefa2022}.\n\n\t\\subsection{Technology Adoption Elasticity and Innovation}')

open('main_body.tex', 'w', encoding='utf-8').write(content)
print(f'\nApplied {applied}/{applied+len(failed)} cuts')
if failed:
    print('Failed:', failed)
print(f'Size: {original_len} → {len(content)} bytes (Δ{len(content)-original_len:+,})')
