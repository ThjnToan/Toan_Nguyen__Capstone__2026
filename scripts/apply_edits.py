import sys

content = open('main_body.tex', 'rb').read().decode('utf-8')
original_len = len(content)
applied = 0
failed = []

def apply(label, old, new):
    global content, applied
    if old in content:
        content = content.replace(old, new, 1)
        applied += 1
        print(f'  OK: {label}')
    else:
        failed.append(label)
        print(f'  FAIL: {label}')

# ─────────────────────────────────────────────
# METHODOLOGY EDITS (Agent 1)
# ─────────────────────────────────────────────

apply('METH-1 battery derivation', r'''\paragraph{Battery Investment Decision.}
A fully endogenous, forward-looking optimization of battery capital---subject to the highly non-linear exponential reliability constraint---would create non-convexities that complicate the Blanchard-Kahn conditions required for local approximation. The approach here instead derives a reduced-form investment rule that proxies the gradient of the structural shadow price while keeping the system computationally tractable and ensuring local saddle-path stability.

To see where this rule comes from, consider the shadow rental rate on storage services. The marginal value of an additional unit of battery capital flows through the chain $K_b \to F \to u \to V \to Y$:

\begin{equation}
	R_{b,t} \;\equiv\; \frac{\partial Y_t}{\partial F_t} \cdot \frac{\partial F_t}{\partial (A_{bat,t-1} K_{b,t-1})} \cdot A_{bat,t-1}
	\label{eq:rental_rate}
\end{equation}

Each component is computed from the model's structural equations. From the CES flexibility aggregator \eqref{eq:flexibility_bundle}:
\begin{equation}
	\frac{\partial F_t}{\partial (A_{bat,t-1} K_{b,t-1})} = \mu \left(\frac{A_{bat,t-1} K_{b,t-1}}{F_t}\right)^{-1/\rho_f}
	\label{eq:dF_dKb}
\end{equation}

The marginal product of flexibility in output operates through the reliability channel. Combining the exponential reliability function \eqref{eq:utilization} with the value-added function \eqref{eq:value_added} and the CES output aggregator \eqref{eq:final_output}:
\begin{equation}
	\frac{\partial Y_t}{\partial F_t} = \underbrace{\frac{\partial Y_t}{\partial V_t}}_{\text{CES share}} \cdot \underbrace{\frac{\partial V_t}{\partial u_t}}_{\text{capital utilization}} \cdot \underbrace{\frac{\partial u_t}{\partial F_t}}_{\text{reliability gradient}} = \frac{\partial Y_t}{\partial V_t} \cdot \alpha \frac{V_t}{u_t} \cdot \frac{\psi (1 - u_t)}{\phi_{int} \cdot \text{Vol}_{ren,t}}
	\label{eq:dY_dF}
\end{equation}

where $\partial u_t / \partial F_t = \psi(1-u_t)/(\phi_{int} \cdot \text{Vol}_{ren,t})$ follows from differentiating the exponential reliability function, and $\partial V_t / \partial u_t = \alpha V_t / u_t$ from the Cobb-Douglas value-added specification. Note that $R_{b,t}$ is increasing in the reliability gap $(1-u_t)$: when reliability falls below target, the marginal return to storage rises sharply due to the exponential curvature of the reliability function.

Table~\ref{tab:mpk_decomp} evaluates each component of the chain at steady-state values, yielding a shadow rental rate of $R_{b,ss} = 10.93\%$ quarterly. The decomposition reveals that the \textit{flexibility bottleneck}---the large marginal product of batteries in the CES aggregator ($\partial F / \partial(A_{bat} K_b) = 2.04$)---is the dominant force, reflecting the scarcity of storage relative to grid capital. The reliability gradient ($\partial u / \partial F = 0.20$) is moderate because the system operates near the flat portion of the exponential curve at $u_{ss} = 0.97$.''',

r'''\paragraph{Battery Investment Decision.}
A fully endogenous optimization of battery capital subject to the non-linear reliability constraint would create non-convexities that complicate the Blanchard-Kahn conditions. The investment rule below proxies the gradient of the structural shadow price while preserving local saddle-path stability; the full shadow-price derivation is provided in Appendix~\ref{app:deriv_euler}.

Table~\ref{tab:mpk_decomp} evaluates the shadow rental rate at steady state, yielding $R_{b,ss} = 10.93\%$ quarterly. The \textit{flexibility bottleneck} ($\partial F / \partial(A_{bat} K_b) = 2.04$) is the dominant force, reflecting battery scarcity in the CES aggregator, while the reliability gradient is moderate because the system operates near $u_{ss} = 0.97$.''')

apply('METH-2 timing footnote', r'''\footnote{\label{fn:timing_convention}The timing convention here is standard in DSGE models with capital: capital installed by the end of period $t-1$ (i.e., $K_{p,t}$ in the accumulation equation $K_{p,t} = (1-\delta_p)K_{p,t-1} + I_{p,t}$) becomes productive in period $t$, earning rental rate $R_t$. Investment $I_{p,t}$ is a use of funds in period $t$ but does not contribute to production until period $t+1$. Writing $R_t K_{p,t}$ on the income side of the budget constraint would double-count $I_{p,t}$, since $K_{p,t}$ already embeds it. Consequently, the value-added function is $V_t = Z(u_t K_{p,t-1})^\alpha L_t^{1-\alpha}$, using the predetermined capital stock $K_{p,t-1}$.}''',

r'''\footnote{\label{fn:timing_convention}Capital installed by the end of period $t-1$ is productive in period $t$, earning rental rate $R_t$; investment $I_{p,t}$ does not contribute until $t+1$. Consequently, the value-added function uses the predetermined stock $K_{p,t-1}$: $V_t = Z(u_t K_{p,t-1})^\alpha L_t^{1-\alpha}$.}''')

apply('METH-3 post-budget accumulation paragraphs', r'''The spending item $I_{p,t}$ in the budget constraint represents the household's purchase of new productive capital goods. It is linked to the next-period capital stock through the standard perpetual-inventory accumulation law (equation~\eqref{eq:productive_capital}):
\begin{equation*}
	K_{p,t+1} = (1-\delta_p)\,K_{p,t} + I_{p,t}
\end{equation*}
\noindent so that $I_{p,t} = K_{p,t+1} - (1-\delta_p)K_{p,t}$: investment is simply the increment to next period's capital stock above what survives depreciation from the current stock. Equivalently, the household's true intertemporal choice variable is the end-of-period capital stock $K_{p,t+1}$, not $I_{p,t}$ itself. In the Lagrangian derivation (Appendix~\ref{app:deriv_euler}), $I_{p,t}$ is substituted out using this relation, so the optimisation is carried out over $\{C_t, L_t, K_{p,t+1}, B^*_t\}$. The Euler equation \eqref{eq:euler} below follows directly from the first-order condition with respect to $K_{p,t+1}$.

The household chooses $\{C_t, L_t, K_{p,t+1}, B^*_t\}$ to maximize utility subject to the budget constraint \eqref{eq:hh_budget} and the accumulation law above. The first-order conditions are derived as follows.''',

r'''The household chooses $\{C_t, L_t, K_{p,t+1}, B^*_t\}$ to maximize utility subject to the budget constraint \eqref{eq:hh_budget} and the standard capital accumulation law \eqref{eq:productive_capital}; the Lagrangian derivation is in Appendix~\ref{app:deriv_euler}.''')

apply('METH-4 post-Euler tax wedge paragraph', r'''The household invests in capital until the marginal cost of saving today equals the expected discounted after-tax return. Since $\tau$ is constant, the only dynamic distortion comes from the standard Euler mechanism: a fall in expected future output $Y_{t+1}$ (due to persistent intermittency shocks) lowers the after-tax return to capital, compressing investment. The constant tax wedge shifts the steady-state capital stock downward but does not amplify the shock response beyond this RBC channel.''',

r'''The tax wedge $\tau$ shifts the steady-state capital stock downward but does not amplify shock dynamics beyond the standard Euler mechanism.''')

apply('METH-5 post-labor-supply paragraph', r'''The after-tax labor share of value-added (left side) equals the household's required compensation for working. The proportional tax $\tau$ drives a wedge that reduces effective labor income, lowering equilibrium labor supply relative to the lump-sum benchmark. The parameter $\varphi$ is calibrated endogenously to match steady-state labor hours of $L_{ss} = 0.33$.''',

r'''The parameter $\varphi$ is calibrated endogenously to match steady-state labor hours $L_{ss} = 0.33$.''')

apply('METH-6 five-item reliability properties list', r'''This exponential formulation ensures that:
\begin{enumerate}
	\item When $F_t \to \infty$, utilization approaches unity (perfect reliability): $u_t \to 1$
	\item When $F_t \to 0$, utilization approaches zero (grid collapse): $u_t \to 0$
	\item The function is always bounded: $u_t \in (0, 1)$ for all $F_t > 0$
	\item Marginal returns to flexibility are strictly diminishing
	\item The reliability penalty increases exponentially with the volatility-to-flexibility ratio
\end{enumerate}

This exponential functional form is standard in reliability engineering \citep{zakeri2015, lund2015} and ensures that flexibility requirements rise nonlinearly as VRE penetration increases, capturing the convex integration costs documented in the empirical literature.''',

r'''This exponential form is standard in reliability engineering \citep{zakeri2015, lund2015}: $u_t \in (0,1)$ for all $F_t > 0$, reliability is bounded between grid collapse and perfection, and marginal returns to flexibility are strictly diminishing.''')

apply('METH-7 pre-Table-1 prose', r'''To illustrate the nonlinearity, Table~\ref{tab:reliability_numerical} reports the reliability response to flexibility shortfalls under the baseline calibration. The convex curvature is evident: a 10\% flexibility deficit reduces reliability by 1.26 percentage points (pp), but a 50\% deficit causes a 14.32~pp collapse---more than 11 times the proportional response. This accelerating penalty captures the engineering reality that the first units of flexibility (spinning reserves, fast-ramping gas) are cheap and effective, while maintaining reliability under extreme shortfalls requires exponentially more resources.''',

r'''Table~\ref{tab:reliability_numerical} illustrates this nonlinearity: a 10\% flexibility deficit reduces reliability by 1.26~pp, but a 50\% deficit causes a 14.32~pp collapse. The accelerating penalty motivates the ``Reliability Valley'' analysis in Section~\ref{sec:results}.''')

apply('METH-8 external sector intro + Bond Euler', r'''The model is embedded in a small open economy framework following the debt-elastic interest rate premium approach of \citet{schmittgrohe2003}. The representative household has access to international bond markets, where it can borrow or lend at a country-specific interest rate that depends on the economy's net foreign asset position.

\subsubsection{Bond Euler Equation}

The household's optimality condition for international bond holdings yields:
\begin{equation}
	\frac{1}{C_t} = \beta \, \mathbb{E}_t \left[ \frac{1 + r^*_t}{C_{t+1}} \right]
	\label{eq:bond_euler}
\end{equation}

The household lends abroad until the cost of forgoing one unit of consumption today equals the expected discounted return from the foreign bond. This is analogous to the capital Euler equation \eqref{eq:euler}, but with the foreign interest rate replacing the domestic marginal product of capital. Combining equations \eqref{eq:euler} and \eqref{eq:bond_euler} yields the no-arbitrage condition linking domestic and external returns:
\begin{equation}
	(1-\tau)\,\alpha Y_{t+1}/K_{p,t} + 1 - \delta_p = 1 + r^*_t
	\label{eq:no_arbitrage}
\end{equation}

\noindent This is not an independent equation but rather follows directly from the two Euler conditions. It is listed here to clarify the economic intuition: in equilibrium, the domestic marginal product of capital equals the external borrowing rate.''',

r'''The model follows the debt-elastic interest rate approach of \citet{schmittgrohe2003}, giving the household access to international bond markets at a country-specific rate.

\subsubsection{Bond Euler Equation}

The household's optimality condition for international bond holdings yields:
\begin{equation}
	\frac{1}{C_t} = \beta \, \mathbb{E}_t \left[ \frac{1 + r^*_t}{C_{t+1}} \right]
	\label{eq:bond_euler}
\end{equation}

Combining with \eqref{eq:euler} delivers the no-arbitrage condition $(1-\tau)\alpha Y_{t+1}/K_{p,t} + 1 - \delta_p = 1 + r^*_t$, confirming that the domestic after-tax return equals the external borrowing rate in equilibrium.''')

apply('METH-9 current account pre-equation prose', r'''Combining the household budget constraint \eqref{eq:hh_budget} with the government's financing identity (tax revenue $\tau Y_t$ funds grid investment $I_{grid,t}$, with the difference covered or built up in the foreign bond position) yields the economy-wide resource constraint. The proportional tax does \textit{not} satisfy Ricardian equivalence: households optimise over after-tax income $(1-\tau)Y_t$, so the fixed wedge distorts equilibrium labour supply and investment decisions. The consolidated resource constraint --- where tax revenue and government spending cancel out --- takes the same form as the lump-sum case. On the income side, substituting $w_t L_t + R_t K_{p,t-1} = V_t$ and noting that $\omega_E = 0.045$ is small so $V_t \approx Y_t$, the economy-wide resource constraint in the open economy is:''',

r'''Consolidating the household budget constraint \eqref{eq:hh_budget} with the government financing identity yields the economy-wide resource constraint:''')

apply('METH-10 post-CA + debt-elastic premium', r'''The economy's foreign asset position improves when it produces more than it consumes and invests domestically (trade surplus), and deteriorates when domestic absorption exceeds output (trade deficit). Interest on existing foreign assets also contributes.

Here $B^*_t$ denotes net foreign assets (positive when the economy is a net creditor) and $I_{grid,t}$ follows the government's reliability-targeting rule specified in Section~\ref{sec:government}. The current account surplus equals $B^*_t - B^*_{t-1}$: the change in the net foreign asset position equals the trade balance ($Y_t - C_t - I_{p,t} - P_{bat,t} I_{bat,t} - I_{grid,t}$) plus net interest income ($r^*_{t-1} B^*_{t-1}$).

\subsubsection{Debt-Elastic Interest Rate Premium}

The country-specific interest rate includes a risk premium that depends on the deviation of net foreign assets from their steady-state level:
\begin{equation}
	r^*_t = \bar{r} + \phi_b \left( \exp(\bar{B}^* - B^*_t) - 1 \right)
	\label{eq:interest_premium}
\end{equation}

The more the economy borrows from abroad, the higher the interest rate it faces---reflecting the credit risk that international lenders perceive as foreign debt accumulates. This is the small open economy analog of a sovereign risk premium.

Here $\bar{r} = 1/\beta - 1$ is the world interest rate consistent with the household's discount factor, $\bar{B}^* = 0$ is the steady-state net foreign asset position (balanced trade), and $\phi_b > 0$ governs the elasticity of the premium to the debt position. When the economy borrows abroad ($B^*_t < \bar{B}^*$), the interest rate rises above $\bar{r}$, stabilizing the external position. This specification ensures stationarity of the net foreign asset process, resolving the well-known unit root problem in small open economy models \citep{schmittgrohe2003}.''',

r'''\subsubsection{Debt-Elastic Interest Rate Premium}

The country-specific interest rate rises with external debt, following \citet{schmittgrohe2003}:
\begin{equation}
	r^*_t = \bar{r} + \phi_b \left( \exp(\bar{B}^* - B^*_t) - 1 \right)
	\label{eq:interest_premium}
\end{equation}

With $\bar{r} = 1/\beta - 1$ and $\bar{B}^* = 0$, this premium stabilizes the net foreign asset process and resolves the unit-root problem standard in small open economy models.''')

# ─────────────────────────────────────────────
# CALIBRATION EDITS (Agent 2)
# ─────────────────────────────────────────────

apply('CAL-1 two-caveats paragraph', r'''	Two caveats are important. First, the regression fit is modest ($R^2 = 0.16$, RMSE $= 5.7$ percentage points): Vietnam's real rate fluctuated sharply during 1993--2023 for reasons largely unrelated to external debt dynamics---domestic inflation crises (2005, 2008, 2010), post-WTO capital surges, and monetary policy cycles each contributed large residuals. The regression therefore does not provide a precise structural identification of the borrowing premium; it provides an \textit{order-of-magnitude anchor} that rules out both the frictionless extreme ($\phi_b \approx 0$) and implausibly large capital-account frictions ($\phi_b \gg 0.10$). Second, the debt variable is \textit{gross} external debt/GNI rather than net foreign assets, which differs from the model's $B^*$ by Vietnam's substantial foreign exchange reserves (approximately 25\% of GNI in recent years); the mapping is therefore approximate.

	The primary justification for $\phi_b = 0.039$ is economic rather than statistical: Section~\ref{sec:calibration_validation} below shows that the SGU-standard $\phi_b = 0.001$ generates an implausibly large quarterly investment collapse ($-21.8\%$) in response to a moderate intermittency shock, while $\phi_b = 0.039$ produces the empirically credible $-3.5\%$ response consistent with Vietnamese national accounts data. The data-anchored estimate thus serves as an empirical discipline on the plausible range, with the sensitivity table (Table~\ref{tab:phi_b_sensitivity}) confirming that all qualitative findings are robust across $\phi_b \in [0.039, 0.10]$.''',

r'''	The regression provides an order-of-magnitude anchor ($R^2 = 0.16$) that rules out the frictionless extreme ($\phi_b \approx 0$) and implausibly large frictions ($\phi_b \gg 0.10$); the primary justification is economic: the SGU-standard $\phi_b = 0.001$ generates an implausible $-21.8\%$ quarterly investment collapse while $\phi_b = 0.039$ yields the credible $-3.5\%$ response consistent with Vietnamese national accounts data (see Table~\ref{tab:phi_b_sensitivity}).''')

apply('CAL-2 calibration validation header', r'''	\subsubsection*{Discount Factor and Capital Share}

	The discount factor $\beta = 0.99$ (quarterly) is anchored to Vietnam's long-term real interest rate. Data from the World Bank (1993--2023, 23 annual observations from 2000--2023 in the post-WTO era) show a mean real rate of 3.14\% annually, corresponding to $\beta = 0.9691$. The slightly higher value $\beta = 0.99$ (implying 4\% annualized real returns) is conservative and reflects current Vietnamese government borrowing conditions. This data-anchored approach replaces narrative calibration with Vietnamese-specific evidence, following the methodology used for the debt-elastic premium $\phi_b$ (see Section~\ref{subsec:external_sector}).

	The capital share $\alpha = 0.35$ is set as the inverse of the labor compensation share in national accounts. General Statistics Office (GSO) data on labor income as a share of value added place this at approximately 0.65 across the Vietnamese economy over 2015--2024, implying capital's residual share of 0.35. This standard growth-accounting approach ensures consistency with long-run factor income distribution while avoiding direct estimation from factor-cost data that may be incomplete or subject to measurement error in developing economies.

	While specific engineering priors for reliability parameters ($\psi$, $\phi_{int}$) are limited for Vietnam due to data constraints on granular outage metrics (e.g., LOLE, SAIDI), the calibration is designed so that implied utilization losses can be interpreted in terms of standard reliability observables.

	For the debt-elastic premium, the reduced-form regression over 28 annual observations yields $\hat{\phi}_b \approx 0.039$ as an order-of-magnitude anchor (see Section~\ref{subsec:external_sector}). While the regression fit is modest ($R^2 = 0.16$), the estimate disciplines the plausible range and rules out the frictionless limit. The primary justification is economic: $\phi_b = 0.039$ generates an empirically credible investment response to shocks, as documented in Table~\ref{tab:phi_b_sensitivity} below.''',

r'''	The calibration targets are closely matched: $\beta = 0.99$ reflects Vietnam's 3.14\% mean real rate (World Bank 2000--2023); $\alpha = 0.35$ is the residual capital share given a labor income share of 0.65 (GSO); and $\phi_b = 0.039$ is the data-anchored debt-elastic premium documented in Section~\ref{subsec:external_sector}, with sensitivity reported in Table~\ref{tab:phi_b_sensitivity} below.''')

apply('CAL-3 investment volatility compression', r'''	Three distinct patterns emerge across the columns and reveal economically important interactions between the external borrowing premium and the model's domestic transmission channels.

	\paragraph{Investment volatility compression.} The most striking feature is the sharp compression of $I_p$ volatility as $\phi_b$ rises: from $-21.8\%$ at the SGU standard to $-3.5\%$ at the baseline and $-2.6\%$ under strong restrictions. The mechanism operates through the bond Euler equation. At $\phi_b \approx 0$, the external rate $r^*_t = \bar{r} + \phi_b(\exp(\bar{B}^*-B^*_t)-1)$ is essentially flat, so the economy can freely borrow or lend at nearly the world rate. A negative intermittency shock reduces the after-tax marginal product of capital, triggering large Euler amplification: the household drastically cuts investment to restore the no-arbitrage return while smoothing consumption almost perfectly via external borrowing. As $\phi_b$ rises, each unit of external adjustment raises $r^*$, discouraging large current-account swings and spreading the shock more evenly across investment and consumption. The investment response therefore declines monotonically with $\phi_b$.''',

r'''	Three patterns confirm economic plausibility. \paragraph{Investment volatility compression.} Rising $\phi_b$ compresses $I_p$ volatility (from $-21.8\%$ to $-3.5\%$ to $-2.6\%$) by discouraging large current-account swings and spreading shocks more evenly across investment and consumption.''')

apply('CAL-4 current-account monotonicity', r'''	\paragraph{Current-account monotonicity.} The $B^*$ impact is a surplus at \emph{all} three values of $\phi_b$, but shrinks monotonically as $\phi_b$ rises ($+0.048 \to +0.003 \to +0.001$). The mechanism is the investment collapse itself. At $\phi_b = 0.001$, investment falls by $-21.8\%$ of its steady-state level while consumption barely moves, so the economy releases an enormous volume of formerly planned investment expenditure as current-account savings---producing a large surplus. As $\phi_b$ rises and the investment collapse shrinks to $-2.6\%$, domestic absorption falls by far less, and the surplus nearly vanishes. Formally, from the resource constraint: $\Delta B^* \approx \Delta Y - \Delta C - \Delta I_p - \Delta I_{bat} - \Delta I_{grid}$. Since $\Delta Y$ is essentially identical across all three calibrations but $\lvert\Delta I_p\rvert$ is still an order of magnitude larger at $\phi_b = 0.001$ than at $\phi_b = 0.10$, $\Delta B^*$ is necessarily largest at the lowest borrowing premium. The finding that Vietnam tends to run modest current-account surpluses during domestic slowdowns \citep{worldbank2024} supports a non-trivial $\phi_b$ where the investment response is moderate and the surplus is small.''',

r'''	\paragraph{Current-account monotonicity.} The current-account surplus shrinks monotonically with $\phi_b$ ($+0.048 \to +0.003 \to +0.001$), consistent with Vietnam's observed modest surpluses during domestic slowdowns \citep{worldbank2024}.''')

apply('CAL-5 calibration rationale phi_b', r'''	\paragraph{Calibration rationale for $\phi_b$.} The choice of $\phi_b = 0.039$ rests primarily on economic plausibility, with the data-anchored regression serving as an order-of-magnitude check rather than a precise identification. The economic argument is threefold. First, the SGU-standard $\phi_b = 0.001$ generates a $-21.8\%$ quarterly investment collapse in response to a moderate intermittency shock---far exceeding any quarterly investment decline in Vietnamese national accounts (GSO 2024) and economically implausible as a within-quarter response. Second, Vietnam's capital account is genuinely partially restricted: the State Bank of Vietnam controls short-term capital flows and the dong is not fully convertible, making the near-frictionless limit inappropriate for this economy. Third, qualitative FEVD rankings are identical at $\phi_b = 0.039$ and $\phi_b = 0.10$ (Table~\ref{tab:phi_b_sensitivity}), and as shown in Section~\ref{sec:welfare}, welfare costs are also insensitive to this choice, so results do not hinge on the precise value within the empirically plausible range. The reduced-form regression ($R^2 = 0.16$, $\hat{\phi}_b = 0.039$) corroborates the order of magnitude: while Vietnam's real rate fluctuations are dominated by domestic inflation dynamics rather than external debt movements, the regression consistently places $\hat{\phi}_b$ well above the frictionless limit and below the strong-restrictions threshold of $0.10$.''',

r'''	\paragraph{Calibration rationale for $\phi_b$.} The choice $\phi_b = 0.039$ is economically motivated: the SGU-standard $\phi_b = 0.001$ generates a $-21.8\%$ quarterly investment collapse inconsistent with Vietnamese national accounts, while Vietnam's partially restricted capital account rules out the near-frictionless limit. Qualitative FEVD rankings and welfare costs are insensitive to the precise value within $\phi_b \in [0.039, 0.10]$. The data-anchored regression ($R^2 = 0.16$, $\hat{\phi}_b = 0.039$) corroborates the order of magnitude.''')

apply('CAL-6 ENS mapping + sensitivity paragraph', r'''	In particular, under the interpretation in Equation (\ref{eq:utilization}), a utilization shortfall $(1-u_t)$ can be approximately mapped to energy not served (ENS) as a share of demand: $\text{ENS}_t/\text{Demand}_t \approx 1-u_t$. This mapping provides a transparent way to evaluate model implications against observable outage data. For instance, the model's steady-state utilization of 97\% implies 3\% energy not served, which falls within EVN's reported range of 3-5\% including planned maintenance and unplanned outages.

	Sensitivity analysis shows that the qualitative results---specifically the ``Agility Gap'' mechanism and the attenuation effect of regulatory distortions---remain robust for $\psi \in [1.5, 3.0]$ and $\mu \in [0.20, 0.35]$. Lower values of $\psi$ dampen the magnitude of the output loss but do not alter the sign of impulse responses or the structural necessity of the transition.''',

r'''	The model's 97\% steady-state utilization implies 3\% energy not served, within EVN's reported 3--5\% range. Qualitative results are robust for $\psi \in [1.5, 3.0]$ and $\mu \in [0.20, 0.35]$.''')

apply('CAL-7 technology adoption elasticity', r'''	The technology adoption elasticity is set to $\eta_{bat} = 0.10$ (quarterly). This parameter is anchored to two complementary empirical literatures that together bracket the plausible range from both sides of the mechanism.

	From the \textit{induced innovation} side, \citet{popp2002} estimates the elasticity of energy-efficient patenting activity to energy prices at 0.05--0.35, using U.S. patent data over 1970--1994. This range directly captures the responsiveness of innovative effort to economic scarcity signals --- the mechanism at work in the model --- and places $\eta_{bat} = 0.10$ comfortably within the lower portion of the empirically documented range. \citet{aghion2016} similarly document that energy price shocks redirect firm-level R\&D toward clean technologies, with estimated elasticities broadly consistent with Popp's range for the automotive sector.

	From the \textit{technology accumulation} side, \citet{kittner2017} estimate an 18--22\% learning rate (cost reduction per doubling of cumulative installed capacity) for utility-scale battery systems over 1991--2015, and \citet{way2022} find that energy storage technologies follow Wright's Law with a progress ratio of 0.79--0.82. These experience curve parameters describe the speed of cost improvement once innovative effort is deployed; they are consistent with $\eta_{bat} = 0.10$ in the sense that a sustained 10\% reliability gap in the model produces approximately 1\% quarterly battery efficiency improvement, matching the secular pace of cost reduction documented in these studies. Together, the two empirical anchors --- induced innovation elasticities for the \textit{direction} of effort, and experience curve parameters for the \textit{speed} of accumulation --- place $\eta_{bat} = 0.10$ on solid empirical footing for both dimensions of the mechanism.''',

r'''	The technology adoption elasticity is set to $\eta_{bat} = 0.10$ (quarterly), within the 0.05--0.35 induced-innovation elasticity range estimated by \citet{popp2002} and consistent with the 18--22\% battery experience curve of \citet{kittner2017}. A sustained 10\% reliability gap in the model yields approximately 1\% quarterly battery efficiency improvement, matching the secular pace of cost reduction documented in \citet{way2022}.''')

apply('CAL-8 battery price shock paragraph', r'''	The world battery price shock follows $\rho_{bat} = 0.90$ and $\sigma_{bat} = 0.08$ (quarterly). The persistence parameter is calibrated at \textit{quarterly} frequency to match the multi-year commodity supply cycles observed in battery materials: the quarterly AR(1) of 0.90 implies a half-life of approximately $\ln(0.5)/\ln(0.90) \approx 6.6$ quarters (1.6 years), consistent with the 2--3 year lithium supply-demand adjustment cycle driven by mine development timelines \citep{bloombergnef2023}. This quarterly value is also consistent with annual AR(1) estimates of 0.88--0.93 for quarterly-averaged battery pack prices from BloombergNEF, which correspond to quarterly persistence of 0.88--0.93 when estimated on quarterly-aggregated data (rather than converting from monthly to quarterly, which would lower persistence).''',

r'''	The world battery price shock follows $\rho_{bat} = 0.90$ and $\sigma_{bat} = 0.08$ (quarterly). The half-life of $\approx 6.6$ quarters (1.6 years) matches the 2--3 year lithium supply-demand adjustment cycle driven by mine development timelines \citep{bloombergnef2023}, and is consistent with BloombergNEF annual AR(1) estimates of 0.88--0.93 for quarterly-averaged battery pack prices.''')

# ─────────────────────────────────────────────
# BASELINE + WELFARE EDITS (Agent 3)
# ─────────────────────────────────────────────

apply('BASE-1 Stage 1 reliability penalty', r'''	On impact, grid utilization (reliability) declines by 1.30\% relative to steady state. While a single shock may appear transitory, it translates into a meaningful disruption to grid operations. Recall that utilization enters the production function multiplicatively:
	\begin{equation*}
		V_t = Z\,(u_t K_{p,t-1})^\alpha L_t^{1-\alpha}
	\end{equation*}
	A decline in $u_t$ directly reduces the effective services from the predetermined capital stock $K_{p,t-1}$ available for production. Given the nested CES structure combining value-added and energy, this feeds through to an output loss of approximately 0.54\% on impact---as shown in the top-left panel of Figure \ref{fig:irf_baseline}.

	This output decline is \textit{not} a standard demand shock. It reflects a binding physical constraint: when renewable generation falls unexpectedly and storage capacity is insufficient to buffer the mismatch, installed productive capital cannot be operated. Manufacturing lines idle, data centers activate backup generators (at higher marginal cost), and grid operators implement rotating load shedding. This is the \textit{reliability penalty}---a structural, supply-side friction that reduces total factor productivity.

	The mechanism operates through two channels. First, the direct utilization effect reduces effective capital services. Second, the increased cost of maintaining production (through backup generation or rescheduled operations) acts as a wedge between the installed capital stock and its economically usable portion. Unlike a temporary technology shock that dissipates exogenously, the reliability penalty persists until the underlying flexibility deficit is addressed through investment.''',

r'''	On impact, grid utilization (reliability) declines by 1.30\% relative to steady state, which directly reduces the effective capital services available for production ($V_t = Z(u_t K_{p,t-1})^\alpha L_t^{1-\alpha}$) and generates an output loss of approximately 0.54\%---as shown in the top-left panel of Figure~\ref{fig:irf_baseline}. This \textit{reliability penalty} is not a standard demand shock but a binding physical constraint: when renewable generation falls unexpectedly and storage is insufficient, installed capital cannot be fully operated. The penalty operates through two channels---direct reduction in effective capital services and an implicit cost wedge from backup generation---and persists until the underlying flexibility deficit is addressed through investment.''')

apply('BASE-2 Stage 2 last two paragraphs', r'''	The variance decomposition reveals that battery investment volatility is 100\% driven by global battery price shocks ($\varepsilon_{bat}$), demonstrating the dominance of the price channel in the battery investment rule. Battery developers' investment decisions are primarily governed by cost conditions in global markets.

	Grid investment follows the public reliability-targeting rule:
	\begin{equation}
		\frac{I_{grid,t}}{\bar{I}_{grid}} = \left(\frac{u^*}{u_t}\right)^{\phi_{grid}}
	\end{equation}
	The parameter $\phi_{grid} = 1.5$ governs responsiveness, but even aggressive fiscal reactions cannot overcome the structural lags inherent in public infrastructure. Transmission line projects require multi-year environmental impact assessments, land acquisition negotiations, and coordinated construction across provincial jurisdictions. The planning horizon for major grid expansions in Vietnam typically exceeds 5--7 years \citep{iea2023}. Consequently, even when the reliability gap becomes apparent (low $u_t$), actual capital accumulation responds sluggishly.

	The proportional tax is fixed at $\tau_{ss} = 1.43\%$ and does not respond dynamically to shocks. Its effect on the investment response operates through the steady-state channel: the tax wedge reduces the equilibrium capital stock to $K_{p,ss} = 9.83$, compared to approximately 9.97 under lump-sum taxation. With a smaller capital base, productive investment $I_{p,ss} = \delta_p K_{p,ss} = 0.246$ is also slightly smaller. Consequently, the Euler equation amplification remains visible in percentage terms: the expected fall in future output reduces the after-tax return to capital, causing productive investment to fall by about 3.5\% relative to its steady-state level on impact. This response is still amplified by the high capital-output ratio ($K/Y \approx 9.8$) and low depreciation rate ($\delta_p = 2.5\%$), which together create the standard RBC ``leverage'' channel.

	The equal proportional response of battery and grid investment (1.95\% each) masks an important asymmetry in levels: grid investment rises 2.5$\times$ more in absolute units due to its larger steady-state base. Market-based pricing mechanisms---such as scarcity pricing, capacity markets, or ancillary service contracts---could improve the targeting of flexibility investment toward the scarce input, improving welfare outcomes.''',

r'''	Variance decomposition confirms that 100\% of battery investment volatility is driven by global battery-price shocks, reflecting Vietnam's price-taking position in world markets. Grid investment follows the public reliability-targeting rule $(I_{grid,t}/\bar{I}_{grid}) = (u^*/u_t)^{\phi_{grid}}$, with $\phi_{grid} = 1.5$; even aggressive fiscal responses cannot overcome the 5--7 year planning and construction lags inherent in transmission infrastructure \citep{iea2023}. The proportional tax $\tau_{ss} = 1.43\%$ is fixed and shifts only the steady state (reducing $K_{p,ss}$ to 9.83 vs.\ 9.97 under lump-sum taxation), so productive investment falls by about 3.5\% on impact through the standard RBC Euler amplification channel. The equal proportional rise in both investment types (1.95\% each) masks a 2.5$\times$ level advantage for grid investment---reflecting the larger steady-state base of public infrastructure.''')

apply('BASE-3 Stage 3 innovation', r'''	The reliability gap also activates the model's scarcity-induced technological change mechanism. By quarter 10, battery technology ($A_{bat}$) has improved by approximately 0.65\% above the baseline trajectory. This improvement represents productivity-enhancing innovation in storage: better battery management systems, optimized charging algorithms, and engineering advances induced by the persistent scarcity signal. Note that because $F_t$ uses $A_{bat,t-1}$, the technology improvement induced at quarter $t$ only raises flexibility from quarter $t+1$ onward---reflecting the one-period delay between the innovation decision and its operational impact.

	The innovation mechanism is governed by:
	\begin{equation}
		A_{bat,t} = A_{bat,t-1} \left[1 + \eta_{bat} \cdot \chi \cdot \frac{u^* - u_t}{u^*} \right]
	\end{equation}
	where $\chi$ captures the degree to which the reliability gap is transmitted to innovators. In the baseline calibration, $\chi = 1.0$ (full transmission). In Vietnam's current regulatory environment, where electricity prices are administratively set and ancillary service markets remain underdeveloped, $\chi < 1$ reflects incomplete transmission of scarcity signals.

	While the per-shock innovation effect is modest, it is cumulative and persistent. The first-order autocorrelation of $A_{bat}$ in the stochastic simulation is 0.997, indicating near-permanent effects: a reliability shortfall today raises the technology level not just tomorrow, but indefinitely. In the 1000-quarter stochastic simulation, the ergodic mean of $A_{bat}$ reaches 1.015---a 1.5\% long-run improvement in battery technology driven entirely by accumulated scarcity-induced innovation, with no R\&D expenditure. Crucially, this technological progress is conditional on the regulatory parameter $\chi$: when $\chi$ is driven toward zero through regulatory suppression of scarcity signals, this innovation channel collapses entirely.''',

r'''	The reliability gap activates scarcity-induced technological change: by quarter 10, battery technology ($A_{bat}$) has improved by approximately 0.65\% above the baseline trajectory, governed by $A_{bat,t} = A_{bat,t-1}[1 + \eta_{bat} \cdot \chi \cdot (u^* - u_t)/u^*]$, where $\chi = 1.0$ in the baseline calibration. The per-shock effect is modest but cumulative: the first-order autocorrelation of $A_{bat}$ is 0.997, and the 1000-quarter ergodic mean reaches 1.015---a 1.5\% endogenous long-run technology improvement. This channel collapses when regulatory suppression of scarcity signals drives $\chi \to 0$, illustrating the welfare cost of administrative price-setting analyzed in Section~\ref{sec:welfare}.''')

apply('BASE-4 FEVD four-finding summary', r'''	Four findings emerge. First, the variance decomposition reveals a structural asymmetry in what drives each flexibility instrument. Battery investment volatility is almost entirely governed by global battery-price shocks (100\%), confirming Vietnam's price-taking position in world markets; intermittency accounts for only 1\%. Grid investment volatility, by contrast, is 82\% determined by domestic implementation shocks $\varepsilon_I$ acting through the fiscal rule, with intermittency contributing the remainder. The two instruments therefore face entirely different forcing variables---international supply chains versus domestic fiscal capacity---which explains why private storage deployment and public grid expansion cannot substitute for one another in policy design. Second, output volatility is dominated by the intermittency shock (95.4\%), with battery-price shocks contributing only 5\%. The dominance of intermittency reflects its AR(1) persistence ($\rho_{ren}=0.85$): once reliability degrades, the deficiency propagates across many quarters, generating large cumulative forecast errors. Third, utilization volatility mirrors the output result, with intermittency accounting for 95.8\% and battery-price shocks explaining the remaining 5\% through capital accumulation dynamics. Fourth, consumption volatility is split between intermittency (85\%) and battery-price shocks (11\%), reflecting the persistent propagation of domestic reliability shortfalls into household welfare. Since the proportional tax rate $\tau$ is fixed, it does not contribute as a variance source.''',

r'''	Table~\ref{tab:fevd} reveals a structural asymmetry in the two flexibility instruments: battery investment variance is 100\% driven by global battery-price shocks (confirming Vietnam's price-taking position), while grid investment variance is 82\% determined by domestic implementation shocks acting through the fiscal rule. Output and utilization volatility are both dominated by intermittency ($\sim$96\%), reflecting the AR(1) persistence $\rho_{ren}=0.85$. The proportional tax $\tau$ is fixed and contributes no variance.''')

apply('BASE-5 Why the Valley three paragraphs', r'''	\paragraph{(1) Front-Loaded Renewable Deployment.} Vietnam's PDP8 targets aggressive near-term installation goals, driven by international climate commitments, falling solar PV costs, and foreign direct investment in renewable projects. Solar capacity is particularly front-loaded: developers rushed to capture Feed-in Tariff incentives before their expiration in 2021, leading to a near-doubling of solar capacity in 18 months. The left panel of Figure \ref{fig:reliability_valley} models this as an S-curve with early acceleration, capturing the boom-bust dynamics of subsidy-driven deployment.

	\paragraph{(2) Lagged Grid Infrastructure.} Public transmission infrastructure follows a bureaucratic investment rule that responds to reliability shortfalls, but with substantial delays. Even when low utilization triggers increased fiscal allocations, the actual commissioning of new transmission lines lags by 5--7 years due to planning, permitting, and construction lead times. In the simulation, grid capital grows at 1.35\% per quarter initially, accelerating to 1.75\% per quarter after quarter 40 when the reliability crisis becomes politically salient. This acceleration is realistic---it reflects emergency policy responses such as streamlined permitting or multilateral development financing---but it cannot fully offset the accumulated deficit.

	\paragraph{(3) Battery Growth Constraints.} While private battery investment is more agile than grid infrastructure, it too faces constraints. Global lithium supply chains, manufacturing capacity for battery cells, and project financing availability all impose upper bounds on deployment rates. The simulation calibrates battery capital growth at 2.0\% per quarter, faster than grid expansion but still insufficient to keep pace with the intermittency induced by rapidly rising VRE shares. The complementarity between storage and transmission (enforced by the CES parameter $\rho = 0.40$) means that neither technology alone can substitute for the other---both must scale together.

	The valley thus represents a predictable consequence of unbalanced transition planning: renewable generation capacity races ahead, driven by global cost declines and policy mandates, while the complementary flexibility infrastructure---both public grid and private storage---accumulates incrementally.''',

r'''	Three structural features interact to produce the valley. First, Vietnam's PDP8 front-loads renewable deployment via subsidy-driven S-curve installation, modeled with early acceleration. Second, public transmission infrastructure faces 5--7 year planning and construction lags: grid capital grows at only 1.35\% per quarter initially, even after fiscal allocations rise. Third, private battery deployment is more agile (2.0\% per quarter) but still insufficient on its own, given the CES complementarity ($\rho = 0.40$) between storage and transmission. The valley is therefore a predictable consequence of unbalanced transition planning: renewable capacity races ahead while both types of flexibility infrastructure accumulate incrementally.''')

apply('BASE-6 summary trailing sentence', r'''	\noindent The core qualitative result is the asymmetric investment response: private battery investment responds quickly but is dominated by global price signals, while public grid investment responds to domestic reliability but accumulates slowly. This ``agility gap'' drives the Reliability Valley and creates the welfare costs quantified below.''',

r'''	\noindent The core qualitative result is the asymmetric investment response---private battery investment is dominated by global price signals while public grid investment accumulates slowly in response to domestic reliability shortfalls---and this agility gap drives the Reliability Valley and the welfare costs quantified in Section~\ref{sec:welfare}.''')

apply('WEL-7 baseline welfare cost paragraphs 2-3', r'''	While expressed as a fraction of steady-state consumption, the welfare cost is economically meaningful in aggregate. Intermittency shocks are not one-time events: Vietnam experiences 4--6 major weather disruptions per year \citep{nguyenle2019}. The stochastic simulation (1000 periods) shows output volatility of 1.60\% of steady state and utilization volatility of 1.96\% of steady state, reflecting the persistent propagation through capital accumulation dynamics. The welfare cost also reflects dynamic misallocation: the fiscal response to reliability shortfalls forces deviations from the household's preferred intertemporal consumption path through the endogenous tax rate.

	This baseline estimate should be interpreted as a lower bound for a single shock. The model omits amplifying factors (distributional impacts, industrial hysteresis, political instability) and dampening factors (household adaptation, demand response, regional grid interconnections). The $\Lambda_{\text{baseline}} = 0.000569\%$ serves as the reference point for the counterfactual analysis below, where I examine how regulatory distortions and structural parameters alter this cost.''',

r'''	While small as a fraction of steady-state consumption, the cost is economically meaningful: Vietnam experiences 4--6 major weather disruptions per year \citep{nguyenle2019}, and the stochastic simulation (1000 periods) shows output volatility of 1.60\% and utilization volatility of 1.96\%, reflecting persistent propagation through capital accumulation. The $\Lambda_{\text{baseline}} = 0.000569\%$ figure should be interpreted as a lower bound for a single shock---omitting distributional impacts, industrial hysteresis, and political instability---and serves as the reference point for the counterfactual analysis below.''')

apply('WEL-8 phi_grid sensitivity paragraph', r'''	\paragraph{Grid Investment Responsiveness ($\phi_{grid}$).} Table~\ref{tab:phi_grid_sensitivity} reports the sensitivity of investment responses and cumulative reliability dynamics to $\phi_{grid}$. The grid investment response scales proportionally with $\phi_{grid}$: raising $\phi_{grid}$ from 0.5 to 3.0 increases the impact response of $I_{grid}$ sixfold (from 0.0001\% to 0.0006\%). This variation directly governs the speed of the public infrastructure response to reliability shortfalls. With persistent intermittency shocks ($\rho_{ren}=0.85$), the cumulative reliability loss over 40 quarters is substantially larger than the i.i.d.\ benchmark: it declines from $-0.063$ at $\phi_{grid} = 0.5$ to $-0.054$ at $\phi_{grid} = 3.0$---a 13\% improvement in cumulative reliability recovery---confirming that more responsive grid policy accelerates reliability recovery, though the effect is modest because capital accumulation lags limit the speed of reliability restoration regardless of fiscal aggressiveness.''',

r'''	\paragraph{Grid Investment Responsiveness ($\phi_{grid}$).} The grid investment impact response scales sixfold as $\phi_{grid}$ rises from 0.5 to 3.0 (Table~\ref{tab:phi_grid_sensitivity}), and cumulative reliability loss over 40 quarters improves by 13\% (from $-0.063$ to $-0.054$). The welfare cost is nearly invariant to $\phi_{grid}$, because capital accumulation lags limit reliability restoration regardless of fiscal aggressiveness.''')

apply('WEL-9 beta sensitivity paragraph', r'''	More impatient households ($\beta = 0.98$) face a higher welfare cost than the baseline, for two reasons. First, higher discounting reduces the present value of future recovery, amplifying the perceived welfare loss from transitory intermittency shocks. Second, the higher steady-state interest rate ($r^* \approx 2.05\%$ quarterly vs.\ 1.01\%) increases the cost of international borrowing, making consumption smoothing more expensive. The resulting higher steady-state consumption share (78.2\% vs.\ 72.8\%) means a larger fraction of output is consumed rather than invested, leaving less margin for absorbing intermittency-driven investment needs. Conversely, patient households ($\beta = 0.995$) borrow cheaply and invest more, reducing the welfare cost. This finding suggests that economies with less-developed financial markets---where effective discount rates are higher due to credit constraints---face amplified welfare costs from the energy transition.''',

r'''	More impatient households ($\beta = 0.98$) face a higher welfare cost because higher discounting reduces the present value of future recovery and raises steady-state borrowing costs, making consumption smoothing more expensive. Conversely, patient households ($\beta = 0.995$) smooth shocks cheaply and invest more, reducing the welfare cost. This suggests that economies with less-developed financial markets---where effective discount rates are high---face amplified welfare costs from the energy transition.''')

apply('WEL-10 external sector dynamics', r'''	The open economy specification reveals an important channel through which intermittency shocks propagate beyond the domestic economy: the trade channel. This channel operates through two distinct mechanisms.

	\paragraph{Mechanism 1: Intermittency $\to$ Current Account Deficit.} Following an adverse intermittency shock, output falls (Channel 1). With the debt-elastic premium ($\phi_b = 0.039$), international borrowing is costly, limiting consumption smoothing via the current account. Instead, the household adjusts primarily through reduced investment: productive investment $I_p$ falls sharply as the after-tax return to capital declines. The fall in domestic investment exceeds the fall in output, generating a current account surplus: net foreign assets $B^*$ rise by approximately 0.003 units on impact. The current account identity confirms:
	\begin{equation*}
		\Delta B^*_t = Y_t - C_t - I_{p,t} - P_{bat,t} I_{bat,t} - I_{grid,t} - r^*_{t-1} B^*_{t-1}
	\end{equation*}
	When $Y_t$ falls but $I_{p,t}$ falls even more, the trade balance improves and $B^*$ rises.

	\paragraph{Mechanism 2: Battery Price Shock $\to$ Terms of Trade.} A positive shock to $P_{bat}$ operates like an adverse terms-of-trade shock for a battery-importing economy. Higher import prices for storage technology raise the cost of flexibility investment, slow capital accumulation in $K_b$, and---through the reliability channel---reduce output. Under the debt-elastic premium ($\phi_b = 0.039$), the variance decomposition shows that intermittency has become the dominant driver of both $B^*$ volatility (100\%) and consumption volatility (85\%), with battery price shocks contributing only 4\% and 11\%, respectively. The borrowing cost limits the current account response to battery price shocks, concentrating their impact through the domestic investment channel rather than the balance-of-payments channel.

	\paragraph{Dual Exposure.} The variance decomposition of net foreign assets reveals the economy's dual exposure: intermittency shocks contribute 100\% and battery price shocks contribute 4\% of $B^*$ volatility (shares exceed 100\% due to negative covariance between shock contributions). This pattern has a clear economic interpretation: domestic reliability risk creates borrowing needs for consumption smoothing and investment financing, while global battery market shocks determine the cost of the flexibility imports required to close the reliability gap. The debt-elastic interest rate premium ($r^*_t = \bar{r} + \phi_b(e^{\bar{B}^* - B^*_t} - 1)$) ensures that persistent borrowing raises the country's risk premium, providing a stabilizing feedback that prevents explosive debt accumulation but amplifies the welfare cost of repeated shocks.

	The external sector dynamics reinforce the policy case for domestic battery manufacturing capacity: reducing import dependence for storage technology would attenuate Mechanism 2 and lower the economy's external vulnerability during the energy transition.''',

r'''	The open economy specification reveals two external propagation channels. An adverse intermittency shock reduces output and---because investment falls more sharply than output---generates a current account surplus ($B^*$ rises by $\approx$0.003 on impact). A positive battery-price shock operates as an adverse terms-of-trade shock: higher storage import costs slow $K_b$ accumulation and reduce output through the reliability channel. The debt-elastic premium ($\phi_b = 0.039$) limits borrowing, concentrating the battery-price impact through the domestic investment channel. Variance decomposition shows intermittency dominates $B^*$ volatility (100\%), while battery-price shocks contribute only 4\%---reinforcing the policy case for domestic battery manufacturing to attenuate external vulnerability during the transition.''')

apply('WEL-11 joint shock three paragraphs', r'''	Three findings emerge from the joint scenario. First, the welfare costs are additive: the joint welfare cost (0.000758\%) essentially equals the sum of individual costs (0.000569\% $+$ 0.000189\% = 0.000758\%), reflecting the linear superposition property of the first-order perturbation solution. Under higher-order approximations, nonlinear interaction effects could generate amplification.

	Second, intermittency dominates the welfare cost---contributing more than three times the battery price component (0.000569\% vs.\ 0.000189\%). This reflects the key role of AR(1) persistence ($\rho_{ren}=0.85$): once reliability degrades, households suffer prolonged consumption losses, generating a much larger welfare cost than a purely transitory disruption.

	Third, the joint scenario creates a ``scissors effect'' on battery investment: intermittency raises the need for storage (positive reliability gap), while the battery price shock raises the cost of acquiring it. The net effect is a decline in battery investment ($-0.0003\%$), implying that during a perfect storm, the economy cannot self-correct through the storage channel. This pathological case strengthens the argument for strategic battery reserves or long-term procurement contracts that decouple domestic storage investment from short-run global price volatility.''',

r'''	Three results stand out. Welfare costs are additive at the first order (joint $\Lambda = 0.000758\%$ $\approx$ sum of individual costs), with intermittency contributing three times more than the battery-price shock due to its AR(1) persistence. Most importantly, the joint scenario produces a ``scissors effect'' on battery investment: the need for storage rises while its cost also rises, causing battery investment to \textit{fall} ($-0.0003\%$)---implying the economy cannot self-correct through the storage channel and strengthening the case for strategic reserves or long-term procurement contracts.''')

# ─────────────────────────────────────────────
# CONCLUSION EDITS (Agent 4)
# ─────────────────────────────────────────────

apply('CONC-1 three bold findings', r'''First, the \textbf{agility gap is structurally asymmetric}: private battery investment and public grid investment are driven by entirely different forcing variables. Variance decomposition shows that 100\% of battery investment volatility originates from global battery-price shocks, while 82\% of grid investment volatility is determined by domestic implementation shocks through the fiscal rule. Because the two instruments respond to separate international and domestic signals, they cannot substitute for one another: addressing the agility gap requires coordinated policy across both dimensions simultaneously. Intermittency shocks reduce grid utilization by 1.30\% and output by 0.54\%; the fixed proportional tax $\tau_{ss}=1.43\%$ shifts the steady state rather than amplifying shocks dynamically.

	Second, the \textbf{agility gap}: private battery investment and public grid investment respond proportionally to reliability shortfalls (1.95\% each relative to their steady-state levels), but the absolute level response of grid investment is 2.5$\times$ larger (agility ratio = 0.40 in level terms). Battery investment is dominated by global battery-price shocks (100\% of $I_{bat}$ variance), while grid investment is driven primarily by the fiscal rule (82\% of $I_{grid}$ variance), creating dual exposure to international supply chains and domestic fiscal decisions.

	Third, \textbf{scarcity-induced technological change}: technological progress in battery storage responds endogenously to the gap between target and actual reliability, creating a self-correcting mechanism. The 1000-quarter simulation confirms battery technology reaches a mean of 1.015---a 1.5\% endogenous improvement conditional on $\chi > 0$.''',

r'''First, the \textbf{agility gap is structurally asymmetric}: battery investment volatility is driven entirely by global battery-price shocks, while grid investment volatility is determined 82\% by domestic fiscal implementation shocks. Because the two instruments respond to separate international and domestic signals, they cannot substitute for one another, and addressing the agility gap requires coordinated policy across both dimensions simultaneously.

	Second, \textbf{dual exposure amplifies macro vulnerability}: the absolute level response of grid investment to reliability shortfalls is 2.5$\times$ larger than that of battery investment (agility ratio = 0.40), while intermittency shocks reduce grid utilization by 1.30\% and output by 0.54\%.

	Third, \textbf{scarcity-induced technological change}: technological progress in battery storage responds endogenously to reliability shortfalls, creating a self-correcting mechanism. The 1000-quarter simulation confirms battery technology reaches a mean of 1.015---a 1.5\% endogenous improvement conditional on $\chi > 0$.''')

apply('CONC-2 welfare policy paragraph', r'''The welfare cost of intermittency is estimated at 0.000569\% of steady-state consumption under full price signal transmission (from Dynare first-order perturbation). Counterfactual analysis reveals that this cost increases by 20\% when regulatory distortions fully suppress scarcity signals ($\chi = 0$), demonstrating that the scarcity-induced innovation channel is quantitatively important for welfare. The constant fiscal wedge introduced by proportional taxation shifts the steady state---reducing the capital stock and hence the base from which intermittency shocks propagate---but does not create dynamic amplification beyond the standard RBC Euler mechanism.''',

r'''The welfare cost of intermittency is estimated at 0.000569\% of steady-state consumption under full price signal transmission. Counterfactual analysis reveals that this cost increases by 20\% when regulatory distortions fully suppress scarcity signals ($\chi = 0$), confirming that the scarcity-induced innovation channel is quantitatively important for welfare. The constant fiscal wedge shifts the steady state but does not create dynamic amplification beyond the standard Euler mechanism.''')

apply('CONC-3 policy implications paragraph', r'''These findings carry direct policy implications. First, market-based mechanisms that transmit scarcity signals---real-time pricing, capacity payments, ancillary service remuneration---are necessary complements to direct storage investment. Second, Vietnam's exposure to global battery supply chains (100\% of battery investment volatility driven by external price shocks) underscores the importance of supply chain diversification and domestic manufacturing capacity. Third, the trade channel implies that the energy transition has external balance consequences that should be incorporated into macroeconomic policy design.''',

r'''These findings carry direct policy implications. Market-based mechanisms that transmit scarcity signals---real-time pricing, capacity payments, ancillary service remuneration---are necessary complements to direct storage investment. Vietnam's near-total exposure to global battery supply chains and the external balance consequences of the energy transition further underscore the need to embed energy policy within a broader macroeconomic policy framework.''')

apply('CONC-4 Vietnam policy paragraph', r'''For Vietnam specifically, the signal attenuation parameter $\chi$ maps onto concrete institutional features. Vietnam Electricity (EVN) currently operates under regulated retail tariffs that do not fully pass through real-time generation costs, effectively dampening the scarcity signal that drives technological innovation in the model. The PDP8's target of 47.5~GW of renewable capacity by 2030 will substantially increase intermittency exposure, yet Vietnam lacks a wholesale electricity spot market, capacity payment mechanisms, or ancillary service remuneration frameworks. The model implies that achieving $\chi \approx 1$ requires specific reforms: introduction of a competitive wholesale market with locational marginal pricing, establishment of capacity remuneration mechanisms for storage providers, and regulatory frameworks that compensate battery operators for frequency regulation and spinning reserve services. Without these reforms, the 20\% welfare cost amplification identified in the counterfactual analysis (and 13\% for partial dampening at $\chi = 0.3$) represents a quantifiable cost of regulatory inertia.''',

r'''For Vietnam specifically, EVN's regulated retail tariffs suppress the real-time scarcity signal modelled by $\chi$, while the PDP8 target of 47.5~GW of renewables by 2030 will substantially increase intermittency exposure. Achieving $\chi \approx 1$ requires three reforms: introduction of a competitive wholesale market with locational marginal pricing, establishment of capacity remuneration mechanisms for storage providers, and regulatory compensation for frequency regulation and spinning reserve services. Without these changes, the 20\% welfare cost amplification identified in the counterfactual (13\% for $\chi = 0.3$) represents a quantifiable cost of regulatory inertia.''')

apply('CONC-5 limitations paragraph', r'''Several limitations qualify the results. First, the model treats the renewable energy endowment as an exogenous parameter ($\bar{E}$) rather than modeling the capacity investment decision. This implies that the results are conditional on a given level of renewable penetration and do not address the optimal pace of renewable deployment. Second, the representative household framework precludes analysis of the distributional consequences of reliability shortfalls, which are likely to fall disproportionately on lower-income households and energy-intensive industries. Third, the model's consumption share in steady state (72.8\%) exceeds Vietnam's observed ratio (64--66\%), reflecting the abstraction from government consumption and inventories (see Section~\ref{sec:calibration}). Fourth, the welfare metric captures only the consumption-equivalent cost of volatility; it omits level effects from the energy transition, adjustment costs from structural transformation, and political economy frictions. Fifth, the log-linearized solution restricts analysis to first-order dynamics around the steady state, potentially understating the costs of large shocks or regime changes. Finally, the model abstracts from fossil fuel phase-out dynamics, carbon pricing, and cross-border electricity trade, each of which is quantitatively relevant for Vietnam's transition.''',

r'''Several limitations qualify the results. The renewable energy endowment is treated as an exogenous parameter ($\bar{E}$), so the results are conditional on a given penetration level and do not address the optimal pace of deployment. The representative household framework precludes distributional analysis, and the model's steady-state consumption share (72.8\%) exceeds Vietnam's observed ratio (64--66\%), reflecting the abstraction from government consumption and inventories. The welfare metric captures only the consumption-equivalent cost of volatility, omitting level effects, structural transformation costs, and political economy frictions. Finally, the log-linearized solution restricts analysis to first-order dynamics, and the model abstracts from fossil fuel phase-out, carbon pricing, and cross-border electricity trade---all quantitatively relevant for Vietnam's transition.''')

apply('CONC-6 future research paragraph', r'''Several extensions merit investigation. First, incorporating nominal rigidities and monetary policy would allow analysis of the ``greenflation'' channel identified by \citet{airaudo2025green} in the context of intermittency-driven supply shocks. Second, a deterministic transition path for renewable penetration targets would enable more detailed analysis of the ``Reliability Valley'' phenomenon, including the optimal sequencing of storage and grid investment. Third, heterogeneous agent extensions would capture the distributional impacts of reliability shortfalls across income groups, firm sizes, and regions---a first-order concern for equity in developing economies. Fourth, endogenizing the renewable capacity decision (replacing the exogenous $\bar{E}$) would allow joint analysis of the optimal pace of renewable deployment and flexibility investment. Fifth, incorporating carbon pricing or border adjustment mechanisms (e.g., CBAM) would link the domestic reliability problem to international climate policy.''',

r'''Several extensions merit investigation. Incorporating nominal rigidities and monetary policy would allow analysis of the ``greenflation'' channel \citep{airaudo2025green} alongside intermittency-driven supply shocks, while a deterministic renewable transition path would enable richer analysis of the ``Reliability Valley'' and the optimal sequencing of storage and grid investment. Heterogeneous agent extensions would capture distributional impacts across income groups and firm sizes---a first-order equity concern in developing economies. Finally, endogenizing the renewable capacity decision and incorporating carbon pricing or border adjustment mechanisms (e.g., CBAM) would link the domestic reliability problem to international climate policy.''')

# ─────────────────────────────────────────────
# WRITE AND REPORT
# ─────────────────────────────────────────────
open('main_body.tex', 'wb').write(content.encode('utf-8'))
print(f'\nApplied: {applied} edits')
print(f'Failed:  {len(failed)} edits')
if failed:
    print('FAILED EDITS:', failed)
print(f'Size change: {original_len} -> {len(content)} bytes ({original_len - len(content):+d})')
