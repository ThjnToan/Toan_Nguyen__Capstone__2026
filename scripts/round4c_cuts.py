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
        failed.append(label); print(f'  FAIL: {label}')

# FIX gov dangling colon (tab+newline+tab between colon and subsection)
S('gov-fix',
    ', calibrated to balance the steady-state budget:\n\t\n\t\\subsection{Battery Price Dynamics}',
    '\\subsection{Battery Price Dynamics}',
    ' ($\\tau_{ss}=\\delta_g K_{g,ss}/Y_{ss}=0.0143$), ensuring tax revenue covers replacement grid investment.\n\n\t\\subsection{Battery Price Dynamics}')

# 3-5b battery post-eq paragraph
S('3-5b battery post-eq',
    'The exponent $\\phi_{grid}$ governs the elasticity of the investment response to the reliability gap---higher values imply more aggressive battery deployment when the grid is stressed.',
    'consistent with capital replacement at $I_{bat,ss} = \\delta_b K_{b,ss}$.',
    'Higher $\\phi_{grid}$ implies more aggressive deployment when the grid is stressed; the price term $1/P_{bat,t}$ reflects Vietnam\'s price-taking position in global markets.')

# 3-8 greenflation paragraph (uses Blanchard-Gali as anchor)
S('3-8 greenflation',
    'This is the ``greenflation\'\' channel: supply-side volatility from renewables generates simultaneous output losses and energy price increases.',
    'generating a policy-relevant trade-off between output stabilization and energy price stabilization.',
    '')

# 3-9 channel2 long -> 2 sentences
S('3-9 channel2',
    'The reliability decline widens the gap $(u^* - u_t)/u^*$, which accelerates scarcity-induced technological change in storage (Equation~\\ref{eq:abat_lin}). The resulting improvement in $A_{bat,t}$ increases the effective services per unit of battery capital in period $t+1$, gradually restoring flexibility and utilization.',
    'adjustment relies on costly capital accumulation rather than efficient technological improvement.',
    'The reliability decline widens the gap $(u^*-u_t)/u^*$, accelerating scarcity-induced innovation in storage and gradually restoring flexibility; the one-period delay means improvements enter $F_{t+1}$ rather than $F_t$.')

# 3-10 channel4 long -> 2 sentences
S('3-10 channel4',
    'This channel identifies the level distortion introduced by proportional taxation, which is absent in the lump-sum baseline.',
    'proportional taxation reduces the economy\'s structural resilience to intermittency shocks even without providing dynamic amplification.',
    'The fixed tax $\\tau$ imposes a permanent wedge: $K_{p,ss}$ is $\\approx$1.4\\% below the lump-sum benchmark, reducing structural resilience without dynamic amplification.')

# 5-1 core contribution (has tab+newline+tab before it)
S('5-1 core contribution',
    'A core contribution of this capstone is the rigorous calibration of the reliability constraint, which is micro-founded on the physical properties of the Vietnamese grid rather than arbitrary assumptions.\n\n\t',
    '\\subsubsection{Renewable Intermittency}',
    '\\subsubsection{Renewable Intermittency}')

# 5-5 chi paragraph
S('5-5 chi paragraph',
    'The price signal transmission parameter is set to $\\chi = 1.0$ in the baseline (full price transmission). The regulated-market counterfactual uses $\\chi = 0.3$',
    'consistent with Vietnam\'s fixed-tariff electricity market structure prior to the 2023 tariff reform.',
    'The baseline sets $\\chi=1.0$ (full signal transmission); the regulated counterfactual uses $\\chi=0.3$, calibrated to Vietnam\'s documented retail tariff pass-through of 25--35\\% \\citep{ieefa2022}.')

open('main_body.tex', 'w', encoding='utf-8').write(content)
print(f'\nApplied {applied}/{applied+len(failed)} cuts')
if failed: print('Failed:', failed)
print(f'Size: {original_len} -> {len(content)} bytes (saved {original_len-len(content):,} chars)')
