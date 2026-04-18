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
    except ValueError:
        failed.append(label); print(f'  FAIL: {label}')

# Remove redundant FEVD summary paragraph (duplicate of table caption)
S('6-fevd repeat paragraph',
    'Table~\\ref{tab:fevd} reveals a structural asymmetry in the two flexibility instruments: battery investment variance is 100\\% driven by global battery-price shocks (confirming Vietnam\'s price-taking position), while grid investment variance is 82\\% determined by domestic implementation shocks acting through the fiscal rule. Output and utilization volatility are both dominated by intermittency ($\\sim$96\\%), reflecting the AR(1) persistence $\\rho_{ren}=0.85$. The proportional tax $\\tau$ is fixed and contributes no variance.',
    '\n\n\t\\subsection{The Reliability Valley During PDP8 Transition}',
    '\n\n\t\\subsection{The Reliability Valley During PDP8 Transition}')

# Remove "Over time" redundant paragraph in Consumption section
S('6-over-time paragraph',
    'Over time, as battery and grid capital accumulate and reliability improves, all variables converge back toward their steady-state paths. The persistence of the response is governed by the capital accumulation dynamics and the AR(1) persistence of the intermittency shock ($\\rho_{ren} = 0.85$).',
    '\n\n\t\\subsection{Forecast Error Variance Decomposition}',
    '\n\n\t\\subsection{Forecast Error Variance Decomposition}')

# Trim Reliability Valley intro to 1 sentence (cut the second paragraph about S-curve)
S('6-valley intro trim',
    'Figure \\ref{fig:reliability_valley} illustrates the transition dynamics. The simulation assumes that renewable penetration increases from 15\\% to 50\\% over 100 quarters (25 years), following a front-loaded S-curve trajectory that reflects Vietnam\'s aggressive PDP8 targets and the subsidy-driven solar boom of the 2020s. Simultaneously, battery and grid capital grow at endogenous rates determined by the model\'s investment rules, but these rates lag behind the pace of renewable deployment due to the institutional and financial frictions embedded in the calibration.',
    '\n\n\t\\subsubsection{Valley Characteristics}',
    'Figure~\\ref{fig:reliability_valley} illustrates the transition dynamics as renewable penetration rises from 15\\% to 50\\% over 25 years.\n\n\t\\subsubsection{Valley Characteristics}')

# Remove the long Stage 2 paragraph that summarizes what the table already shows
S('6-stage2 prose redundancy',
    'Battery investment volatility is 100\\% driven by global battery price shocks, while grid investment follows the public reliability-targeting rule (equation~\\ref{eq:grid_investment_rule}) with parameter $\\phi_{grid} = 1.5$. The proportional tax $\\tau_{ss} = 1.43\\%$ is fixed: its effect operates through the steady-state channel (reducing $K_{p,ss}$ to 9.83 vs.\\ 9.97 under lump-sum taxation), not through dynamic amplification. The equal proportional response of battery and grid investment (1.95\\% each) masks an important asymmetry in levels: grid investment rises 2.5$\\times$ more in absolute units due to its larger steady-state base.',
    '\n\n\t\\subsubsection{Stage 3: Scarcity-Induced Innovation Activation}',
    '\n\n\t\\subsubsection{Stage 3: Scarcity-Induced Innovation Activation}')

open('main_body.tex', 'w', encoding='utf-8').write(content)
print(f'\nApplied {applied}/{applied+len(failed)} cuts')
if failed: print('Failed:', failed)
print(f'Size: {original_len} -> {len(content)} bytes (saved {original_len-len(content):,} chars)')
