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

# Remove FEVD repeat paragraph (unique anchor: 'reveals a structural asymmetry')
S('6-fevd repeat',
    'Table~\\ref{tab:fevd} reveals a structural asymmetry',
    '$\\rho_{ren}=0.85$. The proportional tax $\\tau$ is fixed and contributes no variance.',
    '')

# Remove "Over time" paragraph
S('6-overtime',
    'Over time, as battery and grid capital accumulate and reliability improves, all variables converge back toward their steady-state paths.',
    'AR(1) persistence of the intermittency shock ($\\rho_{ren} = 0.85$).',
    '')

# Trim reliability valley intro 2nd paragraph
S('6-valley-long-intro',
    'Figure \\ref{fig:reliability_valley} illustrates the transition dynamics. The simulation assumes',
    'the institutional and financial frictions embedded in the calibration.',
    'Figure~\\ref{fig:reliability_valley} illustrates the transition dynamics as renewable penetration rises from 15\\% to 50\\% over 25 years.')

# Remove Stage 2 redundant prose paragraph
S('6-stage2-redundant',
    'Battery investment volatility is 100\\% driven by global battery price shocks, while grid investment follows the public reliability-targeting rule',
    'grid investment rises 2.5$\\times$ more in absolute units due to its larger steady-state base.',
    '')

open('main_body.tex', 'w', encoding='utf-8').write(content)
print(f'\nApplied {applied}/{applied+len(failed)} cuts')
if failed: print('Failed:', failed)
print(f'Size: {original_len} -> {len(content)} bytes (saved {original_len-len(content):,})')
