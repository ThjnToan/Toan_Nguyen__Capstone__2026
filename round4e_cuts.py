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

# 6.1: Remove duplicate SS table reference (first paragraph repeats the second)
S('6.1 duplicate SS ref',
    'Table \\ref{tab:steady_state} reports the model\'s steady state. The calibration targets a baseline reliability rate of 97.0\\%, consistent with Vietnam\'s official grid performance objective under PDP8. The battery-to-grid capital ratio of 0.167 reflects the structural imbalance in Vietnam\'s current investment allocation: approximately \\$2.4 billion earmarked for Battery Energy Storage Systems (BESS) versus \\$14.9 billion for transmission and grid infrastructure over the 2021--2030 period.\n\n\t',
    'Table~\\ref{tab:steady_state} (Appendix~\\ref{app:steadystate})',
    'Table~\\ref{tab:steady_state} (Appendix~\\ref{app:steadystate})')

# 6.2: Cut IRF intro paragraph (verbose description of what the shock is)
S('6.2 IRF intro verbose',
    'Figure \\ref{fig:irf_baseline} presents the impulse response functions to a 5\\% positive intermittency shock. This shock is interpreted as an unexpected deterioration in renewable generation conditions---for example, a week of persistent cloud cover reducing solar output, or a sustained lull in wind speeds during the monsoon transition period. The shock follows an AR(1) process with persistence parameter $\\rho_{ren} = 0.85$, consistent with empirical studies of renewable generation volatility in tropical climates \\citep{nguyenle2019}.',
    '\n\n\t\\begin{figure}[H]',
    'Figure~\\ref{fig:irf_baseline} shows the impulse responses to a 5\\% positive intermittency shock ($\\rho_{ren}=0.85$).\n\n\t\\begin{figure}[H]')

# 6.4.3 Policy Implications - check if verbose and trim
# Look for it
idx = content.find('subsubsection{Policy Implications}')
if idx >= 0:
    print(f"Found Policy Implications at {idx}")
    print(repr(content[idx:idx+400]))

open('main_body.tex', 'w', encoding='utf-8').write(content)
print(f'\nApplied {applied}/{applied+len(failed)} cuts')
if failed: print('Failed:', failed)
print(f'Size: {original_len} -> {len(content)} bytes (saved {original_len-len(content):,})')
