#!/usr/bin/env python3
"""
round2_cuts.py – Apply remaining round-2 cuts to main_body.tex.
Uses exact anchor phrases verified against current file state.
"""

FILEPATH = r"C:\Users\Laptop K1\OneDrive\Documents\MATLAB\main_body.tex"

with open(FILEPATH, "r", encoding="utf-8") as f:
    content = f.read()

original_len = len(content)
results = []

def slice_replace(label, start_anchor, end_anchor, replacement, content, start_offset=0):
    try:
        si = content.index(start_anchor, start_offset)
    except ValueError:
        results.append(f"FAIL [{label}]: start anchor not found: {repr(start_anchor[:70])}")
        return content
    try:
        ei = content.index(end_anchor, si) + len(end_anchor)
    except ValueError:
        results.append(f"FAIL [{label}]: end anchor not found: {repr(end_anchor[:70])}")
        return content
    removed = ei - si
    content = content[:si] + replacement + content[ei:]
    results.append(f"OK   [{label}]: -{removed} chars +{len(replacement)} chars")
    return content

# ===========================================================
# B2 – Cut the Validity derivation block
# Start: \noindent \textbf{Validity...} (includes the \n before it)
# End: \subsubsection{Scarcity-Induced...} (no tab prefix in this file)
# ===========================================================

content = slice_replace(
    "B2",
    "\n\t\\noindent \\textbf{Validity of the linearization.}",
    "\n\\subsubsection{Scarcity-Induced Technological Change}",
    "\n\n\t\\noindent Full derivation in Appendix~\\ref{app:log_lin}.\n\n\\subsubsection{Scarcity-Induced Technological Change}",
    content
)

# ===========================================================
# B4a – Cut valley-of-death paragraph in Transmission Channels
# ===========================================================

content = slice_replace(
    "B4a",
    "\n\tThe transition dynamics during which the innovation channel is not yet mature",
    "\n\t\n\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}",
    "\n\n\t\t\\paragraph{Channel 3: The Fiscal Inertia Lag.}",
    content
)

# ===========================================================
# B4b – Cut Channel 3 overshooting sentences
# ===========================================================

content = slice_replace(
    "B4b",
    " When $\\phi_{grid}$ is small, the government under-responds to reliability crises",
    "\n\t\n\t\\paragraph{Channel 4: The Fiscal Wedge Channel.}",
    "\n\n\t\t\\paragraph{Channel 4: The Fiscal Wedge Channel.}",
    content
)

# ===========================================================
# D1 – Replace Efficiency Wedge subsubsection body
# ===========================================================

content = slice_replace(
    "D1",
    "\n\tThe steady-state results reveal a fundamental trade-off between energy reliability",
    "\n\t\n\t\\subsection{Impulse Response to Renewable Intermittency Shock}",
    "\n\tTable~\\ref{tab:steady_state} shows the steady-state values. The 3\\% reliability shortfall ($1-u_{ss}=0.03$) acts as an efficiency wedge, reducing effective capital services by $\\alpha \\times 3\\% = 1.05$~pp relative to perfect reliability. This structural drag is the starting point from which intermittency shocks propagate.\n\n\t\\subsection{Impulse Response to Renewable Intermittency Shock}",
    content
)

# ===========================================================
# D2 – Policy Counterfactual: trim to 1 paragraph
# Find from the beginning of the scenario list to end of sequencing sentence
# ===========================================================

D2_start = "\n\t\n\tThe three scenarios share the same S-curve renewable deployment path (15\\% $\\to$ 50\\% penetration over 100 quarters) but differ in the policy environment:\n\t\n\t\\begin{enumerate}\n\t\t\\item \\textbf{Baseline}: Current parameters ($\\phi_{\\text{grid}} = 1.5$, $P_{\\text{bat}} = 1.0$).\n\t\t\\item \\textbf{Scenario A --- Accelerated Grid Planning}: The government doubles the aggressiveness of its grid investment response ($\\phi_{\\text{grid}} = 3.0$), representing reforms such as streamlined permitting, pre-approved corridor designation, or dedicated grid investment funds.\n\t\t\\item \\textbf{Scenario B --- Battery Subsidy}: A 20\\% reduction in the effective battery price ($P_{\\text{bat}} = 0.8$), representing direct subsidies, tax credits, or concessional financing for storage deployment.\n\t\\end{enumerate}"

D2_end = "For Vietnam's PDP8, the model implies a clear sequencing: \\textit{first}, deploy battery subsidies to reduce the depth of the valley during the critical 2025--2035 decade; \\textit{second}, pursue grid planning reforms to ensure full recovery within the PDP8 horizon."

D2_repl = "\n\n\tTwo counterfactual scenarios share the baseline S-curve renewable path (15\\%$\\to$50\\% penetration over 100 quarters): (A) Accelerated Grid Planning ($\\phi_{\\text{grid}}=3.0$) and (B) Battery Subsidy ($P_{\\text{bat}}=0.8$, a 20\\% price reduction). Table~\\ref{tab:closing_gap} and Figure~\\ref{fig:closing_gap} summarise results. Battery subsidies reduce valley depth by 38.6\\% (from 15.7~pp to 9.6~pp), while accelerated grid planning shortens duration by 20\\% and achieves recovery by Q81. Neither policy alone eliminates the valley, confirming the complementarity proposition: addressing the agility gap requires coordinated investment in both flexibility instruments. Full counterfactual results are tabulated in Appendix~\\ref{app:solution}."

# Need to keep the table and figure between the list and the end phrase.
# Strategy: replace only the scenario list + three-findings block separately
try:
    si_list = content.index(D2_start)
    # Find the table that comes right after the list
    table_start = content.index("\n\t\\begin{table}[H]\n\t\t\\centering\n\t\t\\caption[Closing the Gap", si_list)
    # Replace just the scenario list with our new paragraph + keep the table
    content = content[:si_list] + "\n\n" + content[table_start:]
    results.append(f"OK   [D2-list]: removed scenario list, keeping table")

    # Now remove the three-findings narrative that comes after the figure
    # Find the three-findings block
    three_findings_start = "Three findings emerge from the counterfactual analysis."
    si_three = content.index(three_findings_start)
    end_three = content.index(D2_end, si_three) + len(D2_end)
    content = content[:si_three] + D2_repl.strip() + content[end_three:]
    results.append("OK   [D2-findings]: replaced three-findings narrative")
except ValueError as e:
    results.append(f"FAIL [D2]: {e}")

# ===========================================================
# E1a – Reliability Curvature: trim long paragraph to 2 sentences
# ===========================================================

E1a_start_phrase = "An important analytical result emerges from the sensitivity analysis: the first-order perturbation solution is invariant to $\\psi$."
E1a_end_phrase = "\n\n\t\\paragraph{Grid Investment Responsiveness ($\\phi_{grid}$).}"

E1a_new = "The first-order perturbation solution is invariant to $\\psi$: since $\\phi_{int}$ is calibrated jointly with $\\psi$ to hit $u_{ss}=0.97$, the linearization coefficient $\\xi = \\psi(1-u_{ss})/(\\phi_{int}\\cdot\\text{Vol}_{ren,ss})$ is unchanged. Higher $\\psi$ matters only for nonlinear curvature in the Reliability Valley (Section~\\ref{subsec:valley})."

try:
    si = content.index(E1a_start_phrase)
    ei = content.index(E1a_end_phrase, si) + len(E1a_end_phrase)
    content = content[:si] + E1a_new + content[ei - len(E1a_end_phrase):]
    results.append("OK   [E1a]: trimmed psi paragraph")
except ValueError as e:
    results.append(f"FAIL [E1a]: {e}")

# ===========================================================
# E1e – Cut paragraph after beta table
# ===========================================================

E1e_start_phrase = "\n\n\tMore impatient households ($\\beta = 0.98$) face a higher welfare cost because higher discounting reduces the present value of future recovery and raises steady-state borrowing costs, making consumption smoothing more expensive. Conversely, patient households ($\\beta = 0.995$) smooth shocks cheaply and invest more, reducing the welfare cost. This suggests that economies with less-developed financial markets---where effective discount rates are high---face amplified welfare costs from the energy transition."

E1e_end_phrase = "\n\n\t\t\\subsection{External Sector Dynamics}"

try:
    si = content.index(E1e_start_phrase)
    ei = content.index(E1e_end_phrase, si) + len(E1e_end_phrase)
    content = content[:si] + content[ei - len(E1e_end_phrase):]
    results.append("OK   [E1e]: removed beta after-table paragraph")
except ValueError as e:
    results.append(f"FAIL [E1e]: {e}")

# ===========================================================
# Write result
# ===========================================================

with open(FILEPATH, "w", encoding="utf-8", newline="\n") as f:
    f.write(content)

print(f"\n{'='*60}")
print(f"main_body.tex: {original_len} -> {len(content)} chars ({len(content)-original_len:+d})")
print(f"{'='*60}")
for r in results:
    print(r)
print(f"{'='*60}")
print("Done.")
