TEX = r"C:\Users\Laptop K1\OneDrive\Documents\MATLAB\tex\main_body.tex"
ROOT = r"C:\Users\Laptop K1\OneDrive\Documents\MATLAB\main_body.tex"
import shutil

with open(TEX, 'r', encoding='utf-8') as f:
    c = f.read()

def sr(content, start_anchor, end_anchor, replacement, label):
    s = content.find(start_anchor)
    if s < 0:
        print(f"  FAIL: {label} -- start not found")
        return content
    e = content.find(end_anchor, s)
    if e < 0:
        print(f"  FAIL: {label} -- end not found")
        return content
    e += len(end_anchor)
    print(f"  OK: {label} (removed {e-s} chars, replaced with {len(replacement)} chars)")
    return content[:s] + replacement + content[e:]

# FIX 1: SITC equation - missing \begin{equation}
old1_start = r'The law of motion is:_{\text{inherited knowledge}} \times \left(1 + \underbrace{\eta_{bat}}_{\text{innovation elasticity}} \cdot \underbrace{\chi}_{\substack{\text{signal}\\\text{transmission}}} \cdot \underbrace{\frac{u^* - u_t}{u^*}}_{\text{reliability gap}}\right)'
old1_end = r'\end{equation}'
new1 = 'The law of motion is:\n\t\\begin{equation}\n\t\tA_{bat,t} = A_{bat,t-1} \\times \\left(1 + \\eta_{bat} \\cdot \\chi \\cdot \\frac{u^* - u_t}{u^*}\\right)\n\t\t\\label{eq:learning_by_doing}\n\t\\end{equation}'
c = sr(c, old1_start, old1_end, new1, 'FIX1 SITC equation')

print(f"After FIX1: {len(c)} chars")

# FIX 2: Resource constraint - missing \begin{equation}
old2_start = 'Output equals domestic absorption plus net capital outflows: + P_{bat,t} I_{bat,t} + I_{grid,t} + '
old2_end = r'\end{equation}'
new2 = ('Output equals domestic absorption plus net capital outflows:\n'
        '\t\\begin{equation}\n'
        '\t\tY_t = C_t + I_{p,t} + P_{bat,t} I_{bat,t} + I_{grid,t} + '
        r'\underbrace{(B^*_t - (1+r^*_{t-1})B^*_{t-1})}_{\text{net capital outflow}}'
        '\n\t\t\\label{eq:goods_market_clearing}\n\t\\end{equation}')
c = sr(c, old2_start, old2_end, new2, 'FIX2 Resource constraint equation')

print(f"After FIX2: {len(c)} chars")

# FIX 3: Reliability linearization - stray } then \end{equation}
old3_start = r'reliability elasticity $\xi \equiv \frac{1-u_{ss}}{u_{ss}} \cdot \psi \frac{F_{ss}}{\phi_{int} \cdot \text{Vol}_{ren,ss}}$:'
old3_end = r'\end{equation}'
new3 = ('reliability elasticity $\\xi \\equiv \\frac{1-u_{ss}}{u_{ss}} \\cdot \\psi \\frac{F_{ss}}{\\phi_{int} \\cdot '
        r'\text{Vol}_{ren,ss}}$:'
        '\n\t\\begin{equation}\n'
        '\t\t\\hat{u}_t = -\\xi \\cdot \\hat{\\varepsilon}_{ren,t}\n'
        '\t\t\\label{eq:reliability_lin}\n'
        '\t\\end{equation}')
c = sr(c, old3_start, old3_end, new3, 'FIX3 Reliability linearization')

print(f"After FIX3: {len(c)} chars")

with open(TEX, 'w', encoding='utf-8') as f:
    f.write(c)
print("Written to tex/main_body.tex")

shutil.copy(TEX, ROOT)
print("Copied to root main_body.tex")
