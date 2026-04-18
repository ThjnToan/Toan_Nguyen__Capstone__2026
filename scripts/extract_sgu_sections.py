from pypdf import PdfReader
from pathlib import Path
import re

pdf = Path('sgu_w9270.pdf')
reader = PdfReader(str(pdf))

all_pages = []
for i, page in enumerate(reader.pages, start=1):
    txt = page.extract_text() or ''
    all_pages.append((i, txt))

keywords = [
    'calibrat', 'debt-elastic', 'premium', 'second moments', 'table 3',
    'current account', 'trade balance', 'net foreign asset', 'small open economy',
    'p(d)', 'psi', 'stationarity', 'mendoza'
]

hits = []
for p, txt in all_pages:
    low = txt.lower()
    if any(k in low for k in keywords):
        # Save full page text for manual checking
        hits.append((p, txt))

out = []
out.append(f'Total pages: {len(all_pages)}')
out.append(f'Pages with keyword hits: {len(hits)}')
out.append('')

for p, txt in hits:
    out.append(f'===== PAGE {p} =====')
    # Keep the page text, normalized lightly
    cleaned = re.sub(r'\n{3,}', '\n\n', txt)
    out.append(cleaned)
    out.append('')

Path('sgu_w9270_key_pages.txt').write_text('\n'.join(out), encoding='utf-8')

print(f'Wrote sgu_w9270_key_pages.txt with {len(hits)} matching pages')
