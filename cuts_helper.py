import sys
TEX = r"tex/main_body.tex"

def sr(c, s, e, n, lab):
    try:
        i = c.index(s)
        j = c.index(e, i) + len(e)
        print("  OK:", lab)
        return c[:i] + n + c[j:]
    except ValueError as ex:
        print("  FAIL:", lab, "-", ex)
        return c

with open(TEX, "r", encoding="utf-8") as f:
    c = f.read()

print("Lines before:", c.count(chr(10)))

# Save backup
with open(TEX + ".bak", "w", encoding="utf-8", newline=chr(10)) as f:
    f.write(c)
print("Backup saved")