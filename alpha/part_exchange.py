"""Law of the part exchange in W+W- production (R-121). Rules: part_exchange_rules_2026-11-26.md.
dsigma/dcos ~ |t|^(-n) * ((1+cos)/2)^k, |t| ~ 1 + b^2 - 2 b cos. Run: python -P -E -X utf8 part_exchange.py lep_ww_data.txt"""
import sys
from math import sqrt
MW = 80.37
lines = [l.split() for l in open(sys.argv[1], encoding='utf-8') if l.strip() and not l.startswith('#')]
i = lines.index(['ANGULAR']) + 1; ang = []
while i < len(lines):
    ang.append((float(lines[i][0]), list(map(float, lines[i + 1])), list(map(float, lines[i + 2])))); i += 3

def binavg(f, k, n=400):
    lo = -1 + 0.2 * k
    return sum(f(lo + 0.2 * (j + 0.5) / n) for j in range(n)) / n

print(f'data A_FB: ' + ', '.join(f'{E:.0f}: {(sum(v[5:])-sum(v[:5]))/sum(v):.3f}' for E, v, e in ang) + '\n')
for nexp, nlab in ((0.5, 'n=1/2 (field linear)'), (1.0, 'n=1 (field squared = potential linear)'), (2.0, 'n=2 (potential squared)')):
    for kexp in (0, 1, 2):
        tot = 0.0; afb = []
        for E, v, e in ang:
            b = sqrt(1 - (2 * MW / E)**2)
            f = lambda c: (1 + b * b - 2 * b * c)**(-nexp) * ((1 + c) / 2)**kexp
            m = [binavg(f, j) for j in range(10)]
            w = [1 / x**2 for x in e]
            norm = sum(vi * mi * wi for vi, mi, wi in zip(v, m, w)) / sum(mi * mi * wi for mi, wi in zip(m, w))
            tot += sum(((vi - norm * mi) / ei)**2 for vi, mi, ei in zip(v, m, e))
            afb.append((sum(m[5:]) - sum(m[:5])) / sum(m))
        print(f'{nlab:42s} k={kexp}: chi2 = {tot:8.1f} / 36  A_FB {", ".join(f"{a:.3f}" for a in afb)} -> {"PASS" if tot <= 58.6 else "FAIL"}')
