"""How the coupling g enters the production excess (k): all branches (26.11.2026).
Rules: coupling_k_rules_2026-11-26.md (pushed before run). Record: Р-109.
Run: python -I -X utf8 coupling_k.py
"""
from fractions import Fraction as F
from itertools import product
from math import pi

kappa = 2 / (3 * pi)
c3, s3, c5, s5 = 0.255, 0.018, 0.257, 0.037
rows = []

def judge(name, cn):
    if any(v is None for v in cn.values()):
        rows.append((99, f'{name}: b <= 0 -> excluded', False)); return
    p3, p5 = (cn[3] - c3) / s3, (cn[5] - c5) / s5
    ok = abs(p3) <= 2 and abs(p5) <= 2
    rows.append((abs(p3) + abs(p5), f'{name}: c(3) = {cn[3]:.4f} ({p3:+.1f}s), c(5) = {cn[5]:.4f} ({p5:+.1f}s)' + ('  PASS' if ok else ''), ok))

for (kn, k), V, (Sn, s), Pn in product((('K1', 1.0), ('K2', 0.5)), (3, 6, 8), (('S1', F(1)), ('S1/3', F(1, 3))), ('P=0', 'P=A/V')):
    A = 3 * V * s**2
    PA = F(0) if Pn == 'P=0' else F(1, V)
    cn = {}
    for n in (3, 5):
        d = A * (1 - n * PA)
        cn[n] = None if d <= 0 else k / (kappa * float(d))
    judge(f'{kn} V{V} {Sn:4s} {Pn:5s}', cn)

fs = {'f6a 4/6': (F(4, 6), 6), 'f6b 1': (F(1), 6), 'f8a 4/8': (F(4, 8), 8), 'f8b 6/8': (F(6, 8), 8), 'f8c 2/3': (F(2, 3), 8), 'f3 2/3': (F(2, 3), 3)}
for (fn, (f, V)), Pn in product(fs.items(), ('P=0', 'P=A/V')):
    PA = F(0) if Pn == 'P=0' else F(1, V)
    cn = {}
    for n in (3, 5):
        d = 1 - n * PA
        cn[n] = None if d <= 0 else float(f) / (2 * float(d))
    judge(f'K3 {fn:8s} {Pn:5s}', cn)

print(f'branches {len(rows)}, pass {sum(r[2] for r in rows)}')
for _, t, _ in sorted(rows):
    print('  ' + t)
print(f'\nfor reference: pi/12 = {pi/12:.4f}; 1/4 = 0.25')
