"""Weight of the completion loop and the amplitude c (26.11.2026).
Rules: loop_weight_rules_2026-11-26.md (pushed before run). Record: Р-108.
Run: python -I -X utf8 loop_weight.py
"""
from fractions import Fraction as F
from itertools import product
from math import pi, sqrt

# 1. loop vs eight: pairs of partners
for n in (10, 100, 1000, 10000):
    eight = sum(a * (n - a) for a in range(n + 1)) / (n + 1)       # O at every position, cross pairs n_A*n_B
    loop = n * (n - 1) / 2                                         # all pairs in one loop
    print(f'n = {n:5d}: eight {eight/n**2:.5f} n^2, loop {loop/n**2:.5f} n^2, ratio {loop/eight:.4f}')

# 2. what data say about P/A (R-101)
c3, s3, c5, s5 = 0.255, 0.018, 0.257, 0.037
r = c5 / c3; sr = r * sqrt((s3 / c3)**2 + (s5 / c5)**2)
PA = lambda r: (r - 1) / (5 * r - 3)
print(f'\nc(5)/c(3) = {r:.3f} +- {sr:.3f}; P/A = {PA(r):+.3f}; upper (r + 2 sigma = {r+2*sr:.3f}): P/A <= {PA(r+2*sr):.3f}')

# 3. branches
kappa = 2 / (3 * pi)
axes = {'u': 2, 'd': 1, 's': 1, 'c': 2, 'b': 1}
types = {3: 'uds', 4: 'udsc', 5: 'udscb'}
part = {'p1 own axes': lambda k: k, 'p2 missing': lambda k: 3 - k, 'p3 one axis': lambda k: 1}
total = 0; passed = []; best = []
for R, V, (Sn, unit), (pn, pa), (kn, kv) in product((1, 3), (3, 6, 8), (('S1/3', F(1, 3)), ('S1', F(1))),
                                                   part.items(), (('k=1', 1.0), ('k=1/2', 0.5), ('k=1/pi', 1 / pi))):
    total += 1
    A = R * V * unit**2
    P = {n: sum(3 * (pa(axes[t]) * unit)**2 for t in types[n]) for n in (3, 4, 5)}
    d = {n: A - P[n] for n in P}
    if min(d.values()) <= 0:
        continue
    c = {n: kv / (kappa * float(d[n])) for n in (3, 5)}
    p3_, p5_ = (c[3] - c3) / s3, (c[5] - c5) / s5
    ok = abs(p3_) <= 2 and abs(p5_) <= 2
    tag = f'L{1 if R == 3 else 0} V{V} {Sn} {pn:11s} {kn:6s}: A = {str(A):5s} A-P = {float(d[3]):6.2f}/{float(d[5]):6.2f}  c(3) {c[3]:.3f} ({p3_:+.1f}s) c(5) {c[5]:.3f} ({p5_:+.1f}s)'
    best.append((abs(p3_) + abs(p5_), tag + (' PASS' if ok else '')))
    if ok:
        passed.append(tag)
print(f'\nbranches {total}, with b > 0: {len(best)}, pass: {len(passed)}')
for _, t in sorted(best)[:8]:
    print('  ' + t)
