"""Does the part-attraction amplitude c grow with the number of active part types n? (26.11.2026)
Rules: c_by_flavours_rules_2026-11-26.md (pushed before run). Reuses r_excess_shapes.py (regions, fit (b1) for the Z part).
Run: python -I -X utf8 c_by_flavours.py <path to rpp2022-hadronicrpp_page1001.dat>
"""
import io, os, runpy, contextlib
from math import sqrt, log

here = os.path.dirname(os.path.abspath(__file__))
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(os.path.join(here, 'r_excess_shapes.py'), run_name='shapes')
pts, models, solve = g['pts'], g['models'], g['solve']
basis = models['(b1) c/ln(E/0.3128) + Z']
A = [[sum(f(E) * h(E) / d**2 for _, E, e, d, _, _ in pts) for h in basis] for f in basis]
b = [sum(f(E) * e / d**2 for _, E, e, d, _, _ in pts) for f in basis]
cf = solve(A, b)
E0 = 0.3128
groups = {3: ('A1', 'A2', 'A3'), 4: ('B1', 'B2'), 5: ('C1', 'C2', 'C3')}
res = []
for n, names in groups.items():
    vals = []
    for name, E, e, d, _, _ in pts:
        if name in names:
            Z = cf[1] * basis[1](E) + cf[2] * basis[2](E)
            L = log(E / E0)
            vals.append(((e - Z) * L, d * L)); print(f'  {name}: E {E:.2f}, Z part {Z:+.4f}, c_i = {(e-Z)*L:.3f} +- {d*L:.3f}')
    w = [1 / s**2 for _, s in vals]; W = sum(w)
    cm = sum(wi * v for wi, (v, _) in zip(w, vals)) / W; ce = 1 / sqrt(W)
    res.append((n, cm, ce)); print(f'n = {n}: c = {cm:.3f} +- {ce:.3f}')
S = sum(1 / e**2 for _, _, e in res); Sx = sum(n / e**2 for n, _, e in res); Sy = sum(c / e**2 for _, c, e in res)
Sxx = sum(n * n / e**2 for n, _, e in res); Sxy = sum(n * c / e**2 for n, c, e in res)
D = S * Sxx - Sx**2; slope = (S * Sxy - Sx * Sy) / D; ss = sqrt(S / D)
print(f'slope dc/dn = {slope:+.4f} +- {ss:.4f} ({slope/ss:+.1f} sigma) -> ' +
      ('growth CONFIRMED' if slope > 2 * ss else 'growth REFUTED' if slope < -2 * ss else 'data do not decide'))
c3 = res[0][1]
print(f'for scale: relative growth per type {slope/c3*100:+.1f} % (standard theory ~ +8 %)')
