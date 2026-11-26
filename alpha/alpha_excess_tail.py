"""1/alpha far with the decreasing R excess above 40 GeV (26.11.2026). Rules: alpha_excess_tail_rules_2026-11-26.md (pushed before run).
Run: python -I -X utf8 alpha_excess_tail.py <path to rpp2022-hadronicrpp_page1001.dat>   (pure python; reuses r_excess_shapes.py)
"""
import sys, io, os, runpy, contextlib
from math import pi, sqrt, log

here = os.path.dirname(os.path.abspath(__file__))
with contextlib.redirect_stdout(io.StringIO()):
    g = runpy.run_path(os.path.join(here, 'r_excess_shapes.py'), run_name='shapes')
pts, models, solve = g['pts'], g['models'], g['solve']

kappa = 2 / (3 * pi)
E_step = 197.3269804 / 1.016e-19 / 1000          # GeV
E1, Et = 40.0, 2 * 172.57
base, s_base = 136.774, 0.181                     # P-97
print(f'E_step = {E_step:.3e} GeV')
for name in ('(b1) c/ln(E/0.3128) + Z', '(b2) c/ln(E/0.13957) + Z'):
    basis = models[name]
    A = [[sum(f(E) * h(E) / d**2 for _, E, e, d, _, _ in pts) for h in basis] for f in basis]
    b = [sum(f(E) * e / d**2 for _, E, e, d, _, _ in pts) for f in basis]
    c = solve(A, b)
    cov = [solve(A, [1.0 if i == j else 0.0 for i in range(3)]) for j in range(3)]
    sc = sqrt(cov[0][0])
    E0 = 0.3128 if 'b1' in name else 0.13957
    L = lambda E2, Ea: log(log(E2 / E0) / log(Ea / E0))
    per_c = kappa * (11 / 3 * L(E_step, E1) + 4 / 3 * L(E_step, Et))
    tail, s_tail = c[0] * per_c, sc * per_c
    inv = base + tail; sig = sqrt(s_base**2 + s_tail**2); dev = inv - 137.036
    print(f'{name}: c = {c[0]:.4f} +- {sc:.4f}; tail = {tail:+.3f} +- {s_tail:.3f}; 1/alpha far = {base} {tail:+.3f} = {inv:.3f} +- {sig:.3f}'
          f'  vs 137.036 (dev {dev:+.3f}, {dev/sig:+.1f} sigma) -> {"PASS" if abs(dev) <= 2*sig else "FAIL"}')
