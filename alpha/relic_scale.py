"""Angular acoustic scale of the relic in W (R-122). Rules: relic_scale_rules_2026-11-26.md.
1+z = (T0/T)^p; D2: lengths in today's steps, factor (1+z); theta = r_s / D. Units: T0 = 1, c = 1."""
from math import sqrt, log
theta_data = 0.010411

def integral(p, a, b):                          # int_a^b u^(-p) du
    if abs(p - 1) < 1e-12:
        return float('inf') if a == 0 else log(b / a)
    if a == 0 and p > 1: return float('inf')
    return (b**(1 - p) - a**(1 - p)) / (1 - p)

def theta_D2(p, zr1, cs):
    x = zr1**(-1 / p)
    return cs * integral(p, 0.0, x) / integral(p, x, 1.0)
def theta_D1(p, zr1, cs):
    x = zr1**(-1 / p)
    return cs * x / (1 - x)

cs_set = (1 / sqrt(6), 1 / sqrt(3)); z_set = (1000, 1100, 1400)
print(f'data: theta* = {theta_data}\n')
for label, f, counted in (('D2 (today steps, V-324)', theta_D2, True), ('D1 (steps, excluded by V-324)', theta_D1, False)):
    for p, tag in ((1, 'positions'), (2, 'positions'), (3, 'positions'), (4, 'positions'), (1.41, 'SN fit, not counted')):
        vals = [f(p, z, cs) for cs in cs_set for z in z_set]
        lo, hi = min(vals), max(vals)
        hit = lo <= theta_data <= hi
        verdict = ('PASS' if hit else 'FAIL') if (counted and tag == 'positions') else ('in range' if hit else 'out of range')
        rng = 'infinite (sound horizon diverges)' if lo == float('inf') else f'{lo:.4g} .. {hi:.4g}'
        print(f'{label:31s} p = {p:<4} ({tag:19s}): theta = {rng:34s} -> {verdict}')
    print()

def solve(f, zr1=1100, cs=1 / sqrt(3), lo=0.05, hi=3.0):
    g = lambda p: f(p, zr1, cs) - theta_data
    a, b = lo, hi
    if f is theta_D2: b = 0.999999
    for _ in range(200):
        m = (a + b) / 2
        if (g(a) > 0) == (g(m) > 0): a = m
        else: b = m
    return (a + b) / 2
print(f'p needed for theta* (fit, not derivation; cs = c/sqrt3, 1+z_r = 1100): D2 p = {solve(theta_D2):.3f}; D1 p = {solve(theta_D1):.3f}')
print(f'   with cs = c/sqrt6: D2 p = {solve(theta_D2, cs=1/sqrt(6)):.3f}; D1 p = {solve(theta_D1, cs=1/sqrt(6)):.3f}')
print('SN (V-325, D2): p = 1.41 +- 0.05')
