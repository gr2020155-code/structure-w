"""Steep part at the step (rules: step_part_rules_2026-11-25.md). Delta(1/alpha) = 2 * integral f dr from r = 1 step,
f = H_c * p(r) / (4 pi r^2), p(r) = (2/pi)/sqrt(D^2 - r^2)."""
from math import pi, sqrt, log, exp
D = 7.723e-13 / 1.0155e-34
# integral of p(r)/(4 pi r^2) from 1 to D, done in log variable (r = e^u), accurate near both ends
def I(D, n=200000):
    u0, u1 = 0.0, log(D) - 1e-12
    h = (u1 - u0) / n; s = 0.0
    for k in range(n + 1):
        u = u0 + k * h; r = exp(u)
        val = (2 / pi) / sqrt(max(D * D - r * r, 1e-300)) / (4 * pi * r * r) * r
        s += (1 if k in (0, n) else (4 if k % 2 else 2)) * val
    return s * h / 3
base = I(D)
print(f'D = {D:.3e} steps; integral = {base:.4e}; approx 1/(2 pi^2 D) = {1 / (2 * pi**2 * D):.4e}')
for name, Hc in (('T1 all configurations (2/3) D^2', 2 / 3 * D * D), ('T2 one form, whole circle pi*D', pi * D), ('T3 one form, current "-" act, 2 points', 2.0)):
    d = 2 * Hc * base
    print(f'   {name:42s} Delta(1/alpha) = {d:.4g}  {"HIT" if abs(d - 64.3) < 10 else "miss"} (need ~64, or ~107 with W)')
