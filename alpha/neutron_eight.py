"""Neutron by A-52: where the eight goes. Rules: neutron_eight_rules_2026-11-26.md (written before running).

Step 1: which loop can the eight merge with (axes add with sign, at most 3 axes, V-304).
Step 2: sign of the neutron mean-square charge radius for readings C1, C2, C3 (point charges at loop centres).
"""
from fractions import Fraction as F
from math import sqrt

r2_n_data, r2_n_err = -0.1155, 0.0017      # fm^2, PDG
r_p = 0.8408; r2_p = r_p ** 2               # fm, CODATA 2022

# step 1: merging the eight (3 axes, sign -) with each proton loop (signed axes: u = +2, d = -1)
eight = -3
print('step 1: eight (-3 axes) merged with each proton loop')
for name, ax in (('end loop A-P1 (u)', +2), ('middle loop P1-P2 (d)', -1), ('end loop P2-C (u)', +2)):
    total = ax + eight
    allowed = abs(total) <= 3
    print(f'   {name:22}: {ax:+d} {eight:+d} = {total:+d} axes -> charge {F(total, 3)}  {"allowed" if allowed else "IMPOSSIBLE (more than 3 axes, V-304)"}')
print('   -> the eight can only take an end loop: neutron = (d, d, u) along the line')
print()

# step 2: charge radius, loop centres at -L, 0, +L
L2 = r2_p / F(4, 3).__float__()             # from the proton: <r_p^2> = (4/3) L^2
L = sqrt(L2)
print(f'proton (u,d,u): <r^2> = (4/3) L^2 -> L = {L:.3f} fm')
readings = {
    'C1 full merge (d,d,u)':         (F(1, 3), 'no free parameter'),
    'C3 mirror (d,u,d), control':    (F(-2, 3), 'no free parameter, not from A-52'),
}
print(f'{"reading":34} {"<r_n^2>/L^2":>11} {"<r_n^2> fm^2":>13}  sign   ratio n/p')
for k, (c, note) in readings.items():
    r2 = float(c) * L2
    print(f'{k:34} {str(c):>11} {r2:13.4f}  {"-" if r2 < 0 else "+":4}  {r2 / r2_p:+.3f}   ({note})')
# C2: (1/3 - s^2) L^2, s free
s_sign = sqrt(1 / 3)
s_fit = sqrt(1 / 3 - r2_n_data / L2)
print(f'{"C2 merged by axes, eight apart":34} {"1/3 - s^2":>11}  negative only if s > {s_sign:.3f} L; data need s = {s_fit:.3f} L = {s_fit * L:.3f} fm (free parameter)')
print()
print(f'data: <r_n^2> = {r2_n_data} +- {r2_n_err} fm^2, ratio n/p = {r2_n_data / r2_p:+.3f}')
print('magnetic moment (R-78, S3*R2): -2.00 mu_N for (d,d,u) and (d,u,d) alike; data -1.913')
