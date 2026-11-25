"""Neutron: attachment point of the eight. Rules: neutron_attach_rules_2026-11-26.md (written before running)."""
from math import sqrt

r2_n_data = -0.1155            # fm^2, D-60
r2_p = 0.8408 ** 2             # fm^2
L2 = r2_p / (4 / 3)            # proton: <r^2> = (4/3) L^2  (R-79)
print(f'L = {sqrt(L2):.3f} fm;  data <r_n^2> = {r2_n_data} fm^2, ratio n/p = {r2_n_data / r2_p:+.3f}')
print(f'{"T (point on proton)":24} {"G (point of eight)":22} {"r_O^2/L^2":>9} {"<r_n^2>/L^2":>11} {"fm^2":>8} {"sign":>4} {"n/p":>7}')
for T, x in (('T1 end point A', -1.5), ('T2 inner point P1', -0.5)):
    for G, perp in (('G1 centre O\'', 0.0), ('G2 end point A\'', 1.0)):
        rO2 = x * x + perp * perp
        r2 = 2 / 3 + 0 + 2 / 3 - rO2          # in units of L^2
        print(f'{T:24} {G:22} {rO2:9.2f} {r2:11.3f} {r2 * L2:8.3f} {"-" if r2 < 0 else "+":>4} {r2 * L2 / r2_p:+7.3f}')
