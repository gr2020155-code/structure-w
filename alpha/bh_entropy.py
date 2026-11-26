"""Black-hole entropy with step a = 2*pi*l_P. Rules: bh_entropy_rules_2026-11-26.md (pushed before running).
Units: l_P = 1, GM/c^2 = 1 (ratios do not depend on M)."""
from math import pi, e, log
a = 2 * pi
A_W, A_GR = 4 * pi * e**2, 16 * pi
S_BH = A_GR / 4
print(f'area: W throat {A_W:.2f}, GR horizon {A_GR:.2f}, ratio {A_W / A_GR:.4f} (= e^2/4)')
print(f'Bekenstein-Hawking per square step a^2: a^2/4 = {a * a / 4:.4f} (pi^2 = {pi**2:.4f})')
print()
print(f'{"variant":10} {"S per l_P^2":>12} {"vs 1/4":>8} {"S_W/S_BH (same M)":>18}')
for pn, sig in (('P1', a * a), ('P2', pi * a * a / 4)):
    for sn, s in (('S1', log(2)), ('S2', log(2) / 2)):
        dens = s / sig
        print(f'{pn + "*" + sn:10} {dens:12.5f} {dens / 0.25:8.4f} {A_W * dens / S_BH:18.4f}')
print()
print(f's per point needed for 1/4: P1 {0.25 * a * a:.3f} nats ({e ** (0.25 * a * a):.3g} states), '
      f'P2 {0.25 * pi * a * a / 4:.3f} nats ({e ** (0.25 * pi * a * a / 4):.3g} states)')
