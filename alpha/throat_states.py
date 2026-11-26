"""States of a throat point. Rules: throat_states_rules_2026-11-26.md (pushed before running). Units l_P = 1."""
from math import pi, e, log, sqrt
a = 2 * pi
G = 0.915965594177219                      # Catalan constant
SQ = a * a; DISK = pi * a * a / 4; TRI = sqrt(3) / 2 * a * a
side_pair = log(2) / 2
V = [('K1 side per point, square', log(2), SQ), ('K1 side per point, disk', log(2), DISK),
     ('K1 side per pair, square', side_pair, SQ), ('K1 side per pair, disk', side_pair, DISK),
     ('K2 square matching + side', G / pi + side_pair, SQ),
     ('K3 triangular matching + side', 0.4286 + side_pair, TRI)]
for z in (4, 5, 6, 7):
    V += [(f'K4 free, z={z}, square', log(2 * z), SQ), (f'K4 free, z={z}, disk', log(2 * z), DISK)]
ratio_area = e * e / 4
print(f'{"variant":32} {"s, nats":>8} {"per l_P^2":>10} {"vs 1/4":>7} {"S_W/S_BH":>9}')
rs = []
for name, s, sig in V:
    d = s / sig
    r = d * ratio_area / 0.25
    rs.append(r)
    print(f'{name:32} {s:8.4f} {d:10.5f} {d / 0.25:7.3f} {r:9.3f}')
print()
print(f'range S_W/S_BH (same mass): {min(rs):.3f} - {max(rs):.3f};  per-area vs 1/4: '
      f'{min(s / g for _, s, g in V) / 0.25:.3f} - {max(s / g for _, s, g in V) / 0.25:.3f}')
print(f'needed s for 1/4: square {0.25 * SQ:.2f}, disk {0.25 * DISK:.2f}, triangular {0.25 * TRI:.2f} nats')
