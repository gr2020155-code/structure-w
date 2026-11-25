"""G for the slope (rules: slope_G_rules_2026-11-25.md). Needed: k*H/D^2 = 2/3 (loss alpha/(3pi) per e-fold)."""
from math import pi
need = 2 / 3
# H/D^2 = (forms per tact = 1) * (lifetime) * (share of life in the "-" circuit) * (mean held "-" points = pi*D/2) / D^2
branches = {
    'L1 one circuit (pi*D), half of forms "-"': 1 * pi * 0.5 * (pi / 2),
    'L2 two circuits (2*pi*D), one of them "-"': 1 * 2 * pi * 0.5 * (pi / 2),
}
for lname, HD2 in branches.items():
    for kname, k in (('k = 1', 1.0), ('k = 1/3 <cos^2>', 1 / 3), ('k = 2/3 <sin^2>', 2 / 3)):
        val = k * HD2
        print(f'{lname:45s} {kname:16s} k*H/D^2 = {val:.4f}  (need {need:.4f})  ratio {val/need:.3f}  {"HIT" if abs(val/need-1) < 0.01 else "miss"}')
