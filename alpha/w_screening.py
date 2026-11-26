"""Does the W boson enter the running of alpha? All branches - rules: w_screening_rules_2026-11-26.md (R-114)."""
from math import pi, sqrt, log
hbar = 1.054571817e-34; c = 299792458.0; G = 6.67430e-11; MeV = 1.602176634e-13 / c**2
M_W = 80369.2                                    # MeV (PDG 2024)
l_P = sqrt(hbar * G / c**3); a = 2 * pi * l_P
kappa = 2 / (3 * pi); E = -2.1037
D_W = 2 * hbar / (M_W * MeV * c); L = log(D_W / a) + E
base, err, data = 137.354, 0.185, 137.036
print(f'ln(D_W/a) + E = {log(D_W / a):.3f} {E:+.4f} = {L:.3f};  kappa * L = {kappa * L:.3f} per unit weight')
for name, w, note in (('B0 parts already counted', 0, ''), ('B1 own species, w = 1', 1, ''),
                      ('B2 own species, spin-0 weight 1/4', 0.25, ' (bridge, model)'),
                      ('B3 anti-screening -21/4', -21 / 4, ' (bridge, model)')):
    inv = base + kappa * w * L
    dev = (inv - data) / err
    print(f'{name:36s}: contribution {kappa * w * L:+8.3f}  1/alpha = {inv:8.3f}  ({dev:+6.1f} sigma) -> {"PASS" if abs(dev) <= 2 else "FAIL"}{note}')
