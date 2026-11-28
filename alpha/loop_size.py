"""Size of the global loop from B vs repeated patterns on the relic (R-157).
Rules: loop_size_rules_2026-11-28.md (committed before the run, a0d16e0).
usage: python -P -E -X utf8 loop_size.py"""
from math import sqrt, pi, log, exp
from scipy.integrate import quad

SIG, C, G, MPC, TCMB = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
Or = (4 * SIG / C * TCMB**4) / (3 * (0.7354 * 1e5 / MPC) ** 2 * C * C / (8 * pi * G))
Om = 0.2973; OB = 1 - Or - Om; ZS = 1090.0
SLAB, CUBE = 0.56, 0.97                       # D-94: R_i > x * chi_rec
LY_PER_PC = 3.261563777


def E(z):
    u1 = 1 + z
    return sqrt(Or * u1**4 + Om * u1**3 + OB)


def dist(z):                                  # int_0^z dz/E, via x = ln(1+z)
    return quad(lambda x: exp(x) / E(exp(x) - 1), 0, log(1 + z), limit=500)[0]


Ds = dist(ZS); Dh = dist(1e12)
L1 = 1 / sqrt(OB)                             # in units c/H0
print(f'Om = {Om}, OB = {OB:.4f}, Or = {Or:.3e}; D* = {Ds:.4f} c/H0, D_h = {Dh:.4f} c/H0, 1/sqrt(OB) = {L1:.4f} c/H0\n')
for name, L in (('V1: loop of our half = 1/sqrt(B)', L1), ('V2: whole eight = 1/sqrt(B), our half = half', L1 / 2)):
    r = L / Ds; ok = r >= 2 * SLAB
    print(f'{name}: L/D* = {r:.3f}; needed >= {2 * SLAB:.2f} (slab), {2 * CUBE:.2f} (cubic torus)  =>  {"PASS" if ok else "FAIL"}')
    print(f'   light from the start goes around {Dh / L:.1f} times')
    for h in (0.7354, 0.687):
        cH = C / 1e3 / (100 * h)              # Mpc
        print(f'   h = {h}: L = {L * cH / 1e3:.2f} Gpc = {L * cH * LY_PER_PC / 1e3:.1f} Gly; D* = {Ds * cH / 1e3:.2f} Gpc = {Ds * cH * LY_PER_PC / 1e3:.1f} Gly')
    print()
