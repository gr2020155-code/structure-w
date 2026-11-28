"""Where the complement of a part stops: variants R0-R4 x M1-M2 x F1-F3 (R-154).
Rules: complement_boundary_rules_2026-11-28.md (committed before the run, 5d55ae3).
usage: python -P -E -X utf8 complement_boundary.py"""
from math import log, pi, sqrt, inf

HC = 0.1973269804            # GeV fm
MPART, MN = 0.3128, 0.93827208943
R0 = 2 * HC / MPART          # part form size, fm
C = pi / 12
FD, FE = 0.161, 0.017
XY = 19 / (2 * MN * 100)     # <xy> giving Q^2 = 19 GeV^2 at E = 100 GeV (estimate, bridge)
NUC = [(3, 5, 4), (3, 4, 5)]  # (parts, axes in parts, missing axes): p = uud, n = udd (R-153)


def g_of(r):
    return inf if r >= R0 else C / log(R0 / r)


def phi_of(g, mp):
    if g == inf:
        return 1.0
    return g / (1 + g) if mp == 'M1' else min(g, 1.0)


def f_of(phi, fb):
    num = den = 0.0
    for npart, axes, miss in NUC:
        t = phi * miss
        anti, total = {'F1': (t / 2, axes + t), 'F2': (t / 2, npart + t), 'F3': (t, npart + 2 * t)}[fb]
        num += anti / 2; den += total / 2
    return num / den


def q_r(E, k):               # resolution distance at beam energy E: k*hbar/Q
    return k * HC / sqrt(2 * MN * E * XY)


VAR = [('R0 r0 (part size, full)', R0), ('R1 D_N = 2hbar/(m_N c)', 2 * HC / MN), ('R2 2D_N', 4 * HC / MN),
       ('R3a scattering 0.80 fm', 0.80), ('R3b scattering 0.90 fm', 0.90),
       ('R4a hbar/Q at 100 GeV', q_r(100, 1)), ('R4b 2hbar/Q at 100 GeV', q_r(100, 2))]
print(f'r0 = {R0:.3f} fm, c = pi/12 = {C:.4f}; data f = {FD} +- {FE}\n')
print(f'{"variant":26s} {"r, fm":>6s} {"g":>6s} | ' + ' | '.join(f'{m}/{fb}: phi, f, sigma' for m in ('M1', 'M2') for fb in ('F1', 'F2', 'F3')))
passes = []
for name, r in VAR:
    g = g_of(r); cells = []
    for mp in ('M1', 'M2'):
        for fb in ('F1', 'F2', 'F3'):
            ph = phi_of(g, mp); f = f_of(ph, fb); s = (f - FD) / FE
            cells.append(f'{ph:.3f} {f:.3f} {s:+5.1f}{"*" if abs(s) <= 2 else " "}')
            if abs(s) <= 2:
                passes.append(f'{name} {mp}/{fb}')
    print(f'{name:26s} {r:6.3f} {g:6.3f} | ' + ' | '.join(cells))
print('\n* = within 2 sigma (PASS). Passing:', '; '.join(passes) if passes else 'none')

print('\ncomparison (not counted): energy dependence of f, M1/F1')
for k, lab in ((1, 'R4a hbar/Q'), (2, 'R4b 2hbar/Q')):
    print(f'  {lab}: ' + ', '.join(f'E = {E} GeV: f = {f_of(phi_of(g_of(q_r(E, k)), "M1"), "F1"):.3f}' for E in (30, 100, 200)))
print('  R0-R3: f does not depend on E; data: sigma/E constant for nu and nubar at 30-200 GeV (D-93)')
