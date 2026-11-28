"""Neutrino: growth with energy (R-151) and nubar/nu ratio (R-152).
Rules: nu_xsec_rules_2026-11-28.md (committed before the run, 91ae058).
usage: python -P -E -X utf8 nu_xsec.py"""
from math import pi, sqrt

MMU, ME, MW, MPART = 0.1056583755, 0.00051099895069, 80.3625, 0.3128   # GeV; part form mass V-343 p.2


def pstar(M, m1, m2):
    a = (M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2)
    return sqrt(max(a, 0.0)) / (2 * M)


def I2(M, m1, m2):
    """integral d3p1 d3p2 delta(E) delta3(p) for two forms, GeV^2 (same as ppbar.py, n == 2)."""
    if M <= m1 + m2:
        return 0.0
    p = pstar(M, m1, m2)
    return 4 * pi * p * sqrt(p**2 + m1**2) * sqrt(p**2 + m2**2) / M


def vol(D):
    return 4 / 3 * pi * D**3


# ---------- R-151: shape ----------
def sqrt_s(E):            # nu (massless) on a part at rest
    return sqrt(MPART**2 + 2 * E * MPART)


def weight(E, norm, vbr):
    M = sqrt_s(E)
    D = {'V1': 2 / MMU, 'V2': 2 / MW, 'V3': 2 / E}[vbr]
    w = vol(D) / (2 * pi)**3 * I2(M, MMU, MPART)
    if norm == 'N2':
        w /= vol(D) / (2 * pi)**3 * I2(M, 0.0, MPART)
    return w


def decay_ratio(vbr):     # W -> e nu / W -> mu nu
    if vbr == 'V3':
        return None
    De, Dm = (2 / ME, 2 / MMU) if vbr == 'V1' else (2 / MW, 2 / MW)
    return vol(De) * I2(MW, ME, 0) / (vol(Dm) * I2(MW, MMU, 0))


DATA_W, ERR_W = 10.71 / 10.63, 10.71 / 10.63 * sqrt((0.16 / 10.71)**2 + (0.15 / 10.63)**2)
print('R-151: sigma/E(200 GeV) / sigma/E(30 GeV); K1 needs 1 +- 0.04. K2: W->e nu / W->mu nu, data '
      f'{DATA_W:.3f} +- {ERR_W:.3f}')
for norm in ('N1', 'N2'):
    for vbr in ('V1', 'V2', 'V3'):
        r1 = (weight(200, norm, vbr) / 200) / (weight(30, norm, vbr) / 30)
        k1 = abs(r1 - 1) <= 0.04
        r2 = decay_ratio(vbr)
        k2 = r2 is not None and abs(r2 - DATA_W) <= 2 * ERR_W
        r2s = 'n/a (no incoming nu in a decay)' if r2 is None else f'{r2:.4g}'
        print(f'  {norm}+{vbr}: K1 ratio {r1:.4f} {"ok" if k1 else "no"}; K2 {r2s} {"ok" if k2 else "no"}'
              f' -> {"PASS" if k1 and k2 else "FAIL"}')
print('  shape N1+V2, sigma/E relative to 200 GeV: ' + ', '.join(
    f'{E} GeV {(weight(E, "N1", "V2") / E) / (weight(200, "N1", "V2") / 200):.4f}' for E in (1, 3, 10, 30, 100, 200)))
# comparison only: nubar_e + p -> e+ + n at MeV, W weight vs E_e p_e
MPR, MN = 0.93827208943, 0.93956542052
print('  comparison (not counted), nubar p -> e+ n, N1+V2 weight / (E_e p_e), normalised at 4 MeV:')
ref = None
for Emev in (2.0, 3.0, 4.0, 6.0, 8.0):
    E = Emev / 1000; M = sqrt(MPR**2 + 2 * E * MPR)
    w = I2(M, ME, MN)
    Ee = E - (MN - MPR); pe = sqrt(max(Ee**2 - ME**2, 0))
    q = w / (Ee * pe)
    ref = q if Emev == 4.0 else ref
    print(f'    {Emev} MeV: ratio {q:.4g}')

# ---------- R-152: nubar / nu ----------
N = 200_000
cs = [-1 + 2 * (i + 0.5) / N for i in range(N)]          # equal weight of directions (A-31)
BR = {'S1 halves, half angle (A-08)': lambda c: ((1 + c) / 2)**2,
      'S2 halves, full angle': lambda c: (c * c)**2,
      'S3 whole, half angle': lambda c: (1 + c) / 2,
      'S4 whole, full angle (A-107)': lambda c: c * c}
OK_POS = {'S1 halves, half angle (A-08)': True, 'S4 whole, full angle (A-107)': True}
RD, ERD = 0.334 / 0.677, 0.334 / 0.677 * sqrt((0.008 / 0.334)**2 + (0.014 / 0.677)**2)
print(f'\nR-152: nubar/nu, nu (J = 0) share 1; data {RD:.3f} +- {ERD:.3f}')
for name, f in BR.items():
    R = sum(f(c) for c in cs) / N
    pull = (R - RD) / ERD
    if name in OK_POS:
        verdict = 'PASS' if abs(pull) <= 2 else 'FAIL'
    else:
        verdict = 'excluded by positions before the run, shown only'
    back = f(-1.0)
    print(f'  {name:30s}: R = {R:.4f} ({pull:+.1f} sigma), share backwards (y = 1) {back:.2f} -> {verdict}')
fneed = (3 * RD - 1) / (2 + 2 * RD)
fe = 3 * ERD * (2 + 2 * RD) ** -1 + (3 * RD - 1) * 2 * ERD / (2 + 2 * RD)**2
print(f'  anti-part share needed by data (S1 or S4): f = {fneed:.3f} +- {fe:.3f}  (requirement, not a test)')
