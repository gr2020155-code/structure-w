"""Nucleon magnetic moments from loop fragments: all 24 readings of nucleon_moments_rules_2026-11-26.md.

Units: nuclear magneton mu_N = e*hbar/(2 m_p). A part with charge Q, circulation mass M_part, spin weight w
contributes w * Q * (m_p / M_part) mu_N  (V-14, V-12: mu = Q e hbar / (2 M)).
"""
from fractions import Fraction as F

m_p, m_n = 938.272, 939.565                     # MeV (A-14: masses are data)
mu_p_data, mu_n_data = 2.79284734, -1.91304273  # D-50
TOL = 0.10

u, d = F(2, 3), F(-1, 3)

def weights(S, charges, pair=None):
    """Spin weights of the parts along the particle axis."""
    if S == 'S1':
        return [1] * len(charges)
    if S == 'S2':
        return [1 if q > 0 else -1 for q in charges]
    if S == 'S3':                      # like pair coupled to spin 1, third part attached to total 1/2
        a, b = pair
        return [F(2, 3) if i in (a, b) else F(-1, 3) for i in range(len(charges))]

def moment(S, R, charges, M, pair, extra=None):
    """extra: (charge, weight) of an additional whole part (the eight in N1)."""
    k = len(charges) + (1 if extra else 0)
    unit = m_p / M if R == 'R1' else k * m_p / M
    w = weights(S, charges, pair)
    mu = sum(float(wi * qi) for wi, qi in zip(w, charges))
    if extra:
        mu += extra[0] * extra[1]
    return mu * unit

def ok(x, ref):
    return abs(x / ref - 1) <= TOL

# consistency check: the electron is one whole part, charge -1, spin along the axis -> Dirac moment
print('check: electron as one part, S1/R1 and S1/R2 -> |mu| =',
      abs(moment('S1', 'R1', [F(-1)], 1.0, None) / (m_p / 1.0)), 'Dirac units (must be 1)')
print()

proton = [u, d, u]
print(f'{"S":3} {"R":3} {"neutron":9} {"mu_p":>7} {"err_p":>7} {"mu_n":>7} {"err_n":>7} {"ratio":>7}  pass')
passed = []
for S in ('S1', 'S2', 'S3'):
    for R in ('R1', 'R2'):
        mp = moment(S, R, proton, m_p, pair=(0, 2))
        for N in ('N1-E0', 'N1-E+', 'N1-E-', 'N2'):
            if N == 'N2':
                mn = moment(S, R, [d, d, u], m_n, pair=(0, 1))
            else:
                e = {'N1-E0': 0, 'N1-E+': 1, 'N1-E-': -1}[N]
                mn = moment(S, R, proton, m_n, pair=(0, 2), extra=(-1, e))
            good = ok(mp, mu_p_data) and ok(mn, mu_n_data)
            if good:
                passed.append((S, R, N))
            ratio = mn / mp if mp else float('nan')
            print(f'{S:3} {R:3} {N:9} {mp:7.3f} {100*(mp/mu_p_data-1):+6.1f}% {mn:7.3f} '
                  f'{100*(mn/mu_n_data-1) if mn else -100:+6.1f}% {ratio:7.3f}  {"YES" if good else ""}')
print()
print(f'data: mu_p = {mu_p_data}, mu_n = {mu_n_data}, ratio = {mu_n_data/mu_p_data:.4f}')
print(f'passed ({len(passed)} of 24):', passed if passed else 'none')
