"""Which forms the charge induces: sets H1-H7 (26.11.2026). Rules: induced_forms_rules_2026-11-26.md (pushed before run).
Run: python -X utf8 induced_forms.py   (needs scipy)
"""
from math import pi, sqrt, log
from scipy.integrate import quad
from scipy.special import exp1

GAMMA = 0.5772156649015329
hbar = 1.054571817e-34; c = 299792458.0; G = 6.67430e-11
MeV = 1.602176634e-13 / c**2
a = 2 * pi * sqrt(hbar * G / c**3)
kappa = 2 / (3 * pi)
inv_q0sq = 1 / 0.12282**2; s_W = 2 * 0.00013 / 0.12282**3
target = 137.035999206; alpha0 = 1 / target
opal, opal_err = 440e-5, sqrt(58**2 + 43**2 + 30**2) * 1e-5
t1, t2 = 6.07e6, 1.81e6
w3 = lambda x: 6 * x * (1 - x)
E = quad(lambda x: w3(x) * log(sqrt(x * (1 - x)) / 2), 0, 1, points=[0.5])[0] - GAMMA

m_p, m_pi, m_K, m_D, m_B = 938.27209, 139.57039, 493.677, 1869.66, 5279.34
leptons = [(0.51099895, 1), (105.6583755, 1), (1776.86, 1)]
top = [(172570.0, 4 / 3)]
quarks = {
    'K1 own part mass (Rb)': dict(u=m_p / 3, d=m_p / 3, s=1115.683 - 2 * m_p / 3, c=2286.46 - 2 * m_p / 3, b=5619.60 - 2 * m_p / 3),
    'K2 lightest charged whole form': dict(u=m_pi, d=m_pi, s=m_K, c=m_D, b=m_B),
    'K3 share of lightest two-part form': dict(u=m_pi / 2, d=m_pi / 2, s=m_K - m_pi / 2, c=m_D - m_pi / 2, b=m_B - m_pi / 2)}
NcQ2 = dict(u=4 / 3, d=1 / 3, s=1 / 3, c=4 / 3, b=1 / 3)


def species(reading):
    return leptons + [(quarks[reading][k], NcQ2[k]) for k in 'udscb'] + top


def inv_far(sp):
    return inv_q0sq + sum(kappa * w * (log(2 * hbar / (m * MeV * c) / a) + E) for m, w in sp)


def delta_inv(sp, reading, Q2):
    Q = sqrt(Q2); tot = 0.0
    for m, w in sp:
        f = (lambda x: exp1(m / (Q * sqrt(x * (1 - x))))) if reading == 'M1' else (lambda x: 0.5 * log(1 + Q2 * x * (1 - x) / m**2))
        tot += kappa * w * quad(lambda x: w3(x) * f(x), 0, 1, limit=400, points=[0.5])[0]
    return -tot


def opal_W(sp, reading):
    return alpha0 * (delta_inv(sp, reading, t2) - delta_inv(sp, reading, t1))


H5 = {'pi': (139.57039, 1, 0), 'K': (493.677, 1, 0), 'D': (1869.66, 1, 0), 'Ds': (1968.35, 1, 0), 'B': (5279.34, 1, 0), 'Bc': (6274.47, 1, 0),
      'p': (938.27209, 1, 1), 'Sigma+': (1189.37, 1, 1), 'Sigma-': (1197.449, 1, 1), 'Xi-': (1321.71, 1, 1), 'Omega-': (1672.45, 1, 0),
      'Lc+': (2286.46, 1, 1), 'Xic+': (2467.71, 1, 1), 'Xicc++': (3621.6, 4, 1), 'Xib-': (5797.0, 1, 1), 'Omegab-': (6045.2, 1, 1)}  # (mass, Q^2, spin-1/2 flag)
sets = {f'H{i+1} parts {k}': species(k) for i, k in enumerate(quarks)}
sets['H4 whole lightest pi,K,D,B'] = leptons + [(H5[k][0], 1) for k in ('pi', 'K', 'D', 'B')] + top
sets['H5 whole all long-lived'] = leptons + [(m, q) for m, q, _ in H5.values()] + top
sets['H6 parts K1 + pion tail'] = species('K1 own part mass (Rb)') + [(H5['pi'][0], 1), (m_p / 3, -1)]
sets['H7 spin-1/2 only'] = leptons + [(m, q) for m, q, s in H5.values() if s] + top
print(f'C1: 137.036 +- {2*s_W:.2f};  C2: OPAL 440 within 156 in one reading')
for name, sp in sets.items():
    inv = inv_far(sp); dev = inv - target
    o1, o2 = opal_W(sp, 'M1'), opal_W(sp, 'M2')
    c1 = abs(dev) <= 2 * s_W; c2 = min(abs(o1 - opal), abs(o2 - opal)) <= 2 * opal_err
    print(f'{name}: 1/alpha_W(0) = {inv:.3f} (dev {dev:+.3f}, {dev/s_W:+.1f} sigma) {"PASS" if c1 else "FAIL"}; OPAL M1 {o1*1e5:.0f}, M2 {o2*1e5:.0f} {"PASS" if c2 else "FAIL"}; both: {"yes" if c1 and c2 else "no"}')
