"""Decay length of screening from positions: readings L1-L4 (26.11.2026). Rules: decay_length_rules_2026-11-26.md (pushed before run).
lambda = f * R * sqrt(x(1-x)); all screening lengths x f  <=>  k_len = f in alpha_gap.py terms.
Run: python -X utf8 decay_length.py   (needs scipy)
"""
from math import pi, sqrt, log
from scipy.integrate import quad
from scipy.special import exp1
GAMMA = 0.5772156649015329
hbar = 1.054571817e-34; c = 299792458.0; G = 6.67430e-11
MeV = 1.602176634e-13 / c**2
a = 2 * pi * sqrt(hbar * G / c**3)
kappa = 2 / (3 * pi)
inv_q0sq = 1 / 0.12282**2
target = 137.035999206
alpha0 = 1 / target
opal, opal_err = 440e-5, sqrt(58**2 + 43**2 + 30**2) * 1e-5
t1, t2 = 6.07e6, 1.81e6

m_p = 938.27209
base = {'e': (0.51099895, 1, 'l'), 'mu': (105.6583755, 1, 'l'), 'tau': (1776.86, 1, 'l'),
        'u': (m_p / 3, 4 / 3, 'ud'), 'd': (m_p / 3, 1 / 3, 'ud'), 's': (1115.683 - 2 * m_p / 3, 1 / 3, 'q'),
        'c': (2286.46 - 2 * m_p / 3, 4 / 3, 'q'), 'b': (5619.60 - 2 * m_p / 3, 1 / 3, 'q'), 't': (172570.0, 4 / 3, 't')}
w3 = lambda x: 6 * x * (1 - x)
E = quad(lambda x: w3(x) * log(sqrt(x * (1 - x)) / 2), 0, 1, points=[0.5])[0] - GAMMA   # S3


def masses(k_len=1.0, k_ud=1.0, k_q=1.0):
    out = []
    for m, w, g in base.values():
        mm = m / k_len                                   # all screening lengths x k  <=>  all masses / k
        if g == 'ud': mm *= k_ud * k_q
        if g == 'q': mm *= k_q
        out.append((mm, w))
    return out


def inv_far(add=0.0, **k):
    return inv_q0sq + add + sum(kappa * w * (log(2 * hbar / (m * MeV * c) / a) + E) for m, w in masses(**k))


def delta_inv(reading, Q2, **k):
    Q = sqrt(Q2); tot = 0.0
    for m, w in masses(**k):
        if reading == 'M1':
            f = lambda x: exp1(m / (Q * sqrt(x * (1 - x))))
        else:
            f = lambda x: 0.5 * log(1 + Q2 * x * (1 - x) / m**2)
        tot += kappa * w * quad(lambda x: w3(x) * f(x), 0, 1, limit=400, points=[0.5])[0]
    return -tot


def opal_W(reading, **k):
    return alpha0 * (delta_inv(reading, t2, **k) - delta_inv(reading, t1, **k))


s_W = 2 * 0.00013 / 0.12282**3
readings = {'L1 f=1': 1.0, 'L2 f=sqrt3': sqrt(3), 'L3 f=sqrt(3/2)': sqrt(1.5), 'L4 f=sqrt2': sqrt(2)}
print(f'C1: 137.036, pass if |dev| <= {2*s_W:.2f};  C2: OPAL 440 +- {opal_err*1e5:.0f}, pass if within 156 in one reading')
for name, f in readings.items():
    inv = inv_far(k_len=f); dev = inv - target
    o1, o2 = opal_W('M1', k_len=f), opal_W('M2', k_len=f)
    c1 = abs(dev) <= 2 * s_W; c2 = min(abs(o1 - opal), abs(o2 - opal)) <= 2 * opal_err
    print(f'{name}: edge E = {E + log(f):.4f}; 1/alpha_W(0) = {inv:.3f} +- {s_W:.2f} (dev {dev:+.3f}, {dev/s_W:+.1f} sigma) {"PASS" if c1 else "FAIL"};'
          f' OPAL M1 {o1*1e5:.0f}, M2 {o2*1e5:.0f} {"PASS" if c2 else "FAIL"}; both: {"yes" if c1 and c2 else "no"}')
