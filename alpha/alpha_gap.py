"""Which link of the alpha chain gives the missing +0.94 (26.11.2026). Rules: alpha_gap_rules_2026-11-26.md (pushed before run).
Diagnosis, not derivation: one number per link is FITTED to 1/alpha far = 137.036, then OPAL (not used in the fit) is the test.
Run: python -X utf8 alpha_gap.py   (needs scipy)
"""
from math import pi, sqrt, log
from scipy.integrate import quad
from scipy.special import exp1
from scipy.optimize import brentq

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


print(f'base (P-91, S3): 1/alpha far = {inv_far():.3f}, missing {target - inv_far():+.3f}')
cases = []
add = target - inv_far(); cases.append(('Z1 step/q0 (+add to 1/alpha)', f'add = {add:+.3f} (1/q0^2 = {inv_q0sq+add:.2f}, q0 = {1/sqrt(inv_q0sq+add):.5f}; or a / {2.718281828**(add/(kappa*8)):.3f})', dict()))
k = brentq(lambda v: inv_far(k_len=v) - target, 0.5, 5); cases.append(('Z2 all screening lengths x k', f'k = {k:.4f} (ln k = {log(k):.4f})', dict(k_len=k)))
k = brentq(lambda v: inv_far(k_ud=v) - target, 1e-4, 2); cases.append(('Z3 u = d masses x k_ud', f'k_ud = {k:.4f} -> m_ud = {m_p/3*k:.1f} MeV', dict(k_ud=k)))
k = brentq(lambda v: inv_far(k_q=v) - target, 1e-4, 2); cases.append(('Z4 all quark parts x k_q', f'k_q = {k:.4f} -> u = d {m_p/3*k:.1f}, s {base["s"][0]*k:.1f}, c {base["c"][0]*k:.1f}, b {base["b"][0]*k:.1f} MeV', dict(k_q=k)))

print(f'\nOPAL target 440 +- {opal_err*1e5:.0f} (1e-5); compatible if |dev| <= 156 in at least one reading')
print(f'before fit: M1 {opal_W("M1")*1e5:.0f}, M2 {opal_W("M2")*1e5:.0f}')
for name, fit, kw in cases:
    o1, o2 = opal_W('M1', **kw), opal_W('M2', **kw)
    d1, d2 = (o1 - opal) / opal_err, (o2 - opal) / opal_err
    ok = abs(o1 - opal) <= 2 * opal_err or abs(o2 - opal) <= 2 * opal_err
    print(f'{name}: {fit}\n    OPAL: M1 {o1*1e5:.0f} ({d1:+.1f} sigma), M2 {o2*1e5:.0f} ({d2:+.1f} sigma) -> {"compatible" if ok else "EXCLUDED as sole cause"}')
