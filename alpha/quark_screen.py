"""Size of induced forms carrying quark parts: readings K1-K3 (26.11.2026). Rules: quark_screen_rules_2026-11-26.md (pushed before run).
Run: python -X utf8 quark_screen.py   (needs scipy)
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


el = [leptons[0]]
print(f'check of the rules claim: electron alone, OPAL M1 {opal_W(el, "M1")*1e5:.0f}, M2 {opal_W(el, "M2")*1e5:.0f} (1e-5) vs 440 +- 78')
print(f'C1: 137.036 +- {2*s_W:.2f};  C2: OPAL within 156 in one reading')
for name in quarks:
    sp = species(name); inv = inv_far(sp); dev = inv - target
    o1, o2 = opal_W(sp, 'M1'), opal_W(sp, 'M2')
    c1 = abs(dev) <= 2 * s_W; c2 = min(abs(o1 - opal), abs(o2 - opal)) <= 2 * opal_err
    q = quarks[name]
    print(f'{name}: u=d {q["u"]:.1f}, s {q["s"]:.1f}, c {q["c"]:.1f}, b {q["b"]:.1f} MeV -> 1/alpha_W(0) = {inv:.3f} (dev {dev:+.3f}, {dev/s_W:+.1f} sigma) {"PASS" if c1 else "FAIL"};'
          f' OPAL M1 {o1*1e5:.0f}, M2 {o2*1e5:.0f} {"PASS" if c2 else "FAIL"}; both: {"yes" if c1 and c2 else "no"}')
