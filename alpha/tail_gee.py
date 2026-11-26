"""Excess tail above 40 GeV and Gamma_ee of narrow resonances - rules: tail_gee_rules_2026-11-26.md (R-115).
Run: python -P -E -X utf8 tail_gee.py <path to rpp2022-hadronicrpp_page1001.dat>   (pure python; reuses hadron_dispersion.py)"""
import sys, io, os, runpy, contextlib
from math import pi, sqrt, log

here = os.path.dirname(os.path.abspath(__file__))
argv = sys.argv; sys.argv = [os.path.join(here, 'hadron_dispersion.py'), argv[1]]
with contextlib.redirect_stdout(io.StringIO()):
    H = runpy.run_path(sys.argv[0], run_name='hd')
sys.argv = argv
excl, binned, alpha0 = H['excl'], H['binned'], H['alpha0']

# ---- 1. narrow resonances: PDG 2026 (pdgLive PARTIAL WIDTHS)
narrow = {'J/psi': (3096.900, 5.53), 'psi(2S)': (3686.097, 2.33), 'Y(1S)': (9460.40, 1.340),
          'Y(2S)': (10023.4, 0.612), 'Y(3S)': (10355.1, 0.443)}           # MeV, keV
def N1():
    return sum(3 * G * 1e-3 / (alpha0**2 * M) for M, G in narrow.values())

leptons = [0.51099895, 105.6583755, 1776.93]                              # MeV (PDG)
def lep_R(x, m):                                                          # x = sqrt(s) GeV
    x0 = 2 * m / 1000
    if x <= x0: return 0.0
    b = sqrt(1 - (x0 / x)**2); return b * (3 - b * b) / 2

def delta(Q):
    """(1/(3pi)) * int R(s) Q^2/(s(s+Q^2)) ds, Q in GeV; s = x^2, ds = 2x dx"""
    k = lambda x: Q * Q / (x * x * (x * x + Q * Q)) * 2 * x
    tot = 0.0
    for m in leptons:                                                    # log grid in x from threshold
        x0 = 2 * m / 1000; n = 200000; L = log(1e4 / x0)
        for i in range(n):
            x = x0 * __import__('math').exp((i + 0.5) / n * L); dx = x * L / n
            tot += lep_R(x, m) * k(x) * dx
    pts = [(2 * 0.13957, 0.0)] + [(r[0], r[1]) for r in excl] + [(p[0], p[1]) for p in binned(40.0)]
    pts.append((40.0, pts[-1][1]))
    for (xa, Ra), (xb, Rb) in zip(pts, pts[1:]):
        tot += (xb - xa) / 2 * (Ra * k(xa) + Rb * k(xb))
    tot += 11 / 3 * log(1 + Q * Q / 1600.0)                              # above 40 GeV: int Q^2/(s(s+Q^2)) ds = ln(1+Q^2/s0)
    for M, G in narrow.values():
        Mg = M / 1000; tot += 9 * pi * G * 1e-3 / (alpha0**2 * M) * Q * Q / (Mg * Mg + Q * Q)
    return tot / (3 * pi)

n1 = N1(); n2 = 0.0
print('1. narrow resonances (PDG 2026):')
for name, (M, G) in narrow.items():
    inv_M = 1 / alpha0 - delta(M / 1000); f = (inv_M * alpha0)**2
    v1 = 3 * G * 1e-3 / (alpha0**2 * M); n2 += v1 * f
    print(f'   {name:8s} M {M:9.3f} MeV  Gamma_ee {G:6.3f} keV  N1 {v1:.4f}   1/alpha(M) = {inv_M:.2f}  factor {f:.4f}  N2 {v1*f:.4f}')
print(f'   N1 = {n1:.4f} (R-97: 0.1500);  N2 = {n2:.4f};  N2 - N1 = {n2 - n1:+.4f}')

# ---- 2. tail
kappa = 2 / (3 * pi); c = pi / 12
E_step = 197.3269804 / 1.016e-19 / 1000                                  # GeV, hbar*c/a (as R-100)
E_P = 2 * pi * E_step
E1, mt = 40.0, 172.57
q = 0.12289; base = 136.774 + (1 / q**2 - 66.291); err = 0.130
L = lambda E2, Ea, E0: log(log(E2 / E0) / log(Ea / E0))
def tail(E0s, t_start=2 * mt, E2=E_step):
    w = {'u': 4 / 3, 'd': 1 / 3, 's': 1 / 3, 'c': 4 / 3, 'b': 1 / 3}
    tot = sum(wi * L(E2, E1, E0s[k]) for k, wi in w.items()) + (4 / 3 * L(E2, t_start, E0s['t']) if E2 > t_start else 0.0)
    return kappa * c * tot
one = {k: 0.3128 for k in 'udscbt'}
own = {'u': 0.3128, 'd': 0.3128, 's': 0.4902, 'c': 1.6610, 'b': 4.9941, 't': 172.57}
E58 = 57.8; eps58 = c / log(E58 / 0.3128)
branches = [('H1 one E0, E_step, t from 2m_t', tail(one), True),
            ('H2 own E0 per species', tail(own), False),
            ('H3 t from m_t', tail(one, t_start=mt), True),
            ('H4 upper end m_P', tail(one, E2=E_P), True),
            ('H5 frozen above 58 GeV', tail(one, E2=E58) + kappa * eps58 * (11 / 3 * log(E_step / E58) + 4 / 3 * log(E_step / (2 * mt))), False),
            ('H6 no excess above 58 GeV', tail(one, E2=E58), False)]
print(f'\n2. tail (base without tail {base:.3f}, error {err})')
for name, t, allowed in branches:
    inv = base + t; dev = (inv - 137.036) / err
    print(f'   {name:30s}: tail {t:+.3f}  1/alpha = {inv:.3f}  ({dev:+.1f} sigma) -> {"PASS" if abs(dev) <= 2 else "FAIL"}'
          f'{"" if allowed else "   [excluded by positions]"}')
h1 = branches[0][1]; sysm = max(abs(branches[2][1] - h1), abs(branches[3][1] - h1), abs(n2 - n1))
print(f'\nMAIN: H1 + N1 = {base + h1:.3f} +- {err} (syst +- {sysm:.3f});  with N2: {base + h1 + n2 - n1:.3f}')
