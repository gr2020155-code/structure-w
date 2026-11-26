"""Magnetic anomaly of electron and muon: a = alpha_eff/(2 pi), alpha_eff at the particle's own scale.
Rules: anomaly_rules_2026-11-26.md (R-116).
Run: python -P -E -X utf8 anomaly.py <path to rpp2022-hadronicrpp_page1001.dat>   (pure python; reuses hadron_dispersion.py)"""
import sys, io, os, runpy, contextlib
from math import pi, sqrt, log, exp

here = os.path.dirname(os.path.abspath(__file__))
argv = sys.argv; sys.argv = [os.path.join(here, 'hadron_dispersion.py'), argv[1]]
with contextlib.redirect_stdout(io.StringIO()):
    H = runpy.run_path(sys.argv[0], run_name='hd')
sys.argv = argv
excl, binned = H['excl'], H['binned']

inv0 = 137.035999177; alpha0 = 1 / inv0
a_e, a_mu = 0.00115965218062, 0.00116592072
rho_data = a_mu / a_e - 1
m_e, m_mu, m_tau = 0.51099895, 105.6583755, 1776.93          # MeV
kappa = 2 / (3 * pi); E = -2.1037
print(f'data: a_e = {a_e}, a_mu = {a_mu}, a_mu/a_e - 1 = {rho_data:.7f}')

# 1. electron
need = 2 * pi * a_e
print(f'\n1. electron: 2 pi a_e = 1/{1/need:.4f}; alpha(0) = 1/{inv0:.4f}; needed relative change {need/alpha0 - 1:+.5f}'
      f'  (kappa*alpha = {kappa*alpha0:.5f}, noticed before the rules, not an argument)')

# 2a. distance picture: 1/alpha(r) = 1/alpha0 - kappa * sum max(0, ln(D_f/r) + E), D_f = 2 hbar/(m c); r in units hbar/(MeV c)
def inv_r(r):
    tot = 0.0
    for m in (m_e, m_mu, m_tau):
        x = log(2 / m / r) + E
        if x > 0: tot += x
    return inv0 - kappa * tot
# 2b. energy picture (V-337): 1/alpha(Q) = 1/alpha0 - (1/(3pi)) int R(s) Q^2/(s(s+Q^2)) ds
narrow = {'J/psi': (3096.900, 5.53), 'psi(2S)': (3686.097, 2.33), 'Y(1S)': (9460.40, 1.340), 'Y(2S)': (10023.4, 0.612), 'Y(3S)': (10355.1, 0.443)}
def lep_R(x, m):
    x0 = 2 * m / 1000
    if x <= x0: return 0.0
    b = sqrt(1 - (x0 / x)**2); return b * (3 - b * b) / 2
def delta(Q):                                                    # Q in GeV
    k = lambda x: Q * Q / (x * x * (x * x + Q * Q)) * 2 * x
    tot = 0.0
    for m in (m_e, m_mu, m_tau):
        x0 = 2 * m / 1000; n = 200000; L = log(1e4 / x0)
        for i in range(n):
            x = x0 * exp((i + 0.5) / n * L); tot += lep_R(x, m) * k(x) * x * L / n
    pts = [(2 * 0.13957, 0.0)] + [(r[0], r[1]) for r in excl] + [(p[0], p[1]) for p in binned(40.0)]
    pts.append((40.0, pts[-1][1]))
    for (xa, Ra), (xb, Rb) in zip(pts, pts[1:]):
        tot += (xb - xa) / 2 * (Ra * k(xa) + Rb * k(xb))
    tot += 11 / 3 * log(1 + Q * Q / 1600.0)
    for M, G in narrow.values():
        Mg = M / 1000; tot += 9 * pi * G * 1e-3 / (alpha0**2 * M) * Q * Q / (Mg * Mg + Q * Q)
    return tot / (3 * pi)
def inv_Q(Q): return inv0 - delta(Q)

print('\n2. a_mu/a_e - 1 = alpha_eff(mu)/alpha_eff(e) - 1 by branches of the scale')
res = []
for name, f in (('X1 distance r = R = hbar/(mc)', lambda m: inv_r(1 / m)),
                ('X2 distance r = D = 2hbar/(mc)', lambda m: inv_r(2 / m)),
                ('P1 energy Q = mc^2', lambda m: inv_Q(m / 1000)),
                ('P2 energy Q = mc^2/2', lambda m: inv_Q(m / 2000))):
    ie, im = f(m_e), f(m_mu); rho = ie / im - 1; dev = rho / rho_data - 1
    res.append((name, rho, abs(dev) <= 0.25))
    print(f'   {name:32s}: 1/alpha_eff(e) = {ie:.4f}, 1/alpha_eff(mu) = {im:.4f}, rho = {rho:.5f} ({dev*100:+.1f} % vs data) -> {"AGREES" if abs(dev) <= 0.25 else "NO"}')
n = sum(ok for _, _, ok in res)
print(f'\nagreeing branches: {n} of {len(res)}' + (' -> data do not choose the branch' if n > 1 else ''))
print(f'electron in the same branches: alpha_eff(e)/alpha0 - 1 = ' + ', '.join(f'{nm.split()[0]} {inv0/ f(m_e) - 1:+.5f}' for nm, f in
      (('X1', lambda m: inv_r(1 / m)), ('X2', lambda m: inv_r(2 / m)), ('P1', lambda m: inv_Q(m / 1000)), ('P2', lambda m: inv_Q(m / 2000))))
      + f'   (needed {need/alpha0 - 1:+.5f})')
