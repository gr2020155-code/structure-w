"""Relic under A-105 (squared rate = sum of shares) with the light share (R-133). Rules: relic_light_rules_2026-11-27.md.
g^2 = g0^2 [Or u^-4 + Om u^-3 + OB], u = k/k0 = 1/(1+z).
usage: python -P -E -X utf8 relic_light.py <desi_dr2_bao_data.txt> <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, log, pi
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

bao, snd = sys.argv[1:3]
SIG, C, G, MPC, T0 = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
u_gamma = 4 * SIG / C * T0**4                                  # J/m^3
def Omega_gamma(H0):                                           # H0 in km/s/Mpc
    H = H0 * 1e3 / MPC; return u_gamma / (3 * H * H * C * C / (8 * pi * G))
NU = 1 + 0.2271 * 3.046
theta_data = 0.010411

# --- H0 from Pantheon+ (Cepheid calibrators + Hubble flow), shape from BAO G2 (Om set below, iterate once)
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda nm: np.array([float(r[hdr.index(nm)]) for r in data])
z, m, cal, hf, mu_c = col('zHD'), col('m_b_corr'), col('IS_CALIBRATOR') == 1, col('USED_IN_SH0ES_HF') == 1, col('CEPH_DIST')
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
use = cal | hf; Ci = np.linalg.inv(C_all[np.ix_(use, use)])
def measure_H0(Om):
    E = lambda t: sqrt(Om * (1 + t)**3 + 1 - Om)
    Dc = np.array([quad(lambda t: 1 / E(t), 0, zz)[0] for zz in z[use]])
    y = np.where(cal[use], m[use] - mu_c[use], m[use] - 5 * np.log10((1 + z[use]) * Dc * C / 1e3) - 25)
    X = np.stack([np.ones(use.sum()), np.where(cal[use], 0.0, -1.0)], 1)   # y = M - 5 lg H0 * [HF]
    F = X.T @ Ci @ X; b = np.linalg.solve(F, X.T @ Ci @ y); cov = np.linalg.inv(F)
    lgH = b[1] / 5; return 10**lgH, 10**lgH * log(10) / 5 * sqrt(cov[1, 1]), b[0]

# --- BAO with the light share
L = open(bao).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cb = np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]); Cbi = np.linalg.inv(Cb)
def gfun(g0, Om, Or):
    return lambda zz: g0 * np.sqrt(Or * (1 + zz)**4 + Om * (1 + zz)**3 + 1 - Om - Or)
def chi2_bao(g0, Om, Or):
    g = gfun(g0, Om, Or); out = []
    for zz, k in zip(zb, kb):
        dm = quad(lambda t: 1 / g(t), 0, zz)[0]; dh = 1 / g(zz)
        out.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    r = vb - np.array(out); return float(r @ Cbi @ r)
def fit_bao(Or):
    r = minimize(lambda q: chi2_bao(q[0], q[1], Or), [0.0339, 0.297], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-6))
    return r.x, r.fun

def lnint(g, a, b, w=lambda zz: 1.0):
    return quad(lambda x: w(np.exp(x) - 1) * np.exp(x) / g(np.exp(x) - 1), log(1 + a), log(1 + b), limit=400)[0]
def theta(g, zs1, cs=None, csf=None):
    D = lnint(g, 0, zs1 - 1)
    rs = cs * lnint(g, zs1 - 1, 1e9) if csf is None else lnint(g, zs1 - 1, 1e9, csf)
    return rs / D

(q0, _) = fit_bao(0.0)
H0, eH0, Mabs = measure_H0(q0[1])
H0, eH0, Mabs = measure_H0(fit_bao(Omega_gamma(H0) * NU)[0][1])
h = H0 / 100
print(f'H0 from Pantheon+ (Cepheid calibrators + Hubble flow): {H0:.2f} +- {eH0:.2f} km/s/Mpc, M = {Mabs:.3f}')
print(f'light: u_gamma = {u_gamma:.4e} J/m^3, Omega_gamma h^2 = {Omega_gamma(H0) * h * h:.4e}')
for tag, fac in (('nu0 photons only', 1.0), ('nu1 + neutrinos (official factor, bridge)', NU)):
    Or = Omega_gamma(H0) * fac
    (g0, Om), ch = fit_bao(Or); g = gfun(g0, Om, Or)
    th = [theta(g, zs1, cs) for zs1 in (1000, 1100, 1400) for cs in (1 / sqrt(6), 1 / sqrt(3))]
    lo, hi = min(th), max(th)
    print(f'\n{tag}: Or = {Or:.3e}, 1+z_eq = {Om / Or:.0f}; BAO refit Om = {Om:.4f}, chi2 = {ch:.1f} / 11 (prob {chi2dist.sf(ch, 11):.2f})')
    print(f'   relic theta = {lo:.5f} .. {hi:.5f} (data {theta_data}) -> {"PASS" if lo <= theta_data <= hi else "FAIL"}')
    Ob = 0.0224 / h**2
    csf = lambda zz: 1 / sqrt(3 * (1 + 0.75 * Ob / Omega_gamma(H0) / (1 + zz)))
    tb = [theta(g, zs1, csf=csf) for zs1 in (1000, 1060, 1090, 1100, 1400)]
    print(f'   for comparison, baryon loading (Ob h^2 = 0.0224, bridge): theta at 1+z* = 1000, 1060, 1090, 1100, 1400: '
          + ', '.join(f'{t:.5f}' for t in tb))
    lo_, hi_ = 1e-6, 1e-2                                     # Or giving exact theta at middle cs, 1+z* = 1100
    csm = sqrt(1 / sqrt(6) / sqrt(3))
    for _ in range(60):
        mid = sqrt(lo_ * hi_); gm = gfun(g0, Om, mid)
        if theta(gm, 1100, csm) > theta_data: lo_ = mid
        else: hi_ = mid
    print(f'   Or for exact theta at middle cs (c/{1/csm:.3f}) and 1+z* = 1100: {mid:.3e} (this branch {Or:.3e}, ratio {mid / Or:.2f})')
