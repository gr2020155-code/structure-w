"""Constant share of the rate from the age: test T0*g0 = 1 on raw BAO + SN (R-155).
Rules: age_share_rules_2026-11-28.md (committed before the run, 58276ac). Kinematics as in R-132..R-134 (V-374, V-324).
usage: python -P -E -X utf8 age_share.py <desi_dr2_bao_data.txt> <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, pi, log, exp
from scipy.integrate import quad
from scipy.optimize import minimize, brentq
from scipy.stats import chi2 as chi2dist

bao, snd = sys.argv[1:3]
h = 0.7354
SIG, C, G, MPC, TCMB = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
Or = (4 * SIG / C * TCMB**4) / (3 * (h * 1e5 / MPC) ** 2 * C * C / (8 * pi * G))

L = open(bao).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cbi = np.linalg.inv(np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]))
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
zs_all = np.array([float(r[hdr.index('zHD')]) for r in data]); ms_all = np.array([float(r[hdr.index('m_b_corr')]) for r in data])
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
sel = zs_all > 0.01; zs, ms = zs_all[sel], ms_all[sel]
Csi = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(zs); w1 = Csi @ one; s11 = one @ w1
REF = 1403.3


def E(zz, Om):
    u1 = 1 + zz
    return sqrt(Or * u1**4 + Om * u1**3 + 1 - Or - Om)


def T0g0(Om):                      # age in units 1/g0: int_0^inf dz / ((1+z) E), with x = ln(1+z)
    return quad(lambda x: 1 / E(exp(x) - 1, Om), 0, log(1e12), limit=500)[0]


def chi2_bao(g0, Om):
    g = lambda zz: g0 * E(zz, Om); out = []
    for zz, k in zip(zb, kb):
        dm = quad(lambda t: 1 / g(t), 0, zz)[0]; dh = 1 / g(zz)
        out.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    r = vb - np.array(out); return float(r @ Cbi @ r)


def chi2_sn(Om):
    grid = np.concatenate([[0], np.sort(zs)]); D = np.zeros_like(grid)
    for j in range(1, grid.size):
        D[j] = D[j - 1] + quad(lambda t: 1 / E(t, Om), grid[j - 1], grid[j])[0]
    r = ms - 5 * np.log10((1 + zs) * np.interp(zs, grid, D)); M = (r @ w1) / s11; e = r - M
    return float(e @ Csi @ e) - REF


def bao_fit(Om):
    r = minimize(lambda q: chi2_bao(q[0], Om), [0.034], method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-7))
    return r.fun, r.x[0]


print(f'h = {h}; Or = {Or:.3e}')
Om1 = brentq(lambda om: T0g0(om) - 1.0, 0.05, 0.6, xtol=1e-10)
cb, g0 = bao_fit(Om1); p = chi2dist.sf(cb, 12); sn = chi2_sn(Om1)
print(f'\nA1 (T0*g0 = 1): Om = {Om1:.4f}, OB = {1 - Or - Om1:.4f}')
print(f'   BAO chi2 = {cb:.1f} / 12, probability {p:.3g} -> {"ok" if p >= 0.05 else "no"}; SN d chi2 = {sn:+.1f} -> '
      f'{"ok" if sn <= 11 else "no"}  =>  {"PASS" if p >= 0.05 and sn <= 11 else "FAIL"}')

# comparison: free Om
best = minimize(lambda q: chi2_bao(q[0], q[1]), [0.034, 0.3], method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-7))
g0f, Omf = best.x
print(f'\ncomparison, Om free (BAO): Om = {Omf:.4f}, OB = {1 - Or - Omf:.4f}, BAO chi2 = {best.fun:.1f} / 11, SN d chi2 = {chi2_sn(Omf):+.1f}, '
      f'T0*g0 = {T0g0(Omf):.4f}')
grid = np.linspace(0.20, 0.40, 41)
sng = [chi2_sn(om) for om in grid]; jm = int(np.argmin(sng))
print(f'comparison, Om free (SN only): Om ~ {grid[jm]:.3f}, SN d chi2 = {sng[jm]:+.1f}, T0*g0 = {T0g0(grid[jm]):.4f}')
GYR = 1e9 * 365.25 * 86400
age = T0g0(Om1) / (h * 1e5 / MPC) / GYR
print(f'comparison: A1 age at h = {h}: {age:.2f} Gyr (oldest stars ~13.5); at free BAO Om: {T0g0(Omf) / (h * 1e5 / MPC) / GYR:.2f} Gyr')
