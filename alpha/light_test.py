"""Test of V-385: atomic lines and BAO distances without luminosity (R-131). Rules: light_test_rules_2026-11-26.md.
usage: python -P -E -X utf8 light_test.py <desi_dr2_bao_data.txt>"""
import sys
import numpy as np
from math import sqrt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import chi2 as chi2dist

# test 1: lines (prediction 0)
print('test 1: atomic lines, prediction 0')
for name, v, e in (('d alpha/alpha z = 1.15 (ESPRESSO)', 1.3, sqrt(1.3**2 + 0.4**2)),
                   ('d alpha/alpha z = 0.6-2.4 (29 systems)', -0.5, sqrt(0.5**2 + 0.4**2))):
    print(f'   {name}: {v:+.1f} +- {e:.2f} ppm -> {abs(v)/e:.1f} sigma from 0 -> {"PASS" if abs(v)/e <= 2 else "FAIL"}')
print('   |d mu/mu| < 5e-6 (3 sigma), z = 2.0-4.2: 0 inside -> PASS')

# test 2: BAO
lines = open(sys.argv[1]).read().split('\n')
i0, i1 = lines.index('MEAN'), lines.index('COV')
rows = [l.split() for l in lines[i0 + 1:i1] if l.strip()]
z = np.array([float(r[0]) for r in rows]); val = np.array([float(r[1]) for r in rows]); kind = [r[2] for r in rows]
C = np.array([[float(x) for x in l.split()] for l in lines[i1 + 1:] if l.strip()]); Ci = np.linalg.inv(C)

def pred(DM, DH, zz, k):
    return {'DM_over_rs': DM, 'DH_over_rs': DH, 'DV_over_rs': (zz * DM**2 * DH) ** (1 / 3)}[k]
def vec(DMf, DHf):
    return np.array([pred(DMf(zz), DHf(zz), zz, k) for zz, k in zip(z, kind)])
def fit_rd(v):                                   # r_d enters as 1/r_d: linear, solve exactly
    b = (v @ Ci @ val) / (v @ Ci @ v); r = val - b * v
    return float(r @ Ci @ r), 1 / b
def law(p):                                      # k ~ T^p, D2, units c T0
    DM = lambda zz: np.log1p(zz) if abs(p - 1) < 1e-12 else ((1 + zz) ** ((p - 1) / p) - 1) / (p - 1)
    DH = lambda zz: (1 + zz) ** (-1 / p) / p     # dD/dz, D = int (1+z) dT, T = (1+z)^(-1/p)
    return DM, DH

DM, DH = law(2 / 3)
ch, rd = fit_rd(vec(DM, DH)); lim = chi2dist.ppf(0.95, 12)
print(f'\ntest 2: BAO, k ~ T^(2/3), r_d free: chi2 = {ch:.1f} for 12 dof (pass <= {lim:.1f}), probability {chi2dist.sf(ch, 12):.2e}'
      f' -> {"PASS" if ch <= lim else "FAIL"}')
print(f'   fitted r_d = {rd:.4f} cT0')
rs = [cs * 3 * ((zd) ** -1.5) ** (1 / 3) for zd in (1000, 1400) for cs in (1 / sqrt(6), 1 / sqrt(3))]
print(f'   r_d from the relic calculation: {min(rs):.4f} .. {max(rs):.4f} cT0')
print(f'   ratios (no r_d): DM(2.33)/DM(0.51) model {DM(2.33)/DM(0.51):.3f}, data {val[12]/val[1]:.3f};'
      f' DH(0.51)/DH(2.33) model {DH(0.51)/DH(2.33):.3f}, data {val[2]/val[11]:.3f}')
b = 1 / rd; v = vec(DM, DH) * b
for zz, k, d_, m_, e_ in zip(z, kind, val, v, np.sqrt(np.diag(C))):
    print(f'   z = {zz:5.3f} {k:11s}: data {d_:7.3f} +- {e_:.3f}, model {m_:7.3f} ({(d_ - m_)/e_:+6.1f} sigma)')

print('\nfor comparison (not counted):')
P = np.linspace(0.3, 4, 3701); cs_ = [fit_rd(vec(*law(p)))[0] for p in P]; j = int(np.argmin(cs_)); ok = P[np.array(cs_) - cs_[j] <= 1]
print(f'   k ~ T^p, p free: p = {P[j]:.3f} (1 sigma {ok.min():.3f} .. {ok.max():.3f}), chi2 = {cs_[j]:.1f} (11 dof)  [SN wanted 1.41 +- 0.05]')
def lcdm(Om):
    E = lambda zz: sqrt(Om * (1 + zz)**3 + 1 - Om)
    DMf = lambda zz: quad(lambda u: 1 / E(u), 0, zz)[0]; DHf = lambda zz: 1 / E(zz)
    return fit_rd(vec(DMf, DHf))[0]
r = minimize_scalar(lcdm, bounds=(0.1, 0.6), method='bounded')
print(f'   flat LCDM (official), Om free: Om = {r.x:.3f}, chi2 = {r.fun:.1f} (11 dof)')
