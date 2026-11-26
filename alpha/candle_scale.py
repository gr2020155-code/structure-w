"""Candle scale checks on Pantheon+SH0ES (R-137). Rules: candle_scale_rules_2026-11-27.md.
usage: python -P -E -X utf8 candle_scale.py <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, log
from scipy.integrate import quad

snd = sys.argv[1]; OM = 0.297; CKM = 299792.458
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda nm: np.array([float(r[hdr.index(nm)]) for r in data])
z, m, mu_c = col('zHD'), col('m_b_corr'), col('CEPH_DIST')
cal, hf = col('IS_CALIBRATOR') == 1, col('USED_IN_SH0ES_HF') == 1
mass, c, x1, surv = col('HOST_LOGMASS'), col('c'), col('x1'), col('IDSURVEY').astype(int)
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
E = lambda t: sqrt(OM * (1 + t)**3 + 1 - OM)
Dc = np.array([quad(lambda t: 1 / E(t), 0, zz)[0] if h_ else 0.0 for zz, h_ in zip(z, hf)])

def fit(scal, shf):
    use = scal | shf; Ci = np.linalg.inv(C_all[np.ix_(use, use)]); cu = scal[use]
    y = np.where(cu, m[use] - mu_c[use], m[use] - 5 * np.log10((1 + z[use]) * np.where(cu, 1.0, Dc[use]) * CKM / 100) - 25)
    X = np.stack([np.ones(use.sum()), np.where(cu, 0.0, -1.0)], 1)          # y = M + 5 lg h * (-1 for flow)
    Fm = X.T @ Ci @ X; b = np.linalg.solve(Fm, X.T @ Ci @ y); cov = np.linalg.inv(Fm)
    h = 10**(b[1] / 5); return h, h * log(10) / 5 * sqrt(cov[1, 1]), b[0], sqrt(cov[0, 0]), int(scal.sum()), int(shf.sum())

def show(tag, r):
    h, sh, M, sM, nc, nh = r
    print(f'   {tag}: h = {h:.4f} +- {sh:.4f} (vs 0.685: {(h - 0.685) / sh:+.1f} sigma), M = {M:.3f} +- {sM:.3f}; n_cal = {nc}, n_flow = {nh}')
    return r
def pair(name, a, b):
    print(f'== {name}')
    ra = show('class A', fit(cal & a, hf & a)); rb = show('class B', fit(cal & b, hf & b))
    d = (ra[0] - rb[0]) / sqrt(ra[1]**2 + rb[1]**2); near = any(abs(r[0] - 0.685) <= 2 * r[1] for r in (ra, rb))
    print(f'   difference A - B: {ra[0] - rb[0]:+.4f} ({d:+.1f} sigma); dM-equivalent {5 * np.log10(ra[0] / rb[0]):+.3f} mag'
          f' -> {"EXPLAINS" if near and abs(d) > 2 else "does not explain"}')

print(f'needed candle shift for 0.735 -> 0.685: dM = {5 * np.log10(0.735 / 0.685):.3f} mag')
print('== all'); show('all calibrators + all flow', fit(cal, hf))
ok = mass > 0
pair('C1 host mass: A = log M < 10, B = >= 10', ok & (mass < 10), ok & (mass >= 10))
pair('C2 SN colour: A = c < 0, B = c >= 0', c < 0, c >= 0)
pair('C3 SN stretch: A = x1 < 0, B = x1 >= 0', x1 < 0, x1 >= 0)
print('== C4 surveys: flow only from surveys that have calibrators')
show('matched surveys', fit(cal, hf & np.isin(surv, np.unique(surv[cal]))))
print('== C5 calibrator distance: A = mu < 32.33 (near), B = far; each with all flow')
ra = show('near calibrators', fit(cal & (mu_c < 32.33), hf)); rb = show('far calibrators', fit(cal & (mu_c >= 32.33), hf))
d = (ra[0] - rb[0]) / sqrt(ra[1]**2 + rb[1]**2)
print(f'   difference: {ra[0] - rb[0]:+.4f} ({d:+.1f} sigma); M near - far = {ra[2] - rb[2]:+.3f} mag'
      f' -> {"EXPLAINS" if any(abs(r[0] - 0.685) <= 2 * r[1] for r in (ra, rb)) and abs(d) > 2 else "does not explain"}')
print('== C6 flow by z, all calibrators')
for lo, hi in ((0.023, 0.04), (0.04, 0.07), (0.07, 0.15)):
    show(f'z {lo}-{hi}', fit(cal, hf & (z >= lo) & (z < hi)))
