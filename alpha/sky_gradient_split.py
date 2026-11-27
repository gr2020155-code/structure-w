"""Observation OUTSIDE the R-140 rules (disclosed): where does the zCMB dipole A sit in z?
Dipole A fitted with the flow restricted to z bins (calibrators always in).
usage: python -P -E -X utf8 sky_gradient_split.py <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, log, degrees

snd = sys.argv[1]; OM, CKM, K = 0.297, 299792.458, 5 / log(10)
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda nm: np.array([float(r[hdr.index(nm)]) for r in data])
m, mu_c = col('m_b_corr'), col('CEPH_DIST'); cal = col('IS_CALIBRATOR') == 1
ra, dec = np.radians(col('RA')), np.radians(col('DEC'))
nv = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], 1)
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
zg = np.concatenate([[0], np.geomspace(1e-5, 2.6, 30000)]); Eg = np.sqrt(OM * (1 + zg)**3 + 1 - OM)
Dg = np.concatenate([[0], np.cumsum((1 / Eg[1:] + 1 / Eg[:-1]) / 2 * np.diff(zg))])
for ztag in ('zCMB', 'zHD'):
    z = col(ztag); print(f'== {ztag}')
    for lo, hi in ((0.01, 0.03), (0.03, 0.06), (0.06, 0.15), (0.15, 2.5)):
        idx = np.where(cal | ((z > lo) & (z <= hi)))[0]; Ci = np.linalg.inv(C_all[np.ix_(idx, idx)]); c_ = cal[idx]; fl = (~c_).astype(float)
        D = np.interp(z[idx], zg, Dg)
        y = np.where(c_, m[idx] - mu_c[idx], m[idx] - 5 * np.log10((1 + z[idx]) * np.where(c_, 1, D) * CKM / 100) - 25)
        def gls(cols):
            X = np.stack(cols, 1); Fm = X.T @ Ci @ X; b = np.linalg.solve(Fm, X.T @ Ci @ y); r = y - X @ b
            return b, np.linalg.inv(Fm), float(r @ Ci @ r)
        base = [np.ones(len(idx)), -fl]; _, _, X0 = gls(base)
        b, cov, X1 = gls(base + [-K * fl * nv[idx, j] for j in range(3)]); v = b[2:]; amp = np.linalg.norm(v)
        s = sqrt((v / amp) @ cov[2:, 2:] @ (v / amp))
        vel = amp * CKM * np.median(z[idx][~c_])
        print(f'   z {lo}-{hi} ({int((~c_).sum())} flow SNe): dipole {amp:.4f} +- {s:.4f}, RA {degrees(np.arctan2(v[1], v[0])) % 360:.0f}, '
              f'Dec {degrees(np.arcsin(v[2] / amp)):+.0f}, d chi2 = {X0 - X1:.1f}; as a velocity ~ {vel:.0f} km/s at median z')
