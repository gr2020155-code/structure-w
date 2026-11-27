"""Density gradient centre -> edge: SNe by sky direction (R-140). Rules: sky_gradient_rules_2026-11-27.md.
usage: python -P -E -X utf8 sky_gradient.py <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, log, radians, degrees, acos
from scipy.integrate import quad
from scipy.stats import chi2 as chi2dist

snd = sys.argv[1]; OM, CKM, K = 0.297, 299792.458, 5 / log(10)
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda nm: np.array([float(r[hdr.index(nm)]) for r in data])
zHD, zCMB, m, mu_c = col('zHD'), col('zCMB'), col('m_b_corr'), col('CEPH_DIST')
ra, dec = np.radians(col('RA')), np.radians(col('DEC')); cal = col('IS_CALIBRATOR') == 1
nv = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], 1)
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
E = lambda t: sqrt(OM * (1 + t)**3 + 1 - OM)
zg = np.concatenate([[0], np.geomspace(1e-5, 2.6, 30000)]); Eg = np.sqrt(OM * (1 + zg)**3 + 1 - OM)
cum = lambda y: np.concatenate([[0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(zg))])
Dg = cum(1 / Eg); Wg = np.where(zg > 0, cum(Dg / Eg) / np.where(Dg > 0, Dg, 1), 0)

# CMB dipole direction (memory): l = 264, b = 48 -> equatorial
def gal2eq(l, b):
    l, b = radians(l), radians(b)
    xg = np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])
    Rm = np.array([[-0.0548755604, 0.4941094279, -0.8676661490], [-0.8734370902, -0.4448296300, -0.1980763734], [-0.4838350155, 0.7469822445, 0.4559837762]])
    return Rm @ xg
cmb = gal2eq(264.0, 48.0)

def run(z, ztag):
    use = cal | (z > 0.01); idx = np.where(use)[0]; Ci = np.linalg.inv(C_all[np.ix_(idx, idx)])
    c_ = cal[idx]; zz = z[idx]; D = np.interp(zz, zg, Dg); w = np.interp(zz, zg, Wg); fl = (~c_).astype(float)
    y = np.where(c_, m[idx] - mu_c[idx], m[idx] - 5 * np.log10((1 + zz) * np.where(c_, 1, D) * CKM / 100) - 25)
    base = [np.ones(len(idx)), -fl]                                   # M, 5 lg h
    def gls(cols):
        X = np.stack(cols, 1); Fm = X.T @ Ci @ X; b = np.linalg.solve(Fm, X.T @ Ci @ y); r = y - X @ b
        return b, np.linalg.inv(Fm), float(r @ Ci @ r)
    b0, _, X0 = gls(base)
    print(f'\n=== {ztag}: {len(idx)} SNe ({c_.sum()} calibrators); base h = {10**(b0[1] / 5):.4f}, chi2 = {X0:.1f}')
    for tag, wt, unit in (('A constant tempo dipole d', np.ones(len(idx)), 'fraction'),
                          ('B linear gradient G (per c/H0)', w, 'fraction per c/H0')):
        cols = base + [-K * fl * wt * nv[idx, j] for j in range(3)]
        b, cov, X1 = gls(cols); v = b[2:]; cv3 = cov[2:, 2:]; amp = float(np.linalg.norm(v))
        J = v / amp; samp = float(sqrt(J @ cv3 @ J)); dchi = X0 - X1; p = chi2dist.sf(dchi, 3)
        raD = degrees(np.arctan2(v[1], v[0])) % 360; decD = degrees(np.arcsin(v[2] / amp))
        ang = degrees(acos(np.clip(J @ cmb, -1, 1)))
        lim = sqrt(chi2dist.ppf(0.95, 3)) * sqrt(np.max(np.linalg.eigvalsh(cv3)))
        print(f'   {tag}: amplitude = {amp:.4f} +- {samp:.4f} {unit}; toward RA = {raD:.0f}, Dec = {decD:+.0f} '
              f'(angle to CMB dipole {ang:.0f} deg); d chi2 = {dchi:.1f} (3 dof, p = {p:.3f}) -> '
              f'{"FOUND" if dchi >= 7.81 else "not found"}; 95% bound on any component direction ~ {lim:.4f}')
        print(f'      vector (x, y, z) = ({v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f}), sigmas ({sqrt(cv3[0, 0]):.4f}, {sqrt(cv3[1, 1]):.4f}, {sqrt(cv3[2, 2]):.4f})')
    b, cov, X1 = gls(base + [K * fl * w]); Q, sQ = b[2], sqrt(cov[2, 2]); dchi = X0 - X1
    print(f'   C isotropic outward fall Q (per c/H0): Q = {Q:+.4f} +- {sQ:.4f}; d chi2 = {dchi:.1f} -> {"FOUND" if dchi >= 3.84 else "not found"}; '
          f'local h = {10**(b[1] / 5):.4f}')
    hi = cal[idx] | (zz > 0.15)
    print(f'      (flow SNe with z > 0.15: {int((~c_ & (zz > 0.15)).sum())}, mean w there = {w[~c_ & (zz > 0.15)].mean():.3f} c/H0)')

print(f'CMB dipole direction (memory, l=264 b=48): RA = {degrees(np.arctan2(cmb[1], cmb[0])) % 360:.0f}, Dec = {degrees(np.arcsin(cmb[2])):+.0f}')
run(zHD, 'zHD (with peculiar-velocity model, bridge)')
run(zCMB, 'zCMB (CMB frame only)')
