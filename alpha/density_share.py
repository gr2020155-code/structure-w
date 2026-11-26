"""Link-density law from the share of 3D-distinguishable points (R-123). Rules: density_share_rules_2026-11-26.md.
usage: python -P -E -X utf8 density_share.py <dir with PantheonPlusSH0ES.dat and STATSYS.cov>"""
import sys
import numpy as np
from math import sqrt, log

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n')
hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z)
REF = 1403.3; THR = 11.0; theta_data = 0.010411

def chi2(dl):
    r = m - 5 * np.log10(dl); M = (one @ Ci @ r) / (one @ Ci @ one)
    return float((r - M) @ Ci @ (r - M))
def dl_D2(p):
    dd = np.log1p(z) if abs(p - 1) < 1e-12 else ((1 + z) ** ((p - 1) / p) - 1) / (p - 1)
    return dd * (1 + z)

def I(p, a, b):                                   # int_a^b u^-p du
    if abs(p - 1) < 1e-12: return log(b / a)
    return (b**(1 - p) - a**(1 - p)) / (1 - p)
def theta(p, zr1, cs, zmax1=None):
    x = zr1**(-1 / p)
    x0 = 0.0 if p < 1 else zmax1**(-1 / p)        # finite start (A-100) needed only for p >= 1
    return cs * I(p, x0, x) / I(p, x, 1.0) if (p < 1 or x0 > 0) else float('inf')

branches = [('F1 G1 (s~T, k~s)', 1.0), ('F1 G2 (s~T, k~s^1/3)', 1 / 3), ('F2 G1 (torus front, k~s)', 2.0),
            ('F2 G2 (torus front, k~s^1/3)', 2 / 3), ('F3 G1 (ball front, k~s)', 3.0), ('F3 G2 (ball front, k~s^1/3)', 1.0)]
print(f'SN reference (data-only, R-90) chi2 = {REF}; threshold d chi2 <= {THR}; relic theta* = {theta_data}\n')
for name, p in branches:
    dx = chi2(dl_D2(p)) - REF
    sn = dx <= THR
    if p < 1:
        vals = [theta(p, zr, cs) for zr in (1000, 1100, 1400) for cs in (1 / sqrt(6), 1 / sqrt(3))]
        rel = min(vals) <= theta_data <= max(vals); trng = f'{min(vals):.3g} .. {max(vals):.3g}'
    else:
        out = []
        for zm in (1e4, 1e8, 1e20):
            vals = [theta(p, zr, cs, zm) for zr in (1000, 1100, 1400) for cs in (1 / sqrt(6), 1 / sqrt(3))]
            out.append((min(vals), max(vals)))
        rel = any(lo <= theta_data <= hi for lo, hi in out)
        trng = '; '.join(f'zmax {zm:.0e}: {lo:.3g}..{hi:.3g}' for zm, (lo, hi) in zip((1e4, 1e8, 1e20), out))
    print(f'{name:32s} p = {p:.3f}: SN d chi2 = {dx:+7.1f} ({"ok" if sn else "no"}); relic theta = {trng} ({"ok" if rel else "no"})'
          f' -> {"PASS" if sn and rel else "FAIL"}')
