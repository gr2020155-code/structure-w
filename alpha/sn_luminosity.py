"""SN luminosity depends on link density at the explosion site (R-128). Rules: sn_luminosity_rules_2026-11-26.md.
L_e/L_0 = (1+z)^(-q); k ~ T^(2/3), D2. usage: python -P -E -X utf8 sn_luminosity.py <dir with Pantheon+ files>"""
import sys
import numpy as np

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z)
REF = 1403.3

def fit(mod):
    r = m - mod; M = (one @ Ci @ r) / (one @ Ci @ one); x = r - M
    return float(x @ Ci @ x), x
base = 5 * np.log10((1 + z) * 3 * (1 - (1 + z) ** -0.5))
lum = lambda q: 2.5 * q * np.log10(1 + z)
b1 = (z > 0.6) & (z <= 1.0); b2 = z > 1.0
print(f'reference chi2 = {REF}; branches threshold 11, free q threshold 10')
for q, why in ((0, 'no dependence'), (1/3, 'linear density'), (2/3, 'area density'), (1, 'k (A-36 reading)'), (2, 's^(2/3)'), (3, 's')):
    ch, x = fit(base + lum(q))
    print(f'   q = {q:.3f} ({why:17s}): d chi2 = {ch - REF:+7.1f} -> {"PASS" if ch - REF <= 11 else "FAIL"};'
          f' mean residual z 0.6-1: {x[b1].mean():+.3f}, z > 1: {x[b2].mean():+.3f}')
Q = np.linspace(0, 2, 2001); cs = np.array([fit(base + lum(q))[0] for q in Q]); i = int(np.argmin(cs))
ok = Q[cs - cs[i] <= 1]; x = fit(base + lum(Q[i]))[1]
print(f'free q (observation): q = {Q[i]:.3f} (1 sigma {ok.min():.3f} .. {ok.max():.3f}), d chi2 = {cs[i] - REF:+.1f} -> {"PASS" if cs[i] - REF <= 10 else "FAIL"}')
print(f'   mean residual z 0.6-1: {x[b1].mean():+.3f}, z > 1: {x[b2].mean():+.3f}')
for lo, hi in ((0.01, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 1.0), (1.0, 2.3)):
    s = (z > lo) & (z <= hi); print(f'   z {lo}-{hi}: N = {s.sum():4d}, mean residual {x[s].mean():+.3f}')
