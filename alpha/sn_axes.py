"""SN luminosity as the mean over three tube axes (R-130). Rules: sn_axes_rules_2026-11-26.md.
L/L0 = (1 + 2 x0 T/T0)/(1 + 2 x0), T/T0 = (1+z)^-3/2; k ~ T^(2/3), D2. usage: python -P -E -X utf8 sn_axes.py <dir with Pantheon+ files>"""
import sys
import numpy as np

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z); w1 = Ci @ one; s11 = one @ w1
REF = 1403.3
y = m - 5 * np.log10((1 + z) * 3 * (1 - (1 + z) ** -0.5)); t = (1 + z) ** -1.5
def chi2(x0):
    r = y + 2.5 * np.log10((1 + 2 * x0 * t) / (1 + 2 * x0)); M = (r @ w1) / s11; e = r - M
    return float(e @ Ci @ e)
b0 = chi2(1.0) - REF
print(f'B0 x0 = 1 (now = horn torus), a = 2/3: d chi2 = {b0:+.1f} (threshold 11) -> {"PASS" if b0 <= 11 else "FAIL"}')
X = np.linspace(0.3, 3.0, 2701); cs = np.array([chi2(x) for x in X]); i = int(np.argmin(cs)); ok = X[cs - cs[i] <= 1]
a = lambda x: 2 * x / (1 + 2 * x)
print(f'B1 x0 free: x0 = {X[i]:.3f} (1 sigma {ok.min():.3f} .. {ok.max():.3f}), a = {a(X[i]):.3f}; d chi2 = {cs[i] - REF:+.1f} (threshold 10)'
      f' -> {"PASS" if cs[i] - REF <= 10 else "FAIL"}')
print(f'   x0 = 1 lies {np.sqrt(max(chi2(1.0) - cs[i], 0)):.1f} sigma from the best')
print('B2 geometric mean (L ~ k, q = 1): d chi2 = +63.1 (R-128) -> FAIL')
