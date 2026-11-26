"""Supernova shape test with a data-only yardstick. Rules: sn_empirical_rules_2026-11-26.md (pushed as f23b106 before running).

usage: python sn_empirical.py <dir with Pantheon+SH0ES.dat and STATSYS.cov>
"""
import sys
import numpy as np
from math import erfc, sqrt

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n')
hdr = rows[0].split()
data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01
z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)])
N = len(z)
x = np.log1p(z)

def gls(y, X):                     # min (y - X b)^T Ci (y - X b)
    A = X.T @ Ci @ X
    b = np.linalg.solve(A, X.T @ Ci @ y)
    r = y - X @ b
    return float(r @ Ci @ r)

def pval(chi2, dof):               # Wilson-Hilferty, upper tail
    t = ((chi2 / dof) ** (1 / 3) - (1 - 2 / (9 * dof))) / sqrt(2 / (9 * dof))
    return 0.5 * erfc(t / sqrt(2))

print(f'SNe: {N}')
y0 = m - 5 * np.log10(x)
prev, K = None, None
print('data-only yardstick: mu = M + 5 lg x + b1 x + ... + bK x^K, x = ln(1+z)')
chis = {}
for k in range(1, 6):
    X = np.column_stack([np.ones(N)] + [x**j for j in range(1, k + 1)])
    chis[k] = gls(y0, X)
    print(f'  K = {k}: chi2 = {chis[k]:.1f}')
for k in range(1, 5):
    if chis[k] - chis[k + 1] < 1:
        K = k; break
K = K or 5
ref = chis[K]
thr = 9 + K
print(f'  chosen K = {K}, chi2 = {ref:.1f}, dof {N - K - 1}, p = {pval(ref, N - K - 1):.3f};  pass threshold d chi2 <= {thr}')
print()

def dl(p, D):
    if p == 'exp':
        return {'D1': (1 + z) * np.log1p(z), 'D2': z * (1 + z)}[D]
    if D == 'D1':
        dd = 1 - (1 + z) ** (-1 / p)
    else:
        dd = np.log1p(z) if abs(p - 1) < 1e-9 else ((1 + z) ** ((p - 1) / p) - 1) / (p - 1)
    return dd * (1 + z)

def chi2_branch(dlv):
    return gls(m - 5 * np.log10(dlv), np.ones((N, 1)))

print(f'{"branch":14} {"chi2":>8} {"d chi2":>8} {"p(abs)":>7}  pass')
for p in (1, 2, 3, 4, 'exp'):
    for D in ('D1', 'D2'):
        x2 = chi2_branch(dl(p, D))
        name = f'k~T^{p} {D}' if p != 'exp' else f'k~exp {D}'
        print(f'{name:14} {x2:8.1f} {x2 - ref:+8.1f} {pval(x2, N - 1):7.3f}  {"YES" if x2 - ref <= thr else "no"}')
print('(W1 = T^1 D1, W2 = T^1 D2 = exp D1, W3 = exp D2)')
print()
ps = np.arange(0.5, 4.0001, 0.01)
xs = np.array([chi2_branch(dl(p, 'D2')) for p in ps])
i = int(np.argmin(xs)); ok = ps[xs - xs[i] <= 1]
print(f'observation (fit): free p in D2 = {ps[i]:.2f} ({ok.min():.2f}-{ok.max():.2f}), d chi2 vs yardstick = {xs[i] - ref:+.1f} '
      f'(threshold for a 1-parameter fit: 9 + {K} - 1 = {thr - 1})')
