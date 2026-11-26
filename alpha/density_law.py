"""Law of link-density change k ~ T^p vs Pantheon+. Rules: density_law_rules_2026-11-26.md (pushed as 1726cd6 before running).

usage: python density_law.py <dir with Pantheon+SH0ES.dat and STATSYS.cov>
"""
import sys
import numpy as np

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
one = np.ones_like(z)

def chi2(dl):
    r = m - 5 * np.log10(dl)
    M = (one @ Ci @ r) / (one @ Ci @ one)
    return float((r - M) @ Ci @ (r - M)), r - M

def dl(p, D):                     # in units of c*T0
    if D == 'D1':
        dd = 1 - (1 + z) ** (-1 / p)
    else:
        dd = np.log1p(z) if abs(p - 1) < 1e-12 else ((1 + z) ** ((p - 1) / p) - 1) / (p - 1)
    return dd * (1 + z)

def lcdm(Om):
    zz = np.linspace(0, z.max(), 20001)
    f = 1 / np.sqrt(Om * (1 + zz)**3 + 1 - Om)
    chi = np.concatenate([[0], np.cumsum((f[1:] + f[:-1]) / 2 * np.diff(zz))])
    return (1 + z) * np.interp(z, zz, chi)

ref = min(chi2(lcdm(Om))[0] for Om in np.arange(0, 1.0001, 0.002))
print(f'reference flat LCDM chi2 = {ref:.1f}')
tU = 13.5
H_loc = {'SH0ES': 73.04, 'CCHP': 70.39}
bins = [(0.01, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 1.0), (1.0, 3.0)]
print()
print(f'{"variant":12} {"d chi2":>8}  {"T0 = p/H0, Gyr (SH0ES, CCHP)":>30}  shape age  PASS   mean residual per z bin (diag-weighted not used; plain mean)')
for p, src in ((1, 'F d=1'), (2, 'F d=2, V d=1'), (3, 'F d=3, V d=2'), (4, 'V d=3')):
    for D in ('D1', 'D2'):
        x2, res = chi2(dl(p, D))
        T0 = {k: p * 977.79 / h for k, h in H_loc.items()}
        sh = x2 - ref <= 9
        ag = all(t >= tU for t in T0.values())
        br = ''.join(f'{res[(z > a) & (z <= b)].mean():+8.3f}' for a, b in bins)
        print(f'p={p} {D:7} {x2 - ref:8.1f}  {T0["SH0ES"]:14.1f} {T0["CCHP"]:6.1f}       '
              f'{"ok" if sh else "no":5} {"ok" if ag else "no":4} {"YES" if sh and ag else "no":5}{br}   ({src})')
print()
print('observation (fit, not a derivation): free p')
for D in ('D1', 'D2'):
    ps = np.arange(0.5, 6.0001, 0.01)
    xs = np.array([chi2(dl(p, D))[0] for p in ps])
    i = int(np.argmin(xs)); ok = ps[xs - xs[i] <= 1]
    print(f'  {D}: best p = {ps[i]:.2f} (1 sigma {ok.min():.2f}-{ok.max():.2f}), d chi2 vs reference = {xs[i] - ref:+.1f}, '
          f'T0 = p/H0 = {ps[i] * 977.79 / 73.04:.1f}-{ps[i] * 977.79 / 70.39:.1f} Gyr')
