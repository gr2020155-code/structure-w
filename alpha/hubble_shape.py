"""Hubble diagram shape vs Pantheon+. Rules: hubble_shape_rules_2026-11-26.md (pushed as 972aac1 before running).

usage: python hubble_shape.py <dir with Pantheon+SH0ES.dat and STATSYS.cov>
"""
import sys
import numpy as np

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n')
hdr = rows[0].split()
data = [r.split() for r in rows[1:] if r.strip()]
iz, im, ie = hdr.index('zHD'), hdr.index('m_b_corr'), hdr.index('m_b_corr_err_DIAG')
z_all = np.array([float(r[iz]) for r in data])
m_all = np.array([float(r[im]) for r in data])
e_all = np.array([float(r[ie]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov')
n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01
z, m, e = z_all[sel], m_all[sel], e_all[sel]
C = C_all[np.ix_(sel, sel)]
Ci = np.linalg.inv(C)
Cdi = np.diag(1 / e**2)
print(f'SNe used: {sel.sum()} of {n} (zHD > 0.01), z up to {z.max():.2f}')

def dl_shape(name, z, Om=None):        # d_L in units of c/H0
    if name == 'W1': return z
    if name == 'W2': return (1 + z) * np.log1p(z)
    if name == 'W3': return z * (1 + z)
    zz = np.linspace(0, z.max(), 20001)   # flat LCDM
    integ = 1 / np.sqrt(Om * (1 + zz)**3 + 1 - Om)
    chi = np.concatenate([[0], np.cumsum((integ[1:] + integ[:-1]) / 2 * np.diff(zz))])
    return (1 + z) * np.interp(z, zz, chi)

def chi2(mu_shape, W):                    # free offset M, analytic
    r = m - mu_shape
    one = np.ones_like(r)
    M = (one @ W @ r) / (one @ W @ one)
    res = r - M
    return float(res @ W @ res), res

def mu(name, Om=None):
    return 5 * np.log10(dl_shape(name, z, Om))

# reference: flat LCDM, Om free
best = None
for Om in np.arange(0.0, 1.0001, 0.002):
    x2, _ = chi2(mu('L', Om), Ci)
    if best is None or x2 < best[0]: best = (x2, Om)
ref, Om_b = best
ref_d, res_ref = chi2(mu('L', Om_b), Cdi)
print(f'reference flat LCDM: Om = {Om_b:.3f}, chi2 full = {ref:.1f}, diag = {ref_d:.1f}  (dof {sel.sum() - 2})')
print()
bins = [(0.01, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 1.0), (1.0, 3.0)]
print(f'{"branch":6} {"chi2 full":>10} {"d chi2":>8} {"chi2 diag":>10} {"d diag":>8}  pass   mean residual (data - model), mag, in z bins')
print(' ' * 49 + ''.join(f'{f"{a}-{b}":>11}' for a, b in bins))
for name in ('W1', 'W2', 'W3', 'L'):
    s = mu(name, Om_b)
    x2, res = chi2(s, Ci)
    xd, resd = chi2(s, Cdi)
    br = []
    for a, b in bins:
        k = (z > a) & (z <= b)
        w = 1 / e[k]**2
        br.append(f'{np.sum(resd[k] * w) / np.sum(w):+11.3f}')
    tag = 'ref' if name == 'L' else ('YES' if x2 - ref <= 9 else 'no')
    print(f'{name:6} {x2:10.1f} {x2 - ref:8.1f} {xd:10.1f} {xd - ref_d:8.1f}  {tag:5}' + ''.join(br))
print()
print('number per bin:', [int(((z > a) & (z <= b)).sum()) for a, b in bins])
