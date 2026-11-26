"""Check of Pantheon+ processing (R-127). Rules: sn_processing_rules_2026-11-26.md.
usage: python -P -E -X utf8 sn_processing.py <dir with Pantheon+ files>"""
import sys
import numpy as np

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda name: np.array([float(r[hdr.index(name)]) for r in data])
cv = np.loadtxt(f'{d}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
sel = col('zHD') > 0.01                                  # same set of SN in every variant
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)])
g = {k: col(k)[sel] for k in ('zHD', 'zCMB', 'm_b_corr', 'mB', 'x1', 'c', 'biasCor_m_b')}
A0, B0, K0 = 0.1475, 3.0922, -0.0524                     # decomposition of m_b_corr (rules)
print(f'decomposition check: max |m_b_corr - (mB + a x1 - b c - biasCor + const)| = '
      f'{np.max(np.abs(g["m_b_corr"] - (g["mB"] + A0*g["x1"] - B0*g["c"] - g["biasCor_m_b"] + K0))):.1e}')

def gls(y, X):
    """min chi2 of y = X b + noise with covariance C."""
    XtC = X.T @ Ci; b = np.linalg.solve(XtC @ X, XtC @ y); r = y - X @ b
    return float(r @ Ci @ r), r
def dl_p(z, p):
    dd = np.log1p(z) if abs(p - 1) < 1e-12 else ((1 + z) ** ((p - 1) / p) - 1) / (p - 1)
    return dd * (1 + z)

variants = {
    'P0 authors':                 dict(bias=True,  free=False, z='zHD'),
    'P1 no selection bias corr.': dict(bias=False, free=False, z='zHD'),
    'P2 alpha, beta free':        dict(bias=True,  free=True,  z='zHD'),
    'P3 zCMB':                    dict(bias=True,  free=False, z='zCMB'),
    'P4 all three':               dict(bias=False, free=True,  z='zCMB'),
}
P = np.linspace(0.3, 3.0, 271)
for name, v in variants.items():
    z = g[v['z']]; one = np.ones_like(z)
    if v['free']:
        y = g['mB'] - (g['biasCor_m_b'] if v['bias'] else 0.0); nuis = [one, -g['x1'], g['c']]
    else:
        y = g['mB'] + A0*g['x1'] - B0*g['c'] - (g['biasCor_m_b'] if v['bias'] else 0.0); nuis = [one]
    x = np.log1p(z)
    ref, _ = gls(y - 5*np.log10(x), np.stack(nuis + [x, x*x], 1))
    Xn = np.stack(nuis, 1)
    c23, r23 = gls(y - 5*np.log10(dl_p(z, 2/3)), Xn); c1, _ = gls(y - 5*np.log10(dl_p(z, 1.0)), Xn)
    cp = [gls(y - 5*np.log10(dl_p(z, p)), Xn)[0] for p in P]; ip = int(np.argmin(cp))
    okp = P[np.array(cp) - cp[ip] <= 1]
    b1 = (z > 0.6) & (z <= 1.0); b2 = z > 1.0
    print(f'\n{name}: reference chi2 = {ref:.1f}')
    print(f'   k~T^(2/3): d chi2 = {c23 - ref:+.1f} (threshold 11) -> {"PASS" if c23 - ref <= 11 else "FAIL"};'
          f' mean residual z 0.6-1: {r23[b1].mean():+.3f}, z > 1: {r23[b2].mean():+.3f} mag (data minus model, + = fainter)')
    print(f'   k~T: d chi2 = {c1 - ref:+.1f};  free p (observation): p = {P[ip]:.2f} ({okp.min():.2f} .. {okp.max():.2f}), d chi2 = {cp[ip] - ref:+.1f}')
