"""SN graph near to far and function search (R-129). Rules: sn_shape_rules_2026-11-26.md.
dm(z) = m_b_corr - 5 lg dL - M, dL = (1+z) 3 (1 - (1+z)^-1/2) (k ~ T^(2/3), D2).
usage: python -P -E -X utf8 sn_shape.py <dir with Pantheon+ files> <output png>"""
import sys
import numpy as np
from scipy.optimize import minimize

d, png = sys.argv[1], sys.argv[2]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z); w1 = Ci @ one; s11 = one @ w1
REF = 1403.3
y = m - 5 * np.log10((1 + z) * 3 * (1 - (1 + z) ** -0.5))

def chi2(f):                                   # M free
    r = y - f; M = (r @ w1) / s11; x = r - M
    return float(x @ Ci @ x), M
lgL = lambda L: -2.5 * np.log10(np.clip(L, 1e-9, None))
zz = 1 + z
fam = {
    'F1 L ~ k^q':                    (1, lambda p, u: 2.5 * p[0] * np.log10(u), [0.8]),
    'F2 L = 1 - a(1 - k/k0)':        (1, lambda p, u: lgL(1 - p[0] * (1 - 1 / u)), [0.5]),
    'F3 L = 1 - a(1 - T/T0)':        (1, lambda p, u: lgL(1 - p[0] * (1 - u ** -1.5)), [0.5]),
    'F4 L = 1 - a(1 - s/s0)':        (1, lambda p, u: lgL(1 - p[0] * (1 - u ** -3.0)), [0.4]),
    'F5 dm = a(1 - (1+z)^-b)':       (2, lambda p, u: p[0] * (1 - u ** -p[1]), [0.8, 1.0]),
    'F6 L = 1 - a(1 - (k/k0)^b)':    (2, lambda p, u: lgL(1 - p[0] * (1 - u ** -p[1])), [0.5, 1.0]),
}
best = {}
print(f'reference chi2 = {REF}; threshold 10 (1 number), 9 (2 numbers)')
for name, (k, f, p0) in fam.items():
    obj = lambda p: chi2(f(p, zz))[0]
    r = minimize(obj, p0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-4, maxiter=4000))
    for start in ([0.3] * k, [1.5] * k, [0.9] * k):      # restarts
        r2 = minimize(obj, start, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-4, maxiter=4000))
        if r2.fun < r.fun: r = r2
    dx = r.fun - REF; thr = 11 - k
    # 1-sigma of the first number (profile over the rest)
    a0 = r.x[0]; grid = np.linspace(a0 * 0.7, a0 * 1.3, 121)
    prof = []
    for a in grid:
        if k == 1: prof.append(obj([a]))
        else: prof.append(minimize(lambda q: obj([a, q[0]]), [r.x[1]], method='Nelder-Mead').fun)
    okk = grid[np.array(prof) - r.fun <= 1]
    best[name] = (k, f, r.x)
    extra = '' if k == 1 else f', b = {r.x[1]:.3f}'
    print(f'   {name:28s}: a = {r.x[0]:.3f} (1 sigma {okk.min():.3f} .. {okk.max():.3f}){extra}; d chi2 = {dx:+.1f} -> {"PASS" if dx <= thr else "FAIL"}')

# binned values: one free number per bin (no separate M)
edges = np.array([0.01, 0.02, 0.035, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0, 1.3, 2.3])
idx = np.digitize(z, edges) - 1; X = np.zeros((z.size, edges.size - 1)); X[np.arange(z.size), idx] = 1
F = X.T @ Ci @ X; cov_b = np.linalg.inv(F); vb = cov_b @ (X.T @ Ci @ y); eb = np.sqrt(np.diag(cov_b))
zb = np.array([np.exp(np.mean(np.log(z[idx == i]))) for i in range(edges.size - 1)]); nb = np.bincount(idx)
# zero of the scale from the best 1-number family
order = sorted(fam, key=lambda nm: chi2(best[nm][1](best[nm][2], zz))[0] - REF)
nm1 = [nm for nm in order if best[nm][0] == 1][0]
M0 = chi2(best[nm1][1](best[nm1][2], zz))[1]
print(f'\nbinned dm (zero from M of {nm1.split()[0]}), z = geometric mean in bin:')
for i in range(edges.size - 1):
    L = 10 ** (-0.4 * (vb[i] - M0))
    print(f'   z {edges[i]:.3f}-{edges[i+1]:.3f}: N = {nb[i]:4d}, z = {zb[i]:.3f}, k/k0 = {1/(1+zb[i]):.3f}, dm = {vb[i]-M0:+.3f} +- {eb[i]:.3f}, L/L0 = {L:.3f}')

# graph
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK, INK2, MUTED, GRID, AXIS, SURF = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#6250d6']
show = [nm1] + [nm for nm in order if nm != nm1][:3]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURF)
zl = np.exp(np.linspace(np.log(0.01), np.log(2.3), 300)); ul = 1 + zl
for a in ax:
    a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.6); a.tick_params(colors=MUTED, labelsize=9)
    for s in a.spines.values(): s.set_color(AXIS)
    a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
for j, nm in enumerate(show):
    k, f, p = best[nm]; Mj = chi2(f(p, zz))[1]
    dmc = f(p, ul) + (Mj - M0)
    lab = f'{nm.split()[0]}: ' + nm.split(' ', 1)[1] + f'  (Δχ² {chi2(f(p, zz))[0] - REF:+.1f})'
    ax[0].plot(zl, dmc, color=SER[j], lw=2, label=lab)
    ax[1].plot(1 / ul, 10 ** (-0.4 * dmc), color=SER[j], lw=2)
ax[0].errorbar(zb, vb - M0, yerr=eb, fmt='o', ms=6, color=INK, ecolor=INK2, elinewidth=1.2, capsize=0,
               mec=SURF, mew=1.5, label='Pantheon+, среднее по интервалу', zorder=5)
L_b = 10 ** (-0.4 * (vb - M0)); eL = L_b * 0.4 * np.log(10) * eb
ax[1].errorbar(1 / (1 + zb), L_b, yerr=eL, fmt='o', ms=6, color=INK, ecolor=INK2, elinewidth=1.2, capsize=0, mec=SURF, mew=1.5, zorder=5)
ax[0].axhline(0, color=AXIS, lw=1); ax[1].axhline(1, color=AXIS, lw=1)
ax[0].set_xscale('log'); ax[0].set_xlabel('z (от близких к дальним)', color=INK2)
ax[0].set_ylabel('Δm: насколько тусклее закона k ∝ T^(2/3), зв. вел.', color=INK2)
ax[0].set_title('Сверхновые: добавочное потускнение', color=INK, loc='left', fontsize=12)
ax[1].set_xlabel('плотность связей в месте взрыва k/k₀ = 1/(1+z)', color=INK2)
ax[1].set_ylabel('нужная светимость L/L₀', color=INK2)
ax[1].set_title('То же как светимость от плотности', color=INK, loc='left', fontsize=12)
ax[1].invert_xaxis()
ax[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc='upper left')
fig.tight_layout(); fig.savefig(png, dpi=150, facecolor=SURF)
print(f'\ngraph: {png}')
