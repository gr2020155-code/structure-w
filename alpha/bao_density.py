"""Density law k(T) from the BAO ruler (R-132). Rules: bao_density_rules_2026-11-26.md.
g = d ln k / dT = 1/D_H (units 1/r_d), u = k/k0 = 1/(1+z).
usage: python -P -E -X utf8 bao_density.py <desi_dr2_bao_data.txt> <dir with Pantheon+ files> <output png>"""
import sys
import numpy as np
from math import sqrt, log
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

bao, snd, png = sys.argv[1:4]
L = open(bao).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rows = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rows]); vb = np.array([float(r[1]) for r in rows]); kb = [r[2] for r in rows]
Cb = np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]); Cbi = np.linalg.inv(Cb)

fams = {   # name: (n params, g(u, params), start)
    'G1 power k~T^p':        (2, lambda u, q: q[0] * u ** (-1 / q[1]), [0.035, 1.0]),
    'G2 g2 = A/u^3 + B':     (2, lambda u, q: np.sqrt(np.maximum(q[0] * u ** -3 + q[1], 1e-12)), [3e-4, 8e-4]),
    'G3 g = A/u^1.5 + B':    (2, lambda u, q: q[0] * u ** -1.5 + q[1], [0.015, 0.02]),
    'G4 g2 = A/u^2 + B':     (2, lambda u, q: np.sqrt(np.maximum(q[0] * u ** -2 + q[1], 1e-12)), [6e-4, 5e-4]),
    'G5 g2 = A/u^n + B':     (3, lambda u, q: np.sqrt(np.maximum(q[0] * u ** -q[2] + q[1], 1e-12)), [3e-4, 8e-4, 3.0]),
}
def DM(g, q, z):
    return quad(lambda t: 1 / g(1 / (1 + t), q), 0, z, limit=200)[0]
def bao_vec(g, q):
    out = []
    for z, k in zip(zb, kb):
        dm, dh = DM(g, q, z), 1 / g(1 / (1 + z), q)
        out.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (z * dm * dm * dh) ** (1 / 3)}[k])
    return np.array(out)
def chi2_bao(g, q):
    try:
        r = vb - bao_vec(g, q)
    except Exception:
        return 1e12
    return float(r @ Cbi @ r) if np.all(np.isfinite(r)) else 1e12

# SN
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
zs_all = np.array([float(r[hdr.index('zHD')]) for r in data]); ms_all = np.array([float(r[hdr.index('m_b_corr')]) for r in data])
c = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = zs_all > 0.01; zs, ms = zs_all[sel], ms_all[sel]
Csi = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(zs); w1 = Csi @ one; s11 = one @ w1
REF = 1403.3
def chi2_sn(g, q):
    grid = np.concatenate([[0], np.sort(zs)]); D = np.zeros_like(grid)
    for j in range(1, grid.size):
        D[j] = D[j - 1] + quad(lambda t: 1 / g(1 / (1 + t), q), grid[j - 1], grid[j])[0]
    Dz = np.interp(zs, grid, D); r = ms - 5 * np.log10((1 + zs) * Dz); M = (r @ w1) / s11; e = r - M
    return float(e @ Csi @ e)

def lnint(g, q, a, b):           # int_a^b dz/g with x = ln(1+z)
    return quad(lambda x: np.exp(x) / g(np.exp(-x), q), log(1 + a), log(1 + b), limit=400)[0]
def theta_rng(g, q):
    th = []
    for zr1 in (1000, 1100, 1400):
        D = lnint(g, q, 0, zr1 - 1); I = lnint(g, q, zr1 - 1, 1e9)
        th += [cs * I / D for cs in (1 / sqrt(6), 1 / sqrt(3))]
    return min(th), max(th)
def T_of_z(g, q, z):              # T0 - T(z) and T0 (integral of dz/((1+z) g))
    f = lambda x: 1 / g(np.exp(-x), q)
    return quad(f, 0, log(1 + z), limit=200)[0], quad(f, 0, log(1e9), limit=400)[0]

theta_data = 0.010411
print('BAO: chi2, dof, probability; SN d chi2 (no new numbers, threshold 11); relic theta range (data 0.010411)')
res = {}
for name, (k, g, q0) in fams.items():
    best = None
    for s in (1.0, 0.5, 2.0):
        r = minimize(lambda q: chi2_bao(g, q), np.array(q0) * np.array([s, s] + [1.0] * (k - 2)),
                     method='Nelder-Mead', options=dict(xatol=1e-9, fatol=1e-6, maxiter=6000))
        if best is None or r.fun < best.fun: best = r
    q = best.x; dof = 13 - k; p = chi2dist.sf(best.fun, dof)
    ok_bao = p >= 0.05
    sn = chi2_sn(g, q) - REF
    tr = theta_rng(g, q)
    gt = np.array([g(u, q) for u in np.linspace(0.02, 1, 400)])
    mono = bool(np.all(np.diff(gt) < 0))             # g decreases as k grows (R2)
    ug = np.linspace(0.02, 1, 4000)
    ugv = ug * np.array([g(u, q) for u in ug]); dd = np.diff(ugv)
    conv = ug[1:][dd > 0]; z_conv = (1 / conv.min() - 1) if conv.size else None
    dT, T0 = T_of_z(g, q, 0); T0g0 = T0 * g(1.0, q)
    res[name] = (k, g, q, best.fun)
    print(f'\n{name}: params = {np.array2string(q, precision=5)}')
    print(f'   BAO chi2 = {best.fun:.1f} / {dof} dof, probability {p:.3g} -> {"PASS" if ok_bao else "FAIL"}')
    print(f'   SN d chi2 = {sn:+.1f} -> {"PASS" if sn <= 11 else "FAIL"};  relic theta = {tr[0]:.4g} .. {tr[1]:.4g} -> '
          f'{"PASS" if tr[0] <= theta_data <= tr[1] else "FAIL"}')
    print(f'   g decreasing everywhere (R2): {mono};  growth of k convex for z < {z_conv:.3f}' if z_conv is not None else
          f'   g decreasing everywhere (R2): {mono};  growth of k never convex')
    print(f'   T0 * g0 = {T0g0:.3f} (age of the whole in units 1/g0); g0 = {g(1.0, q):.5f} / r_d')

# graph: g vs k/k0 with data points from D_H; k/k0 vs T/T0
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
INK, INK2, MUTED, GRID, AXIS, SURF = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
SER = ['#2a78d6', '#eb6834', '#1baf7a', '#6250d6']
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURF)
for a in ax:
    a.set_facecolor(SURF); a.grid(True, color=GRID, lw=0.6); a.tick_params(colors=MUTED, labelsize=9)
    for s in a.spines.values(): s.set_color(AXIS)
    a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
isel = [j for j, k in enumerate(kb) if k == 'DH_over_rs']
ud = 1 / (1 + zb[isel]); gd = 1 / vb[isel]; eg = np.sqrt(np.diag(Cb))[isel] / vb[isel] ** 2
uu = np.linspace(0.25, 1.0, 300)
show = ['G2 g2 = A/u^3 + B', 'G3 g = A/u^1.5 + B', 'G4 g2 = A/u^2 + B', 'G1 power k~T^p']
for j, nm in enumerate(show):
    k, g, q, ch = res[nm]
    ax[0].plot(uu, [g(u, q) for u in uu], color=SER[j], lw=2, label=f'{nm.split()[0]}: {nm.split(" ", 1)[1]}  (χ² {ch:.0f} / {13 - k})')
ax[0].errorbar(ud, gd, yerr=eg, fmt='o', ms=6, color=INK, ecolor=INK2, elinewidth=1.2, capsize=0, mec=SURF, mew=1.5,
               label='DESI DR2: g = 1/(D_H/r_d)', zorder=5)
ax[0].invert_xaxis(); ax[0].set_xlabel('плотность связей k/k₀ (сейчас = 1)', color=INK2)
ax[0].set_ylabel('g = d ln k / dT  (на r_d)', color=INK2)
ax[0].set_title('Относительный темп роста плотности по линейке', color=INK, loc='left', fontsize=12)
ax[0].legend(frameon=False, fontsize=8, labelcolor=INK2, loc='upper left')
k, g, q, ch = res['G2 g2 = A/u^3 + B']
zz = np.exp(np.linspace(0, log(1e4), 400)) - 1
_, T0 = T_of_z(g, q, 0)
Tz = np.array([T0 - T_of_z(g, q, z)[0] for z in zz]) / T0
ax[1].plot(Tz, 1 / (1 + zz), color=SER[0], lw=2, label='Г2 по линейке')
ax[1].plot(Tz, Tz ** (2 / 3), color=MUTED, lw=1.5, ls='--', label='k ∝ T^(2/3) (В-376)')
Td = np.array([T0 - T_of_z(g, q, z)[0] for z in zb[isel]]) / T0
ax[1].plot(Td, 1 / (1 + zb[isel]), 'o', ms=6, color=INK, mec=SURF, mew=1.5, label='такты точек DESI (по Г2)', zorder=5)
ax[1].set_xlabel('такт T/T₀ (сейчас = 1)', color=INK2); ax[1].set_ylabel('плотность связей k/k₀', color=INK2)
ax[1].set_title('Закон плотности k(T)', color=INK, loc='left', fontsize=12)
ax[1].legend(frameon=False, fontsize=8, labelcolor=INK2, loc='upper left')
fig.tight_layout(); fig.savefig(png, dpi=150, facecolor=SURF)
print(f'\ngraph: {png}')
