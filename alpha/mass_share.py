"""Mass share = local eights (Ob, from D/H) + pre-material links (Op u^-n) (R-134). Rules: mass_share_rules_2026-11-27.md.
g^2 = g0^2 [Or u^-4 + Ob u^-3 + Op u^-n + OB].
usage: python -P -E -X utf8 mass_share.py <desi_dr2_bao_data.txt> <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, pi
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2dist

bao, snd = sys.argv[1:3]
h = 0.7354                                           # own H0 (R-133)
SIG, C, G, MPC, T0 = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
Or = (4 * SIG / C * T0**4) / (3 * (h * 1e5 / MPC) ** 2 * C * C / (8 * pi * G))

L = open(bao).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cbi = np.linalg.inv(np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]))
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
zs_all = np.array([float(r[hdr.index('zHD')]) for r in data]); ms_all = np.array([float(r[hdr.index('m_b_corr')]) for r in data])
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
sel = zs_all > 0.01; zs, ms = zs_all[sel], ms_all[sel]
Csi = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(zs); w1 = Csi @ one; s11 = one @ w1
REF = 1403.3

def E2(zz, Ob, Op, nn):
    u1 = 1 + zz; return Or * u1**4 + Ob * u1**3 + Op * u1**nn + 1 - Or - Ob - Op
def chi2_bao(g0, Ob, Op, nn):
    if np.any([E2(zz, Ob, Op, nn) <= 0 for zz in (0, 0.5, 1, 2.4)]): return 1e9
    g = lambda zz: g0 * sqrt(E2(zz, Ob, Op, nn)); out = []
    for zz, k in zip(zb, kb):
        dm = quad(lambda t: 1 / g(t), 0, zz)[0]; dh = 1 / g(zz)
        out.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    r = vb - np.array(out); return float(r @ Cbi @ r)
def chi2_sn(Ob, Op, nn):
    grid = np.concatenate([[0], np.sort(zs)]); D = np.zeros_like(grid)
    for j in range(1, grid.size):
        D[j] = D[j - 1] + quad(lambda t: 1 / sqrt(E2(t, Ob, Op, nn)), grid[j - 1], grid[j])[0]
    r = ms - 5 * np.log10((1 + zs) * np.interp(zs, grid, D)); M = (r @ w1) / s11; e = r - M
    return float(e @ Csi @ e) - REF

print(f'h = {h}; Or (photons) = {Or:.3e}')
for btag, obh2 in (('b1 (100 Ob h2 = 2.166, computed rate)', 0.02166), ('b2 (100 Ob h2 = 2.235, measured rate)', 0.02235)):
    Ob = obh2 / h**2
    print(f'\n{btag}: Ob = {Ob:.4f}')
    r = minimize(lambda q: chi2_bao(q[0], Ob, 0.0, 3), [0.034], method='Nelder-Mead')
    sn = chi2_sn(Ob, 0.0, 3)
    print(f'   P0  eights only (Op = 0): BAO chi2 = {r.fun:.1f} / 12, prob {chi2dist.sf(r.fun, 12):.2e}; SN d chi2 = {sn:+.1f} -> '
          f'{"PASS" if chi2dist.sf(r.fun, 12) >= 0.05 and sn <= 11 else "FAIL"}')
    for nn, tag in ((3, 'P3  n = 3 (fixed count, ~1/s)'), (2, 'P2  n = 2'), (1, 'P1  n = 1'), (0, "P0' n = 0 (merges with B)")):
        best = None
        for op0 in (0.05, 0.25, 0.6):
            rr_ = minimize(lambda q: chi2_bao(q[0], Ob, q[1], nn), [0.034, op0], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-6))
            if best is None or rr_.fun < best.fun: best = rr_
        g0, Op = best.x; p = chi2dist.sf(best.fun, 11); sn = chi2_sn(Ob, Op, nn)
        print(f'   {tag}: Op = {Op:.4f}, Op/Ob = {Op / Ob:.2f}; BAO chi2 = {best.fun:.1f} / 11, prob {p:.2e}; SN d chi2 = {sn:+.1f} -> '
              f'{"PASS" if p >= 0.05 and sn <= 11 else "FAIL"}')
    best = None
    for st in ([0.034, 0.25, 3.0], [0.034, 0.4, 2.5], [0.034, 0.2, 3.3]):
        rr_ = minimize(lambda q: chi2_bao(q[0], Ob, q[1], q[2]), st, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-6, maxiter=4000))
        if best is None or rr_.fun < best.fun: best = rr_
    g0, Op, nn = best.x
    print(f'   observation, n free: n = {nn:.2f}, Op = {Op:.4f}, Op/Ob = {Op / Ob:.2f}; BAO chi2 = {best.fun:.1f} / 10; SN d chi2 = {chi2_sn(Ob, Op, nn):+.1f}')
