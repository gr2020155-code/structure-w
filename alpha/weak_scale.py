"""Weak scale lambda_F: branches (b) number of points and (c) heaviest form (R-144).
Rules: weak_scale_rules_2026-11-27.md (committed before the run, 670a0dc).
usage: python -P -E -X utf8 weak_scale.py"""
import numpy as np
from fractions import Fraction as Fr
from math import log, sqrt

GF, MP, ME = 1.1663788e-5, 1.220890e19, 0.51099895e-3          # GeV^-2, GeV, GeV
MT_DIR, MT_DIR_E = 172.60, 0.27
MT_MS, MT_MS_LO = 162.5, 1.5
MT_POLE, MT_POLE_E = 172.1, 0.6
MW, MZ = 80.3625, 91.1876
MF = GF ** -0.5; V = (sqrt(2) * GF) ** -0.5
N = lambda m: MP / m
Ne, NF = N(ME), N(MF)
print(f'M_F = {MF:.2f} GeV, v = {V:.2f} GeV; N_e = {Ne:.4e}, N_F = {NF:.4e}, ln N_e = {log(Ne):.4f}, ln N_F = {log(NF):.4f}')

# ---- Test B
ps = sorted({Fr(i, j) for j in range(1, 5) for i in range(1, j + 1)})
ratios = sorted({pf / pe for pe in ps for pf in ps if pf / pe < 1})
print(f'\n=== Test B: exponents p = i/j (j <= 4): {len(ps)} values, {len(ratios)} ratios p_F/p_e < 1')
def best(lnT):
    r = min(ratios, key=lambda q: abs(float(q) * log(Ne) - lnT)); return r, float(r) * log(Ne) - lnT
for tag, m, counted in (('M_F (counted)', MF, True), ('v', V, False), ('M_W', MW, False), ('M_Z', MZ, False)):
    lnT = log(N(m)); r, d = best(lnT)
    print(f'  {tag:14s}: ln N/ln N_e = {lnT / log(Ne):.5f}; nearest ratio {r} = {float(r):.5f} -> N model/data = {np.exp(d):.3f} '
          f'({"PASS" if abs(d) <= 0.02 else "FAIL"}{"" if counted else ", not counted"})')
rng = np.random.default_rng(1); T = rng.uniform(0.5, 1.0, 200000) * log(Ne)
rv = np.array([float(q) for q in ratios]) * log(Ne)
hit = np.min(np.abs(T[:, None] - rv[None, :]), axis=1) <= 0.02
print(f'  random targets ln N in [0.5, 1]*ln N_e hitting some ratio within 2 %: {hit.mean() * 100:.1f} %')

# ---- Test C
S = [('1 saturation', 1.0), ('2/3 |Q_t|', 2 / 3), ('1/2 two ends', 0.5), ('4/9 Q_t^2', 4 / 9), ('1/3 one axis', 1 / 3)]
lo, hi = GF * (MT_MS - MT_MS_LO) ** 2, GF * (MT_DIR + MT_DIR_E) ** 2
print(f'\n=== Test C: beta_w = G_F m_t^2: direct {GF * MT_DIR**2:.4f}, pole {GF * MT_POLE**2:.4f}, MS-bar {GF * MT_MS**2:.4f}; window [{lo:.4f}, {hi:.4f}]')
inside = [n for n, s in S if lo <= s <= hi]
for n, s in S:
    mt = MF * sqrt(s)
    print(f'  {n:14s} s = {s:.4f}: {"inside" if lo <= s <= hi else "outside"}; m_t = M_F*sqrt(s) = {mt:.2f} GeV '
          f'(direct {100 * (mt / MT_DIR - 1):+.1f} %, MS-bar {100 * (mt / MT_MS - 1):+.1f} %)')
print(f'  -> {len(inside)} inside: {inside} -> ' + ('PASS (not rejected)' if len(inside) == 1 else 'FAIL' if not inside else 'undecided'))
X = np.exp(rng.uniform(log(0.25), 0.0, 200000)); w = log(hi / lo)
c = np.zeros_like(X, bool)
for _, s in S: c |= (np.log(X) >= log(s) - w / 2) & (np.log(X) <= log(s) + w / 2)
print(f'  window width ln = {w:.4f}; random values (log-uniform on [0.25, 1]) hitting some share: {c.mean() * 100:.1f} %')
if len(inside) == 1:
    s = dict(S)[inside[0]]
    print(f'  consequence: M_F = m_t/sqrt(s) = {MT_DIR / sqrt(s):.1f} GeV (direct), {MT_MS / sqrt(s):.1f} (MS-bar) vs {MF:.1f}; '
          f'lambda_F then needs m_t/m_e = {MT_DIR / ME:.4e} from W (not available)')
