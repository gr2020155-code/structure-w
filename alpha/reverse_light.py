"""Relic as reverse light along the eight (R-125). Rules: reverse_light_rules_2026-11-26.md.
k ~ T^(2/3), units c = T0 = 1, horn torus r = R = cT0. usage: python -P -E -X utf8 reverse_light.py <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, pi, cos

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n')
hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z)
REF = 1403.3; theta_data = 0.010411; FLAT = 702.0

def chi2(dl):
    r = m - 5 * np.log10(dl); M = (one @ Ci @ r) / (one @ Ci @ one)
    return float((r - M) @ Ci @ (r - M))
D = lambda zz: 3 * (1 - (1 + zz) ** -0.5)
relic = [(cs * 3 * (zr1 ** -1.5) ** (1 / 3), D(zr1 - 1)) for zr1 in (1000, 1100, 1400) for cs in (1 / sqrt(6), 1 / sqrt(3))]
def rng(Sfun):
    th = [rs / Sfun(Dr) if Sfun(Dr) > 0 else 0.0 for rs, Dr in relic]
    return min(th), max(th)
inside = lambda t: t[0] <= theta_data <= t[1]

# O1: zero does nothing
t1 = rng(lambda x: x)
print(f'O1 zero does nothing: theta = {t1[0]:.4g} .. {t1[1]:.4g} -> {"PASS" if inside(t1) else "FAIL"}')
print(f'   needed: S/D = {t1[0]/theta_data:.3f} .. {t1[1]/theta_data:.3f} (spreading must be larger than straight)')

# O2: zero focuses (converging thin lens at D_L, power P > 0)
best = 0.0; hit = False
for DL in np.linspace(0.05, 2.9, 60):
    for P in np.linspace(0.01, 20, 400):
        t = rng(lambda x: abs(DL + (1 - P * DL) * (x - DL)))
        hit |= inside(t)
        best = max(best, max(abs(DL + (1 - P * DL) * (Dr - DL)) / Dr for _, Dr in relic))
print(f'O2 zero focuses: max S/D over grid = {best:.3f} (needs >= {t1[0]/theta_data:.3f}) -> {"PASS" if hit else "FAIL"}')

# O3: horn torus surface, light along the meridian (the eight); J'' + K J = 0, K = cos v / (1 + cos v)
K = lambda v: cos(v) / (1 + cos(v))
def jacobi(v0, sgn, smax, h=2e-4):
    """J(s) along v = v0 + sgn*s; returns s grid, J (inf after the path hits 'zero' v = +-pi)."""
    N = int(smax / h) + 1; s = np.linspace(0, smax, N); h = s[1]
    J = np.full(N, np.inf); y, yp = 0.0, 1.0; J[0] = 0.0
    f = lambda si, y: -K(v0 + sgn * si) * y
    for i in range(N - 1):
        si = s[i]
        if pi - abs(v0 + sgn * (si + h)) < 1e-4: break   # reaches zero: spreading diverges
        k1 = (yp, f(si, y)); k2 = (yp + h/2*k1[1], f(si + h/2, y + h/2*k1[0]))
        k3 = (yp + h/2*k2[1], f(si + h/2, y + h/2*k2[0])); k4 = (yp + h*k3[1], f(si + h, y + h*k3[0]))
        y += h/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); yp += h/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        J[i + 1] = y
    return s, J
def Sfun(s, J):
    return lambda x: float(np.interp(x, s, J)) if np.isfinite(np.interp(x, s, J)) else np.inf

# check: axis family J1 = 1 + cos v solves the equation
v = np.linspace(0, 3, 7); print(f'   check J1 = 1+cos v: residual = {max(abs(-np.cos(v) + np.cos(v)/(1+np.cos(v))*(1+np.cos(v)))):.1e}')

print('\nO3 horn torus, r = R = cT0; v0 = our place (0 outer equator, pi = zero)')
print('   v0    | reverse (towards zero) theta    | direct theta                  | SN direct d chi2')
smax = 2.93
res = []
for v0 in np.linspace(0, pi - 0.02, 40):
    sr, Jr = jacobi(v0, +1, smax); sd, Jd = jacobi(v0, -1, smax)
    tr = rng(Sfun(sr, Jr)); td = rng(Sfun(sd, Jd))
    Sd = Sfun(sd, Jd); dl = (1 + z) * np.array([Sd(x) for x in D(z)])
    dx = chi2(dl) - REF if np.all(np.isfinite(dl)) and np.all(dl > 0) else np.inf
    res.append((v0, tr, td, dx))
    print(f'   {v0:5.3f} | {tr[0]:.4g} .. {tr[1]:.4g} {"IN " if inside(tr) else "   "}'.ljust(45) +
          f'| {td[0]:.4g} .. {td[1]:.4g} {"IN " if inside(td) else "   "}'.ljust(32) + f'| {dx:+.1f}')
rev_ok = [(v0, dx) for v0, tr, td, dx in res if inside(tr)]
both = [(v0, dx) for v0, dx in rev_ok if dx <= FLAT]
dir_ok = [v0 for v0, tr, td, dx in res if inside(td)]
print(f'\nreverse relic in range at v0 = {[round(a, 3) for a, _ in rev_ok]}')
print(f'   of these, SN not worse than flat (<= +{FLAT}): {[(round(a, 3), round(b, 1)) for a, b in both]}; threshold 10: '
      f'{[round(a, 3) for a, b in rev_ok if b <= 10]}')
print(f'direct relic in range at v0 = {[round(a, 3) for a in dir_ok]}')
print(f'O3 -> {"PASS" if both else "FAIL"}')

# Addition after the first run (disclosed in R-125): finer grid of the same v0 near the zero, direct light only.
print('\nfiner grid near the zero, direct light (added after the first run; one fitted number v0)')
fine = []
for v0 in np.arange(2.55, 2.905, 0.005):
    sd, Jd = jacobi(v0, -1, smax); Sd = Sfun(sd, Jd)
    td = rng(Sd); dx = chi2((1 + z) * np.array([Sd(x) for x in D(z)])) - REF
    fine.append((v0, td, dx))
i = int(np.argmin([f[2] for f in fine])); vb, tb, db = fine[i]
one_s = [f[0] for f in fine if f[2] - db <= 1]
in_rel = [f[0] for f in fine if inside(f[1])]
print(f'   SN best v0 = {vb:.3f} (1 sigma {min(one_s):.3f} .. {max(one_s):.3f}); SN d chi2 = {db:+.1f} (threshold 10)')
print(f'   relic theta at this v0 = {tb[0]:.4g} .. {tb[1]:.4g} -> {"in range" if inside(tb) else "out of range"}')
print(f'   v0 range the relic alone accepts: {min(in_rel):.3f} .. {max(in_rel):.3f}' if in_rel else '   relic accepts no v0 here')
print(f'   distance from us to the zero at SN best: {pi - vb:.3f} cT0')
print(f'   direct observation -> {"PASS" if db <= 10 and inside(tb) else "FAIL"}')
