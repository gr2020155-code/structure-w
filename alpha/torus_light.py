"""Light spreading with curvature from the horn-torus shape (R-124). Rules: torus_light_rules_2026-11-26.md.
k ~ T^(2/3), units c = T0 = 1. usage: python -P -E -X utf8 torus_light.py <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n')
hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
iz, im = hdr.index('zHD'), hdr.index('m_b_corr')
z_all = np.array([float(r[iz]) for r in data]); m_all = np.array([float(r[im]) for r in data])
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z)
REF = 1403.3; theta_data = 0.010411

def chi2(dl):
    r = m - 5 * np.log10(dl); M = (one @ Ci @ r) / (one @ Ci @ one)
    return float((r - M) @ Ci @ (r - M))
D = lambda zz: 3 * (1 - (1 + zz) ** -0.5)                  # int (1+z) dT with 1+z = T^(-2/3)
def S(Dv, geo, rho):
    if geo == '0': return Dv
    return rho * (np.sin(Dv / rho) if geo == '+' else np.sinh(Dv / rho))
def theta_rng(geo, rho):
    out = []
    for zr1 in (1000, 1100, 1400):
        x = zr1 ** -1.5                                     # T_r
        for cs in (1 / sqrt(6), 1 / sqrt(3)):
            rs = cs * 3 * x ** (1 / 3)                      # int_0^x T^(-2/3) dT
            Dr = D(zr1 - 1)
            Sv = S(Dr, geo, rho)
            out.append(rs / Sv if Sv > 0 else float('inf'))
    return min(out), max(out)

print(f'reference chi2 = {REF}; relic theta* = {theta_data}\n')
for geo, rho, tag in (('0', None, 'flat'), ('+', 1.0, 'rho = cT0'), ('+', 2.0, 'rho = 2cT0'), ('-', 1.0, 'rho = cT0'), ('-', 2.0, 'rho = 2cT0')):
    if geo == '+' and D(1400 - 1) / rho >= np.pi:
        trng = (float('nan'), float('nan'))
    else:
        trng = theta_rng(geo, rho)
    dx = chi2((1 + z) * S(D(z), geo, rho)) - REF
    ok = dx <= 11 and trng[0] <= theta_data <= trng[1]
    print(f'G{geo} {tag:11s}: SN d chi2 = {dx:+7.1f}; relic theta = {trng[0]:.4g} .. {trng[1]:.4g} -> {"PASS" if ok else "FAIL"}')

print('\nobservation (one fitted number): G- with free rho fitted on SN, then relic predicted')
grid = np.exp(np.linspace(np.log(0.1), np.log(20), 2000))
vals = [chi2((1 + z) * S(D(z), '-', r)) for r in grid]
i = int(np.argmin(vals)); best = grid[i]; dx = vals[i] - REF
okr = grid[np.array(vals) - vals[i] <= 1]
tr = theta_rng('-', best)
print(f'   best rho = {best:.3f} cT0 (1 sigma {okr.min():.3f} .. {okr.max():.3f}); SN d chi2 = {dx:+.1f} (threshold 10)')
print(f'   relic theta at this rho = {tr[0]:.4g} .. {tr[1]:.4g} (data {theta_data}) -> {"PASS" if dx <= 10 and tr[0] <= theta_data <= tr[1] else "FAIL"}')
need = [r for r in grid if theta_rng('-', r)[0] <= theta_data <= theta_rng('-', r)[1]]
print(f'   rho range that the relic alone accepts: {min(need):.2f} .. {max(need):.2f} cT0' if need else '   no rho fits the relic')
