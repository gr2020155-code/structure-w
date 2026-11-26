"""Light in the volume of the horn torus (R-126). Rules: torus_volume_rules_2026-11-26.md.
Space T1: ds^2 = dsigma^2 + rho^2 dOmega^2, rho = r(1+cos v), r = R = cT0 = 1. k ~ T^(2/3), D2.
usage: python -P -E -X utf8 torus_volume.py <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, pi

d = sys.argv[1]
rows = open(f'{d}/PantheonPlusSH0ES.dat').read().split('\n')
hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda name: np.array([float(r[hdr.index(name)]) for r in data])
z_all, m_all, ra_all, de_all = col('zHD'), col('m_b_corr'), np.radians(col('RA')), np.radians(col('DEC'))
c = np.loadtxt(f'{d}/STATSYS.cov'); n = int(c[0]); C_all = c[1:].reshape(n, n)
sel = z_all > 0.01; z, m = z_all[sel], m_all[sel]
nhat = np.stack([np.cos(de_all) * np.cos(ra_all), np.cos(de_all) * np.sin(ra_all), np.sin(de_all)], 1)[sel]
Ci = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(z); w1 = Ci @ one; s11 = one @ w1
REF = 1403.3; THR = 8.0; theta_data = 0.010411
D = lambda zz: 3 * (1 - (1 + zz) ** -0.5)
Dsn = D(z)
relic = [(cs * 3 * (zr1 ** -1.5) ** (1 / 3), D(zr1 - 1)) for zr1 in (1000, 1100, 1400) for cs in (1 / sqrt(6), 1 / sqrt(3))]

def chi2_many(R):                       # rows of residuals r = m - 5 log10 dl; M free
    M = (R @ w1) / s11; X = R - M[:, None]
    return np.einsum('ij,ij->i', X @ Ci, X)

B = np.linspace(0, pi, 73)              # beta0: 0 = towards the zero, pi = away
H = 5e-4; SMAX = 2.93; NS = int(SMAX / H); EVERY = 10
sgrid = np.arange(0, NS + 1, EVERY) * H
def rays(v0):
    """S(s, beta0) = sqrt|J1 J2| on sgrid; inf once the ray reaches the zero (rho -> 0)."""
    y = np.zeros((6, B.size)); y[0] = v0; y[1] = B; y[3] = 1.0; y[5] = 1.0   # v, beta, J1, J1', J2, J2'
    dead = np.zeros(B.size, bool); out = np.zeros((sgrid.size, B.size)); clairaut0 = (1 + np.cos(v0)) * np.sin(B)
    def f(y):
        v, b, j1, p1, j2, p2 = y
        cv, sv, sb, cb = np.cos(v), np.sin(v), np.sin(b), np.cos(b)
        rho = np.maximum(1 + cv, 1e-12)
        K1 = cv / rho; K2 = cv * cv / rho**2
        return np.array([cb, sv * sb / rho, p1, -K1 * j1, p2, -(cb * cb * K1 + sb * sb * K2) * j2])
    for i in range(NS + 1):
        if i % EVERY == 0:
            out[i // EVERY] = np.where(dead, np.inf, np.sqrt(np.abs(y[2] * y[4])))
        if i == NS: break
        k1 = f(y); k2 = f(y + H / 2 * k1); k3 = f(y + H / 2 * k2); k4 = f(y + H * k3)
        y = y + H / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        dead |= (1 + np.cos(y[0])) < 1e-5
    drift = np.max(np.abs(np.where(dead, 0, (1 + np.cos(y[0])) * np.sin(y[1]) - clairaut0)))
    return out, drift

# axes on the sky (Fibonacci sphere)
NA = 600; k = np.arange(NA) + 0.5; ph = np.arccos(1 - 2 * k / NA); th = pi * (1 + 5**0.5) * k
AX = np.stack([np.sin(ph) * np.cos(th), np.sin(ph) * np.sin(th), np.cos(ph)], 1)
BETA = np.arccos(np.clip(AX @ nhat.T, -1, 1))            # (NA, Nsn)
bi = np.clip(np.searchsorted(B, BETA) - 1, 0, B.size - 2); bt = (BETA - B[bi]) / (B[1] - B[0])

def relic_theta(Ssb):
    """sky-average theta (weight sin beta) for each (cs, z_r), and max/min over directions."""
    avg, ratio = [], []
    for rs, Dr in relic:
        Sb = np.array([np.interp(Dr, sgrid, Ssb[:, j]) if np.all(np.isfinite(Ssb[:, j][sgrid <= Dr + H])) else np.inf for j in range(B.size)])
        tb = np.where(np.isfinite(Sb) & (Sb > 0), rs / np.maximum(Sb, 1e-300), 0.0)
        wts = np.sin(B); avg.append(np.trapezoid(tb * wts, B) / np.trapezoid(wts, B))
        pos = tb[tb > 0]; ratio.append(pos.max() / pos.min() if pos.size and np.all(tb > 0) else np.inf)
    return min(avg), max(avg), max(ratio)

print(f'SN reference chi2 = {REF}; threshold d chi2 <= {THR} (3 fitted numbers); relic theta* = {theta_data}')
print('   v0    | SN best d chi2 (over axis) | relic sky-avg theta       | relic max/min over sky | Clairaut drift')
V0 = np.linspace(0, pi - 0.02, 64)
res = []
for v0 in V0:
    Ssb, drift = rays(v0)
    # S at each SN distance for each beta0, then interpolate in beta for each axis
    Tsn = np.array([np.interp(Dsn, sgrid, np.where(np.isfinite(Ssb[:, j]), Ssb[:, j], 1e30)) for j in range(B.size)])   # (NB, Nsn)
    Tsn = np.where(Tsn > 1e29, np.inf, Tsn)
    idx = np.arange(z.size)[None, :]
    S_ax = Tsn[bi, idx] * (1 - bt) + Tsn[bi + 1, idx] * bt
    with np.errstate(divide='ignore', invalid='ignore'):
        R = m[None, :] - 5 * np.log10(np.maximum((1 + z)[None, :] * S_ax, 1e-30))
    bad = ~np.all(np.isfinite(R), 1)
    R[bad] = 0.0
    ch = chi2_many(R) - REF; ch[bad] = np.inf
    ia = int(np.argmin(ch)); lo, hi, rat = relic_theta(Ssb)
    res.append((v0, ch[ia], AX[ia], lo, hi, rat))
    print(f'   {v0:5.3f} | {ch[ia]:+10.1f}'.ljust(38) + f'| {lo:.4g} .. {hi:.4g} {"IN" if lo <= theta_data <= hi else "  "}'.ljust(28)
          + f'| {rat:8.3g}'.ljust(25) + f'| {drift:.1e}')

i = int(np.argmin([r[1] for r in res])); v0b, chb, axb, lob, hib, ratb = res[i]
dec = np.degrees(np.arcsin(axb[2])); ra = np.degrees(np.arctan2(axb[1], axb[0])) % 360
print(f'\nSN best: v0 = {v0b:.3f} (distance to the zero {pi - v0b:.3f} cT0), axis to the zero RA = {ra:.1f}, DEC = {dec:.1f};'
      f' d chi2 = {chb:+.1f} (threshold {THR})')
print(f'   relic at this place: sky-avg theta = {lob:.4g} .. {hib:.4g} (data {theta_data}) -> {"in range" if lob <= theta_data <= hib else "out of range"};'
      f' max/min over sky = {ratb:.3g}')
print(f'   v0 within d chi2 + 1 of best: {[round(r[0], 3) for r in res if r[1] - chb <= 1]}')
print(f'   v0 the relic alone accepts: {[round(r[0], 3) for r in res if r[3] <= theta_data <= r[4]]}')
print(f'T1 -> {"PASS" if chb <= THR and lob <= theta_data <= hib else "FAIL"}')
