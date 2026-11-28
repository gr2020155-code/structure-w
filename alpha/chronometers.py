"""Rate of link growth from cosmic chronometers: g(z) = -1/(1+z) dz/dT, ages in ticks directly (R-158).
Rules: chronometers_rules_2026-11-28.md (committed before the run, 649437d). Data: cc_moresco_data.txt (D-95).
usage: python -P -E -X utf8 chronometers.py <cc_moresco_data.txt>"""
import sys
import numpy as np
from math import sqrt, log, exp
from scipy.integrate import quad
from scipy.stats import chi2 as chi2dist

blocks, cur = {}, None
for line in open(sys.argv[1], encoding='utf-8'):
    s = line.strip()
    if not s or s.startswith('#'):
        continue
    key = s.split()[0]
    if key in ('BC03', 'M11', 'SYS'):
        cur = key; blocks[cur] = []; continue
    blocks[cur].append(s)
bc = np.array([[float(x) for x in l.split(',')[:3]] for l in blocks['BC03']])
m11 = np.array([[float(x) for x in l.split(',')[:3]] for l in blocks['M11']])
sy = np.array([[float(x) for x in l.split()] for l in blocks['SYS']])

Or, OmL = 4.573e-5, 0.2973
REFS = [('Cepheids 0.7354', 0.7354, 0.011), ('red giants 0.710', 0.710, 0.011), ('far rate 0.685', 0.685, 0.010)]


def E(z, Om):
    return np.sqrt(Or * (1 + z) ** 4 + Om * (1 + z) ** 3 + 1 - Or - Om)


def cov(d, sps_col=4, stlib=False, diag_only=False):
    z, H, e = d.T
    C = np.diag(e ** 2)
    if diag_only:
        return C
    cols = [1, sps_col] + ([2] if stlib else [])
    for c in cols:
        eta = np.interp(z, sy[:, 0], sy[:, c]) / 100
        C = C + np.outer(H * eta, H * eta)
    return C


def fit_h(d, C, Om):
    z, H, _ = d.T; Ci = np.linalg.inv(C); m = 100 * E(z, Om)
    F = m @ Ci @ m; h = (m @ Ci @ H) / F; r = H - h * m
    return h, 1 / sqrt(F), float(r @ Ci @ r)


def T0g0(Om):
    return quad(lambda x: 1 / E(exp(x) - 1, Om), 0, log(1e12), limit=500)[0]


def report(name, d, C, main=False):
    h, sh, c2 = fit_h(d, C, OmL); p = chi2dist.sf(c2, len(d) - 1)
    dist = [(lab, (h - v) / sqrt(sh ** 2 + s ** 2)) for lab, v, s in REFS]
    print(f'{name}: h = {h:.4f} +- {sh:.4f}; chi2 = {c2:.2f} / {len(d) - 1}, probability {p:.3f}'
          + (f' -> shape {"PASS" if p >= 0.05 else "FAIL"}' if main else ''))
    print('   distance: ' + '; '.join(f'{lab}: {s:+.2f} sigma' for lab, s in dist))
    if main:
        sC, sF = abs(dist[0][1]), abs(dist[2][1])
        dec = 'far rate' if (sC > 2 and sF < 1) else 'Cepheids' if (sF > 2 and sC < 1) else 'none'
        print(f'   clocks decide: {"no" if dec == "none" else "yes, for " + dec}')
    return h


print(f'Om = {OmL} (shape from the BAO ruler, R-156), Or = {Or:.3e}; only h free\n')
hm = report('MAIN  BC03, cov = diag + SPS(mod_ooo) + IMF', bc, cov(bc), main=True)
print('\ncomparison (not counted):')
report('M11, same cov', m11, cov(m11))
report('BC03, SPS full (mod)', bc, cov(bc, sps_col=3))
report('BC03, + stellar library', bc, cov(bc, stlib=True))
report('BC03, diagonal only', bc, cov(bc, diag_only=True))
grid = np.linspace(0.05, 0.80, 301)
res = [(fit_h(bc, cov(bc), om)[2], om) for om in grid]; c2b, omb = min(res)
hb, shb, _ = fit_h(bc, cov(bc), omb)
ok = [om for c2, om in res if c2 - c2b <= 1]
print(f'BC03 free Om and h: Om = {omb:.3f} (1 sigma {min(ok):.2f}-{max(ok):.2f}), h = {hb:.4f} +- {shb:.4f} (at that Om), chi2 = {c2b:.2f} / 13')
tg = T0g0(OmL); GYR_H = 977.79222 / 100      # 1/(100 km/s/Mpc) in Gyr
print(f'age: T0*g0 = {tg:.4f}; at main h = {hm:.4f}: {tg * GYR_H / hm:.2f} Gyr (oldest stars ~13.5)')
