"""Z decay weights and the spin-0 partner of W (26.11.2026).
Rules: z_weight_w_partner_rules_2026-11-26.md (pushed before run). Record: Р-107.
Run: python -I -X utf8 z_weight_w_partner.py
"""
from math import log, sqrt

MW, MZ = 80.3625, 91.1876
c, dc = 0.251, 0.017
L = log(MZ / 0.3128); eps, deps = c / L, dc / L
data = {'inv': (20.000, 0.055), 'l': (3.3658, 0.0023), 'had': (69.911, 0.056), 'Rb': (0.21629, 0.00066), 'Rc': (0.1721, 0.0030)}
kinds = {'nu': (+0.5, 0.0, 1), 'l': (-0.5, -1.0, 1), 'u': (+0.5, 2/3, 3), 'd': (-0.5, -1/3, 3)}

def fractions(s2, e):
    w = {k: n * (T**2 + (T - 2 * Q * s2)**2) for k, (T, Q, n) in kinds.items()}
    up, dn = w['u'] * (1 + e), w['d'] * (1 + e)
    had = 2 * up + 3 * dn
    T = 3 * w['nu'] + 3 * w['l'] + had
    return {'inv': 100 * 3 * w['nu'] / T, 'l': 100 * w['l'] / T, 'had': 100 * had / T, 'Rb': dn / had, 'Rc': up / had}

def report(name, s2):
    r = fractions(s2, eps)
    r2 = fractions(s2, eps + deps)
    out, ok = [], True
    for k in data:
        th = abs(r2[k] - r[k]); s = sqrt(data[k][1]**2 + th**2); p = (r[k] - data[k][0]) / s
        ok &= abs(p) <= 2
        out.append(f'{k} {r[k]:.4f} ({p:+.1f}s)')
    print(f'  {name:26s} s2 = {s2:.4f}: ' + ', '.join(out) + f' -> {"PASS" if ok else "fail"}')

print(f'Part 1. Z weights (eps(MZ) = {eps:.4f} +- {deps:.4f})')
report('ZA from masses', 1 - (MW / MZ)**2)
lo, hi = 0.15, 0.30                                # ZB: fit s2 to invisible fraction (monotonic)
for _ in range(60):
    mid = (lo + hi) / 2
    if fractions(mid, eps)['inv'] > data['inv'][0]: hi = mid
    else: lo = mid
report('ZB fitted to invisible', (lo + hi) / 2)
for name, s2 in (('ZC 1/4', 0.25), ('ZC 2/9', 2 / 9), ('ZC 1/3', 1 / 3)):
    report(name, s2)

print('\nPart 2. Spin-0 partner of W')
pairs = [('rho/pi', 775.3, 139.6), ('K*/K', 891.7, 493.7), ('D*/D', 2010.3, 1869.7),
         ('B*/B', 5324.7, 5279.3), ('J/psi/eta_c', 3096.9, 2983.9), ('Upsilon/eta_b', 9460.3, 9398.7)]
for n, m1, m0 in pairs:
    print(f'  {n:14s} spin 0 lighter: {m0 < m1}; Delta = {m1-m0:6.1f} MeV, Delta*M1 = {(m1-m0)*m1/1e6:.3f} GeV^2')
heavy = pairs[2:]
xs = [log(m1) for _, m1, _ in heavy]; ys = [log(m1 - m0) for _, m1, m0 in heavy]
n = len(xs); mx, my = sum(xs) / n, sum(ys) / n
p = -sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx)**2 for x in xs)
A = my + p * mx
res = [y - (A - p * x) for x, y in zip(xs, ys)]; sres = sqrt(sum(r * r for r in res) / (n - 2))
k1 = [(m1 - m0) * m1 for _, m1, m0 in heavy]
for M, name in ((MW * 1000, 'W'), (MZ * 1000, 'Z')):
    d_free = M ** (-p) * 2.718281828 ** A
    d1 = [k / M for k in k1]
    print(f'  {name}: p free = {p:.2f} -> Delta = {d_free:.2f} MeV (x/÷ {2.718281828**sres:.1f}); p = 1 -> Delta = {min(d1):.1f}-{max(d1):.1f} MeV;'
          f' partner at {M/1000 - max(max(d1), d_free*2.718281828**sres)/1000:.4f}-{M/1000 - min(min(d1), d_free/2.718281828**sres)/1000:.4f} GeV')
print('  LEP: charged spin-0 (tau nu, c s) excluded below 80 GeV (type II), 72.5 GeV (type I)')
