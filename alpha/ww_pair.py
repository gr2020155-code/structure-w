"""W+W- pair production: induced form (A) vs exchange of parts (B) - rules: ww_pair_rules_2026-11-26.md (R-120).
Run: python -P -E -X utf8 ww_pair.py lep_ww_data.txt"""
import sys
from math import sqrt
MW = 80.37
lines = [l.split() for l in open(sys.argv[1], encoding='utf-8') if l.strip() and not l.startswith('#')]
i = lines.index(['TOTAL']) + 1; total = []
while lines[i] != ['ANGULAR']:
    total.append(tuple(map(float, lines[i]))); i += 1
i += 1; ang = []
while i < len(lines):
    E = float(lines[i][0]); v = list(map(float, lines[i + 1])); e = list(map(float, lines[i + 2])); ang.append((E, v, e)); i += 3

print('1. total: sigma_WW / sigma_mumu(point), sigma_mumu = 86.8 nb / s')
for rs, s, ds in total:
    smm = 86.8e3 / rs**2                  # pb
    print(f'   {rs:6.1f} GeV: sigma_WW = {s:5.2f} +- {ds:.2f} pb, sigma_mumu = {smm:.3f} pb, ratio = {s/smm:5.2f}')

shapes = {
    'A0 induced form, flat': lambda c, b: 1.0,
    'A1 induced form, 1+cos^2': lambda c, b: 1 + c * c,
    'B1 exchange, potential 1/r -> 1/q^2': lambda c, b: 1 / (1 + b * b - 2 * b * c)**2,
    'B2 exchange, field 1/r^2 -> 1/q': lambda c, b: 1 / (1 + b * b - 2 * b * c),
}
def binavg(f, b, k, n=400):
    lo = -1 + 0.2 * k
    return sum(f(lo + 0.2 * (j + 0.5) / n, b) for j in range(n)) / n

print('\n2. data: forward-backward asymmetry of W- (cos > 0 vs < 0), and last/first bin')
for E, v, e in ang:
    F, B = sum(v[5:]), sum(v[:5]); dF2 = sum(x * x for x in e[5:]); dB2 = sum(x * x for x in e[:5])
    A = (F - B) / (F + B); dA = 2 * sqrt(B * B * dF2 + F * F * dB2) / (F + B)**2
    print(f'   {E:6.2f} GeV: A_FB = {A:.3f} +- {dA:.3f} ({A/dA:.0f} sigma from 0); bin10/bin1 = {v[9]/v[0]:.1f}')

print('\n3. shape chi2 (normalisation free in each range, 9 dof each, 36 total; pass <= 58.6)')
for name, f in shapes.items():
    tot = 0.0; parts = []
    for E, v, e in ang:
        b = sqrt(1 - (2 * MW / E)**2)
        m = [binavg(f, b, k) for k in range(10)]
        w = [1 / x**2 for x in e]
        norm = sum(vi * mi * wi for vi, mi, wi in zip(v, m, w)) / sum(mi * mi * wi for mi, wi in zip(m, w))
        chi = sum(((vi - norm * mi) / ei)**2 for vi, mi, ei in zip(v, m, e))
        F, B = sum(m[5:]), sum(m[:5])
        parts.append(f'{E:.0f}: chi2 {chi:6.1f}, A_FB {(F-B)/(F+B):.3f}, b10/b1 {m[9]/m[0]:.1f}')
        tot += chi
    print(f'   {name:38s}: chi2 = {tot:7.1f} / 36 -> {"PASS" if tot <= 58.6 else "FAIL"}')
    for p in parts: print('        ' + p)
