"""Is there an excess in Bhabha scattering? Effective fraction f (26.11.2026).
Rules: bhabha_excess_rules_2026-11-26.md (pushed before run). Record: Р-112.
Run: python -I -X utf8 bhabha_excess.py
"""
from math import cos, sin, pi, sqrt, acos

a = 1 / 137.036
d = a - 3 * a / (4 * pi)

def parts(t):
    s_ = 1 + cos(t)**2
    tt = 2 * (1 + cos(t / 2)**4) / sin(t / 2)**4
    it = -4 * cos(t / 2)**4 / sin(t / 2)**2
    return s_, tt, it

t0 = pi / 2
s_, tt, it = parts(t0)
print(f'check at 90 deg: s {s_:.3f}, t {tt:.3f}, int {it:.3f}, sum {s_+tt+it:.3f} vs full {((3+cos(t0)**2)/(1-cos(t0)))**2:.3f}')

def f_eff(lo, hi, N=200000):
    S = T = I = 0.0
    h = (hi - lo) / N
    for i in range(N):
        t = lo + (i + 0.5) * h
        s_, tt, it = parts(t); w = sin(t) * h
        S += s_ * w; T += tt * w; I += it * w
    return (S + I / 2) / (S + T + I), S / (S + T + I), I / (S + T + I)

data = {'CMD-3': ((1.0, pi - 1.0), 1.0017, 0.0063), 'BaBar': ((acos(0.7), acos(-0.7)), 1.0040, 0.0111)}
for name, ((lo, hi), v, s) in data.items():
    f, fs, fi = f_eff(lo, hi)
    pred = 1 + d * (1 - f)
    print(f'{name}: theta {lo:.3f}-{hi:.3f} rad; share s {fs:.4f}, interference {fi:+.4f}, f = {f:+.4f}; '
          f'W = {pred:.4f} ({(pred - v) / s:+.2f} sigma), QED 1.0000 ({(1 - v) / s:+.2f} sigma)')
