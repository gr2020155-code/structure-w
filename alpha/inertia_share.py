"""How is inertia shared between the core and the global eight? All branches (26.11.2026).
Rules: inertia_share_rules_2026-11-26.md (pushed before run). Record: Р-105.
Run: python -I -X utf8 inertia_share.py
"""
from fractions import Fraction as F
from itertools import product

N = 10**22   # tacts per orbit, order of magnitude (B-245 p.4); only matters that N >> 1
share = {
    'S1 per participant':   lambda ng: F(1, 1 + ng),
    'S2 per exchange end':  lambda ng: F(1, 2),
    'S3 per point (3 vs 8)': lambda ng: F(3, 3 + 8 * ng),
    'S4 per tact':          lambda ng: F(N, N + 1),
}
passed = []
for name, f in share.items():
    v = {ng: f(ng) for ng in (1, 3)}
    ok = all(x == F(1, 2) for x in v.values())
    passed += [name] if ok else []
    show = ', '.join(f'n_g={ng}: f={float(x):.4g}' for ng, x in v.items())
    print(f'{name:24s} {show:40s} -> {"PASS" if ok else "fail"}')
print(f'pass: {passed}\n')

def spins(units):            # sum of units each 1/2 (orbit share), all couplings
    tot = set()
    for signs in product((F(1, 2), F(-1, 2)), repeat=units):
        tot.add(abs(sum(signs)))
    return sorted(tot)

def show(s):
    return '/'.join(str(x) for x in s)

# (name, units under S2, spin under b2 (B-27), data spin (QM reading of angular distributions except e, gamma))
table = [
    ('e, mu, tau (one eight)',          [1],  '1/2',     '1/2'),
    ('nu (eight at a point, B-65)',     [1],  '1/2',     '1/2'),
    ('p, n (three parts, B-312)',       [3],  '1/2/3/2', '1/2'),
    ('Omega- (three parts)',            [3],  '1/2/3/2', '3/2'),
    ('pi (two parts)',                  [2],  '0/1',     '0'),
    ('rho (two parts)',                 [2],  '0/1',     '1'),
    ('W, Z read 1: one loop with own inertia (B-294 p.1)', ['loop'], '1', '1'),
    ('W, Z read 2: two parts (k, 3-k) (B-305 p.3)',        [2],      '0/1', '1'),
    ('gamma (one loop, no inertia)',    ['photon'], '1', '1'),
]
for name, u, b2, data in table:
    if u == ['loop']:
        s2 = '1/2/3/2 (loop 1 + orbit 1/2)'
        ok = 'NO'
    elif u == ['photon']:
        s2 = '1 (no lag, no orbit share)'
        ok = 'ok'
    else:
        s = spins(u[0]); s2 = show(s)
        ok = 'ok' if data in [str(x) for x in s] else 'NO'
    print(f'{name:52s} S2: {s2:30s} b2: {b2:9s} data: {data:4s} S2 {ok}')
