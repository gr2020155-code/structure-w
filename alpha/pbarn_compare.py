"""Comparison of the fixed prediction R-149 with p-bar n data at rest (Klempt, Batty, Richard, Phys. Rep. 413 (2005) 197,
hep-ex/0501020, Table 10; data of Bettini et al. 1967, liquid D2, spectator cut 250 MeV/c). No parameter is changed.
usage: python -E -X utf8 pbarn_compare.py"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pbarn as P
w = {(k, j): s for s, (k, j, _) in zip(P.S, P.SETS)}
prong = lambda c: sum(v for (k, j), v in w.items() if 2 * k + 1 == c)
DATA = {1: (16.4, 0.5), 3: (59.7, 1.2), 5: (23.4, 0.7), 7: (0.39, 0.07)}
tot = sum(v[0] for v in DATA.values())
mean = sum(c * v[0] for c, v in DATA.items()) / tot
print('\n=== comparison with data (Table 10) ===')
print('charged prongs, %: ' + ', '.join(f'{c}: W {100 * prong(c):.1f} / data {v[0]} +- {v[1]}' for c, v in DATA.items()))
print(f'<n+-> from prongs: data {mean:.3f}  (<n+> {(mean - 1) / 2:.3f}, <n-> {(mean + 1) / 2:.3f});  W {P.npl + P.nmi:.3f}')
EXC = [('pi- pi0', 0, 1, 0.40, 0.04), ('pi- 2pi0', 0, 2, 0.68, 0.07), ('pi- 4pi0', 0, 4, 1.32, 0.20),
       ('2pi- pi+', 1, 0, 1.57, 0.21), ('2pi- pi+ pi0', 1, 1, 21.8, 2.2), ('2pi- pi+ 2pi0', 1, 2, 7.0, 1.1),
       ('3pi- 2pi+', 2, 0, 5.15, 0.47), ('3pi- 2pi+ pi0', 2, 1, 15.1, 1.0)]
print('exclusive channels, % (W / data):')
for n, k, j, d, e in EXC:
    print(f'  {n:16s} W {100 * w.get((k, j), 0):5.2f}  data {d:5.2f} +- {e:.2f}')
print(f'  pi- 3pi0         W {100 * w.get((0, 3), 0):5.2f}  data: not listed')
