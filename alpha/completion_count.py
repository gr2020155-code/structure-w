"""How many completion forms does the medium give a part? Two readings, all branches (26.11.2026).
Rules: completion_count_rules_2026-11-26.md (pushed before run, commit 5298072). Record: Р-102.
Run: python -I -X utf8 completion_count.py
"""
from fractions import Fraction as F
from math import pi

kappa = 2 / (3 * pi)
types = {3: ['u', 'd', 's'], 4: ['u', 'd', 's', 'c'], 5: ['u', 'd', 's', 'c', 'b']}
axes = {'u': 2, 'd': 1, 's': 1, 'c': 2, 'b': 1}           # k axes of 3 (В-304, В-305)
data = {3: (0.255, 0.018), 5: (0.257, 0.037)}            # Р-101

A_br = {
    'a1 transfer, ordered pairs (6)': lambda n: 6 * F(1, 9),
    'a2 transfer, unordered (3)':     lambda n: 3 * F(1, 9),
    'a3 transfer, 9-1 (8)':           lambda n: 8 * F(1, 9),
    'b1 add missing, one form':       lambda n: sum(3 * F(3 - axes[t], 3) ** 2 for t in types[n]),
    'b2 add missing, axis by axis':   lambda n: sum((3 - axes[t]) * 3 * F(1, 9) for t in types[n]),
}
A_br = {k: A_br[k] for k in sorted(A_br)}
P_br = {
    'p1 as charge':      lambda n: sum(3 * F(axes[t], 3) ** 2 for t in types[n]),
    'p2 missing share':  lambda n: sum(3 * F(3 - axes[t], 3) ** 2 for t in types[n]),
    'p3 equal 1/3':      lambda n: len(types[n]) * F(1, 3),
}
k_br = {'k1 = 1': 1.0, 'k2 = 1/2': 0.5, 'k3 = 1/pi (bridge)': 1 / pi}

total = excl = passed = 0
best = None
for an, A in A_br.items():
    for pn, P in P_br.items():
        d = {n: A(n) - P(n) for n in (3, 4, 5)}
        for kn, k in k_br.items():
            total += 1
            if min(d.values()) <= 0:
                excl += 1
                print(f'{an:32s} {pn:17s} {kn:20s} A-P = {", ".join(str(d[n]) for n in (3,4,5))} -> excluded (b <= 0)')
                continue
            c = {n: k / (kappa * float(d[n])) for n in (3, 5)}
            pulls = {n: (c[n] - data[n][0]) / data[n][1] for n in (3, 5)}
            ok = all(abs(p) <= 2 for p in pulls.values())
            passed += ok
            need = c[3] / data[3][0]
            if best is None or abs(pulls[3]) < abs(best[0]):
                best = (pulls[3], an, pn, kn, c[3])
            print(f'{an:32s} {pn:17s} {kn:20s} A-P = {", ".join(str(d[n]) for n in (3,4,5))}  '
                  f'c(3) = {c[3]:.3f} ({pulls[3]:+.0f} sigma), c(5) = {c[5]:.3f} ({pulls[5]:+.0f} sigma), x{need:.1f} of data -> {"PASS" if ok else "fail"}')
print(f'\nbranches {total}: excluded by sign {excl}, pass {passed}, fail {total - excl - passed}')
print(f'closest: {best[1]} / {best[2]} / {best[3]}: c(3) = {best[4]:.3f} ({best[0]:+.0f} sigma)')
for kn, k in k_br.items():
    print(f'needed A - P at n = 3 for {kn}: {k / (kappa * data[3][0]):.1f}')
