"""How the returned push combines with the orbit (R-119). Rules: orbit_sum_rules_2026-11-26.md.
a_e = a0 (1 - c a0), a0 = alpha/(2 pi); branches give c."""
from math import pi, sqrt, log
alpha = 1 / 137.035999177; a0 = alpha / (2 * pi); a_e = 0.00115965218062
c_need = (1 - a_e / a0) / a0
print(f'a0 = {a0:.10f}; a_e = {a_e}; needed c = {c_need:.4f}  (delta = {c_need*a0:.6f})\n')

def c_of(a): return (1 - a / a0) / a0             # exact c from a branch's a
rows = [('K1 linear', c_of(a0)),
        ('K2 repeat, replaces (A-64, A-69)', c_of(a0 / (1 + a0))),
        ("K2' repeat, adds", c_of(a0 * (1 + a0))),
        ('K3 by squares (A-39)', c_of(sqrt(1 + 2 * a0) - 1)),
        ('K4 repeat at both ends (A-97)', c_of(a0 / (1 + a0)**2))]
for name, c in rows:
    dev = c / c_need - 1
    print(f'{name:34s}: c = {c:+.4f}  ({dev*100:+7.1f} %)  -> {"PASS" if abs(dev) <= 0.05 else "FAIL"}')

# K5: act absent from the orbit while away; tau = (r/R - 1)/pi turns
n = 200000
tau_M2 = sum(((1 / (1 - ((i + 0.5) / n) / 2)) - 1) / pi for i in range(n)) / n      # p = 2/r^2 on [1,2]
print(f'\nK5 conservation (A-28): M1 <tau> = infinite (p = 1/r^2, <r> diverges) -> excluded')
d5 = 2 * pi * tau_M2
print(f'K5 conservation (A-28), M2: <tau> = {tau_M2:.4f} turn (exact (2 ln2 - 1)/pi = {(2*log(2)-1)/pi:.4f}); delta = 2 pi <tau> = {d5:.4f}'
      f' vs allowed <= {0.05*c_need*a0:.6f} -> {"excluded" if d5 > 0.05 * c_need * a0 else "allowed"}')
