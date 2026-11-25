"""Quarks as fragments of direction (rules: quark_fragments_rules_2026-11-25.md)."""
from math import comb
from fractions import Fraction as F
data = {0: ('neutrino', F(0), 1, F(0)), 1: ('d-quark', F(1, 3), 3, F(1, 3)), 2: ('u-quark', F(2, 3), 3, F(4, 3)), 3: ('electron', F(1), 1, F(1))}
print(' k  fragment   |Q| rule  data   variants rule  data   screening P1 rule  data   P2 (lego)')
tot1 = tot_d = 0
for k in range(4):
    Q = F(k, 3); n = comb(3, k); s1 = n * Q * Q; s2 = (n * Q) ** 2
    name, Qd, nd, sd = data[k]; tot1 += s1; tot_d += sd
    print(f' {k}  {name:9s}  {str(Q):4s}  {str(Qd):4s}   {n:2d}            {nd:2d}     {str(s1):4s}              {str(sd):4s}   {s2}')
print(f'per generation: P1 {tot1} (data {tot_d})')
print('R = sum over quark flavours of N_c*Q^2: u,d,s ->', 3*(F(4,9)+F(1,9)+F(1,9)), '; +c ->', 3*(F(4,9)*2+F(1,9)*2), '; +b ->', 3*(F(4,9)*2+F(1,9)*3))
