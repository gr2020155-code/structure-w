"""Slope via the split of the figure eight (rules: slope_split_rules_2026-11-25.md).
Share of unordered pairs of the n diameter points that O splits between loops, averaged over all positions of O."""
from math import pi
from fractions import Fraction
for n in (10, 100, 1000, 10000):
    pairs = n * (n - 1) // 2
    # O between point k and k+1: loop A has k points, loop B has n-k, k = 1..n-1
    s = sum(Fraction(k * (n - k), pairs) for k in range(1, n))
    mean = s / (n - 1)
    print(f'n = {n:6d}: mean share of split pairs = {float(mean):.6f}  (exact {mean if n <= 100 else "..."})')
print(f'limit 1/3; slope (2/pi)*(1/3) = {2/(3*pi):.6f} = 2/(3pi) = {2/(3*pi):.6f}')
