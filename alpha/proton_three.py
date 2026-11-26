"""Why mu_p = 3. Rules: proton_three_rules_2026-11-26.md (pushed as 1c301a4 before running)."""
from fractions import Fraction as F

m_p, m_n = 938.272, 939.565
mu_p, mu_n = 2.79284734, -1.91304273
u, d = F(2, 3), F(-1, 3)
W_PAIR, W_THIRD = F(2, 3), F(-1, 3)              # V-312

print('step 1: decomposition of the data')
print(f'  mu_p + mu_n: model 1 (any weights summing to 1, equal masses), data {mu_p + mu_n:.4f} -> factor {mu_p + mu_n:.3f}')
print(f'  mu_p - mu_n: model 5,                                          data {mu_p - mu_n:.4f} -> factor {(mu_p - mu_n) / 5:.3f}')
print()
# mass share of the loops along the line: (end, middle, end)
K = {'K0 loops equal': (F(1, 3), F(1, 3), F(1, 3)),
     'K1 by points':   (F(3, 8), F(1, 4), F(3, 8)),
     'K2 by eights':   (F(1, 4), F(1, 2), F(1, 4))}
print(f'{"reading":16} {"mu_p":>7} {"err":>7} {"mu_n":>7} {"err":>7} {"sum":>7} {"diff":>7}  pass (+-2 %)')
for name, (e1, mid, e2) in K.items():
    # proton: end u (pair), middle d (third), end u (pair)
    p = float(W_PAIR * u / e1 + W_THIRD * d / mid + W_PAIR * u / e2)
    # neutron: end d (merged, pair), middle d (pair), end u (third); masses are shares of m_n
    n = float(W_PAIR * d / e1 + W_PAIR * d / mid + W_THIRD * u / e2) * m_p / m_n
    ep, en = p / mu_p - 1, n / mu_n - 1
    ok = abs(ep) <= .02 and abs(en) <= .02
    print(f'{name:16} {p:7.3f} {ep:+7.1%} {n:7.3f} {en:+7.1%} {p + n:7.3f} {p - n:7.3f}  {"YES" if ok else "no"}')
print(f'{"data":16} {mu_p:7.3f} {"":7} {mu_n:7.3f} {"":7} {mu_p + mu_n:7.3f} {mu_p - mu_n:7.3f}')
