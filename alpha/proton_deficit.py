"""Nucleon moment deficit: branches (a) spin share in the closure, (b) motion of parts.
Rules: proton_deficit_rules_2026-11-26.md (pushed as cf5442a before running)."""
from fractions import Fraction as F

m_p, m_n = 938.272, 939.565
mu_p, mu_n = 2.79284734, -1.91304273
u, d = F(2, 3), F(-1, 3)
# parts: (charge, weight V-312, role, position); mass m/3 -> unit 3
P = [(u, F(2, 3), 'pair', 'end'), (d, F(-1, 3), 'third', 'mid'), (u, F(2, 3), 'pair', 'end')]
N = [(d, F(2, 3), 'pair', 'end'), (d, F(2, 3), 'pair', 'mid'), (u, F(-1, 3), 'third', 'end')]
BR = {'A1/B1 all parts': lambda r, p: True,
      'A2 third part':   lambda r, p: r == 'third',
      'A3 pair':         lambda r, p: r == 'pair',
      'A4/B3 middle':    lambda r, p: p == 'mid',
      'B2 end loops':    lambda r, p: p == 'end'}

def mu(parts, sel, x, scale):
    return sum(3 * float(w * q) * (x if sel(r, p) else 1) for q, w, r, p in parts) * scale

def solve(sel):                                   # mu_p linear in x
    a, b = mu(P, sel, 0, 1), mu(P, sel, 1, 1)
    return (mu_p - a) / (b - a)

def eta(parts, sel, x):                           # share of the whole's spin taken from the selected parts
    return float(sum(w for q, w, r, p in parts if sel(r, p))) * (1 - x)

print(f'{"branch":17} {"factor x":>8} {"eta (p)":>8} {"mu_n":>8} {"err":>7} {"n/p":>7}  pass (+-2 %)')
passed = []
for name, sel in BR.items():
    x = solve(sel)
    n = mu(N, sel, x, m_p / m_n)
    e = n / mu_n - 1
    ok = abs(e) <= .02
    if ok: passed.append((name, sel))
    print(f'{name:17} {x:8.4f} {eta(P, sel, x):8.4f} {n:8.3f} {e:+7.1%} {n / mu_p:7.3f}  {"YES" if ok else "no"}')
print(f'{"data":17} {"":8} {"":8} {mu_n:8.3f} {"":7} {mu_n / mu_p:7.3f}')
print()
print('values without fitting, for passed branches (both within +-2 %):')
for name, sel in passed:
    w_sel = float(sum(w for q, w, r, p in P if sel(r, p)))
    for label, x in (('f = 8/9', 8 / 9), ('eta = 1/4', 1 - 0.25 / w_sel), ('eta = 1/3', 1 - (1 / 3) / w_sel)):
        p_, n_ = mu(P, sel, x, 1), mu(N, sel, x, m_p / m_n)
        ok = abs(p_ / mu_p - 1) <= .02 and abs(n_ / mu_n - 1) <= .02
        print(f'  {name:15} {label:10} x = {x:.4f}: mu_p {p_:.3f} ({p_ / mu_p - 1:+.1%}), mu_n {n_:.3f} ({n_ / mu_n - 1:+.1%})  {"YES" if ok else "no"}')
print()
print('observation (seen before the rules, NOT counted): eta = 1/9 for passed branches')
for name, sel in passed:
    w_sel = float(sum(w for q, w, r, p in P if sel(r, p)))
    x = 1 - (1 / 9) / w_sel
    p_, n_ = mu(P, sel, x, 1), mu(N, sel, x, m_p / m_n)
    print(f'  {name:15} x = {x:.4f}: mu_p {p_:.3f} ({p_ / mu_p - 1:+.1%}), mu_n {n_:.3f} ({n_ / mu_n - 1:+.1%})')
