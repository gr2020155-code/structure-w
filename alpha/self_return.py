"""Negative correction to the electron anomaly: return of the self-coupling - rules: self_return_rules_2026-11-26.md (R-117)."""
from math import pi, log, exp, floor
alpha = 1 / 137.035999177; kappa = 2 / (3 * pi); E = -2.1037
delta_data = 0.001513
n = 400000

def reach(model):                       # r/R samples by quantiles of p(r) (deterministic midpoints)
    out = []
    for i in range(n):
        u = (i + 0.5) / n                # fraction of pushes reaching beyond r: M1 R/r -> r = 1/u ; M2 on [1,2] with p = 2/r^2
        out.append(1 / u if model == 'M1' else 1 / (1 - u / 2))
    return out

print(f'needed delta = {delta_data} (all orders); kappa*alpha = {kappa*alpha:.6f}')
rows = []
for model in ('M1', 'M2'):
    rs = reach(model)
    # A: alternation, phase phi uniform: P(odd) = average over r of measure of phi with floor(phi + tau) odd
    def p_odd(tau):
        k, f = floor(tau), tau - floor(tau)   # phi in [0,1): phi + tau in [tau, tau+1): first part (1-f) gives floor = k, rest f gives k+1
        return (1 - f) * (k % 2) + f * ((k + 1) % 2)
    P = sum(p_odd((r - 1) / pi) for r in rs) / n
    rows.append((f'A1 alternation, reversed -> 0, {model}', P))
    rows.append((f'A2 alternation, reversed -> -1, {model}', 2 * P))
    for legs in (1, 2):
        for bname, rb in (('sharp D = 2R', 2.0), ('edge D*e^E', 2.0 * exp(E))):
            L = 0.0 if rb <= 1 else sum(log(min(r, rb)) for r in rs) / n
            rows.append((f'B{legs} screening, {bname}, {model}', legs * kappa * alpha * L))
print()
for name, d in rows:
    dev = d / delta_data - 1
    ok = abs(dev) <= 0.05
    print(f'{name:40s}: delta = {d:.6f}  ({dev*100:+8.1f} %)  c2 = {-d*pi/(2*alpha):+.4f}  -> {"PASS" if ok else "FAIL"}')
print('\nc2 = coefficient of (alpha/pi)^2 in a: c2 = -delta*pi/(2*alpha); official second order -0.3285; all orders from data -0.3257')
