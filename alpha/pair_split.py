"""Edge as a pair of loops (rules: pair_split_rules_2026-11-25.md).
w(rho) = sum over splits x of weight(x) * exp(-rho / sqrt(x(1-x))), rho = r/R, R = hbar/(mc), D = 2R."""
from math import exp, log, sqrt
src = open('electron_edge.py', encoding='utf-8').read().split("print('Data side")[0]
exec(src)                                   # brings uehling(rho, n), GAMMA

def mix(rho, weight, n=4000):
    h = 1 / n; s = 0; norm = 0
    for k in range(n + 1):
        x = k * h
        c = 1 if k in (0, n) else (4 if k % 2 else 2)
        if 0 < x < 1:
            s += c * weight(x) * exp(-rho / sqrt(x * (1 - x))); norm += c * weight(x)
    return s / norm

print('1. Identity check: data profile (Uehling) vs mixture with weight 6x(1-x)')
worst = max(abs(uehling(r, 20000)[1] - mix(r, lambda x: x * (1 - x), 20000)) for r in (0.001, 0.01, 0.1, 0.3, 1, 3))
print(f'   max difference over rho = 0.001..3: {worst:.1e}')

def E_of(weight, n=200000):
    # ln(r_eff/D) = E_w[ln(R sqrt(x(1-x)))/D] - gamma
    h = 1 / n; s = 0; norm = 0
    for k in range(1, n):
        x = k * h; w = weight(x); s += w * 0.5 * log(x * (1 - x)); norm += w
    return s / norm - log(2) - GAMMA

def half(f):
    lo, hi = 1e-4, 5.0
    for _ in range(60):
        m = sqrt(lo * hi)
        lo, hi = (m, hi) if f(m) > 0.5 else (lo, m)
    return lo / 2                          # in units of D (rho = r/R, r/D = rho/2)

print('\n2. Variants of the split x (target: ln(r_eff/D) = -2.1037, half at 0.150 D, tail length R/2)')
S1E = log(0.5) - log(2) - GAMMA
print(f'   S1 symmetric eight x = 1/2 (A-74): ln(r_eff/D) = {S1E:+.4f}; half at {log(2)/2/2:.4f} D; tail e^(-2r/R), length R/2')
for name, wgt in (('S2 any split, weight 1', lambda x: 1.0), ('S3 weight x(1-x) (product rule)', lambda x: x * (1 - x))):
    E = E_of(wgt)
    print(f'   {name}: ln(r_eff/D) = {E:+.4f}; half at {half(lambda rho: mix(rho, wgt)):.4f} D; tail length R/2 (from x near 1/2)')
print(f'   data: {-(5/6 + GAMMA + log(2)):+.4f}')
kappa = 2 / (3 * pi)
print('\n3. In 1/alpha for three leptons the edge error is 3*kappa*(E - target):')
for name, E in (('S1', S1E), ('S2', E_of(lambda x: 1.0))):
    print(f'   {name}: {3 * kappa * (E - (-(5/6 + GAMMA + log(2)))):+.3f}')
