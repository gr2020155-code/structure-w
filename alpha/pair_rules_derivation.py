"""Weight of the split and decay length from positions (rules: pair_rules_derivation_2026-11-25.md).
Diameter of n steps, junction O splits it into loops of n_A and n_B = n - n_A steps (diameters).
Weight: pairs n_A*n_B (A-31) or 1.  Decay: lambda = 1/sqrt(<k^2>) (A-39, squares) or 1/<k> (linear),
time share of a loop = its length share (A-84 + one step per tact), k = 1/radius."""
from math import log, sqrt, exp
GAMMA = 0.5772156649015329
TARGET = -(5 / 6 + GAMMA + log(2)); HALF_DATA = 0.1502
n = 2000
D = float(n)

def lam(nA, nB, mode):
    rA, rB = nA / 2, nB / 2                        # loop radii in steps
    xA, xB = nA / n, nB / n                         # time shares = length shares
    if mode == 'squares':
        return 1 / sqrt(xA / rA ** 2 + xB / rB ** 2)
    return 1 / (xA / rA + xB / rB)

def profile(weight, mode):
    items = [(weight(nA, n - nA), lam(nA, n - nA, mode)) for nA in range(1, n)]
    W = sum(w for w, _ in items)
    E = sum(w * log(l / D) for w, l in items) / W - GAMMA
    w_at = lambda r: sum(w * exp(-r / l) for w, l in items) / W
    lo, hi = 1e-4 * D, 5 * D
    for _ in range(60):
        m = sqrt(lo * hi); lo, hi = (m, hi) if w_at(m) > 0.5 else (lo, m)
    return E, lo / D, max(l for _, l in items) / D

print(f'n = {n} steps; data: ln(r_eff/D) = {TARGET:+.4f}, half at {HALF_DATA} D, longest decay length D/4 = {0.25}')
for wname, wf in (('pairs nA*nB (A-31)', lambda a, b: a * b), ('weight 1', lambda a, b: 1.0)):
    for mode in ('squares', 'linear'):
        E, h, lmax = profile(wf, mode)
        ok = abs(E - TARGET) < 0.05 and abs(h - HALF_DATA) < 0.03
        print(f'   {wname:20s} x {mode:7s}: ln(r_eff/D) = {E:+.4f}  half at {h:.4f} D  longest decay length {lmax:.4f} D  {"HIT" if ok else "miss"}')
