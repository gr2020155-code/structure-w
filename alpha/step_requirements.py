"""What the step needs (no W rule here, only bookkeeping with the derived edge and slope).
Variant 1: extra charged forms below 5e-20 m (q0 = 1). Variant 2: q0 < 1 with known forms only."""
from math import pi, log, sqrt, exp
hbar = 1.054571817e-34; c = 299792458; MeV = 1.602176634e-13 / c**2
a = 1.0155e-34; k = 2 / (3 * pi); E = -2.1037            # edge ln(r_eff/D)
inv_far = 137.035999177
D = lambda m: 2 * hbar / (m * MeV * c)
def fermions(light=None):
    forms = [('e', 0.51099895, 1, 1), ('mu', 105.6583755, 1, 1), ('tau', 1776.86, 1, 1),
             ('u', 2.16, 3, 4/9), ('d', 4.70, 3, 1/9), ('s', 93.5, 3, 1/9),
             ('c', 1273.0, 3, 4/9), ('b', 4183.0, 3, 1/9), ('t', 172570.0, 3, 4/9)]
    tot = 0
    for name, m, Nc, Q2 in forms:
        if light and name in ('u', 'd', 's'): m = max(m, light)
        tot += k * Nc * Q2 * (log(D(m) / a) + E)
    return tot
for lab, light in (('current quark masses', None), ('light quarks as 300 MeV (hadron estimate)', 300.0)):
    F = fermions(light)
    for wlab, W in (('no W boson', 0.0), ('W boson, standard sign', -42.80)):
        need = inv_far - F - W            # what the step must supply including the 1 of q0 = 1
        q0 = 1 / sqrt(inv_far - F - W) if inv_far - F - W > 0 else float('nan')
        print(f'{lab:42s} {wlab:22s}: fermions {F:6.2f}; 1/alpha at the step {inv_far - F - W:6.2f} -> q0 = 1/{1/q0:5.2f};'
              f' with q0 = 1 missing {inv_far - F - W - 1:6.2f}')
# variant 1: a log-uniform tower of unit-charge forms between 5e-20 m and the step, density rho per e-fold of size
L = log(5e-20 / a)
for need in (64.3, 107.1):
    rho = 2 * need / (k * L * L)
    print(f'tower needed for {need}: {rho:.2f} unit-charge forms per e-fold of size below 5e-20 m ({rho * L:.0f} forms in all, mass ratio between neighbours {exp(1/rho):.1f})')
