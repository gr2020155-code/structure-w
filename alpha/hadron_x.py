"""X rule for wholes made of cross and tau-form fragments (R-145).
Rules: hadron_x_rules_2026-11-27.md (committed before the run, 46351ac).
usage: python -P -E -X utf8 hadron_x.py"""
from math import comb

ME, IA = 0.51099895069, 137.035999177
U = 1.5 * IA
X = lambda m: (m / ME - 1) / U          # m in MeV -> X
M = lambda x: (1 + U * x) * ME          # X -> m in MeV
XLEP = {1: 0, 2: 1, 3: 17}
K = {'u': 2, 'd': 1, 'c': 2, 's': 1, 't': 2, 'b': 1}
G = {'u': 1, 'd': 1, 'c': 2, 's': 2, 't': 3, 'b': 3}
BASE = {'p': 938.27208943, 'omega': 782.66}
MAIN = [('Omega-', 'sss', 1672.43, 'p'), ('phi', 'ss', 1019.460, 'omega'),
        ('J/psi', 'cc', 3096.900, 'omega'), ('Upsilon', 'bb', 9460.40, 'omega')]
SIDE = [('Lambda', 'uds', 1115.683, 'p'), ('Sigma0', 'uds', 1192.642, 'p'), ('Xi-', 'dss', 1321.70, 'p'),
        ('Xi0', 'uss', 1314.82, 'p'), ('eta_c', 'cc', 2984.09, 'omega'), ('eta_b', 'bb', 9398.7, 'omega')]
F = {'A1': lambda k, g: XLEP[g], 'A2': lambda k, g: k / 3 * XLEP[g], 'A3': lambda k, g: (k / 3) ** 2 * XLEP[g],
     'A4': lambda k, g: comb(3, k) * (k / 3) ** 2 * XLEP[g], 'M1': lambda k, g: XLEP[g], 'M2': lambda k, g: k / 3 * XLEP[g]}

def pred(br, frag, base):
    xb = X(BASE[base]); f = [F[br](K[q], G[q]) for q in frag]
    if br[0] == 'A': return M(xb + sum(f))
    p = 1.0
    for v in f: p *= (1 + v) ** (1 / len(f))
    return M(xb * p)

print('X of data: ' + ', '.join(f'{n} {X(m):.3f}' for n, _, m, _ in MAIN) + f'; p {X(BASE["p"]):.3f}, omega {X(BASE["omega"]):.3f}')
for br in F:
    res = [(n, pred(br, fr, b) / m - 1) for n, fr, m, b in MAIN]
    ok = all(abs(d) <= 0.01 for _, d in res)
    side = ', '.join(f'{n} {100 * (pred(br, fr, b) / m - 1):+.1f}' for n, fr, m, b in SIDE)
    print(f'\n{br}: ' + ', '.join(f'{n} {100 * d:+.1f} %' for n, d in res) + f' -> {"PASS" if ok else "FAIL"}')
    print(f'    side (not counted, %): {side}')

xs_b, xs_m = (X(1672.43) - X(BASE['p'])) / 3, (X(1019.460) - X(BASE['omega'])) / 2
ms_b, ms_m = (X(1672.43) / X(BASE['p'])) ** (1 / 3), (X(1019.460) / X(BASE['omega'])) ** 0.5
print(f'\nDiagnostic (not counted): s term additive: baryons {xs_b:.3f}, mesons {xs_m:.3f}, ratio {xs_b / xs_m:.2f}; '
      f'multiplicative: baryons {ms_b:.4f}, mesons {ms_m:.4f}')
print(f'  additive s from Omega-p, applied to Lambda/Xi: Lambda {100 * (M(X(BASE["p"]) + xs_b) / 1115.683 - 1):+.1f} %, '
      f'Xi- {100 * (M(X(BASE["p"]) + 2 * xs_b) / 1321.70 - 1):+.1f} %')
print(f'  mass ratios: J/psi/phi {3096.900 / 1019.460:.4f}, Upsilon/J/psi {9460.40 / 3096.900:.4f}, J/psi/eta_c {3096.900 / 2984.09:.4f}')
