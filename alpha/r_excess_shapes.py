"""Three directions for the R excess, tested by the energy shape of eps (26.11.2026). Rules: r_excess_shapes_rules_2026-11-26.md (pushed before run).
Run: python -I -X utf8 r_excess_shapes.py <path to rpp2022-hadronicrpp_page1001.dat>   (pure python)
"""
import sys
from math import sqrt, log, exp, lgamma

rows = []
for line in open(sys.argv[1], encoding='latin-1'):
    if line.startswith('*') or not line.strip():
        continue
    q = line.split("'"); num = q[0].split()
    if len(num) < 7 or len(q) < 4:
        continue
    try:
        x, R, st, norm = float(num[0]), float(num[3]), abs(float(num[4])), float(num[6])
    except ValueError:
        continue
    tail = q[4].split() if len(q) > 4 else []
    ptp = abs(float(tail[0])) if tail else 0.0
    rows.append(dict(x=x, R=R, e=sqrt(st**2 + (R * ptp / 100)**2), norm=norm, exp=(q[1].strip(), q[3].strip())))

regions = [('A1', 2.0, 2.6, 2.0), ('A2', 2.6, 3.05, 2.0), ('A3', 3.15, 3.65, 2.0), ('B1', 5.0, 7.0, 10 / 3), ('B2', 7.0, 9.3, 10 / 3),
           ('C1', 11.2, 20.0, 11 / 3), ('C2', 20.0, 30.0, 11 / 3), ('C3', 30.0, 40.0, 11 / 3), ('D', 50.0, 64.0, 11 / 3)]
pts = []
for name, lo, hi, want in regions:
    sel = [r for r in rows if lo <= r['x'] <= hi and 'EXCLSUM' not in r['exp'][0]]
    exps = {}
    for r in sel:
        exps.setdefault(r['exp'], []).append(r)
    comb = []
    for k, p in exps.items():
        w = [1 / a['e']**2 for a in p]; W = sum(w)
        Rm = sum(wi * a['R'] for wi, a in zip(w, p)) / W; xm = sum(wi * a['x'] for wi, a in zip(w, p)) / W
        nm = sum(a['norm'] for a in p) / len(p)
        comb.append((Rm, sqrt(1 / W + (Rm * nm / 100)**2), xm))
    w = [1 / c[1]**2 for c in comb]; W = sum(w)
    Rm = sum(wi * c[0] for wi, c in zip(w, comb)) / W; err = 1 / sqrt(W); xm = sum(wi * c[2] for wi, c in zip(w, comb)) / W
    chi2 = sum(((c[0] - Rm) / c[1])**2 for c in comb); n = len(comb) - 1
    if n > 0 and chi2 > n:
        err *= sqrt(chi2 / n)
    pts.append((name, xm, Rm / want - 1, err / want, len(sel), len(comb)))
    print(f'{name} {lo}-{hi} GeV: {len(sel)} pts, {len(comb)} exp, E ~ {xm:.2f} GeV, eps = {Rm/want-1:+.4f} +- {err/want:.4f}')

MZ = 91.1876
V1 = lambda E: E * E / (MZ * MZ - E * E)
V2 = lambda E: V1(E)**2
models = {
    '(a) const + Z': [lambda E: 1.0, V1, V2],
    '(b1) c/ln(E/0.3128) + Z': [lambda E: 1 / log(E / 0.3128), V1, V2],
    '(b2) c/ln(E/0.13957) + Z': [lambda E: 1 / log(E / 0.13957), V1, V2],
    '(v) Z only': [V1, V2]}


def solve(A, b):
    n = len(b); M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i])); M[i], M[p] = M[p], M[i]
        for r in range(n):
            if r != i:
                f = M[r][i] / M[i][i]
                M[r] = [a - f * c for a, c in zip(M[r], M[i])]
    return [M[i][n] / M[i][i] for i in range(n)]


def chi2_prob(chi2, dof):            # upper regularized gamma Q(dof/2, chi2/2) by series for P
    a, x = dof / 2, chi2 / 2
    s = term = 1 / a; k = 1
    while term > 1e-14 * s:
        term *= x / (a + k); s += term; k += 1
    P = s * exp(-x + a * log(x) - lgamma(a))
    return 1 - P


print()
best = {}
for name, basis in models.items():
    A = [[sum(f(E) * g(E) / d**2 for _, E, e, d, _, _ in pts) for g in basis] for f in basis]
    b = [sum(f(E) * e / d**2 for _, E, e, d, _, _ in pts) for f in basis]
    c = solve(A, b)
    model = lambda E: sum(ci * f(E) for ci, f in zip(c, basis))
    chi2 = sum(((e - model(E)) / d)**2 for _, E, e, d, _, _ in pts); dof = len(pts) - len(basis)
    pr = chi2_prob(chi2, dof)
    best[name] = chi2
    print(f'{name}: params ' + ', '.join(f'{v:+.4f}' for v in c) + f'; chi2 = {chi2:.2f}/{dof}, prob = {pr:.3f} -> {"EXCLUDED" if pr < 0.01 else "acceptable"}')
    print('     residuals (sigma): ' + ' '.join(f'{n}:{(e-model(E))/d:+.1f}' for n, E, e, d, _, _ in pts))
    print('     model eps at 3, 8, 14, 31 GeV: ' + ', '.join(f'{model(E):+.3f}' for E in (3.1, 7.7, 13.7, 31.0)) +
          (f'; strong-like part only at 31 GeV: {c[0]*basis[0](31.0):+.3f}' if len(basis) == 3 else ''))
d1 = best['(a) const + Z'] - best['(b1) c/ln(E/0.3128) + Z']; d2 = best['(a) const + Z'] - best['(b2) c/ln(E/0.13957) + Z']
print(f'\nchi2(a) - chi2(b1) = {d1:+.2f}; chi2(a) - chi2(b2) = {d2:+.2f}  (threshold 4) -> ' +
      ('data prefer (b)' if max(d1, d2) > 4 else 'data prefer (a)' if min(d1, d2) < -4 else 'data do not choose between (a) and (b)'))
