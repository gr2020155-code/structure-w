"""Masses from the form: a predeclared, finite search - written BEFORE running (16.11.2026).

Forms (author's positions, PC 4, A-07, A-23, A-52); features read from them:
            E eights  Pt points  Lp loops  D dims  X crosses  S shared loops  T extra point  O opposite cross
  e           1         3          2        1        0          0               0              0
  mu          2         5          4        2        1          0               0              0
  tau         3         7          6        3        2          0               0              0
  p           2         4          3        1        0          1               0              0
  n           3         6          5        2        1          1               0              0
  pi+-        2         6          4        2        1          0               1              0
  pi0         2         5          4        2        1          0               0              1
Family (ties masses to alpha, as W requires and as Barut's precedent suggests):
  m/m_e = 1 + (1/alpha) * (a*g + b*h),  g, h - two different functions from the list G below.
a, b are fixed by mu and p (two equations, two unknowns). Predicted: tau, n, pi+-, pi0.
Success: all four within 1 %. The number of tried pairs is reported: a hit is meaningful only if it is rare.
"""
import itertools, numpy as np
alpha = 1 / 137.035999084
F = {  # E, Pt, Lp, D, X, S, T, O
 'e':   (1,3,2,1,0,0,0,0), 'mu': (2,5,4,2,1,0,0,0), 'tau': (3,7,6,3,2,0,0,0),
 'p':   (2,4,3,1,0,1,0,0), 'n':  (3,6,5,2,1,1,0,0), 'pi+': (2,6,4,2,1,0,1,0), 'pi0': (2,5,4,2,1,0,0,1)}
mass = {'e': 1, 'mu': 206.7682830, 'tau': 3477.23, 'p': 1836.15267343, 'n': 1838.68366173, 'pi+': 273.13211, 'pi0': 264.14251}
names = ['E', 'Pt', 'Lp', 'D', 'X', 'S', 'T', 'O']
def base(f):
    E, Pt, Lp, D, X, S, T, O = f
    b = dict(zip(names, f))
    b['Sd'] = sum(k ** 4 for k in range(D))           # Barut-type sum over added dimensions
    b['E-1'] = E - 1; b['D-1'] = D - 1; b['Pt-3'] = Pt - 3; b['Lp-2'] = Lp - 2
    return b
G = {}
keys = list(base(F['e']).keys())
for k in keys:
    for pw in (1, 2, 3, 4):
        G[f'{k}^{pw}'] = (lambda f, k=k, pw=pw: base(f)[k] ** pw)
for k1, k2 in itertools.combinations(keys, 2):
    G[f'{k1}*{k2}'] = (lambda f, k1=k1, k2=k2: base(f)[k1] * base(f)[k2])
pred_names = ['tau', 'n', 'pi+', 'pi0']
tried = hits = 0; best = []
for (gn, g), (hn, h) in itertools.combinations(G.items(), 2):
    A = np.array([[g(F['mu']), h(F['mu'])], [g(F['p']), h(F['p'])]], float)
    y = np.array([(mass['mu'] - 1) * alpha, (mass['p'] - 1) * alpha])
    if abs(np.linalg.det(A)) < 1e-9: continue
    a, b = np.linalg.solve(A, y)
    # electron itself must give m = m_e: a*g(e) + b*h(e) = 0
    if abs(a * g(F['e']) + b * h(F['e'])) > 1e-9: continue
    tried += 1
    errs = [(1 + (a * g(F[x]) + b * h(F[x])) / alpha) / mass[x] - 1 for x in pred_names]
    worst = max(abs(e) for e in errs)
    best.append((worst, gn, hn, a, b, errs))
    if worst < 0.01: hits += 1
best.sort()
print(f'candidate functions {len(G)}; pairs valid (electron gives m_e) {tried}; all four within 1 %: {hits}')
for w, gn, hn, a, b, errs in best[:10]:
    print(f'  worst {w*100:6.2f} %  F = {a:+.4f}*{gn} {b:+.4f}*{hn}   tau {errs[0]*100:+.2f} %  n {errs[1]*100:+.2f} %  pi+ {errs[2]*100:+.2f} %  pi0 {errs[3]*100:+.2f} %')
