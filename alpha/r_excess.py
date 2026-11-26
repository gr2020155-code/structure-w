"""Excess of R over the W weights sum NcQ^2: per region, experiments combined with their normalisation errors (26.11.2026).
Rules: r_excess_rules_2026-11-26.md (pushed before run). Data only (D-68).
Run: python -I -X utf8 r_excess.py <path to rpp2022-hadronicrpp_page1001.dat>
"""
import sys
from math import sqrt, log

rows = []
for line in open(sys.argv[1], encoding='latin-1'):
    if line.startswith('*') or not line.strip():
        continue
    q = line.split("'")
    num = q[0].split()
    if len(num) < 7 or len(q) < 4:
        continue
    try:
        x, R, st, norm = float(num[0]), float(num[3]), abs(float(num[4])), float(num[6])
    except ValueError:
        continue
    tail = q[4].split() if len(q) > 4 else []
    ptp = abs(float(tail[0])) if tail else 0.0
    rows.append(dict(x=x, R=R, e=sqrt(st**2 + (R * ptp / 100)**2), norm=norm, exp=(q[1].strip(), q[3].strip())))

regions = [('A', 2.0, 3.7, 2.0, [(3.05, 3.15), (3.65, 3.72)]), ('B', 5.0, 9.3, 10 / 3, []),
           ('C1', 11.2, 20.0, 11 / 3, []), ('C2', 20.0, 40.0, 11 / 3, [])]
res = []
for name, lo, hi, want, cut in regions:
    sel = [r for r in rows if lo <= r['x'] <= hi and 'EXCLSUM' not in r['exp'][0] and not any(a <= r['x'] <= b for a, b in cut)]
    exps = {}
    for r in sel:
        exps.setdefault(r['exp'], []).append(r)
    comb = []
    for k, pts in exps.items():
        w = [1 / p['e']**2 for p in pts]; W = sum(w)
        Rm = sum(wi * p['R'] for wi, p in zip(w, pts)) / W
        xm = sum(wi * p['x'] for wi, p in zip(w, pts)) / W
        nm = sum(p['norm'] for p in pts) / len(pts)
        comb.append((Rm, sqrt(1 / W + (Rm * nm / 100)**2), xm, k, len(pts)))
    w = [1 / c[1]**2 for c in comb]; W = sum(w)
    Rm = sum(wi * c[0] for wi, c in zip(w, comb)) / W; err = 1 / sqrt(W)
    xm = sum(wi * c[2] for wi, c in zip(w, comb)) / W
    chi2 = sum(((c[0] - Rm) / c[1])**2 for c in comb); n = len(comb) - 1
    scale = sqrt(chi2 / n) if n > 0 and chi2 > n else 1.0
    err *= scale
    eps, de = Rm / want - 1, err / want
    res.append((name, xm, eps, de))
    print(f'{name} {lo}-{hi} GeV: {len(sel)} points, {len(comb)} experiments, chi2/n = {chi2:.1f}/{n} (error x{scale:.2f}); '
          f'R = {Rm:.3f} +- {err:.3f} at sqrt(s) ~ {xm:.1f} GeV vs {want:.3f}: eps = {eps:+.4f} +- {de:.4f} ({eps/de:+.1f} sigma)')
    for c in sorted(comb, key=lambda c: c[2]):
        print(f'      {c[3][0]:16s} {c[3][1]:16s} n={c[4]:3d}  sqrt(s)~{c[2]:6.2f}  R = {c[0]:.3f} +- {c[1]:.3f}')

# straight line eps = e0 + b ln(x/10)
S = Sx = Sy = Sxx = Sxy = 0.0
for _, x, e, d in res:
    w = 1 / d**2; t = log(x / 10)
    S += w; Sx += w * t; Sy += w * e; Sxx += w * t * t; Sxy += w * t * e
Dn = S * Sxx - Sx**2
b = (S * Sxy - Sx * Sy) / Dn; e0 = (Sxx * Sy - Sx * Sxy) / Dn; sb = sqrt(S / Dn); se0 = sqrt(Sxx / Dn)
mean = Sy / S
chi_const = sum(((e - mean) / d)**2 for _, _, e, d in res)
chi_line = sum(((e - e0 - b * log(x / 10)) / d)**2 for _, x, e, d in res)
print(f'\nconstant eps: {mean:+.4f} +- {1/sqrt(S):.4f}, chi2 = {chi_const:.2f}/3')
print(f'line eps = e0 + b ln(sqrt(s)/10 GeV): e0 = {e0:+.4f} +- {se0:.4f}, b = {b:+.4f} +- {sb:.4f} ({b/sb:+.1f} sigma), chi2 = {chi_line:.2f}/2 -> '
      f'{"energy dependence FOUND" if abs(b) > 2 * sb else "no significant energy dependence"}')
print('reality: ' + ', '.join(f'{n} {e/d:+.1f} sigma' for n, _, e, d in res) + f' -> {"REAL (>3 sigma somewhere)" if any(abs(e/d) > 3 for _, _, e, d in res) else "not established"}')
