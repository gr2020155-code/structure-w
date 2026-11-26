"""Hadronic part of 1/alpha from e+e- -> hadrons data R(s) via the W link 'screening <-> pair production' (26.11.2026).
Rules: hadron_dispersion_rules_2026-11-26.md (pushed before run).
Run: python -I -X utf8 hadron_dispersion.py <path to rpp2022-hadronicrpp_page1001.dat>    (pure python)
"""
import sys
from math import pi, sqrt, log

path = sys.argv[1]
alpha0 = 1 / 137.035999206
inv_W_parts, s_W = 136.097, 0.14                     # P-91
m_p = 938.27209
parts = [(m_p / 3, 4 / 3), (m_p / 3, 1 / 3), (1115.683 - 2 * m_p / 3, 1 / 3), (2286.46 - 2 * m_p / 3, 4 / 3), (5619.60 - 2 * m_p / 3, 1 / 3)]  # MeV, NcQ2
narrow = {'J/psi': (3096.9, 5.53, 0.10), 'psi(2S)': (3686.1, 2.33, 0.04), 'Y(1S)': (9460.3, 1.340, 0.018),
          'Y(2S)': (10023.3, 0.612, 0.011), 'Y(3S)': (10355.2, 0.443, 0.008)}             # MeV, keV, keV

# ---- read data
excl, incl = [], []
for line in open(path, encoding='latin-1'):
    if line.startswith('*') or not line.strip():
        continue
    q = line.split("'")
    num = q[0].split()
    if len(num) < 7 or len(q) < 4:
        continue
    try:
        x, R, st, sysp = float(num[0]), float(num[3]), float(num[4]), float(num[6])
    except ValueError:
        continue
    name = q[1].strip()
    tail = q[4].split() if len(q) > 4 else []
    ptp = abs(float(tail[0])) if tail else 0.0
    rec = (x, R, abs(st), sqrt(sysp**2 + ptp**2))
    if 'EXCLSUM' in name:
        excl.append(rec)
    elif x >= 2.0:
        incl.append(rec)
excl.sort(); incl.sort()
print(f'read: EXCLSUM {len(excl)} points {excl[0][0]}-{excl[-1][0]} GeV; inclusive >= 2 GeV: {len(incl)} points')


def bins_of(lo, hi, w):
    out = []; b = lo
    while b < hi - 1e-9:
        out.append((b, min(b + w, hi))); b += w
    return out


def binned(xmax):
    edges = bins_of(2.0, 5.0, 0.05) + bins_of(5.0, 11.0, 0.2) + bins_of(11.0, 40.0, 1.0)
    pts = []
    for lo, hi in edges:
        if lo >= xmax:
            break
        sel = [r for r in incl if lo <= r[0] < hi and r[0] <= xmax]
        if not sel:
            continue
        ws = [1 / r[2]**2 for r in sel]; W = sum(ws)
        pts.append((sum(w * r[0] for w, r in zip(ws, sel)) / W, sum(w * r[1] for w, r in zip(ws, sel)) / W,
                    1 / sqrt(W), sum(w * r[3] for w, r in zip(ws, sel)) / W))
    return pts


def trap(points, x_end=None):
    """int R ds/s = int 2R/x dx over points (x, R, stat, sys%): value, stat err, sys err"""
    val = st2 = sy = 0.0
    wts = [0.0] * len(points)
    for i in range(len(points) - 1):
        x0, x1 = points[i][0], points[i + 1][0]; h = x1 - x0
        wts[i] += h / 2 * 2 / x0; wts[i + 1] += h / 2 * 2 / x1
    for w, (x, R, st, sp) in zip(wts, points):
        val += w * R; st2 += (w * st)**2; sy += w * R * sp / 100
    return val, sqrt(st2), sy


def data_integral(xmax):
    pts = [(2 * 0.13957, 0.0, 0.0, 0.0)] + excl + [p for p in binned(xmax)]
    if pts[-1][0] < xmax:                       # flat to xmax from last bin
        pts.append((xmax,) + pts[-1][1:])
    return trap(pts)


def parts_integral(xmax, n=200000):
    """int R_parts ds/s from threshold to xmax (GeV); R_parts = sum NcQ2 * beta(3-beta^2)/2"""
    tot = 0.0
    for m, w in parts:
        x0 = 2 * m / 1000
        if x0 >= xmax:
            continue
        # substitution x = x0*cosh(u) keeps the threshold smooth
        from math import acosh, cosh, sinh
        U = acosh(xmax / x0); du = U / n; s = 0.0
        for k in range(n):
            u = (k + 0.5) * du; x = x0 * cosh(u); b = sqrt(1 - (x0 / x)**2)
            s += b * (3 - b * b) / 2 * 2 / x * x0 * sinh(u) * du
        tot += w * s
    return tot


def narrow_part():
    v = e = 0.0
    for M, G, dG in narrow.values():
        v += 3 * G * 1e-3 / (alpha0**2 * M); e += 3 * dG * 1e-3 / (alpha0**2 * M)
    return v, e


nv, ne = narrow_part()
print(f'narrow resonances: +{nv:.4f} +- {ne:.4f} in 1/alpha (' + ', '.join(f'{k} {3*G*1e-3/(alpha0**2*M):.4f}' for k, (M, G, _) in narrow.items()) + ')')
for xmax in (40.0, 30.0, 20.0):
    D, Dst, Dsy = data_integral(xmax); P = parts_integral(xmax)
    corr = (D - P) / (3 * pi)
    err_d = sqrt(Dst**2 + Dsy**2) / (3 * pi)
    inv = inv_W_parts + corr + nv
    sig = sqrt(s_W**2 + err_d**2 + ne**2)
    tag = 'MAIN' if xmax == 40.0 else 'info'
    print(f'[{tag}] s_max = ({xmax:.0f} GeV)^2: int R_data ds/s = {D:.4f} (stat {Dst:.4f}, sys {Dsy:.4f}); int R_parts ds/s = {P:.4f};'
          f' correction (data - parts)/(3pi) = {corr:+.4f} +- {err_d:.4f}')
    print(f'       1/alpha(0) = {inv_W_parts} + {corr:+.4f} + {nv:.4f} = {inv:.3f} +- {sig:.3f}  vs 137.036  (dev {inv-137.036:+.3f}, {(inv-137.036)/sig:+.1f} sigma)'
          f' -> {"PASS" if abs(inv-137.036) <= 2*sig else "FAIL"}' if tag == 'MAIN' else f'       1/alpha(0) = {inv:.3f} +- {sig:.3f}')

print('\nlink check: weighted mean R of inclusive points vs sum NcQ2 of active parts')
for lo, hi, want in ((5.0, 9.3, 10 / 3), (12.0, 40.0, 11 / 3)):
    sel = [r for r in incl if lo <= r[0] <= hi]
    ws = [1 / r[2]**2 for r in sel]; W = sum(ws)
    Rm = sum(w * r[1] for w, r in zip(ws, sel)) / W; st = 1 / sqrt(W); sp = sum(w * r[3] for w, r in zip(ws, sel)) / W
    err = sqrt(st**2 + (Rm * sp / 100)**2)
    print(f'  {lo}-{hi} GeV: {len(sel)} points, R = {Rm:.3f} +- {err:.3f} (stat {st:.3f}, sys {Rm*sp/100:.3f})  vs {want:.3f}: {100*(Rm/want-1):+.1f} %, {(Rm-want)/err:+.1f} sigma')
