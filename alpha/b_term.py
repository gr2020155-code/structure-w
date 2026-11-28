"""Is the term B constant or tied to the tick count? g^2 = g0^2 [Or u^-4 + Om u^-3] + C T^-n, n = 0, 1, 2 (R-156).
Rules: b_term_rules_2026-11-28.md (committed before the run, bf6cd09). Kinematics as in R-132..R-134, R-155 (V-374, V-324, V-387).
usage: python -P -E -X utf8 b_term.py <desi_dr2_bao_data.txt> <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, pi, log
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2 as chi2dist

bao, snd = sys.argv[1:3]
h = 0.7354
SIG, C_L, G, MPC, TCMB = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
Or = (4 * SIG / C_L * TCMB**4) / (3 * (h * 1e5 / MPC) ** 2 * C_L * C_L / (8 * pi * G))

L = open(bao).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cbi = np.linalg.inv(np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]))
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
zs_all = np.array([float(r[hdr.index('zHD')]) for r in data]); ms_all = np.array([float(r[hdr.index('m_b_corr')]) for r in data])
cv = np.loadtxt(f'{snd}/STATSYS.cov'); nn = int(cv[0]); C_all = cv[1:].reshape(nn, nn)
sel = zs_all > 0.01; zs, ms = zs_all[sel], ms_all[sel]
Csi = np.linalg.inv(C_all[np.ix_(sel, sel)]); one = np.ones_like(zs); w1 = Csi @ one; s11 = one @ w1
REF = 1403.3
X0 = log(1e-8)
ZG = np.linspace(0, 2.5, 25001)


def history(Om, C, n):
    """T(x), x = ln u, from u = 1e-8 to 1; returns dense solution (units g0 = 1)."""
    u0 = 1e-8
    if n == 2:
        T0i = u0**2 * sqrt((0.25 - C) / Or)
    else:
        T0i = u0**2 / (2 * sqrt(Or))
    def rhs(x, T):
        return [1 / sqrt(Or * np.exp(-4 * x) + Om * np.exp(-3 * x) + C * T[0] ** (-n))]
    return solve_ivp(rhs, (X0, 0.0), [T0i], method='LSODA', rtol=1e-10, atol=1e-30, dense_output=True)


def solve_C(Om, n):
    """C such that g = g0 at u = 1. None if no solution."""
    if n == 0:
        return 1 - Or - Om
    def f(C):
        s = history(Om, C, n); T0 = s.y[0, -1]
        return Or + Om + C * T0 ** (-n) - 1
    hi = 0.25 - 1e-9 if n == 2 else 5.0
    if f(hi) < 0:
        return None
    return brentq(f, 1e-12, hi, xtol=1e-13)


def model(Om, n):
    C = solve_C(Om, n)
    if C is None:
        return None
    s = history(Om, C, n)
    T0 = s.y[0, -1]
    x = -np.log1p(ZG); T = s.sol(x)[0]
    E = np.sqrt(Or * (1 + ZG) ** 4 + Om * (1 + ZG) ** 3 + C * T ** (-n))
    D = cumulative_trapezoid(1 / E, ZG, initial=0)
    return dict(C=C, T0=T0, E=E, D=D, OB=C * T0 ** (-n))


def chi2_bao(mo):
    f = []
    for zz, k in zip(zb, kb):
        dm = np.interp(zz, ZG, mo['D']); dh = 1 / np.interp(zz, ZG, mo['E'])
        f.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    f = np.array(f); s = (f @ Cbi @ vb) / (f @ Cbi @ f)       # s = 1/g0 in units of r_d, best analytically
    r = vb - s * f
    return float(r @ Cbi @ r), 1 / s


def chi2_sn(mo):
    r = ms - 5 * np.log10((1 + zs) * np.interp(zs, ZG, mo['D'])); M = (r @ w1) / s11; e = r - M
    return float(e @ Csi @ e) - REF


def cb_of(Om, n):
    mo = model(Om, n)
    return (1e9, None) if mo is None else (chi2_bao(mo)[0], mo)


print(f'h = {h}; Or = {Or:.3e}\n')
GYR = 1e9 * 365.25 * 86400
for n in (0, 1, 2):
    grid = np.linspace(0.05, 0.60, 56)
    vals = [cb_of(om, n)[0] for om in grid]
    ok = [om for om, v in zip(grid, vals) if v < 1e8]
    if not ok:
        print(f'n = {n}: no solution for Om in 0.05-0.60 (C < 1/4 cannot give g = g0 today)  =>  FAIL\n')
        continue
    j = int(np.argmin(vals))
    lo, hi = grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]
    res = minimize_scalar(lambda om: cb_of(om, n)[0], bounds=(lo, hi), method='bounded', options=dict(xatol=1e-6))
    Om = res.x; cb, mo = cb_of(Om, n); g0r = chi2_bao(mo)[1]
    p = chi2dist.sf(cb, 11); sn = chi2_sn(mo)
    verdict = 'PASS' if p >= 0.05 and sn <= 11 else 'FAIL'
    print(f'n = {n}: Om = {Om:.4f}, OB today = {mo["OB"]:.4f}, C = {mo["C"]:.4f}  (solutions for Om in {min(ok):.2f}-{max(ok):.2f})')
    print(f'   BAO chi2 = {cb:.1f} / 11, probability {p:.3g} -> {"ok" if p >= 0.05 else "no"}; '
          f'SN d chi2 = {sn:+.1f} -> {"ok" if sn <= 11 else "no"}  =>  {verdict}')
    sng = [(chi2_sn(m2), om) for om in np.linspace(0.10, 0.50, 81) if (m2 := model(om, n)) is not None]
    sbest = min(sng)
    age = mo['T0'] / (h * 1e5 / MPC) / GYR
    print(f'   comparison: T0*g0 = {mo["T0"]:.4f}, w today = {-1 + n / (3 * mo["T0"]):+.3f}, age at h = {h}: {age:.2f} Gyr; '
          f'SN only: best Om ~ {sbest[1]:.3f}, d chi2 = {sbest[0]:+.1f}\n')
