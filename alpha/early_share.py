"""Early light-like share checked with primordial helium (R-139). Rules: early_share_rules_2026-11-27.md.
usage: python -P -E -X utf8 early_share.py <desi_dr2_bao_data.txt>"""
import sys
import numpy as np
from math import sqrt, pi, exp, log
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar, brentq

KB, HP, ME, MH, SIGT = 1.380649e-23, 6.62607015e-34, 9.1093837015e-31, 1.6735575e-27, 6.6524587321e-29
SIG, C, G, MPC, T0 = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
EV = 1.602176634e-19; BH = 13.605693 * EV; LAM_A = 121.567e-9; EA = HP * C / LAM_A; L2S = 8.2246
NU = 1 + 0.2271 * 3.046; THETA = 0.010411; H100 = 1e5 / MPC; KN = NU / 0.2271
rho_gam = 4 * SIG / C**3 * T0**4
def Om_gamma(h): return rho_gam / (3 * (h * H100)**2 / (8 * pi * G))

L = open(sys.argv[1]).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cbi = np.linalg.inv(np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]))
def bao_chi2(Om, Or):
    E = lambda t: sqrt(Or * (1 + t)**4 + Om * (1 + t)**3 + 1 - Om - Or); v = []
    for zz, k in zip(zb, kb):
        dm = quad(lambda t: 1 / E(t), 0, zz)[0]; dh = 1 / E(zz)
        v.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    v = np.array(v); b = (v @ Cbi @ vb) / (v @ Cbi @ v); r = vb - b * v; return float(r @ Cbi @ r)
def fit_Om(Or):
    return minimize_scalar(lambda o: bao_chi2(o, Or), bounds=(0.2, 0.4), method='bounded', options=dict(xatol=1e-6)).x

def theta(h, Obh2, Om, nu, F, f=0.0, Y=0.245):            # as R-135 / R-136
    Og = Om_gamma(h); Or = Og * nu * (1 + f); Ob = Obh2 / h**2; H0 = h * H100
    E = lambda z: sqrt(Or * (1 + z)**4 + Om * (1 + z)**3 + 1 - Om - Or)
    nH0 = (1 - Y) * Obh2 * 3 * H100**2 / (8 * pi * G) / MH
    def rhs(z, x):
        x = min(max(x[0], 1e-12), 1.0); T = T0 * (1 + z); t = T / 1e4; nH = nH0 * (1 + z)**3; H = H0 * E(z)
        aB = F * 1e-19 * 4.309 * t**-0.6166 / (1 + 0.6703 * t**0.53)
        beta = aB * (2 * pi * ME * KB * T / HP**2)**1.5 * exp(-BH / 4 / (KB * T))
        K = LAM_A**3 / (8 * pi * H); n1 = nH * (1 - x)
        Cf = (1 + K * L2S * n1) / (1 + K * (L2S + beta) * n1)
        return [(x * x * nH * aB - beta * (1 - x) * exp(-EA / (KB * T))) * Cf / (H * (1 + z))]
    z_start = 1700.0; T = T0 * (1 + z_start); nH = nH0 * (1 + z_start)**3
    S = (2 * pi * ME * KB * T / HP**2)**1.5 * exp(-BH / (KB * T)) / nH; x0 = (-S + sqrt(S * S + 4 * S)) / 2
    sol = solve_ivp(rhs, (z_start, 200.0), [x0], method='Radau', rtol=1e-7, atol=1e-10, dense_output=True)
    xe = lambda z: float(sol.sol(z)[0]) if z >= 200 else float(sol.sol(200.0)[0])
    dtau = lambda z: xe(z) * nH0 * (1 + z)**3 * SIGT * C / ((1 + z) * H0 * E(z))
    zs = np.linspace(200, 1700, 3001); dt = np.array([dtau(z) for z in zs])
    tau = np.concatenate([[0], np.cumsum((dt[1:] + dt[:-1]) / 2 * np.diff(zs))]) + quad(dtau, 0, 200)[0]
    zstar = float(np.interp(1.0, tau, zs))
    R = lambda z: 0.75 * Ob / Og / (1 + z)
    rs = quad(lambda x_: np.exp(x_) / sqrt(3 * (1 + R(np.exp(x_) - 1))) / E(np.exp(x_) - 1), log(1 + zstar), log(1e9), limit=400)[0]
    return rs / quad(lambda z: 1 / E(z), 0, zstar, limit=400)[0]

def Yp(obh2, dN): return 0.24696 * (273.75 * obh2 / 6.129)**0.039 * ((3 + dN) / 3)**0.163
DATA = (('Aver 2021', 0.2453, 0.0034), ('Hsyu 2020', 0.2436, 0.0040), ('EMPRESS 2022', 0.2370, 0.0034))

F = 1.14
for btag, obh2_0 in (('b1', 0.02166), ('b2', 0.02235)):
    y0 = Yp(obh2_0, 0.0)
    print(f'\n{btag}: Obh2 (dN = 0) = {obh2_0}, Yp without share = {y0:.4f}; ' +
          '; '.join(f'{t}: {(y0 - v) / s:+.1f} sigma' for t, v, s in DATA))
    for t, v, s in DATA:
        dmax = brentq(lambda d: Yp(obh2_0 * ((3 + d) / 3)**(0.396 / 1.597), d) - (v + 2 * s), -2.9, 5)
        print(f'   max dN allowed by {t} (2 sigma): {dmax:+.2f}  (f = {dmax / KN:+.4f})')
    for h in (0.7354, 0.710):
        Om = fit_Om(Om_gamma(h) * NU); f = 0.14; obh2 = obh2_0
        for it in range(6):
            dN = KN * f; obh2 = obh2_0 * ((3 + dN) / 3)**(0.396 / 1.597)
            fa = theta(h, obh2, Om, NU, F, 0.0) - THETA
            f_new = brentq(lambda x: theta(h, obh2, Om, NU, F, x) - THETA, 0.0, 1.5, xtol=1e-5) if fa > 0 else 0.0
            if abs(f_new - f) < 1e-4: f = f_new; break
            f = f_new
        dN = KN * f; y = Yp(obh2, dN)
        res = [(t, (y - v) / s) for t, v, s in DATA]; npass = sum(abs(d) <= 2 for _, d in res)
        print(f'   h = {h}: f = {f:.4f} (dN = {dN:.2f}), Obh2 = {obh2:.5f}, Yp = {y:.4f}; ' +
              '; '.join(f'{t}: {d:+.1f} sigma' for t, d in res) + f' -> {"PASS" if npass >= 2 else "FAIL"} ({npass}/3)')
