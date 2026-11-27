"""All masses together at recombination (R-143). Rules: all_masses_rules_2026-11-27.md.
usage: python -P -E -X utf8 all_masses.py <desi_dr2_bao_data.txt>"""
import sys
import numpy as np
from math import sqrt, pi, exp, log
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar, brentq

KB, HP, ME0, MH, SIGT0 = 1.380649e-23, 6.62607015e-34, 9.1093837015e-31, 1.6735575e-27, 6.6524587321e-29
SIG, C, G, MPC, T0 = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
EV = 1.602176634e-19; BH0 = 13.605693 * EV; LAM_A0 = 121.567e-9; L2S0 = 8.2246
NU = 1 + 0.2271 * 3.046; THETA = 0.010411; H100 = 1e5 / MPC
rho_gam = 4 * SIG / C**3 * T0**4
def Om_gamma(h): return rho_gam / (3 * (h * H100)**2 / (8 * pi * G))

L = open(sys.argv[1]).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cbi = np.linalg.inv(np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]))
def bao_v(Om, Or):
    E = lambda t: sqrt(Or * (1 + t)**4 + Om * (1 + t)**3 + 1 - Om - Or); v = []
    for zz, k in zip(zb, kb):
        dm = quad(lambda t: 1 / E(t), 0, zz)[0]; dh = 1 / E(zz)
        v.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    return np.array(v)
def bao_chi2(Om, Or):
    v = bao_v(Om, Or); b = (v @ Cbi @ vb) / (v @ Cbi @ v); r = vb - b * v; return float(r @ Cbi @ r), b
def fit_Om(Or):
    r = minimize_scalar(lambda o: bao_chi2(o, Or)[0], bounds=(0.2, 0.4), method='bounded', options=dict(xatol=1e-6))
    return r.x, bao_chi2(r.x, Or)[1], r.fun
def b_sigma(Or):
    Om, b, ch = fit_Om(Or)
    prof = lambda bb: minimize_scalar(lambda o: float((lambda r: r @ Cbi @ r)(vb - bb * bao_v(o, Or))), bounds=(0.2, 0.4),
                                      method='bounded', options=dict(xatol=1e-6)).fun - ch - 1
    return b, (brentq(prof, b, b * 1.05) - brentq(prof, b * 0.95, b)) / 2, Om

def model(h, Obh2, Om, lam=1.0, links=False, zt=200.0, nu=NU, F=1.14, Y=0.245):
    le = lam
    ME = ME0 * le; BH = BH0 * le; LAM_A = LAM_A0 / le; EA = HP * C / LAM_A
    L2S = L2S0 * le; SIGT = SIGT0 / le**2; sc = le; arate = 1 / le**2
    Og = Om_gamma(h); Or = Og * nu; Ob = Obh2 / h**2; Op = Om - Ob; H0 = h * H100
    s = lambda z: 1 + (lam - 1) / (1 + exp(-min(max((z - zt) / 20.0, -700), 700)))
    E = lambda z: sqrt(Or * (1 + z)**4 + (Ob * s(z) + Op * (s(z) if links else 1.0)) * (1 + z)**3 + 1 - Om - Or)
    nH0 = (1 - Y) * Obh2 * 3 * H100**2 / (8 * pi * G) / MH
    def rhs(z, x):
        x = min(max(x[0], 1e-12), 1.0); T = T0 * (1 + z); t = T / sc / 1e4; nH = nH0 * (1 + z)**3; H = H0 * E(z)
        aB = arate * F * 1e-19 * 4.309 * t**-0.6166 / (1 + 0.6703 * t**0.53)
        beta = aB * (2 * pi * ME * KB * T / HP**2)**1.5 * exp(-BH / 4 / (KB * T))
        K = LAM_A**3 / (8 * pi * H); n1 = nH * (1 - x)
        Cf = (1 + K * L2S * n1) / (1 + K * (L2S + beta) * n1)
        return [(x * x * nH * aB - beta * (1 - x) * exp(-EA / (KB * T))) * Cf / (H * (1 + z))]
    z_start = 1700.0 * sc; T = T0 * (1 + z_start); nH = nH0 * (1 + z_start)**3
    S = (2 * pi * ME * KB * T / HP**2)**1.5 * exp(-BH / (KB * T)) / nH; x0 = (-S + sqrt(S * S + 4 * S)) / 2
    zend = 200.0 * sc
    sol = solve_ivp(rhs, (z_start, zend), [x0], method='Radau', rtol=1e-7, atol=1e-10, dense_output=True)
    xe = lambda z: float(sol.sol(z)[0]) if z >= zend else float(sol.sol(zend)[0])
    R = lambda z: 0.75 * Ob * s(z) / Og / (1 + z)
    dtau = lambda z: xe(z) * nH0 * (1 + z)**3 * SIGT * C / ((1 + z) * H0 * E(z))
    zs = np.linspace(zend, z_start, 3001); dt = np.array([dtau(z) for z in zs]); dd = dt / np.array([R(z) for z in zs])
    cum = lambda a: np.concatenate([[0], np.cumsum((a[1:] + a[:-1]) / 2 * np.diff(zs))])
    tau = cum(dt) + quad(dtau, 0, zend, limit=200)[0]; taud = cum(dd) + quad(lambda z: dtau(z) / R(z), 0, zend, limit=200)[0]
    zstar = float(np.interp(1.0, tau, zs)); zdrag = float(np.interp(1.0, taud, zs))
    rsf = lambda za: quad(lambda x_: np.exp(x_) / sqrt(3 * (1 + R(np.exp(x_) - 1))) / E(np.exp(x_) - 1), log(1 + za), log(1e9), limit=400)[0]
    DM = quad(lambda z: 1 / E(z), 0, zstar, limit=400, points=[zt])[0]
    return dict(theta=rsf(zstar) / DM, zstar=zstar, zdrag=zdrag, rsd=rsf(zdrag))

OBH2, eh = 0.02235, 0.0105
chk = model(0.7354, OBH2, fit_Om(Om_gamma(0.7354) * NU)[0])
print(f'check (lam = 1, as R-142): theta = {chk["theta"]:.6f}, z* = {chk["zstar"]:.1f}, z_drag = {chk["zdrag"]:.1f}')
for h in (0.7354, 0.710, 0.685):
    b, sb, Om = b_sigma(Om_gamma(h) * NU)
    print(f'\n=== h = {h}: BAO Om = {Om:.4f}, b = {b:.4f} +- {sb:.4f}')
    for tag, links in (('P1 links unchanged', False), ('P2 links scale too', True)):
        fun = lambda x: model(h, OBH2, Om, lam=x, links=links)['theta'] - THETA
        fa, fb = fun(0.90), fun(1.10)
        if fa * fb > 0:
            print(f'   {tag}: no root in 0.90..1.10 (theta - data {fa:+.5f} .. {fb:+.5f}) -> FAIL'); continue
        x = brentq(fun, 0.90, 1.10, xtol=1e-5); mm = model(h, OBH2, Om, lam=x, links=links); bm = 1 / mm['rsd']
        Omh = fit_Om(Om_gamma(h + eh) * NU)[0]; sh = abs(1 / model(h + eh, OBH2, Omh, lam=x, links=links)['rsd'] - bm)
        dev = (bm - b) / sqrt(sb**2 + sh**2)
        print(f'   {tag}: all masses {(x - 1) * 100:+.2f} %; z* = {mm["zstar"]:.1f}, z_drag = {mm["zdrag"]:.1f}, '
              f'r_d = {mm["rsd"] * 299792.458 / (100 * h):.2f} Mpc; b model = {bm:.4f} vs BAO {b:.4f} (+- {sb:.4f}, +- {sh:.4f} from h): '
              f'{dev:+.1f} sigma -> {"PASS" if abs(dev) <= 2 and abs(x - 1) <= 0.10 else "FAIL"}', flush=True)
