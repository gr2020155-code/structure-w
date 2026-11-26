"""Relic sound speed and recombination, exact angle test (R-135). Rules: relic_sound_rules_2026-11-27.md.
usage: python -P -E -X utf8 relic_sound.py <desi_dr2_bao_data.txt>"""
import sys
import numpy as np
from math import sqrt, pi, exp, log
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar, brentq

KB, HP, ME, MH, SIGT = 1.380649e-23, 6.62607015e-34, 9.1093837015e-31, 1.6735575e-27, 6.6524587321e-29
SIG, C, G, MPC, T0 = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
EV = 1.602176634e-19; BH = 13.605693 * EV; LAM_A = 121.567e-9; EA = HP * C / LAM_A; L2S = 8.2246
NU = 1 + 0.2271 * 3.046; THETA = 0.010411
H100 = 1e5 / MPC                                            # 100 km/s/Mpc in 1/s
rho_gam = 4 * SIG / C**3 * T0**4                             # kg/m^3 (energy/c^2)
def Om_gamma(h): return rho_gam / (3 * (h * H100)**2 / (8 * pi * G))

# BAO refit of Om at given Or (as R-133)
L = open(sys.argv[1]).read().split('\n'); i0, i1 = L.index('MEAN'), L.index('COV')
rr = [l.split() for l in L[i0 + 1:i1] if l.strip()]
zb = np.array([float(r[0]) for r in rr]); vb = np.array([float(r[1]) for r in rr]); kb = [r[2] for r in rr]
Cbi = np.linalg.inv(np.array([[float(x) for x in l.split()] for l in L[i1 + 1:] if l.strip()]))
def bao_chi2(Om, Or):                                        # g0 (= 1/r_d scale) enters linearly as 1/g0: solve exactly
    E = lambda t: sqrt(Or * (1 + t)**4 + Om * (1 + t)**3 + 1 - Om - Or); v = []
    for zz, k in zip(zb, kb):
        dm = quad(lambda t: 1 / E(t), 0, zz)[0]; dh = 1 / E(zz)
        v.append({'DM_over_rs': dm, 'DH_over_rs': dh, 'DV_over_rs': (zz * dm * dm * dh) ** (1 / 3)}[k])
    v = np.array(v); b = (v @ Cbi @ vb) / (v @ Cbi @ v); r = vb - b * v; return float(r @ Cbi @ r)
def fit_Om(Or):
    r = minimize_scalar(lambda o: bao_chi2(o, Or), bounds=(0.2, 0.4), method='bounded', options=dict(xatol=1e-6))
    lo = brentq(lambda o: bao_chi2(o, Or) - r.fun - 1, 0.2, r.x); hi = brentq(lambda o: bao_chi2(o, Or) - r.fun - 1, r.x, 0.4)
    return r.x, (hi - lo) / 2

def model(h, Obh2, Om, nu, F, Y=0.245):
    Or = Om_gamma(h) * nu; Og = Om_gamma(h); Ob = Obh2 / h**2
    H0 = h * H100
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
    tau = np.concatenate([[0], np.cumsum((dt[1:] + dt[:-1]) / 2 * np.diff(zs))])
    tau += quad(lambda z: xe(200.0) * nH0 * (1 + z)**3 * SIGT * C / ((1 + z) * H0 * E(z)), 0, 200)[0]
    zstar = float(np.interp(1.0, tau, zs))
    R = lambda z: 0.75 * Ob / Og / (1 + z)
    rs = quad(lambda x_: np.exp(x_) / sqrt(3 * (1 + R(np.exp(x_) - 1))) / E(np.exp(x_) - 1), log(1 + zstar), log(1e9), limit=400)[0]
    DM = quad(lambda z: 1 / E(z), 0, zstar, limit=400)[0]
    return rs / DM, zstar

hs, eh = 0.7354, 0.0105
print(f'h = {hs} +- {eh}; theta* data = {THETA}')
for ntag, nu in (('nu0 photons', 1.0), ('nu1 +neutrinos (bridge)', NU)):
    Om, eOm = fit_Om(Om_gamma(hs) * nu)
    print(f'\n{ntag}: BAO Om = {Om:.4f} +- {eOm:.4f}')
    for btag, obh2, eob in (('b1', 0.02166, sqrt(0.015**2 + 0.011**2) / 100), ('b2', 0.02235, sqrt(0.016**2 + 0.033**2) / 100)):
        for F in (1.00, 1.14):
            th, zst = model(hs, obh2, Om, nu, F)
            d = [model(hs + eh, obh2, fit_Om(Om_gamma(hs + eh) * nu)[0], nu, F)[0] - th,
                 model(hs, obh2, Om + eOm, nu, F)[0] - th,
                 model(hs, obh2 + eob, Om, nu, F)[0] - th]
            sth = sqrt(sum(x * x for x in d)); dev = (th - THETA) / sth
            fh = lambda hh: model(hh, obh2, fit_Om(Om_gamma(hh) * nu)[0], nu, F)[0] - THETA
            fa, fb = fh(0.40), fh(0.95)
            hneed = f'{brentq(fh, 0.40, 0.95, xtol=1e-4):.4f}' if fa * fb < 0 else f'none in 0.40..0.95 (theta - data: {fa:+.5f} .. {fb:+.5f})'
            print(f'   {btag} F = {F:.2f}: z* = {zst:.1f}, theta = {th:.6f} +- {sth:.6f} ({(th / THETA - 1) * 100:+.2f} %, {dev:+.1f} sigma)'
                  f' -> {"PASS" if abs(dev) <= 2 else "FAIL"};  h for exact theta = {hneed}')
th1, _ = model(hs, 0.02235, fit_Om(Om_gamma(hs) * NU)[0], NU, 1.14, Y=0.235); th2, _ = model(hs, 0.02235, fit_Om(Om_gamma(hs) * NU)[0], NU, 1.14, Y=0.255)
print(f'\nhelium check (nu1 b2 F 1.14): Y = 0.235 -> {th1:.6f}, Y = 0.255 -> {th2:.6f}')
