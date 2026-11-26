"""Hubble tension variants inside W (R-136). Rules: hubble_variants_rules_2026-11-27.md.
usage: python -P -E -X utf8 hubble_variants.py <desi_dr2_bao_data.txt> <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, pi, exp, log
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar, minimize, brentq

KB, HP, ME, MH, SIGT = 1.380649e-23, 6.62607015e-34, 9.1093837015e-31, 1.6735575e-27, 6.6524587321e-29
SIG, C, G, MPC, T0 = 5.670374419e-8, 299792458.0, 6.67430e-11, 3.0856775814913673e22, 2.72548
EV = 1.602176634e-19; BH = 13.605693 * EV; LAM_A = 121.567e-9; EA = HP * C / LAM_A; L2S = 8.2246
NU = 1 + 0.2271 * 3.046; THETA = 0.010411; CKM = 299792.458
H100 = 1e5 / MPC
rho_gam = 4 * SIG / C**3 * T0**4
def Om_gamma(h): return rho_gam / (3 * (h * H100)**2 / (8 * pi * G))

# ---------------- BAO (as R-135)
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
def bao_chi2_b(Om, Or, b):
    r = vb - b * bao_v(Om, Or); return float(r @ Cbi @ r)
def bao_chi2(Om, Or):
    v = bao_v(Om, Or); b = (v @ Cbi @ vb) / (v @ Cbi @ v); r = vb - b * v; return float(r @ Cbi @ r), b
def fit_Om(Or):
    r = minimize_scalar(lambda o: bao_chi2(o, Or)[0], bounds=(0.2, 0.4), method='bounded', options=dict(xatol=1e-6))
    return r.x, bao_chi2(r.x, Or)[1], r.fun
def b_sigma(Or):                                            # profile over Om, delta chi2 = 1
    Om, b, ch = fit_Om(Or)
    prof = lambda bb: minimize_scalar(lambda o: bao_chi2_b(o, Or, bb), bounds=(0.2, 0.4), method='bounded',
                                      options=dict(xatol=1e-6)).fun - ch - 1
    lo = brentq(prof, b * 0.95, b); hi = brentq(prof, b, b * 1.05); return b, (hi - lo) / 2, Om

# ---------------- relic model (as R-135) + drag epoch + extra early shares
def model(h, Obh2, Om, nu, F, f=0.0, q=0.0, Y=0.245):
    Og = Om_gamma(h); Or = Og * nu * (1 + f); Ob = Obh2 / h**2; H0 = h * H100
    Ox = q * Or / 1101.0**2
    E = lambda z: sqrt(Or * (1 + z)**4 + Om * (1 + z)**3 + Ox * (1 + z)**6 + 1 - Om - Or - Ox)
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
    R = lambda z: 0.75 * Ob / Og / (1 + z)
    dtau = lambda z: xe(z) * nH0 * (1 + z)**3 * SIGT * C / ((1 + z) * H0 * E(z))
    zs = np.linspace(200, 1700, 3001); dt = np.array([dtau(z) for z in zs]); dd = dt / np.array([R(z) for z in zs])
    cum = lambda a: np.concatenate([[0], np.cumsum((a[1:] + a[:-1]) / 2 * np.diff(zs))])
    tau = cum(dt) + quad(dtau, 0, 200)[0]; taud = cum(dd) + quad(lambda z: dtau(z) / R(z), 0, 200)[0]
    zstar = float(np.interp(1.0, tau, zs)); zdrag = float(np.interp(1.0, taud, zs))
    rsf = lambda za: quad(lambda x_: np.exp(x_) / sqrt(3 * (1 + R(np.exp(x_) - 1))) / E(np.exp(x_) - 1),
                          log(1 + za), log(1e9), limit=400)[0]
    DM = quad(lambda z: 1 / E(z), 0, zstar, limit=400)[0]
    return dict(theta=rsf(zstar) / DM, zstar=zstar, zdrag=zdrag, rsd=rsf(zdrag))

hs, eh, F = 0.7354, 0.0105, 1.14
BR = (('b1', 0.02166, sqrt(0.015**2 + 0.011**2) / 100), ('b2', 0.02235, sqrt(0.016**2 + 0.033**2) / 100))

print('=== V0: mean tempo from the sound ruler (BAO + drag), no relic; nu1, F = 1.14', flush=True)
for btag, obh2, eob in BR:
    def resid(h):
        Om = fit_Om(Om_gamma(h) * NU)[0]; return 1 / model(h, obh2, Om, NU, F)['rsd'] - b_sigma(Om_gamma(h) * NU)[0]
    h0 = brentq(resid, 0.55, 0.85, xtol=1e-5)
    b, sb, Om = b_sigma(Om_gamma(h0) * NU); mm = model(h0, obh2, Om, NU, F)
    dh = 1e-3
    dbdh = (1 / model(h0 + dh, obh2, fit_Om(Om_gamma(h0 + dh) * NU)[0], NU, F)['rsd'] - b_sigma(Om_gamma(h0 + dh) * NU)[0]
            - (1 / mm['rsd'] - b)) / dh
    dbdo = 1 / model(h0, obh2 + eob, Om, NU, F)['rsd'] - 1 / mm['rsd']
    sh = sqrt(sb**2 + dbdo**2) / abs(dbdh)
    print(f'   {btag}: h = {h0:.4f} +- {sh:.4f}; z_drag = {mm["zdrag"]:.1f}, r_d = {mm["rsd"] * CKM / (100 * h0):.2f} Mpc, '
          f'BAO b = {b:.4f} +- {sb:.4f}, Om = {Om:.4f}; vs relic-needed 0.685: {(h0 - 0.685) / sh:+.1f} sigma -> '
          f'{"AGREES" if abs(h0 - 0.685) <= 2 * sh else "DISAGREES"}; vs local 0.7354: {(h0 - hs) / sqrt(sh**2 + eh**2):+.1f} sigma',
          flush=True)

# ---------------- SNe: calibrators (Cepheid distances) + all z > 0.01, full covariance
snd = sys.argv[2]
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda nm: np.array([float(r[hdr.index(nm)]) for r in data])
z, m, cal, mu_c = col('zHD'), col('m_b_corr'), col('IS_CALIBRATOR') == 1, col('CEPH_DIST')
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
use = cal | (z > 0.01); z, m, cal, mu_c = z[use], m[use], cal[use], mu_c[use]
Ci = np.linalg.inv(C_all[np.ix_(use, use)]); one = np.ones(use.sum()); w1 = Ci @ one; s11 = one @ w1
OmS = fit_Om(Om_gamma(hs) * NU)[0]
zg = np.concatenate([[0], np.geomspace(1e-5, 2.6, 40000)]); Eg = np.sqrt(OmS * (1 + zg)**3 + 1 - OmS)
def cumint(y): return np.concatenate([[0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(zg))])
I0 = cumint(1 / Eg)
def chi2_sn(h, eps=0.0, A=0.0, zv=None):
    zt = (1 + z) / (1 + eps) - 1
    I = I0 if A == 0 else cumint(1 / (Eg * (1 + A / (1 + np.exp(np.clip((zg - zv) / (0.1 * zv), -700, 700))))))
    mu = 5 * np.log10((1 + z) * CKM / (100 * h) * np.interp(zt, zg, I)) + 25
    r = m - np.where(cal, mu_c, mu); M = (r @ w1) / s11; e = r - M; return float(e @ Ci @ e)
print(f'\n=== SNe: {use.sum()} (calibrators {cal.sum()}), Om = {OmS:.4f}', flush=True)
b0 = minimize_scalar(lambda hh: chi2_sn(hh), bounds=(0.6, 0.85), method='bounded', options=dict(xatol=1e-6)); X0 = b0.fun
print(f'   base: h = {b0.x:.4f}, chi2_0 = {X0:.1f}')
print('=== V1: local share of redshift, 1+z_obs = (1+eps)(1+z)')
r1 = minimize(lambda p: chi2_sn(p[0], p[1]), [0.73, 0.0], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-4))
c1 = minimize_scalar(lambda e: chi2_sn(0.685, e), bounds=(-0.02, 0.02), method='bounded', options=dict(xatol=1e-7))
print(f'   free: h = {r1.x[0]:.4f}, eps = {r1.x[1]:+.5f}, chi2 - chi2_0 = {r1.fun - X0:+.1f}')
print(f'   h = 0.685: eps = {c1.x:+.5f}, chi2 - chi2_0 = {c1.fun - X0:+.1f} -> {"PASS" if c1.fun - X0 <= 4 else "FAIL"}')
print('=== V2: local zone of other tempo, g = g_far (1 + A S(z)), S step at z_v')
for zv in (0.03, 0.05, 0.1, 0.15, 0.3, 0.5):
    r2 = minimize(lambda p: chi2_sn(p[0], 0.0, p[1], zv), [0.73, 0.0], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-4))
    c2 = minimize_scalar(lambda a: chi2_sn(0.685, 0.0, a, zv), bounds=(-0.3, 0.5), method='bounded', options=dict(xatol=1e-7))
    print(f'   z_v = {zv}: free h_far = {r2.x[0]:.4f}, A = {r2.x[1]:+.4f}, d chi2 = {r2.fun - X0:+.1f};  h_far = 0.685: A = {c2.x:+.4f}, '
          f'd chi2 = {c2.fun - X0:+.1f} -> {"PASS" if c2.fun - X0 <= 4 else "FAIL"}', flush=True)

print('\n=== V3: extra early share at h = 0.7354 (nu1, F = 1.14)')
b, sb, Om = b_sigma(Om_gamma(hs) * NU)
for btag, obh2, eob in BR:
    for tag, key, hi in (('(a) light-like u^-4, Or*(1+f)', 'f', 1.5), ('(b) u^-6, q = ratio to light at 1+z = 1101', 'q', 3.0)):
        fun = lambda x: model(hs, obh2, Om, NU, F, **{key: x})['theta'] - THETA
        fa, fb = fun(0.0), fun(hi)
        if fa * fb > 0:
            print(f'   {btag} {tag}: no root in 0..{hi} (theta - data {fa:+.5f} .. {fb:+.5f})'); continue
        x = brentq(fun, 0.0, hi, xtol=1e-5); mm = model(hs, obh2, Om, NU, F, **{key: x})
        bm = 1 / mm['rsd']
        sh = abs(1 / model(hs + eh, obh2, fit_Om(Om_gamma(hs + eh) * NU)[0], NU, F, **{key: x})['rsd'] - bm)
        dev = (bm - b) / sqrt(sb**2 + sh**2)
        print(f'   {btag} {tag}: {key} = {x:.4f}; z* = {mm["zstar"]:.1f}, z_drag = {mm["zdrag"]:.1f}, '
              f'r_d = {mm["rsd"] * CKM / (100 * hs):.2f} Mpc; b model = {bm:.4f} vs BAO {b:.4f} +- {sb:.4f} (+- {sh:.4f} from h): '
              f'{dev:+.1f} sigma -> {"PASS" if abs(dev) <= 2 else "FAIL"}', flush=True)
