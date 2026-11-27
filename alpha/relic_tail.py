"""Full relic spectrum and tail, two branches of pre-material links (R-141). Rules: relic_tail_rules_2026-11-27.md.
usage: python -P -E -X utf8 relic_tail.py <dir with Planck binned spectra>"""
import sys
import numpy as np
import camb
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2 as chi2dist

d = sys.argv[1]
DAT = {k: np.loadtxt(f'{d}/COM_PowerSpect_CMB-{k}-binned_{v}.txt') for k, v in (('TT', 'R3.01'), ('TE', 'R3.02'), ('EE', 'R3.02'))}
THETA100 = 1.0411; OM = 0.297; OBH2_0 = 0.02235; TAU = 0.054; LMAX = 2600

def obh2_of(dN): return OBH2_0 * ((3 + dN) / 3)**(0.396 / 1.597)
def params(branch, h, dN, ns=0.965):
    ob = obh2_of(dN); p = camb.CAMBparams()
    oc = OM * h * h - ob if branch == 'K' else 1e-5
    p.set_cosmology(H0=100 * h, ombh2=ob, omch2=oc, mnu=0.0, nnu=3.046 + dN, tau=TAU)
    if branch == 'S':
        Ob = ob / h**2; Op = OM - Ob; OB = 1 - OM
        a = np.logspace(-9, 0, 800); p.DarkEnergy = camb.dark_energy.DarkEnergyFluid(cs2=1.0)
        p.DarkEnergy.set_w_a_table(a, -OB / (Op * a**-3 + OB))
    p.InitPower.set_params(As=2.1e-9, ns=ns); p.set_for_lmax(LMAX); p.NonLinear = camb.model.NonLinear_none
    return p
def theta100(branch, h, dN):
    p = params(branch, h, dN); p.WantCls = False
    return camb.get_background(p).get_derived_params()['thetastar']
def spectra(branch, h, dN, ns):
    r = camb.get_results(params(branch, h, dN, ns)); cl = r.get_cmb_power_spectra(CMB_unit='muK', raw_cl=False)['total']
    return cl                                                # columns TT EE BB TE, D_l
def model_bins(cl):
    out = {}
    for k, c in (('TT', 0), ('EE', 1), ('TE', 3)):
        ell = DAT[k][:, 0]; out[k] = np.interp(ell, np.arange(cl.shape[0]), cl[:, c])   # D_l at bin centre
    return out
def chi2_fit(branch, h, dN):
    def at(ns):
        mb = model_bins(spectra(branch, h, dN, ns))
        num = sum(np.sum(mb[k] * DAT[k][:, 1] / DAT[k][:, 2]**2) for k in DAT); den = sum(np.sum((mb[k] / DAT[k][:, 2])**2) for k in DAT)
        s = num / den
        parts = {k: float(np.sum(((DAT[k][:, 1] - s * mb[k]) / DAT[k][:, 2])**2)) for k in DAT}
        tail = DAT['TT'][:, 0] > 1000
        parts['TTtail'] = float(np.sum(((DAT['TT'][tail, 1] - s * mb['TT'][tail]) / DAT['TT'][tail, 2])**2)); parts['ntail'] = int(tail.sum())
        return sum(parts[k] for k in DAT), s, parts
    r = minimize_scalar(lambda ns: at(ns)[0], bounds=(0.85, 1.10), method='bounded', options=dict(xatol=1e-4))
    c, s, parts = at(r.x); return c, r.x, s * 2.1e-9, parts

print('Planck BestFit column, same diagonal errors: chi2 = %.1f' % sum(np.sum(((DAT[k][:, 1] - DAT[k][:, 4]) / DAT[k][:, 2])**2) for k in DAT))
NDOF = 215 - 2
for branch, bname in (('K', 'K: links cluster (as cold matter)'), ('S', 'S: links smooth (cs2 = 1)')):
    print(f'\n=== {bname}')
    hH = brentq(lambda h: theta100(branch, h, 0.0) - THETA100, 0.45, 0.95, xtol=1e-5)
    dF = {}
    for tag, h in (('F1', 0.7354), ('F2', 0.710)):
        fa = theta100(branch, h, 0.0) - THETA100
        dF[tag] = brentq(lambda x: theta100(branch, h, x) - THETA100, 0.0, 4.0, xtol=1e-4) if fa > 0 else None
    subs = [('H', hH, 0.0), ('N0', 0.7354, 0.0)] + [(t, h, dF[t]) for t, h in (('F1', 0.7354), ('F2', 0.710)) if dF[t] is not None]
    for t in ('F1', 'F2'):
        if dF[t] is None: print(f'   {t}: no dN >= 0 gives exact theta')
    for tag, h, dN in subs:
        c, ns, As, parts = chi2_fit(branch, h, dN); p = chi2dist.sf(c, NDOF)
        print(f'   {tag}: h = {h:.4f}, dN = {dN:.3f}, Obh2 = {obh2_of(dN):.5f}, 100 theta* = {theta100(branch, h, dN):.5f}; '
              f'best ns = {ns:.4f}, As = {As:.3e}; chi2 = {c:.1f} / {NDOF} (TT {parts["TT"]:.1f}, TE {parts["TE"]:.1f}, EE {parts["EE"]:.1f}; '
              f'TT tail l>1000 {parts["TTtail"]:.1f} / {parts["ntail"]}), prob {p:.2e} -> {"PASS" if p >= 0.05 else "FAIL"}', flush=True)
