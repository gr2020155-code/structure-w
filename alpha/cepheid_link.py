"""Link 'geometry -> Cepheids' checked with TRGB / JAGB (R-138). Rules: cepheid_link_rules_2026-11-27.md.
usage: python -P -E -X utf8 cepheid_link.py <cchp_trgb_jagb_2024.txt> <dir with Pantheon+ files>"""
import sys
import numpy as np
from math import sqrt, log
from scipy.integrate import quad

OM, CKM = 0.297, 299792.458
L = open(sys.argv[1]).read().split('\n')
def table(tag):
    i = L.index(tag); out = []
    for l in L[i + 1:]:
        if not l.strip() or l.startswith('#') or l.strip() in ('TOT', 'JWST'): break
        p = l.split(); out.append((p[0], p[1], *map(float, p[2:])))
    return out
TOT = {r[0]: r for r in table('TOT')}; JW = table('JWST')

snd = sys.argv[2]
rows = open(f'{snd}/PantheonPlusSH0ES.dat').read().split('\n'); hdr = rows[0].split(); data = [r.split() for r in rows[1:] if r.strip()]
col = lambda nm: np.array([float(r[hdr.index(nm)]) for r in data])
cid = np.array([r[0] for r in data]); z, m, mu_c, raw = col('zHD'), col('m_b_corr'), col('CEPH_DIST'), col('m_b_corr_err_RAW')
cal, hf = col('IS_CALIBRATOR') == 1, col('USED_IN_SH0ES_HF') == 1
cv = np.loadtxt(f'{snd}/STATSYS.cov'); n = int(cv[0]); C_all = cv[1:].reshape(n, n)
E = lambda t: sqrt(OM * (1 + t)**3 + 1 - OM)
Dc = np.array([quad(lambda t: 1 / E(t), 0, zz)[0] if h_ else 0.0 for zz, h_ in zip(z, hf)])

def wmean(d, s):
    w = 1 / np.asarray(s)**2; return float(np.sum(w * d) / np.sum(w)), float(1 / sqrt(np.sum(w)))

print('needed shift: 5 lg(0.735/0.685) = %.3f mag' % (5 * np.log10(0.735 / 0.685)))
# ---- T1 direct comparison, one value per host
hosts = {}
for sn, r in TOT.items():
    if r[4] > 0: hosts.setdefault(r[1], (sn, r))
dT_R22, sR22, dT_P, sP = [], [], [], []
for host, (sn, r) in hosts.items():
    dT_R22.append(r[4] - r[2]); sR22.append(sqrt(r[3]**2 + r[5]**2))
    I = np.where((cid == sn) & cal)[0]
    if len(I) == 0:
        I = np.where(np.isin(cid, [s for s, rr in TOT.items() if rr[1] == host]) & cal)[0]
    if len(I): dT_P.append(mu_c[I[0]] - r[2]); sP.append(sqrt(r[3]**2 + r[5]**2))
a, sa = wmean(dT_R22, sR22); b, sb = wmean(dT_P, sP)
print(f'== T1 Cepheid - TRGB, per host: R22 column: {a:+.3f} +- {sa:.3f} mag ({len(dT_R22)} hosts, {a / sa:+.1f} sigma); '
      f'Pantheon+ CEPH_DIST: {b:+.3f} +- {sb:.3f} ({len(dT_P)} hosts, {b / sb:+.1f} sigma); unweighted R22 {np.mean(dT_R22):+.3f}, rms {np.std(dT_R22):.3f}')
jh = {}
for sn, host, mt, st, mj, sj in JW:
    if mj > 0: jh.setdefault(host, (sn, mt, st, mj, sj))
dJC, sJC, dJT, sJT = [], [], [], []
for host, (sn, mt, st, mj, sj) in jh.items():
    r = TOT[sn]
    if r[4] > 0: dJC.append(r[4] - mj); sJC.append(sqrt(r[5]**2 + sj**2))
    dJT.append(mj - mt); sJT.append(sqrt(sj**2 + st**2))
a, sa = wmean(dJC, sJC); b, sb = wmean(dJT, sJT)
print(f'== T1J Cepheid(R22) - JAGB: {a:+.3f} +- {sa:.3f} mag ({len(dJC)} hosts, {a / sa:+.1f} sigma); JAGB - TRGB(JWST): {b:+.3f} +- {sb:.3f} ({len(dJT)})')

# ---- fits
def fit(rows_cal, mu_cal, rows_hf, extra=None):
    idx = np.concatenate([rows_cal, rows_hf]); nc = len(rows_cal)
    Cm = C_all[np.ix_(idx, idx)].copy()
    if extra is not None:
        k = extra[0]; Cm[np.ix_(k, k)] = 0; Cm[np.ix_(k, np.arange(len(idx)))] = 0; Cm[np.ix_(np.arange(len(idx)), k)] = 0
        Cm[np.ix_(k, k)] = extra[1]
    Ci = np.linalg.inv(Cm)
    y = np.concatenate([m[rows_cal] - mu_cal, m[rows_hf] - 5 * np.log10((1 + z[rows_hf]) * Dc[rows_hf] * CKM / 100) - 25])
    X = np.stack([np.ones(len(idx)), np.r_[np.zeros(nc), -np.ones(len(rows_hf))]], 1)
    Fm = X.T @ Ci @ X; bb = np.linalg.solve(Fm, X.T @ Ci @ y); cov = np.linalg.inv(Fm)
    h = 10**(bb[1] / 5); return h, h * log(10) / 5 * sqrt(cov[1, 1]), bb[0]
def show(tag, r):
    print(f'   {tag}: h = {r[0]:.4f} +- {r[1]:.4f} (vs 0.685: {(r[0] - 0.685) / r[1]:+.1f} sigma), M = {r[2]:.3f}'); return r

HF = np.where(hf)[0]
both = [i for i in np.where(cal)[0] if cid[i] in TOT and TOT[cid[i]][4] > 0]
both = np.array(both); mu_T = np.array([TOT[cid[i]][2] for i in both])
print(f'== T2 same SNe ({len(both)} rows, {len(set(cid[both]))} SNe, {len(set(TOT[c][1] for c in cid[both]))} hosts) + {len(HF)} flow')
rC = show('Cepheids (CEPH_DIST)', fit(both, mu_c[both], HF)); rT = show('TRGB (CCHP)', fit(both, mu_T, HF))
print(f'   TRGB / Cepheid h ratio = {rT[0] / rC[0]:.4f}')

def t3(names, tag):
    new = np.array([i for i in range(len(cid)) if cid[i] in names and not cal[i]])
    rc = np.concatenate([both, new]); mu = np.concatenate([mu_T, [TOT[cid[i]][2] for i in new]])
    k = np.arange(len(both), len(both) + len(new)); B = np.zeros((len(new), len(new)))
    for a_, i in enumerate(new):
        for b_, j in enumerate(new):
            if TOT[cid[i]][1] == TOT[cid[j]][1]: B[a_, b_] += TOT[cid[i]][3]**2
            if cid[i] == cid[j]: B[a_, b_] += 0.10**2
            if a_ == b_: B[a_, b_] += raw[i]**2
    print(f'   added rows: {", ".join(cid[new])}')
    return show(tag, fit(rc, mu, HF, (k, B)))
print('== T3 full TRGB sample')
t3({'2006dd', '1980N', '1981D', '2011iv'}, 'TRGB, without 2007on (as CCHP)')
t3({'2006dd', '1980N', '1981D', '2011iv', '2007on'}, 'TRGB, with 2007on')
print('== T4 sample effect, Cepheids')
show('Cepheids, SNe of T2 only', rC); show('Cepheids, all 77 calibrators', fit(np.where(cal)[0], mu_c[cal], HF))
