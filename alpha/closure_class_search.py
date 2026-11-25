"""Class 1 (junction ON the closure element): one count k for all six cases.  16.11.2026.

Theory: W-Structure (Hrachya Danielyan).  Prior broad search: patterns_search.py / patterns_report_2026-11-16.md.

Everything below the line "FIXED BEFORE FITTING" was written before any prediction was computed:
encoding of the six cases, the candidate definitions of k, the combination laws, the constant modes,
the scores and the null models.  The same function `evaluate` is applied to the data and to the nulls.

Target: c = (mass change) * alpha, in units of m_e, so that  delta m = c / alpha.
Output: closure_class_output.txt (all numbers quoted in the report).
"""
import itertools, math, sys
import numpy as np

ALPHA = 1 / 137.035999084
ME_MEV = 0.51099895
RNG = np.random.default_rng(161126)

M = {'e': 1.0, 'mu': 206.7682830, 'tau': 3477.23, 'p': 1836.15267343, 'pi+': 273.13211,
     'pi0': 264.14251, 'Ds': 1968.35 / ME_MEV}
HOST_OF = {'mu': 'e', 'tau': 'mu', 'p': 'e', 'pi+': 'mu', 'pi0': 'mu', 'Ds': 'tau'}
CASES = ['mu', 'tau', 'p', 'pi+', 'pi0', 'Ds']
C_OBS = np.array([(M[k] - M[HOST_OF[k]]) * ALPHA for k in CASES])

# =============================================================================================
# FIXED BEFORE FITTING
# =============================================================================================
# 1. Encoding of each junction (host H, added part X, closure element E).
#    eH   eights of the host passing through E (before the junction)
#    lH   loops of the host passing through / touching E (before)
#    eP   eights passing through E after the junction (a point is not an eight)
#    lP   loops passing through / touching E after the junction (a point = one loop)
#    aX   loops of the added part: eight -> 2, point -> 1
#    aNet loops that are new in the product (an eight sharing one loop with the host -> 1)
#    d    dimension of the closure at E after the junction (1, 2, 3)
#    opp  number of opposite-direction pairs of eights at E after the junction
#
#  mu : H = e,   X = eight,                      E = O (contracted eight)
#  tau: H = mu,  X = eight (perpendicular),      E = O
#  p  : H = e,   X = eight sharing one loop,     E = the common double loop (1D)
#       in e the future common loop touches 2 loops (itself, the other loop at O);
#       in p it touches 3 (A-P1, P1-P2, P2-C).  Both eights pass P1-P2 with the same orientation (-),
#       so no opposite pair.
#  pi+: H = mu,  X = point at O; the point takes the dimension of the cross (author) -> d = 2
#  pi0: H = mu,  X = the second place of the double point (one point); the cross becomes
#       opposite-direction -> opp = 1; no new loop is drawn (aNet = 0), the loops through the double
#       centre are the same four (lP = 4).  [encoding choice; alternative aNet = 1, lP = 5 in section S]
#  Ds : H = tau, X = point at O, d = 3 (candidate form, author)
ENC = {
    #        eH lH eP lP aX aNet d opp
    'mu':  dict(eH=1, lH=2, eP=2, lP=4, aX=2, aNet=2, d=2, opp=0),
    'tau': dict(eH=2, lH=4, eP=3, lP=6, aX=2, aNet=2, d=3, opp=0),
    'p':   dict(eH=1, lH=2, eP=2, lP=3, aX=2, aNet=1, d=1, opp=0),
    'pi+': dict(eH=2, lH=4, eP=2, lP=5, aX=1, aNet=1, d=2, opp=0),
    'pi0': dict(eH=2, lH=4, eP=2, lP=4, aX=1, aNet=0, d=2, opp=1),
    'Ds':  dict(eH=3, lH=6, eP=3, lP=7, aX=1, aNet=1, d=3, opp=0),
}
# forms NOT in the list (blind predictions; encoded before any comparison)
NEW = {
    # point attached at the common double loop of the proton: E 1D, two eights share it, 3 loops touch it
    'p + point at common loop': ('p', dict(eH=2, lH=3, eP=2, lP=4, aX=1, aNet=1, d=1, opp=0)),
    # third eight at O of the muon cross, perpendicular, opposite in direction to both eights of the cross
    'mu + opposite third eight': ('mu', dict(eH=2, lH=4, eP=3, lP=6, aX=2, aNet=2, d=3, opp=2)),
}

# 2. Candidate definitions of ONE count k (each applied to all cases identically)
KDEF = {
    'D0  k = 1 (dimension only)':             lambda s: 1,
    'D1  k = eH (host eights through E)':     lambda s: s['eH'],
    'D2  k = lH (host loops through E)':      lambda s: s['lH'],
    'D3  k = eP (eights through E after)':    lambda s: s['eP'],
    'D4  k = lP (loops through E after)':     lambda s: s['lP'],
    'D5  k = aX (loops of added part)':       lambda s: s['aX'],
    'D6  k = aNet (new loops)':               lambda s: s['aNet'],
    'D7  k = d - 1 (index of new plane)':     lambda s: s['d'] - 1,
    'D8  k = d (planes involved)':            lambda s: s['d'],
    'D9  k = eH * aX/2 (point = half eight)': lambda s: s['eH'] * s['aX'] / 2,
    'D10 k = lH * aX/2':                      lambda s: s['lH'] * s['aX'] / 2,
    'D11 k = eH - opp (direction)':           lambda s: s['eH'] - s['opp'],
}

# 3. Combination law "multiply over planes": delta c = C * k^p(d)
#    A (author): 1D -> k^2, 2D -> k^4, 3D -> the new plane adds its own k^4 (sum over planes; tau - e = 1^4 + 2^4)
#    B (rejected alternative, counted as a trial): 3D -> k^6 (product over three planes)
LAWS = {'A: p = 2,4,4': {1: 2, 2: 4, 3: 4}, 'B: p = 2,4,6': {1: 2, 2: 4, 3: 6}}

# 4. Constant: 'fit_mu' = C fitted on the muon only (predicts the other 5); 'C=3/2' = fixed a priori (predicts 6)
CMODES = ['fit_mu', 'C=3/2']

# 5. Scores (computed over the predicted cases): hits within 1 %, 3 %, 10 %; max relative error.
TOLS = (0.01, 0.03, 0.10)
# 6. Nulls: each c independently x exp(U(-ln2, ln2)) ['x2'], or x (1 + U(-0.1, 0.1)) ['10%'].
#    The statistic is the BEST over all variants (definitions x laws x constant modes), so the
#    look-elsewhere effect over the variants is included in P.
# =============================================================================================

def hvals(enc, kf, law):
    return np.array([float(kf(enc[k])) ** law[enc[k]['d']] for k in CASES])

def build_variants(enc=ENC):
    V = []
    for (kn, kf), (ln, law), cm in itertools.product(KDEF.items(), LAWS.items(), CMODES):
        V.append((kn, ln, cm, hvals(enc, kf, law)))
    return V

def predict(h, cm, c):
    """predicted c for all six cases; mask of predicted (tested) cases"""
    if cm == 'fit_mu':
        C = c[0] / h[0] if h[0] != 0 else np.nan
        mask = np.array([False, True, True, True, True, True])
    else:
        C = 1.5
        mask = np.ones(6, bool)
    return C * h, mask, C

def scores(h, cm, c):
    pr, mask, _ = predict(h, cm, c)
    err = np.abs(pr[mask] / c[mask] - 1)
    err = np.where(np.isfinite(err), err, np.inf)
    return [int((err <= t).sum()) for t in TOLS], err.max()

def evaluate(V, c):
    best_hits = [0, 0, 0]; best_max = np.inf
    for _, _, cm, h in V:
        hits, mx = scores(h, cm, c)
        best_hits = [max(a, b) for a, b in zip(best_hits, hits)]
        best_max = min(best_max, mx)
    return best_hits, best_max

def null(c, kind, rng=RNG):
    if kind == 'x2':
        return c * np.exp(rng.uniform(-math.log(2), math.log(2), size=c.shape))
    return c * (1 + rng.uniform(-0.1, 0.1, size=c.shape))

# secondary (declared in advance): points only.  C_pt fitted on pi+, predicts pi0 and Ds.
PT = ['pi+', 'pi0', 'Ds']
IPT = [CASES.index(k) for k in PT]

def evaluate_pt(Vpt, cpt):
    best_hits = [0, 0, 0]; best_max = np.inf
    for _, _, h in Vpt:
        if h[0] == 0: continue
        pr = cpt[0] / h[0] * h
        err = np.abs(pr[1:] / cpt[1:] - 1)
        best_hits = [max(a, int((err <= t).sum())) for a, t in zip(best_hits, TOLS)]
        best_max = min(best_max, err.max())
    return best_hits, best_max

def main(nrep=20000):
    out = []
    pr = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))
    pr('== 0. Targets c = delta(m/m_e) * alpha ==')
    for k, c in zip(CASES, C_OBS):
        pr(f'  {k:4s} (rel. {HOST_OF[k]:3s}) delta = {M[k]-M[HOST_OF[k]]:10.4f} m_e   c = {c:.5f}   c/1.5 = {c/1.5:.5f}')

    pr('\n== 1. Encoding and k values ==')
    pr('  ' + f'{"definition":40s}' + ''.join(f'{k:>6s}' for k in CASES))
    for kn, kf in KDEF.items():
        pr('  ' + f'{kn:40s}' + ''.join(f'{kf(ENC[k]):6.3g}' for k in CASES))

    V = build_variants()
    pr(f'\n== 2. Predictions: {len(KDEF)} definitions x {len(LAWS)} laws x {len(CMODES)} constant modes = {len(V)} variants ==')
    pr('  predicted c (error %); for fit_mu the muon is the fit point')
    rows = []
    for kn, ln, cm, h in V:
        p, mask, C = predict(h, cm, C_OBS)
        err = p / C_OBS - 1
        hits, mx = scores(h, cm, C_OBS)
        rows.append((kn, ln, cm, C, p, err, hits, mx))
        cells = ' '.join(f'{k} {pv:8.3f}({100*e:+7.1f})' if m else f'{k} {pv:8.3f}(  fit )'
                         for k, pv, e, m in zip(CASES, p, err, mask))
        pr(f'  {kn[:4]} {ln[:1]} {cm:6s} C={C:8.4f} | {cells} | hits1/3/10 {hits} max {100*mx:.1f}%')
    obs_hits, obs_max = evaluate(V, C_OBS)
    pr(f'\n  BEST over variants: hits within 1/3/10 % = {obs_hits}; minimal max-error = {100*obs_max:.1f} %')
    pr('  variants achieving the best hit counts:')
    for t_i, t in enumerate(TOLS):
        who = [f'{r[0][:4]}{r[1][:1]}/{r[2]}' for r in rows if r[6][t_i] == obs_hits[t_i]]
        pr(f'    within {100*t:.0f} %: {obs_hits[t_i]} hits by {len(who)} variants: {", ".join(who[:12])}' + (' ...' if len(who) > 12 else ''))
    # which case is hit
    pr('  hits within 3 % by case:')
    for r in rows:
        p, mask = r[4], predict(np.ones(6), r[2], C_OBS)[1]
        hit = [k for k, e, m in zip(CASES, r[5], mask) if m and abs(e) <= 0.03]
        if hit: pr(f'    {r[0][:4]} {r[1][:1]} {r[2]:6s}: {hit}')

    # what is beyond Barut: among variants that reproduce tau within 1 % (Barut part), hits on p, pi+, pi0, Ds
    inew = [CASES.index(k) for k in ('p', 'pi+', 'pi0', 'Ds')]
    pr('  beyond Barut: variants with tau within 1 % and their hits on (p, pi+, pi0, Ds) within 10 %:')
    for r in rows:
        if abs(r[5][1]) <= 0.01:
            nh = int((np.abs(r[5][inew]) <= 0.10).sum())
            pr(f'    {r[0][:4]} {r[1][:1]} {r[2]:6s}: {nh} of 4; errors ' +
               ', '.join(f'{CASES[i]} {100*r[5][i]:+.0f} %' for i in inew))

    pr(f'\n== 3. Significance: identical procedure on null targets, {nrep} replicates each ==')
    for kind in ('x2', '10%'):
        H = np.zeros((nrep, 3), int); MX = np.zeros(nrep)
        for r in range(nrep):
            hh, mm = evaluate(V, null(C_OBS, kind)); H[r] = hh; MX[r] = mm
        ps = [(H[:, i] >= obs_hits[i]).mean() for i in range(3)]
        pr(f'  null {kind:4s}: P(best hits >= observed) at 1/3/10 %: ' + ', '.join(f'{p:.4f}' for p in ps) +
           f';  P(best max-error <= {100*obs_max:.1f} %) = {(MX <= obs_max + 1e-12).mean():.4f}')
        pr(f'        null distribution of best hits@1%: ' +
           ', '.join(f'{v}: {(H[:,0]==v).mean():.3f}' for v in range(0, 4)))

    # single-variant reference: the Barut-consistent variant D1/A/fit_mu, tau only
    i_t = CASES.index('tau')
    h1 = hvals(ENC, KDEF['D1  k = eH (host eights through E)'], LAWS['A: p = 2,4,4'])
    e_tau = abs(C_OBS[0] / h1[0] * h1[i_t] / C_OBS[i_t] - 1)
    for kind in ('x2', '10%'):
        z = np.array([abs(C_OBS[0] / h1[0] * h1[i_t] / null(C_OBS, kind)[i_t] - 1) for _ in range(nrep)])
        pr(f'  reference (one pre-chosen variant D1/A, tau only, err {100*e_tau:.2f} %): null {kind}: P = {(z <= e_tau).mean():.4f}')

    pr('\n== 4. Secondary (declared in advance): points only, C_pt fitted on pi+, predict pi0 and Ds ==')
    Vpt = [(kn, ln, hvals(ENC, kf, law)[IPT]) for (kn, kf), (ln, law) in itertools.product(KDEF.items(), LAWS.items())]
    cpt = C_OBS[IPT]
    pr(f'  observed ratios: pi0/pi+ = {cpt[1]/cpt[0]:.4f}, Ds/pi+ = {cpt[2]/cpt[0]:.4f}')
    for kn, ln, h in Vpt:
        if h[0] == 0:
            pr(f'  {kn[:4]} {ln[:1]}: h(pi+) = 0, no fit'); continue
        r = h / h[0]
        pr(f'  {kn[:4]} {ln[:1]}: predicted ratios pi0/pi+ = {r[1]:.4f} ({100*(r[1]/(cpt[1]/cpt[0])-1):+.1f} %), '
           f'Ds/pi+ = {r[2]:.4f} ({100*(r[2]/(cpt[2]/cpt[0])-1):+.1f} %)')
    ohp, omp = evaluate_pt(Vpt, cpt)
    pr(f'  BEST: hits within 1/3/10 % = {ohp}; minimal max-error {100*omp:.1f} %  ({len(Vpt)} variants)')
    for kind in ('x2', '10%'):
        H = np.zeros((nrep, 3), int); MX = np.zeros(nrep)
        for r in range(nrep):
            hh, mm = evaluate_pt(Vpt, null(cpt, kind)); H[r] = hh; MX[r] = mm
        ps = [(H[:, i] >= ohp[i]).mean() for i in range(3)]
        pr(f'  null {kind:4s}: P(best hits >= observed) at 1/3/10 %: ' + ', '.join(f'{p:.4f}' for p in ps) +
           f';  P(best max-error <= {100*omp:.1f} %) = {(MX <= omp + 1e-12).mean():.4f}')

    pr('\n== 5. Eights only (mu, tau, p): is there one definition giving tau AND p? ==')
    for kn, ln, cm, C, p, err, hits, mx in rows:
        e2 = [err[CASES.index('tau')], err[CASES.index('p')]]
        if max(abs(x) for x in e2) < 0.10:
            pr(f'  {kn[:4]} {ln[:1]} {cm}: tau {100*e2[0]:+.2f} %, p {100*e2[1]:+.2f} %')
    best_tp = min(max(abs(r[5][CASES.index('tau')]), abs(r[5][CASES.index('p')])) for r in rows)
    pr(f'  best max(|err tau|, |err p|) over all variants: {100*best_tp:.1f} %')
    pr('  required k (law A, fit on mu with k_mu = 1): tau k^4 = c_tau/c_mu = '
       f'{C_OBS[1]/C_OBS[0]:.3f} -> k = {(C_OBS[1]/C_OBS[0])**0.25:.4f};  p k^2 = c_p/c_mu = {C_OBS[2]/C_OBS[0]:.3f} -> k = {(C_OBS[2]/C_OBS[0])**0.5:.4f}')
    pr('  required k for points (law A, C from mu): ' + ', '.join(
        f'{k}: k^{LAWS["A: p = 2,4,4"][ENC[k]["d"]]} = {C_OBS[CASES.index(k)]/C_OBS[0]:.4f} -> k = '
        f'{(C_OBS[CASES.index(k)]/C_OBS[0])**(1/LAWS["A: p = 2,4,4"][ENC[k]["d"]]):.4f}' for k in ('pi+', 'pi0', 'Ds')))

    pr('\n== S. Sensitivity: alternative encodings ==')
    alts = {
        'pi0: aNet=1, lP=5 (double point drawn as a new point with a loop)': ('pi0', dict(aNet=1, lP=5)),
        'p: lH=1 (only the loop itself before)': ('p', dict(lH=1)),
        'p: lP=4 (passages: common loop twice + two adjacent loops)': ('p', dict(lP=4)),
        'p: opp=1 (loop orientations alternate +,-,+)': ('p', dict(opp=1)),
        'pi+: d=3 (point adds an axis)': ('pi+', dict(d=3)),
    }
    for nm, (case, ch) in alts.items():
        enc = {k: dict(v) for k, v in ENC.items()}; enc[case].update(ch)
        Va = build_variants(enc)
        hh, mm = evaluate(Va, C_OBS)
        pr(f'  {nm}: best hits 1/3/10 % = {hh}, minimal max-error {100*mm:.1f} %')

    pr('\n== 6. Formal predictions for forms not in the list (written before any lookup; NOT tests) ==')
    pr('  law A, C fitted on mu; value = predicted mass change in m_e and resulting mass in MeV')
    for nm, (host, s) in NEW.items():
        pr(f'  {nm} (host {host}, {M[host]*ME_MEV:.2f} MeV):')
        for kn, kf in KDEF.items():
            hmu = float(kf(ENC['mu'])) ** LAWS['A: p = 2,4,4'][2]
            if hmu == 0: continue
            C = C_OBS[0] / hmu
            dc = C * float(kf(s)) ** LAWS['A: p = 2,4,4'][s['d']]
            dm = dc / ALPHA
            pr(f'    {kn:40s} k={kf(s):5.3g}  delta = {dm:9.1f} m_e  -> m = {(M[host]+dm)*ME_MEV:9.1f} MeV')
    return out

if __name__ == '__main__':
    nrep = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    out = main(nrep)
    with open('/home/claude/alpha/closure_class_output.txt', 'w') as fh:
        fh.write('\n'.join(out) + '\n')
