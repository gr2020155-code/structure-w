"""Pattern search: masses of W-forms from the STATES of their points (16.11.2026).

Written BEFORE the fitting runs.  All choices below (encoding, features, coefficient set,
families, nulls, thresholds) are fixed here; the search and the null tests use the SAME code.

Theory: W-Structure (Hrachya Danielyan).  Forms = author's positions (see report, section 1).
Rule already found by the author ("multiply over planes"):
    m/m_e = 1 + (3/2)(1/alpha) X,  X_mu = 1^4 = 1,  X_tau = 1^4 + 2^4 = 17,  X_p = 3^2 = 9.

Outputs: patterns_search_output.txt (all numbers quoted in the report).
"""
import itertools, math, sys, time, json
from fractions import Fraction
import numpy as np
from scipy.spatial import cKDTree

ALPHA = 1 / 137.035999084
ME_MEV = 0.51099895
RNG = np.random.default_rng(20261116)

# ---------------------------------------------------------------------------------------------
# 0. Data (m / m_e).  Fitting targets: mu, tau, p, n, pi+, pi0.  Held out: D_s and all others.
# ---------------------------------------------------------------------------------------------
MASS = {'e': 1.0, 'mu': 206.7682830, 'tau': 3477.23, 'p': 1836.15267343, 'n': 1838.68366173,
        'pi+': 273.13211, 'pi0': 264.14251}
HELD_OUT_MASS = {'Ds': 1968.35 / ME_MEV}
TARGETS6 = ['mu', 'tau', 'p', 'n', 'pi+', 'pi0']
TARGETS3 = ['mu', 'tau', 'p']

# ---------------------------------------------------------------------------------------------
# 1. Encoding: every point of a form with its STATE.
#    n    = passages of loops through the point per full cycle (a shared double loop counts 2,
#           an attached point counts 1 per exchange with its host + 1 per loop through it)
#    d    = dimensionality of closure at the point (number of independent axes of its loops/links)
#    role = end | closure | shared | attached | double
#    rel  = relation of the eights closing at the point: 'opp' / 'same' / None
#    carries = number of eights threaded through an attached point
#    host_d  = dimensionality of the junction an attached point is attached to
# Junction elements (where eights meet or something is attached) carry the same state variables.
# ---------------------------------------------------------------------------------------------
def P(n, d, role, rel=None, carries=0, host_d=None):
    return dict(n=n, d=d, role=role, rel=rel, carries=carries, host_d=host_d)

def J(kind, d, loops, eights, pas, rel=None, host_d=None, carries=0):
    return dict(kind=kind, d=d, loops=loops, eights=eights, pas=pas, rel=rel, host_d=host_d, carries=carries)

END = P(1, 1, 'end')

def form_e():
    return dict(points=[END, P(2, 1, 'closure', 'opp'), END], loops=2, eights=1, junctions=[])

def form_mu():
    return dict(points=[END] * 4 + [P(4, 2, 'closure', 'same')], loops=4, eights=2,
                junctions=[J('cross', 2, 4, 2, 4, 'same')])

def form_tau():
    return dict(points=[END] * 6 + [P(6, 3, 'closure', 'same')], loops=6, eights=3,
                junctions=[J('cross', 3, 6, 3, 6, 'same')])

def form_p():
    return dict(points=[END, P(3, 1, 'shared', 'opp'), P(3, 1, 'shared', 'opp'), END], loops=3, eights=2,
                junctions=[J('shared', 1, 3, 2, 3, 'opp')])

def add_point(f, host_d):
    """an attached point (neutrino-like): one exchange per cycle with its host, takes host's dimension"""
    g = {k: (list(v) if isinstance(v, list) else v) for k, v in f.items()}
    g['points'] = g['points'] + [P(1, host_d, 'attached', host_d=host_d)]
    g['junctions'] = g['junctions'] + [J('attach', host_d, 0, 0, 1, host_d=host_d)]
    return g

def add_point_with_eight(f, host_d):
    """attached point N with an electron-type eight threaded through it, in a perpendicular plane:
       N: exchange with host + one loop passage -> n = 2; d = host_d + 1 (max 3); the eight brings 3 points."""
    g = {k: (list(v) if isinstance(v, list) else v) for k, v in f.items()}
    dN = min(3, host_d + 1)
    g['points'] = g['points'] + [P(2, dN, 'attached', carries=1, host_d=host_d), END, P(2, 1, 'closure', 'opp'), END]
    g['loops'] += 2; g['eights'] += 1
    g['junctions'] = g['junctions'] + [J('attach', dN, 1, 1, 2, host_d=host_d, carries=1)]
    return g

def form_pi0():
    # cross of OPPOSITE-direction eights; the centre is a double point: two coincident places,
    # each passed by one eight's two loops (n = 2), lying in the plane of the cross (d = 2)
    return dict(points=[END] * 4 + [P(2, 2, 'double', 'opp'), P(2, 2, 'double', 'opp')], loops=4, eights=2,
                junctions=[J('cross', 2, 4, 2, 2, 'opp')])

FORMS = {
    'e': form_e(), 'mu': form_mu(), 'tau': form_tau(), 'p': form_p(),
    'n': add_point_with_eight(form_p(), 1),          # proton + eight through a point, perpendicular
    'pi+': add_point(form_mu(), 2),                  # same-direction cross + point
    'pi0': form_pi0(),                               # opposite-direction cross, double point
}
# forms NOT used in fitting (blind predictions); built by the same constructors
NEW_FORMS = {
    'tau+pt (D_s?)': add_point(form_tau(), 3),
    'e+pt': add_point(form_e(), 1),
    'p+pt': add_point(form_p(), 1),
    'pi0+pt': add_point(form_pi0(), 2),
    'e+8 (via pt)': add_point_with_eight(form_e(), 1),
    'mu+8 (via pt)': add_point_with_eight(form_mu(), 2),
    'tau+8 (via pt)': add_point_with_eight(form_tau(), 3),
}

# ---------------------------------------------------------------------------------------------
# 2. Base features (integers, read off the point states).  Each is taken as EXCESS over the
#    electron (f(form) - f(e)), so the electron gives exactly m_e for every rule.
# ---------------------------------------------------------------------------------------------
def barut(k, q):
    return sum(j ** q for j in range(1, k + 1))

def raw_features(f):
    pts, js = f['points'], f['junctions']
    cross = [j for j in js if j['kind'] == 'cross']
    shared = [j for j in js if j['kind'] == 'shared']
    att = [j for j in js if j['kind'] == 'attach']
    Jc = sum(j['eights'] - 1 for j in cross)
    F = {
        'Npts': len(pts), 'Nloops': f['loops'], 'Neights': f['eights'],
        'D': max(p['d'] for p in pts),
        'Jc': Jc,                                            # order of crossing (eights added at a centre)
        'Barut4': barut(Jc, 4), 'Barut2': barut(Jc, 2),
        'Js': len(shared), 'Ksh': sum(j['loops'] for j in shared),
        'Natt': len(att), 'AttD': sum(j['d'] for j in att), 'AttHostD': sum(j['host_d'] for j in att),
        'AttCarry': sum(j['carries'] for j in att),
        'Ndbl': sum(1 for p in pts if p['role'] == 'double'),
        'Nopp': sum(1 for j in cross if j['rel'] == 'opp'),
        'Pmax': max(p['n'] for p in pts),
        'SumN2': sum(p['n'] ** 2 for p in pts),
        'SumN2D': sum(p['n'] ** (2 * p['d']) for p in pts),
        'Planes': sum(math.comb(j['d'], 2) for j in js if j['kind'] in ('cross',)),
        # junction powers, "multiply over planes": contribution k^(2d) of each junction element
        'JPloops': sum(j['loops'] ** (2 * j['d']) for j in js),
        'JPeights': sum((j['eights'] - 1) ** (2 * j['d']) for j in js if j['eights'] > 0),
        'JPpas': sum(j['pas'] ** (2 * j['d']) for j in js),
        # the author's rule written as ONE feature: cross -> sum over planes of j^4 (j = 1..d-1),
        # shared -> loops^2, attached point -> 0 (variant a) or 1^(2d) = 1 (variant b)
        'AuthA': sum(barut(j['d'] - 1, 4) for j in cross) + sum(j['loops'] ** 2 for j in shared),
        'AuthB': sum(barut(j['d'] - 1, 4) for j in cross) + sum(j['loops'] ** 2 for j in shared) + len(att),
    }
    return F

BASE_NAMES = list(raw_features(FORMS['e']).keys())

def base_vector(f):
    r, r0 = raw_features(f), raw_features(FORMS['e'])
    return np.array([r[k] - r0[k] for k in BASE_NAMES], float)

def feature_table(forms):
    """F = {b^p, p = 1..4} U {b_i * b_j, i < j} over base features (excess over e)."""
    B = {k: base_vector(f) for k, f in forms.items()}
    names, cols = [], []
    nb = len(BASE_NAMES)
    for i in range(nb):
        for p in (1, 2, 3, 4):
            names.append(f'{BASE_NAMES[i]}^{p}' if p > 1 else BASE_NAMES[i])
            cols.append({k: B[k][i] ** p for k in forms})
    for i, j in itertools.combinations(range(nb), 2):
        names.append(f'{BASE_NAMES[i]}*{BASE_NAMES[j]}')
        cols.append({k: B[k][i] * B[k][j] for k in forms})
    return names, cols

# coefficient set: rationals p/q, 1 <= p, q <= 6, both signs; powers of 1/alpha: s in {0, 1/2, 1}
RATS = sorted({Fraction(p, q) for p in range(1, 7) for q in range(1, 7)})
S_POW = [0.0, 0.5, 1.0]

def build_terms(forms_fit, forms_extra):
    allforms = dict(forms_fit); allforms.update(forms_extra)
    fnames, fcols = feature_table(allforms)
    keep = [i for i, c in enumerate(fcols) if any(c[k] != 0 for k in forms_fit)]
    terms, vec, vec_extra = [], [], []
    for i in keep:
        v = np.array([fcols[i][k] for k in forms_fit])
        ve = np.array([fcols[i][k] for k in forms_extra])
        for r in RATS:
            for s in S_POW:
                for sg in (1, -1):
                    c = sg * float(r) * ALPHA ** (-s)
                    terms.append((fnames[i], sg * r, s)); vec.append(c * v); vec_extra.append(c * ve)
    return terms, np.array(vec), np.array(vec_extra), len(keep)

def term_str(t):
    name, r, s = t
    a = {0.0: '', 0.5: '/sqrt(a)', 1.0: '/a'}[s]
    return f'({r}){a}*[{name}]'

# ---------------------------------------------------------------------------------------------
# 3. Searches.  Score = max over targets of |relative error of m|.
# ---------------------------------------------------------------------------------------------
def search_additive(V, t, topk=0, L=None):
    """m = 1 + L + T1 + T2 (T2 may be empty; L a fixed vector, default 0).  Chebyshev NN in error space."""
    L = np.zeros(V.shape[1]) if L is None else L
    Y = np.vstack([V, np.zeros(V.shape[1])]) / t        # T2 / t
    tree = cKDTree(Y)
    Q = (t - 1 - L - V) / t
    Q = np.vstack([Q, (t - 1 - L) / t])
    k = 1 if not topk else 4
    d, idx = tree.query(Q, k=k, p=np.inf)
    if not topk:
        return d.min()
    res = []
    for i in range(Q.shape[0]):
        for j in range(k):
            res.append((d[i, j], i, idx[i, j]))
    res.sort()
    return res

def search_multiplicative(V, t, topk=0):
    """m = (1 + T1)(1 + T2) (T2 may be empty), log space Chebyshev (= relative error to first order)."""
    ok = (1 + V > 0.05).all(axis=1)
    Vk = V[ok]
    L = np.log(1 + Vk)
    Y = np.vstack([L, np.zeros(V.shape[1])])
    tree = cKDTree(Y)
    Q = np.vstack([np.log(t) - L, np.log(t)])
    k = 1 if not topk else 4
    d, idx = tree.query(Q, k=k, p=np.inf)
    if not topk:
        return np.expm1(d.min())                        # upper bound of relative error
    orig = np.flatnonzero(ok)
    res = []
    for i in range(Q.shape[0]):
        for j in range(k):
            ii = orig[i] if i < len(orig) else -1
            jj = orig[idx[i, j]] if idx[i, j] < len(orig) else -1
            res.append((d[i, j], ii, jj))
    res.sort()
    return res

def ladder_search(V, t):
    """m = 1 + T (one term, all forms by the same feature).  Best max error."""
    E = np.abs((1 + V) / t - 1).max(axis=1)
    return E.min(), E

# Family "author's freedom": m = 1 + (C/alpha) X_f, C in RATS, each form picks X from its own
# dimension-tied set (2D: k^4; 3D: sums or products over 2 or 3 planes of a^4; 1D: k^2).
def author_sets(kmax=3):
    ks = range(1, kmax + 1)
    S2 = sorted({k ** 4 for k in ks})
    S3 = set()
    for a, b in itertools.product(ks, repeat=2):
        S3 |= {a ** 4 + b ** 4, (a * b) ** 4}
    for a, b, c in itertools.product(ks, repeat=3):
        S3 |= {a ** 4 + b ** 4 + c ** 4, (a * b * c) ** 4}
    S3 |= {k ** 6 for k in ks}
    S1 = sorted({k ** 2 for k in range(1, kmax + 3)})
    return {'mu': np.array(S2, float), 'tau': np.array(sorted(S3), float), 'p': np.array(S1, float)}

AUTH = author_sets()

def author_freedom_best(t3):
    """t3 = masses of (mu, tau, p).  Best over C and per-form X of max relative error."""
    best = np.inf; arg = None
    for r in RATS:
        c = float(r) / ALPHA
        errs, picks = [], []
        for name, tm in zip(TARGETS3, t3):
            m = 1 + c * AUTH[name]
            e = np.abs(m / tm - 1); i = e.argmin()
            errs.append(e[i]); picks.append(AUTH[name][i])
        w = max(errs)
        if w < best:
            best, arg = w, (r, picks)
    return best, arg

def n_author_formulas():
    return len(RATS) * len(AUTH['mu']) * len(AUTH['tau']) * len(AUTH['p'])

# Family "per-point state rule" (F5): every point contributes v = h(state) by ONE function h,
# consistent with "multiply over planes" variants (q = 2d etc.).  Attached points may get their own
# scale (two tiers) - "an attached point's contribution depends on the state it acquires".
KDEFS = {'n': lambda p: p['n'], 'n-1': lambda p: p['n'] - 1, 'n-2': lambda p: max(p['n'] - 2, 0),
         'n/2': lambda p: p['n'] / 2, 'd-1': lambda p: p['d'] - 1}
QDEFS = {'2d': lambda p: 2 * p['d'], '2': lambda p: 2, '4': lambda p: 4, 'd': lambda p: p['d'],
         'd+1': lambda p: p['d'] + 1, '2(d-1)': lambda p: 2 * (p['d'] - 1)}
AGGS = ('pow', 'sum')
FILTERS = ('all', 'nonend')
COMPS = ('add', 'mul')
CVALS = [(sg * r, s) for r in RATS for s in S_POW for sg in (1, -1)]
CNUM = np.array([float(r) * ALPHA ** (-s) for r, s in CVALS])

def point_value(p, kd, qd, agg):
    k, q = KDEFS[kd](p), QDEFS[qd](p)
    if agg == 'pow':
        return float(k) ** q if k != 0 else 0.0
    return float(sum(j ** q for j in range(1, int(math.floor(k)) + 1)))

def perpoint_struct_arrays(forms):
    """for each structure: per form lists of (main values, attached values)."""
    structs = []
    for kd, qd, agg, flt, comp in itertools.product(KDEFS, QDEFS, AGGS, FILTERS, COMPS):
        per = {}
        for name, f in forms.items():
            pts = [p for p in f['points'] if flt == 'all' or p['role'] != 'end']
            main = [point_value(p, kd, qd, agg) for p in pts if p['role'] != 'attached']
            att = [point_value(p, kd, qd, agg) for p in pts if p['role'] == 'attached']
            per[name] = (main, att)
        structs.append(((kd, qd, agg, flt, comp), per))
    return structs

def perpoint_predictions(structs, names, c1=CNUM, c2=None):
    """returns array [struct, c1, c2, form] of predicted m/m_e (two tiers; c2=None -> one tier)."""
    e_ref = None
    out = []
    for (kd, qd, agg, flt, comp), per in structs:
        e_main = per['e'][0] if 'e' in per else None
        rows = []
        for nm in names:
            main, att = per[nm]
            if comp == 'add':
                Sm = sum(main) - sum(e_main); Sa = sum(att)
                if c2 is None:
                    rows.append(1 + c1 * (Sm + Sa))
                else:
                    rows.append(1 + c1[:, None] * Sm + c2[None, :] * Sa)
            else:
                with np.errstate(all='ignore'):
                    Pm = np.prod([1 + c1 * v for v in main], axis=0) / np.prod([1 + c1 * v for v in e_main], axis=0)
                    bad = np.zeros_like(c1, bool)
                    for v in main + e_main:
                        bad |= (1 + c1 * v) <= 0
                    Pm = np.where(bad, np.nan, Pm)
                    if c2 is None:
                        Pa = np.prod([1 + c1 * v for v in att], axis=0) if att else np.ones_like(c1)
                        for v in att: Pa = np.where(1 + c1 * v <= 0, np.nan, Pa)
                        rows.append(Pm * Pa)
                    else:
                        Pa = np.prod([1 + c2 * v for v in att], axis=0) if att else np.ones_like(c2)
                        for v in att: Pa = np.where(1 + c2 * v <= 0, np.nan, Pa)
                        rows.append(Pm[:, None] * Pa[None, :])
        out.append(np.stack(rows, axis=-1))
    return np.array(out)

C2NUM = np.concatenate([[0.0], CNUM])
C2VALS = [(0, 0.0)] + CVALS

def perpoint_best(PRED, t):
    E = np.abs(PRED / t - 1).max(axis=-1)
    E = np.where(np.isnan(E), np.inf, E)
    return E.min(), E

# ---------------------------------------------------------------------------------------------
# 4. Nulls
# ---------------------------------------------------------------------------------------------
def null_targets(t, kind):
    if kind == 'x2':      # log-uniform within a factor 2 of each true mass
        return t * np.exp(RNG.uniform(-math.log(2), math.log(2), size=t.shape))
    if kind == '10%':
        return t * np.exp(RNG.uniform(-0.1, 0.1, size=t.shape))
    if kind == '3%':
        return t * np.exp(RNG.uniform(-0.03, 0.03, size=t.shape))
    if kind == 'perm':    # masses shuffled among the forms
        return RNG.permutation(t)
    raise ValueError(kind)

KNOWN_MEV = {  # well-established hadrons/leptons, PDG central values, MeV
    'pi0': 134.9768, 'pi+': 139.57039, 'K+': 493.677, 'K0': 497.611, 'eta': 547.862, 'rho(770)': 775.26,
    'omega(782)': 782.66, "K*(892)+": 891.67, 'K*(892)0': 895.55, 'p': 938.27209, 'n': 939.56542,
    "eta'(958)": 957.78, 'f0(980)': 990.0, 'a0(980)': 980.0, 'phi(1020)': 1019.461, 'Lambda': 1115.683,
    'h1(1170)': 1166.0, 'Sigma+': 1189.37, 'Sigma0': 1192.642, 'Sigma-': 1197.449, 'b1(1235)': 1229.5,
    'a1(1260)': 1230.0, 'Delta(1232)': 1232.0, 'K1(1270)': 1253.0, 'f2(1270)': 1275.5, 'f1(1285)': 1281.9,
    'eta(1295)': 1294.0, 'a2(1320)': 1318.2, 'Xi0': 1314.86, 'Xi-': 1321.71, 'Sigma(1385)': 1384.0,
    'Lambda(1405)': 1405.1, 'K1(1400)': 1403.0, 'N(1440)': 1440.0, 'N(1520)': 1515.0, 'Lambda(1520)': 1519.0,
    'f2p(1525)': 1517.4, 'Xi(1530)': 1531.8, 'Omega-': 1672.45, 'rho3(1690)': 1688.8, 'tau': 1776.86,
    'D0': 1864.84, 'D+': 1869.66, 'Ds+': 1968.35, 'D*0': 2006.85, 'D*+': 2010.26, 'Ds*': 2112.2,
    'Lambda_c': 2286.46, 'Sigma_c(2455)': 2453.97, 'Xi_c+': 2467.71, 'Omega_c': 2695.2, 'eta_c': 2983.9,
    'J/psi': 3096.9, 'chi_c0': 3414.71, 'chi_c1': 3510.67, 'h_c': 3525.38, 'chi_c2': 3556.17,
    'psi(2S)': 3686.10, 'psi(3770)': 3773.7,
}

def nearest_known(m_mev):
    best = min(KNOWN_MEV.items(), key=lambda kv: abs(kv[1] / m_mev - 1))
    return best[0], best[1], best[1] / m_mev - 1

def chance_hit(m_mev, tol):
    """fraction of the log-interval [m/1.5, m*1.5] covered by +-tol windows around known masses."""
    lo, hi = m_mev / 1.5, m_mev * 1.5
    xs = np.exp(np.linspace(math.log(lo), math.log(hi), 20001))
    km = np.array(list(KNOWN_MEV.values()))
    hit = (np.abs(xs[:, None] / km[None, :] - 1) <= tol).any(axis=1)
    return hit.mean()

# ---------------------------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------------------------
def diff_null(t6):
    """keep mu, tau, p; randomize the three DIFFERENCES n-p, pi+-mu, pi0-mu (x/÷2, log-uniform)"""
    m = dict(zip(TARGETS6, t6))
    u = np.exp(RNG.uniform(-math.log(2), math.log(2), size=3))
    out = dict(m)
    out['n'] = m['p'] + (m['n'] - m['p']) * u[0]
    out['pi+'] = m['mu'] + (m['pi+'] - m['mu']) * u[1]
    out['pi0'] = m['mu'] + (m['pi0'] - m['mu']) * u[2]
    return np.array([out[k] for k in TARGETS6])

def main(n_null=200):
    out = []
    pr = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)), sys.stdout.flush())
    trials = {}

    pr('== 1. ENCODING (states of points) ==')
    for name, f in list(FORMS.items()) + list(NEW_FORMS.items()):
        pts = ', '.join(f"{p['role'][:3]}(n{p['n']},d{p['d']}" + (f",{p['rel']}" if p['rel'] else '') +
                        (f",carries{p['carries']}" if p['carries'] else '') + ')' for p in f['points'])
        pr(f'  {name:15s} loops {f["loops"]}, eights {f["eights"]}: {pts}')
    pr('\n  base features (excess over e):')
    pr('  ' + ' '.join(f'{k:>8s}' for k in ['form'] + BASE_NAMES))
    for name, f in list(FORMS.items()) + list(NEW_FORMS.items()):
        pr('  ' + f'{name[:8]:>8s} ' + ' '.join(f'{int(x):8d}' for x in base_vector(f)))

    pr('\n== 2. TARGETS: X = (m/m_e - 1) * alpha * 2/3 ==')
    for k in TARGETS6:
        pr(f'  {k:4s} m/m_e = {MASS[k]:.6f}   X = {(MASS[k]-1)*ALPHA*2/3:.5f}')
    t6 = np.array([MASS[k] for k in TARGETS6])
    t3 = np.array([MASS[k] for k in TARGETS3])

    # ---- 3a. author's ladder
    pr('\n== 3a. Author ladder m = 1 + (3/2)(1/a) X, X = 1, 17, 9 ==')
    for k, X in (('mu', 1), ('tau', 17), ('p', 9)):
        m = 1 + 1.5 / ALPHA * X
        pr(f'  {k:4s} X={X:2d}  m = {m:.3f}  data {MASS[k]:.3f}  err {100*(m/MASS[k]-1):+.3f} %')
    obs_auth = max(abs((1 + 1.5 / ALPHA * X) / MASS[k] - 1) for k, X in (('mu', 1), ('tau', 17), ('p', 9)))
    b, arg = author_freedom_best(t3)
    trials['F1 author freedom'] = n_author_formulas()
    pr(f'  F1 "author freedom": C in {len(RATS)} rationals, X_mu in {AUTH["mu"].tolist()}, '
       f'X_tau in {len(AUTH["tau"])} values {AUTH["tau"].astype(int).tolist()[:12]}..., X_p in {AUTH["p"].tolist()}')
    pr(f'  formulas: {n_author_formulas()}; best on data: max err {100*b:.3f} % with C={arg[0]}, X={[int(x) for x in arg[1]]}')
    for kind in ('x2', '10%', 'perm'):
        if kind == 'perm':
            vals = np.array([author_freedom_best(t3[list(pm)])[0] for pm in itertools.permutations(range(3)) if pm != (0, 1, 2)])
        else:
            vals = np.array([author_freedom_best(null_targets(t3, kind))[0] for _ in range(20000)])
        pr(f'  null {kind:5s}: replicates {len(vals)}, P(best <= {100*obs_auth:.2f} %) = {(vals <= obs_auth + 1e-12).mean():.4f}; '
           f'median best {100*np.median(vals):.2f} %')
    # loose: any form may take X from the union of all sets
    U = np.unique(np.concatenate(list(AUTH.values())))
    def loose(t):
        best = np.inf
        for r in RATS:
            c = float(r) / ALPHA
            w = max(np.abs((1 + c * U) / tm - 1).min() for tm in t)
            best = min(best, w)
        return best
    vals = np.array([loose(null_targets(t3, 'x2')) for _ in range(5000)])
    pr(f'  F1-loose (each form X from the union, {len(U)} values): data {100*loose(t3):.3f} %; null x2: P = {(vals <= obs_auth + 1e-12).mean():.4f}')
    # C = 3/2 taken as given (Barut 1979 / Nambu 1952): tau & p only, and p only
    c = 1.5 / ALPHA
    def fixed_C(tt, tp):
        return max(np.abs((1 + c * AUTH['tau']) / tt - 1).min(), np.abs((1 + c * AUTH['p']) / tp - 1).min())
    obs2 = fixed_C(MASS['tau'], MASS['p'])
    for kind in ('x2', '10%'):
        v = np.array([fixed_C(*null_targets(np.array([MASS['tau'], MASS['p']]), kind)) for _ in range(20000)])
        pr(f'  C = 3/2 given, tau & p: observed {100*obs2:.3f} %; null {kind}: P = {(v <= obs2 + 1e-12).mean():.4f}')
    obsp = np.abs((1 + c * AUTH['p']) / MASS['p'] - 1).min()
    for kind in ('x2', '10%'):
        v = np.array([np.abs((1 + c * AUTH['p']) / null_targets(np.array([MASS['p']]), kind)[0] - 1).min() for _ in range(20000)])
        pr(f'  C = 3/2 given, p only (X_p in k^2): observed {100*obsp:.3f} %; null {kind}: P = {(v <= obsp + 1e-12).mean():.4f}')

    # ---- 3b. feature grammar, one term
    fit = {k: FORMS[k] for k in TARGETS6}
    terms, V6, Vx, nfeat = build_terms(fit, NEW_FORMS)
    pr(f'\n== 3b. F2 feature grammar: {len(BASE_NAMES)} base features -> {nfeat} non-zero features; '
       f'{len(RATS)} rationals x 3 powers of 1/a x 2 signs -> {len(terms)} terms ==')
    idx3 = [TARGETS6.index(k) for k in TARGETS3]
    V3 = V6[:, idx3]
    best1_3, E3 = ladder_search(V3, t3)
    trials['F2 one term (mu,tau,p)'] = len(terms)
    order = np.argsort(E3)
    pr(f'  one term, targets mu,tau,p: best max err {100*best1_3:.3f} %; terms within 1 %: {(E3 < 0.01).sum()}')
    for i in order[:6]:
        pr(f'    {100*E3[i]:6.3f} %  m = 1 + {term_str(terms[i])}')
    for kind in ('10%', 'x2'):
        vals = np.array([ladder_search(V3, null_targets(t3, kind))[0] for _ in range(3000)])
        pr(f'  null {kind}: P(best <= {100*best1_3:.3f} %) = {(vals <= best1_3 + 1e-12).mean():.4f}')
    # without the author's own feature (AuthA/AuthB): does the grammar find it by itself?
    keep = [i for i, tm in enumerate(terms) if 'Auth' not in tm[0]]
    b_noauth = ladder_search(V3[keep], t3)[0]
    pr(f'  same without the features AuthA/AuthB: best {100*b_noauth:.2f} %')
    best1_6, E6 = ladder_search(V6, t6)
    pr(f'  one term, all six targets: best max err {100*best1_6:.2f} %  ({term_str(terms[int(np.argmin(E6))])})')

    # ---- 3b'. F2j: junction-level universal rule X = sum_J h(k_J, d_J), ONE k-definition for all junctions
    JK = {'loops': lambda j: j['loops'], 'eights-1': lambda j: j['eights'] - 1, 'eights': lambda j: j['eights'],
          'pas': lambda j: j['pas'], 'loops-eights': lambda j: j['loops'] - j['eights'], 'loops/2': lambda j: j['loops'] / 2,
          'pas-1': lambda j: j['pas'] - 1, 'd-1': lambda j: j['d'] - 1}
    JH = {'k^(2d)': lambda k, d: k ** (2 * d), 'k^2': lambda k, d: k ** 2, 'k^4': lambda k, d: k ** 4,
          'k^(d+1)': lambda k, d: k ** (d + 1), 'sum j^(2d)': lambda k, d: barut(int(k), 2 * d),
          'sum j^4': lambda k, d: barut(int(k), 4), 'sum j^2': lambda k, d: barut(int(k), 2),
          'C(d,2) k^4': lambda k, d: math.comb(d, 2) * k ** 4, 'k^(2 min(d,2))': lambda k, d: k ** (2 * min(d, 2)),
          'sum_{j<d} (jk)^4': lambda k, d: sum((jj * k) ** 4 for jj in range(1, d))}
    rows = []
    for (kn, kf), (hn, hf) in itertools.product(JK.items(), JH.items()):
        X = np.array([sum(hf(kf(j), j['d']) for j in FORMS[k]['junctions'] if j['kind'] in ('cross', 'shared')) for k in TARGETS3], float)
        for r in RATS:
            m = 1 + float(r) / ALPHA * X
            rows.append((np.abs(m / t3 - 1).max(), kn, hn, r, X))
    rows.sort(key=lambda z: z[0])
    trials['F2j junction rule'] = len(rows)
    pr(f"\n  F2j junction-level rule, one k for all junctions: {len(JK)} k-definitions x {len(JH)} h x {len(RATS)} C = {len(rows)} formulas")
    for w, kn, hn, r, X in rows[:5]:
        pr(f'    {100*w:6.2f} %  k={kn}, h={hn}, C={r}: X(mu,tau,p) = {X.tolist()}')

    # ---- 3c. two free terms, six targets
    n_pairs = len(terms) * (len(terms) + 1) // 2 + len(terms)
    trials['F3 two terms add'] = n_pairs; trials['F3 two terms mul'] = n_pairs
    resA = search_additive(V6, t6, topk=1)
    resM = search_multiplicative(V6, t6, topk=1)
    pr(f'\n== 3c. F3 two free terms on all six targets; nominal formulas: additive {n_pairs}, multiplicative {n_pairs} ==')

    def show(res, kind, L=None, nshow=6):
        L = np.zeros(6) if L is None else L
        seen = set(); shown = 0; best = None
        for d, i, j in res:
            key = tuple(sorted((i, j)))
            if key in seen: continue
            seen.add(key)
            T1 = V6[i] if 0 <= i < len(terms) else np.zeros(6)
            T2 = V6[j] if 0 <= j < len(terms) else np.zeros(6)
            m = 1 + L + T1 + T2 if kind == 'add' else (1 + T1) * (1 + T2)
            err = m / t6 - 1
            if best is None: best = (np.abs(err).max(), i, j)
            s1 = term_str(terms[i]) if 0 <= i < len(terms) else '0'
            s2 = term_str(terms[j]) if 0 <= j < len(terms) else '0'
            f = f'1 + L + {s1} + {s2}' if kind == 'add' else f'(1 + {s1}) * (1 + {s2})'
            dd = (m[TARGETS6.index('n')] - m[TARGETS6.index('p')], m[TARGETS6.index('pi+')] - m[TARGETS6.index('mu')],
                  m[TARGETS6.index('pi0')] - m[TARGETS6.index('mu')])
            pr(f'    max {100*np.abs(err).max():6.3f} %  [' + ' '.join(f'{k} {100*e:+.2f}' for k, e in zip(TARGETS6, err)) +
               f']  n-p {dd[0]:.2f} (2.53), pi+-mu {dd[1]:.1f} (66.4), pi0-mu {dd[2]:.1f} (57.4)   m = {f}')
            shown += 1
            if shown >= nshow: break
        return best
    pr('  additive m = 1 + T1 + T2:')
    bestA = show(resA, 'add')
    pr('  multiplicative m = (1 + T1)(1 + T2):')
    bestM = show(resM, 'mul')
    obsA, obsM = bestA[0], bestM[0]

    # ---- 3d. F4: author ladder for mu, tau, p, plus up to two terms that VANISH on e, mu, tau, p.
    #      Those terms alone must give the differences n-p, pi+-mu, pi0-mu (ladder gives n=p, pi=mu).
    L6 = 1.5 / ALPHA * np.array([raw_features(FORMS[k])['AuthA'] for k in TARGETS6], float)
    Lx = 1.5 / ALPHA * np.array([raw_features(NEW_FORMS[k])['AuthA'] for k in NEW_FORMS], float)
    pr('\n== 3d. F4: author ladder L = (3/2)(1/a)*AuthA, extension terms vanishing on e, mu, tau, p ==')
    eL = (1 + L6) / t6 - 1
    pr('  ladder alone: ' + ' '.join(f'{k} {100*e:+.2f} %' for k, e in zip(TARGETS6, eL)))
    iD = [TARGETS6.index(k) for k in ('n', 'pi+', 'pi0')]
    zero = np.all(V6[:, [TARGETS6.index(k) for k in TARGETS3]] == 0, axis=1)
    VD = V6[zero][:, iD]; VDx = Vx[zero]; termsD = [terms[i] for i in np.flatnonzero(zero)]
    D_obs = np.array([MASS['n'] - MASS['p'], MASS['pi+'] - MASS['mu'], MASS['pi0'] - MASS['mu']])
    pr(f'  extension terms (vanish on e, mu, tau, p): {len(termsD)}; targets n-p {D_obs[0]:.3f}, pi+-mu {D_obs[1]:.3f}, pi0-mu {D_obs[2]:.3f} m_e')
    nD = len(termsD) * (len(termsD) + 1) // 2 + len(termsD)
    trials['F4 ladder + two ext. terms'] = nD
    def searchD(D, topk=0):
        Y = np.vstack([VD, np.zeros(3)]) / D
        tree = cKDTree(Y)
        Q = np.vstack([(D - VD) / D, D / D])
        k = 1 if not topk else 6
        d, idx = tree.query(Q, k=k, p=np.inf)
        if not topk: return d.min()
        res = sorted((d[i, j], i, idx[i, j]) for i in range(Q.shape[0]) for j in range(k))
        return res
    resD = searchD(D_obs, topk=1)
    obsD = resD[0][0]
    pr(f'  formulas: {nD}; best max relative error on the three DIFFERENCES: {100*obsD:.2f} %')
    seen = set(); shown = 0; nearD = []
    for d, i, j in resD:
        key = tuple(sorted((i, j)))
        if key in seen: continue
        seen.add(key)
        T = (VD[i] if i < len(termsD) else 0) + (VD[j] if j < len(termsD) else 0)
        Tx = (VDx[i] if i < len(termsD) else 0) + (VDx[j] if j < len(termsD) else 0)
        if d <= 1.5 * obsD: nearD.append(Tx)
        if shown < 8:
            s1 = term_str(termsD[i]) if i < len(termsD) else '0'; s2 = term_str(termsD[j]) if j < len(termsD) else '0'
            pr(f'    max {100*d:6.2f} %  n-p {T[0]:7.3f}  pi+-mu {T[1]:7.2f}  pi0-mu {T[2]:7.2f}   ext = {s1} + {s2}')
            shown += 1
    bestD = resD[0]
    nearD = np.array(nearD)
    # one extension term per difference type is excluded by construction? report single-term best too
    E1 = np.abs(VD / D_obs - 1).max(axis=1)
    pr(f'  ONE extension term: best {100*E1.min():.2f} %  ({term_str(termsD[int(np.argmin(E1))])})')
    # pions only (n left out): two terms for two numbers
    EP = searchD_p = None
    Yp = np.vstack([VD[:, 1:], np.zeros(2)]) / D_obs[1:]
    treep = cKDTree(Yp)
    dp, _ = treep.query(np.vstack([(D_obs[1:] - VD[:, 1:]) / D_obs[1:], [1, 1]]), k=1, p=np.inf)
    obsP = dp.min()
    pr(f'  pions only (pi+-mu, pi0-mu), up to two terms: best {100*obsP:.3f} %')
    def searchP(D2):
        Yq = np.vstack([VD[:, 1:], np.zeros(2)]) / D2
        tq = cKDTree(Yq)
        dq, _ = tq.query(np.vstack([(D2 - VD[:, 1:]) / D2, [1, 1]]), k=1, p=np.inf)
        return dq.min()
    obsL = obsD

    # ---- 3e. F5: per-point state rule
    allforms = dict(FORMS); allforms.update(NEW_FORMS)
    structs = perpoint_struct_arrays(allforms)
    pred1 = perpoint_predictions(structs, TARGETS6)                      # [S, C, 6]
    pred2 = perpoint_predictions(structs, TARGETS6, c2=C2NUM).astype(np.float32)  # [S, C, C2, 6]
    n1, n2 = pred1.shape[0] * pred1.shape[1], pred2.shape[0] * pred2.shape[1] * pred2.shape[2]
    trials['F5 per-point one tier'] = n1; trials['F5 per-point two tiers'] = n2
    pr(f'\n== 3e. F5 per-point state rule: {len(structs)} structures (k: {list(KDEFS)}, q: {list(QDEFS)}, '
       f'{AGGS}, {FILTERS}, {COMPS}) x {len(CVALS)} coefficients; one tier {n1}, two tiers {n2} formulas ==')
    b1, E1p = perpoint_best(pred1, t6)
    b2, E2p = perpoint_best(pred2, t6)
    def show_pp(E, two, nshow=6):
        flat = np.argsort(E, axis=None)[:nshow]
        res = []
        for fl in flat:
            ix = np.unravel_index(fl, E.shape)
            st = structs[ix[0]][0]
            c1 = CVALS[ix[1]]
            m = (pred2[ix[0], ix[1], ix[2]] if two else pred1[ix[0], ix[1]]).astype(float)
            err = m / t6 - 1
            c2s = f', attached: ({C2VALS[ix[2]][0]})a^-{C2VALS[ix[2]][1]}' if two else ''
            pr(f'    max {100*E[ix]:6.3f} %  [' + ' '.join(f'{k} {100*e:+.2f}' for k, e in zip(TARGETS6, err)) +
               f']  struct k={st[0]} q={st[1]} {st[2]} {st[3]} {st[4]}; C=({c1[0]})a^-{c1[1]}{c2s}')
            res.append(ix)
        return res
    pr(f'  one tier: best {100*b1:.3f} %')
    show_pp(E1p, False)
    pr(f'  two tiers: best {100*b2:.3f} %')
    top2 = show_pp(E2p, True)
    # restriction: rules that reproduce the author's ladder on mu, tau, p (each within 1 %)
    i3 = [TARGETS6.index(k) for k in TARGETS3]
    E3pp = np.abs(pred1[..., i3] / t3 - 1).max(axis=-1)
    E3pp = np.where(np.isnan(E3pp), np.inf, E3pp)
    pr(f'  per-point one-tier rules giving mu, tau, p each within 1 %: {(E3pp < 0.01).sum()} of {n1}; best {100*E3pp.min():.2f} %')
    ix = np.unravel_index(np.argmin(E3pp), E3pp.shape)
    pr(f'    best of those on mu,tau,p: struct {structs[ix[0]][0]}, C=({CVALS[ix[1]][0]})a^-{CVALS[ix[1]][1]}, '
       f'm = {np.round(pred1[ix[0], ix[1]], 1).tolist()} (data {np.round(t6, 1).tolist()})')

    # ---- 4. Significance
    pr(f'\n== 4. Significance: identical searches on null targets ({n_null} replicates per null) ==')
    pr('  nulls: "10%" = each target x exp(U(-0.1, 0.1)); "x2" = x exp(U(-ln2, ln2)); "perm" = masses shuffled among forms;')
    pr('         "diff" = mu, tau, p kept, the differences n-p, pi+-mu, pi0-mu each x exp(U(-ln2, ln2))')
    sig = {}
    for kind in ('10%', 'x2', 'perm', 'diff'):
        acc = {k: [] for k in ('F3add', 'F3mul', 'F4', 'F4pi', 'F5one', 'F5two')}
        t0 = time.time()
        for r in range(n_null):
            tn = diff_null(t6) if kind == 'diff' else null_targets(t6, kind)
            acc['F3add'].append(search_additive(V6, tn)); acc['F3mul'].append(search_multiplicative(V6, tn))
            if kind in ('diff', '10%', 'x2'):
                Dn = np.array([tn[3] - tn[2], tn[4] - tn[0], tn[5] - tn[0]]) if kind == 'diff' else \
                     D_obs * np.exp(RNG.uniform(-(0.1 if kind == '10%' else math.log(2)), (0.1 if kind == '10%' else math.log(2)), 3))
                acc['F4'].append(searchD(np.abs(Dn))); acc['F4pi'].append(searchP(np.abs(Dn[1:])))
            else:
                acc['F4'].append(np.nan); acc['F4pi'].append(np.nan)
            acc['F5one'].append(perpoint_best(pred1, tn)[0]); acc['F5two'].append(perpoint_best(pred2, tn)[0])
        obs = {'F3add': obsA, 'F3mul': obsM, 'F4': obsL, 'F4pi': obsP, 'F5one': b1, 'F5two': b2}
        line = []
        for k, v in acc.items():
            v = np.array(v, float)
            if np.isnan(v).all(): continue
            pv = (v <= obs[k] + 1e-9).mean()
            sig[(kind, k)] = (pv, np.median(v))
            line.append(f'{k}: P={pv:.3f} (obs {100*obs[k]:.2f} %, null median {100*np.median(v):.2f} %)')
        pr(f'  null {kind:4s} [{time.time()-t0:.0f} s]:')
        for l in line: pr('      ' + l)

    pr('\n  nominal numbers of formulas tried:')
    for k, v in trials.items():
        pr(f'    {k:28s} {v:>14,d}')
    pr(f'    {"total":28s} {sum(trials.values()):>14,d}')

    # ---- 5. Blind predictions (stated BEFORE any comparison)
    names_new = list(NEW_FORMS.keys())
    pr('\n== 5. Blind predictions for forms NOT used in fitting (m in MeV) ==')
    preds = {}
    def pp_new(ix, two):
        st = structs[ix[0]]
        c1 = CNUM[ix[1]:ix[1] + 1]
        if two:
            c2 = C2NUM[ix[2]:ix[2] + 1]
            arr = perpoint_predictions([st], names_new, c1=c1, c2=c2)[0, 0, 0]
        else:
            arr = perpoint_predictions([st], names_new, c1=c1)[0, 0]
        return arr
    i, j = bestA[1], bestA[2]
    preds['F3 best additive'] = 1 + (Vx[i] if i < len(terms) else 0) + (Vx[j] if j < len(terms) else 0)
    i, j = bestD[1], bestD[2]
    preds['F4 ladder+ext best'] = 1 + Lx + (VDx[i] if i < len(termsD) else 0) + (VDx[j] if j < len(termsD) else 0)
    preds['F4 ladder alone'] = 1 + Lx
    pr('  ' + f'{"form":16s}' + ''.join(f'{k:>20s}' for k in preds))
    for q, nm in enumerate(names_new):
        pr('  ' + f'{nm:16s}' + ''.join(f'{preds[k][q]*ME_MEV:20.1f}' for k in preds))
    # spread among near-best F4 formulas (difference error <= 1.5 x best)
    pr(f'\n  spread of the EXTENSION part (added to the ladder value) among F4 formulas with error <= 1.5 x best '
       f'({len(nearD)} formulas), in m_e:')
    for q, nm in enumerate(names_new):
        col = nearD[:, q]
        pr(f'    {nm:16s} 10th pct {np.percentile(col, 10):9.1f}  median {np.median(col):9.1f}  90th pct {np.percentile(col, 90):9.1f}')
    ext_ds = nearD[:, 0]
    D_ds = HELD_OUT_MASS['Ds'] - MASS['tau']
    pr(f'  (comparison of the first row is in section 6)')

    # class of rules "fit the pion point, predict the tau point": every single extension term that gives
    # pi+ - mu within 3.2 % (the accuracy of the best pion fits); its value for tau+pt is its prediction of D_s - tau.
    okpi = np.abs(VD[:, 1] / D_obs[1] - 1) <= 0.032
    pred_ds = VDx[okpi, 0]
    pr(f'\n  single extension terms giving pi+ - mu within 3.2 %: {okpi.sum()}; their prediction of (tau+pt) - tau, m_e:')
    pr(f'    10th pct {np.percentile(pred_ds, 10):.1f}, median {np.median(pred_ds):.1f}, 90th pct {np.percentile(pred_ds, 90):.1f}; '
       f'share equal to 0: {(pred_ds == 0).mean():.3f}')
    globals()['_PRED_DS'] = pred_ds
    # ---- 6. Comparison with known particles (AFTER predictions)
    pr('\n== 6. Comparison with known particles (after the predictions above) ==')
    for k in preds:
        for q, nm in enumerate(names_new):
            mm = preds[k][q] * ME_MEV
            if mm <= 50: continue
            kn, km, rel = nearest_known(mm)
            pr(f'  {k:18s} {nm:16s} {mm:9.1f} MeV -> nearest {kn} {km} MeV ({100*rel:+.2f} %); '
               f'chance that a random value hits some listed particle within {abs(100*rel):.2f} %: {chance_hit(mm, abs(rel)):.2f}')
    pr(f'\n  D_s - tau (held out) = {D_ds:.1f} m_e. F4 near-best formulas: extension for tau+pt within 10 % of it: '
       f'{(np.abs(ext_ds / D_ds - 1) < 0.10).mean():.3f}, within 1 %: {(np.abs(ext_ds / D_ds - 1) < 0.01).mean():.3f}; '
       f'with 0 for tau+pt: {(ext_ds == 0).mean():.3f}')
    pdz = globals()['_PRED_DS']
    for tol in (0.01, 0.03, 0.10):
        hit = (np.abs(pdz / D_ds - 1) <= tol).mean()
        # reference: share of the same predictions that fall within tol of a RANDOM value log-uniform in [D/2, 2D]
        ref = np.mean([(np.abs(pdz / (D_ds * math.exp(RNG.uniform(-math.log(2), math.log(2)))) - 1) <= tol).mean() for _ in range(400)])
        pr(f'  pion-fitted terms: share predicting D_s - tau = {D_ds:.1f} m_e within {100*tol:.0f} %: {hit:.3f} '
           f'(same share for a random target in [D/2, 2D]: {ref:.3f})')
    pr(f'\n  D_s (held out) {HELD_OUT_MASS["Ds"]*ME_MEV:.2f} MeV; tau + point predictions: ' +
       ', '.join(f'{k}: {preds[k][0]*ME_MEV:.1f} MeV ({100*(preds[k][0]/HELD_OUT_MASS["Ds"]-1):+.1f} %)' for k in preds))
    return out

def ds_class_check():
    """Blind D_s test for the class 'one extension term that reproduces pi+ - mu within tol':
    each such term (vanishing on e, mu, tau, p) predicts (tau+pt) - tau.  Written before looking at its output."""
    fit = {k: FORMS[k] for k in TARGETS6}
    terms, V6, Vx, _ = build_terms(fit, NEW_FORMS)
    zero = np.all(V6[:, [TARGETS6.index(k) for k in TARGETS3]] == 0, axis=1)
    VD = V6[zero]; VDx = Vx[zero]; tD = [terms[i] for i in np.flatnonzero(zero)]
    dpi = MASS['pi+'] - MASS['mu']; dds = HELD_OUT_MASS['Ds'] - MASS['tau']
    lines = [f'D_s - tau = {dds:.1f} m_e (held out); pi+ - mu = {dpi:.2f} m_e']
    for tol in (0.032, 0.05, 0.10):
        ok = np.abs(VD[:, TARGETS6.index('pi+')] / dpi - 1) <= tol
        pdz = VDx[ok, 0]
        lines.append(f'tol {100*tol:.1f} %: {ok.sum()} terms fit pi+ - mu; predictions of (tau+pt)-tau: '
                     + ', '.join(sorted({f'{v:.1f}' for v in pdz}, key=float)[:40]))
        for t2 in (0.01, 0.03, 0.10):
            hit = (np.abs(pdz / dds - 1) <= t2).mean() if len(pdz) else float('nan')
            ref = np.mean([(np.abs(pdz / (dds * math.exp(RNG.uniform(-math.log(2), math.log(2)))) - 1) <= t2).mean()
                           for _ in range(2000)]) if len(pdz) else float('nan')
            lines.append(f'   share within {100*t2:.0f} % of D_s - tau: {hit:.3f}; for a random target in [D/2, 2D]: {ref:.3f}')
        for i in np.flatnonzero(ok)[:12]:
            lines.append(f'     {term_str(tD[i])}: pi+ - mu {VD[i, TARGETS6.index("pi+")]:.2f}, n - p {VD[i, TARGETS6.index("n")]:.2f}, '
                         f'pi0 - mu {VD[i, TARGETS6.index("pi0")]:.2f}, (tau+pt)-tau {VDx[i, 0]:.1f}')
    print('\n'.join(lines))
    with open('/home/claude/alpha/patterns_search_ds_output.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'ds':
    ds_class_check(); sys.exit(0)

if __name__ == '__main__':
    n_null = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    out = main(n_null)
    with open('/home/claude/alpha/patterns_search_output.txt', 'w') as fh:
        fh.write('\n'.join(out) + '\n')
