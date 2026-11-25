#!/usr/bin/env python3
"""
junction_count.py -- счёт конфигураций стыка "(−)-петля лёгкой части -> (+)-конец тяжёлой"
для водорода (протон), позитрония (позитрон) и мюония (антимюон, крест двух восьмёрок).

Определения (R1..R7, Ω_all, Ω, Ω_s, Ω_o, кандидаты Q1..Q5) записаны в
junction_count_report_2026-11-16.md, разделы 1-2, до первого запуска.

Метод: оценка Кнута (случайный спуск по дереву перебора, вес = произведение
ветвлений) по дереву "размещение -> фаза -> акты свободных точек за 4 шага".
Оценка несмещённая; для проверки тот же код работает в режиме точного перебора.
"""
import sys, math, itertools, time, json
from functools import lru_cache
import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

ALPHA_INV = 137.035999084
P_TARGET = 2 * ALPHA_INV  # 274.071998

# ---------------------------------------------------------------- arrangements
def delaunay_periodic(N, d, seed):
    rng = np.random.default_rng(seed)
    L = N ** (1.0 / d)
    pts = rng.random((N, d)) * L
    shifts = list(itertools.product([-1, 0, 1], repeat=d))
    c = shifts.index(tuple([0] * d))
    allp = np.concatenate([pts + np.array(s) * L for s in shifts])
    idx = np.tile(np.arange(N), len(shifts))
    tri = Delaunay(allp)
    edges = set()
    lo, hi = c * N, (c + 1) * N
    for simp in tri.simplices:
        for i, j in itertools.combinations(simp, 2):
            if not (lo <= i < hi or lo <= j < hi):
                continue
            a, b = idx[i], idx[j]
            if a != b:
                edges.add((min(a, b), max(a, b)))
    return N, edges

def lattice(kind, L):
    if kind == 'square':
        n = L * L; f = lambda x, y: (x % L) * L + (y % L)
        E = {tuple(sorted((f(x, y), f(x + dx, y + dy)))) for x in range(L) for y in range(L)
             for dx, dy in [(1, 0), (0, 1)]}
    elif kind == 'tri':
        n = L * L; f = lambda x, y: (x % L) * L + (y % L)
        E = {tuple(sorted((f(x, y), f(x + dx, y + dy)))) for x in range(L) for y in range(L)
             for dx, dy in [(1, 0), (0, 1), (1, -1)]}
    elif kind == 'cubic':
        n = L ** 3; f = lambda x, y, z: ((x % L) * L + (y % L)) * L + (z % L)
        E = {tuple(sorted((f(x, y, z), f(x + a, y + b, z + cc)))) for x in range(L) for y in range(L)
             for z in range(L) for a, b, cc in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]}
    elif kind == 'fcc':  # sites x+y+z even in an LxLxL periodic cube, L even
        sites = [(x, y, z) for x in range(L) for y in range(L) for z in range(L) if (x + y + z) % 2 == 0]
        ix = {s: i for i, s in enumerate(sites)}; n = len(sites)
        vecs = [v for v in itertools.product([-1, 0, 1], repeat=3) if sorted(map(abs, v)) == [0, 1, 1]]
        E = set()
        for s in sites:
            for v in vecs:
                t = tuple((s[k] + v[k]) % L for k in range(3))
                E.add(tuple(sorted((ix[s], ix[t]))))
    E = {e for e in E if e[0] != e[1]}
    return n, E

class Arr:
    def __init__(self, n, edges, name):
        self.n, self.name = n, name
        nb = [[] for _ in range(n)]
        for a, b in edges:
            nb[a].append(b); nb[b].append(a)
        self.nb = [tuple(sorted(x)) for x in nb]
        self.adj = [frozenset(x) for x in nb]
        r = [a for a, b in edges] + [b for a, b in edges]
        cc = [b for a, b in edges] + [a for a, b in edges]
        M = csr_matrix((np.ones(len(r)), (r, cc)), shape=(n, n))
        self.dist = shortest_path(M, unweighted=True, directed=False).astype(np.int16)
        self.mean_deg = 2 * len(edges) / n

# ---------------------------------------------------------------- forms
LIGHT = dict(labels=['eA', 'eO', 'eB'], edges=[('eO', 'eA'), ('eO', 'eB')],
             dirs={'eA': '+', 'eO': '+-', 'eB': '-'}, L=4,
             forced={0: [('eO', 'eA')], 1: [('eO', 'eA')], 2: [('eO', 'eB')], 3: [('eO', 'eB')]})
HEAVY = {
    'proton': dict(labels=['Ap', 'P1', 'P2', 'C'], edges=[('Ap', 'P1'), ('P1', 'P2'), ('P2', 'C')],
                   dirs={'Ap': '+', 'P1': '+-', 'P2': '+-', 'C': '+'}, L=4,
                   forced={0: [('P1', 'Ap'), ('P2', 'C')], 1: [('P1', 'Ap'), ('P2', 'C')],
                           2: [('P1', 'P2')], 3: [('P1', 'P2')]}),
    'positron': dict(labels=['pA', 'pO', 'pB'], edges=[('pO', 'pA'), ('pO', 'pB')],
                     dirs={'pA': '-', 'pO': '+-', 'pB': '+'}, L=4,
                     forced={0: [('pO', 'pA')], 1: [('pO', 'pA')], 2: [('pO', 'pB')], 3: [('pO', 'pB')]}),
    'muon': dict(labels=['Om', 'Ma', 'Mb', 'Mc', 'Md'],
                 edges=[('Om', 'Ma'), ('Om', 'Mb'), ('Om', 'Mc'), ('Om', 'Md')],
                 dirs={'Om': '+-', 'Ma': '+', 'Mb': '-', 'Mc': '+', 'Md': '-'}, L=8,
                 forced={0: [('Om', 'Ma')], 1: [('Om', 'Ma')], 2: [('Om', 'Mb')], 3: [('Om', 'Mb')],
                         4: [('Om', 'Mc')], 5: [('Om', 'Mc')], 6: [('Om', 'Md')], 7: [('Om', 'Md')]}),
}

def allowed_pairs(sysname, version):
    H = HEAVY[sysname]; out = []
    for l in LIGHT['labels']:
        for h in H['labels']:
            dl, dh = LIGHT['dirs'][l], H['dirs'][h]
            if version == 'V1':
                ok = (dl == '-' and dh == '+')
            else:
                ok = (len(dl) == 1 and len(dh) == 1 and dl != dh)
            if ok:
                out.append((l, h))
    return out

def form_info(F):
    """free labels per step, next forced step (offset) and partner per (label, step)."""
    L = F['L']
    forced_lab = {s: {x for p in F['forced'][s] for x in p} for s in range(L)}
    free = {s: [x for x in F['labels'] if x not in forced_lab[s]] for s in range(L)}
    partner = {}
    for s in range(L):
        for a, b in F['forced'][s]:
            partner[(a, s)] = b; partner[(b, s)] = a
    nxt = {}
    for x in F['labels']:
        for s in range(L):
            for k in range(1, L + 1):
                s2 = (s + k) % L
                if x in forced_lab[s2]:
                    nxt[(x, s)] = (k, partner[(x, s2)]); break
    # is step s the second step of a loop (same pair forced at s-1)?
    second = {s: [p for p in F['forced'][s] if p in F['forced'][(s - 1) % L]] for s in range(L)}
    return dict(free=free, nxt=nxt, second=second, forced_lab=forced_lab)

LINFO = form_info(LIGHT)
HINFO = {k: form_info(v) for k, v in HEAVY.items()}

# ---------------------------------------------------------------- Q4 (schedules only)
def q4(sysname, version):
    H = HEAVY[sysname]; hi = HINFO[sysname]
    al = set(allowed_pairs(sysname, version))
    Lh = H['L']; T = 4 * Lh // math.gcd(4, Lh)
    slots = 0; cyc_ok = 0
    for phi in range(Lh):
        for c in range(T // 4):
            any_c = False
            for t4 in range(4):
                t = 4 * c + t4
                lf = LINFO['free'][t % 4]; hf = hi['free'][(t + phi) % Lh]
                if any((l, h) in al for l in lf for h in hf):
                    slots += 1; any_c = True
            cyc_ok += any_c
    return (T * Lh) / slots if slots else math.inf, (T // 4 * Lh) / cyc_ok if cyc_ok else math.inf, T

# ---------------------------------------------------------------- the counting tree
class Counter:
    def __init__(self, arr, sysname, version, tree):
        self.a, self.sys, self.V, self.tree = arr, sysname, version, tree  # tree: 'cons' (Ω) or 'all' (Ω_all)
        self.H = HEAVY[sysname]; self.hi = HINFO[sysname]
        self.allowed = set(allowed_pairs(sysname, version))
        self.contact_types = allowed_pairs(sysname, version)
        self.light = set(LIGHT['labels'])
        self.Lh = self.H['L']
        # BFS growth order of heavy form from each root
        self.grow = {}
        for root in {h for _, h in self.contact_types}:
            order, seen = [], {root}
            frontier = [root]
            while frontier:
                nf = []
                for u in frontier:
                    for a, b in self.H['edges']:
                        for x, y in ((a, b), (b, a)):
                            if x == u and y not in seen:
                                seen.add(y); order.append((y, u)); nf.append(y)
                frontier = nf
            self.grow[root] = order
        self.reach_cache = {}

    # --- schedule helpers in light-step time t (heavy step = t+phi)
    def forced(self, t, phi):
        return LIGHT['forced'][t % 4] + self.H['forced'][(t + phi) % self.Lh]
    def free(self, t, phi):
        return LINFO['free'][t % 4] + self.hi['free'][(t + phi) % self.Lh]
    def next_forced(self, x, t, phi):
        if x in self.light:
            k, p = LINFO['nxt'][(x, t % 4)]
        else:
            k, p = self.hi['nxt'][(x, (t + phi) % self.Lh)]
        return t + k, p

    def reach(self, q, prev, r, target):
        """exists walk of exactly r moves from q ending adjacent to target (≠ target)."""
        key = (q, prev, r, target)
        c = self.reach_cache.get(key)
        if c is not None:
            return c
        a = self.a
        if a.dist[q, target] > r + 1:
            res = False
        elif r == 0:
            res = target in a.adj[q]
        else:
            res = any(self.reach(q2, q if self.tree == 'cons' else None, r - 1, target)
                      for q2 in a.nb[q] if q2 != target and (self.tree != 'cons' or q2 != prev))
        if len(self.reach_cache) > 2_000_000:
            self.reach_cache.clear()
        self.reach_cache[key] = res
        return res

    def prune_ok(self, x, q, t, phi):
        t2, p = self.next_forced(x, t, phi)
        r = t2 - (t + 1)
        return self.a.dist[q, self.home[p]] <= r + 1 and q != self.home[p]

    # --- placement (Knuth level); returns weight or 0
    def place(self, rng):
        a = self.a; w = float(a.n)
        o = int(rng.integers(a.n)); nbo = a.nb[o]
        if len(nbo) < 2:
            return 0.0, None
        A = nbo[rng.integers(len(nbo))]
        rest = [x for x in nbo if x != A]
        B = rest[rng.integers(len(rest))]
        w *= len(nbo) * len(rest)
        home = {'eO': o, 'eA': A, 'eB': B}
        occ = {o: 'eO', A: 'eA', B: 'eB'}
        ct = self.contact_types
        l, h = ct[rng.integers(len(ct))]; w *= len(ct)
        opts = [x for x in a.nb[home[l]] if x not in occ]
        if not opts:
            return 0.0, None
        x = opts[rng.integers(len(opts))]; w *= len(opts)
        home[h] = x; occ[x] = h
        for y, u in self.grow[h]:
            opts = [z for z in a.nb[home[u]] if z not in occ]
            if not opts:
                return 0.0, None
            z = opts[rng.integers(len(opts))]; w *= len(opts)
            home[y] = z; occ[z] = y
        contacts = [(l2, h2) for l2, h2 in ct if home[h2] in a.adj[home[l2]]]
        w /= len(contacts)
        return w, (home, contacts)

    def start(self, home, contacts, phi):
        """t0 positions: a loop in its second step at t0 starts swapped."""
        pos = dict(home)
        for x, y in LINFO['second'][0]:
            pos[x], pos[y] = pos[y], pos[x]
        for x, y in self.hi['second'][phi % self.Lh]:
            pos[x], pos[y] = pos[y], pos[x]
        self.home = home
        self.contact_h = {h for _, h in contacts}
        return pos

    # --- recursive core: mode 'exact' sums all branches, mode 'knuth' follows one random branch
    def run_steps(self, pos, phi, mode, rng, record=None):
        self.phi = phi
        return self._step(0, pos, {x: None for x in pos}, False, False, mode, rng, record)

    def _step(self, t, pos, prev, sh, ot, mode, rng, record):
        if t == 4:
            return self._end(pos, prev, sh, ot)
        phi = self.phi
        F = self.forced(t, phi)
        for x, y in F:
            if pos[y] not in self.a.adj[pos[x]]:
                return (0.0, 0.0, 0.0)
        occ = {p: x for x, p in pos.items()}
        return self._choose(t, pos, prev, sh, ot, mode, rng, F, self.free(t, phi), 0, frozenset(), frozenset(), (), occ, record)

    def _choose(self, t, pos, prev, sh, ot, mode, rng, F, fl, i, claimed, paired, acts, occ, record):
        while i < len(fl) and fl[i] in paired:
            i += 1
        if i == len(fl):
            # R1 for the medium: the remaining places must still split into exchange pairs
            used = set(pos.values()) | {q for _, q, _ in acts}
            if not self.medium_ok(used):
                return (0.0, 0.0, 0.0)
            # apply all acts simultaneously
            new = dict(pos); nprev = {x: None for x in pos}
            for x, y in F:
                new[x], new[y] = pos[y], pos[x]
            sh2, ot2 = sh, ot
            for f, q, lab in acts:
                new[f] = q; nprev[f] = pos[f]
                if lab is not None:
                    new[lab] = pos[f]; nprev[lab] = q
                    sh2 = True
                    if lab not in self.contact_h:
                        ot2 = True
            if record is not None:
                record.append((t, dict(pos), acts, list(F)))
            return self._step(t + 1, new, nprev, sh2, ot2, mode, rng, record)
        f = fl[i]; pf = pos[f]; isl = f in self.light; cons = self.tree == 'cons'
        opts = []
        for q in self.a.nb[pf]:
            if q in claimed:
                continue
            if cons and prev[f] is not None and q == prev[f]:
                continue
            lab = occ.get(q)
            if lab is None:
                if self.prune_ok(f, q, t, self.phi):
                    opts.append((q, None))
            else:
                if lab in paired or lab not in fl or (lab in self.light) == isl or not isl:
                    continue
                if cons and (f, lab) not in self.allowed:
                    continue
                if self.prune_ok(f, q, t, self.phi) and self.prune_ok(lab, pf, t, self.phi):
                    opts.append((q, lab))
        if not opts:
            return (0.0, 0.0, 0.0)
        def go(o):
            q, lab = o
            p2 = paired | {f} | ({lab} if lab else set())
            return self._choose(t, pos, prev, sh, ot, mode, rng, F, fl, i + 1, claimed | {q}, p2,
                                acts + ((f, q, lab),), occ, record)
        if mode == 'exact':
            s = [0.0, 0.0, 0.0]
            for o in opts:
                r = go(o)
                for k in range(3):
                    s[k] += r[k]
            return tuple(s)
        o = opts[rng.integers(len(opts))]
        r = go(o)
        m = len(opts)
        return (r[0] * m, r[1] * m, r[2] * m)

    def medium_ok(self, used, cap=14):
        """R1 for the medium (local test): every pocket of medium places closed off by `used`
        (component of size <= cap) must have a perfect matching; the rest is the open medium.
        Validated against a full blossom matching in medium_check()."""
        nb = self.a.nb; seen = set()
        for u in used:
            for v in nb[u]:
                if v in used or v in seen:
                    continue
                comp = [v]; seen_c = {v}; k = 0; big = False
                while k < len(comp):
                    x = comp[k]; k += 1
                    for y in nb[x]:
                        if y not in used and y not in seen_c:
                            seen_c.add(y); comp.append(y)
                            if len(comp) > cap:
                                big = True; break
                    if big:
                        break
                seen |= seen_c   # a big component's places need no second search
                if big:
                    continue
                if len(comp) % 2 or not self._pm_small(comp):
                    return False
        return True

    def _pm_small(self, comp):
        if not comp:
            return True
        cs = set(comp); x = comp[0]; rest = comp[1:]
        for y in self.a.nb[x]:
            if y in cs:
                if self._pm_small([z for z in rest if z != y]):
                    return True
        return False

    def _end(self, pos, prev, sh, ot):
        # R6 at the horizon: every point free in the last step must be able to reach its next loop partner
        for x in pos:
            if x in self.free(3, self.phi):
                t2, p = self.next_forced(x, 3, self.phi)
                r = t2 - 4
                if not self.reach(pos[x], prev[x] if self.tree == 'cons' else None, r, self.home[p]):
                    return (0.0, 0.0, 0.0)
        return (1.0, 1.0 if sh else 0.0, 1.0 if ot else 0.0)

    # --- one Knuth sample over placement + phase + acts
    def sample(self, rng):
        w, pl = self.place(rng)
        if w == 0.0:
            return 0.0, 0.0, 0.0
        phi = int(rng.integers(self.Lh)); w *= self.Lh
        pos = self.start(pl[0], pl[1], phi)
        r = self.run_steps(pos, phi, 'knuth', rng)
        return w * r[0], w * r[1], w * r[2]

# ---------------------------------------------------------------- driver
def estimate(arr, sysname, version, tree, nsamp, seed, nbatch=20):
    C = Counter(arr, sysname, version, tree)
    rng = np.random.default_rng(seed)
    S = np.zeros((nbatch, 3))
    per = nsamp // nbatch
    for b in range(nbatch):
        for _ in range(per):
            S[b] += C.sample(rng)
    return S  # batch sums of (valid, shared, other)

def ratio_err(num, den):
    """ratio of totals and batch-jackknife standard error."""
    N, D = num.sum(), den.sum()
    if D == 0:
        return math.inf, math.nan
    R = N / D
    k = len(num)
    jk = [(N - num[i]) / (D - den[i]) if D - den[i] > 0 else math.inf for i in range(k)]
    jk = np.array(jk)
    if not np.all(np.isfinite(jk)):
        return R, math.nan
    se = math.sqrt((k - 1) / k * ((jk - jk.mean()) ** 2).sum())
    return R, se

def make_arr(spec):
    kind = spec[0]
    if kind == 'D2':
        n, E = delaunay_periodic(spec[1], 2, spec[2]); return Arr(n, E, f"Delaunay2D N={spec[1]} seed={spec[2]}")
    if kind == 'D3':
        n, E = delaunay_periodic(spec[1], 3, spec[2]); return Arr(n, E, f"Delaunay3D N={spec[1]} seed={spec[2]}")
    n, E = lattice(kind, spec[1]); return Arr(n, E, f"lattice {kind} L={spec[1]}")

def job(args):
    spec, sysname, version, tree, nsamp, seed = args
    arr = ARR_CACHE.get(spec)
    if arr is None:
        arr = make_arr(spec); ARR_CACHE[spec] = arr
    t0 = time.time()
    S = estimate(arr, sysname, version, tree, nsamp, seed)
    return (spec, sysname, version, tree, S.tolist(), arr.mean_deg, time.time() - t0)

ARR_CACHE = {}

# ---------------------------------------------------------------- checks
def exact_vs_knuth(spec=('square', 8), sysname='proton', version='V1', tree='cons', nk=20000, seed=7):
    """fix one placement, sum over phases exactly, compare with Knuth estimate on the same placement."""
    arr = make_arr(spec); rng = np.random.default_rng(seed)
    C = Counter(arr, sysname, version, tree)
    while True:
        w, pl = C.place(rng)
        if w > 0:
            break
    ex = np.zeros(3)
    for phi in range(C.Lh):
        pos = C.start(pl[0], pl[1], phi)
        ex += C.run_steps(pos, phi, 'exact', rng)
    kn = np.zeros(3)
    for _ in range(nk):
        phi = int(rng.integers(C.Lh))
        pos = C.start(pl[0], pl[1], phi)
        kn += np.array(C.run_steps(pos, phi, 'knuth', rng)) * C.Lh
    kn /= nk
    return ex, kn

def medium_check(spec=('D2', 144, 1), sysname='proton', version='V1', nleaf=15, seed=3):
    """for random valid leaves of Ω: does the rest of the places admit a perfect matching at every step?"""
    import networkx as nx
    arr = make_arr(spec); rng = np.random.default_rng(seed)
    C = Counter(arr, sysname, version, 'cons')
    G = nx.Graph(); G.add_nodes_from(range(arr.n))
    for u in range(arr.n):
        for v in arr.nb[u]:
            G.add_edge(u, v)
    ok = bad = 0
    tries = 0
    while ok + bad < nleaf and tries < 200000:
        tries += 1
        w, pl = C.place(rng)
        if w == 0:
            continue
        phi = int(rng.integers(C.Lh))
        pos = C.start(pl[0], pl[1], phi)
        rec = []
        r = C.run_steps(pos, phi, 'knuth', rng, rec)
        if r[0] == 0:
            continue
        good = True
        for t, p, acts, F in rec:
            used = set(p.values()) | {q for _, q, _ in acts}
            H = G.subgraph([v for v in range(arr.n) if v not in used])
            m = nx.max_weight_matching(H, maxcardinality=True)
            if 2 * len(m) != H.number_of_nodes():
                good = False; break
        ok += good; bad += (not good)
    return ok, bad

def summary(path='junction_count_results.json'):
    """universality (same value for proton/positron/muon), size dependence, comparison with 274.07."""
    from scipy.stats import chi2
    T = json.load(open(path)); out = []
    P = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))
    def key(r): return (tuple(r['spec']), r['sys'], r['V'])
    D = {key(r): r for r in T}
    specs = []
    for r in T:
        if tuple(r['spec']) not in specs: specs.append(tuple(r['spec']))
    sets = {'V1': [('proton', 'V1'), ('positron', 'V1'), ('muon', 'V1')],
            'V2': [('proton', 'V1'), ('positron', 'V2'), ('muon', 'V2')]}
    def chi_eq(vals):
        v = np.array([x[0] for x in vals]); e = np.array([x[1] for x in vals])
        if not np.all(np.isfinite(v)) or not np.all(np.isfinite(e)) or np.any(e <= 0):
            return math.nan, 0.0 if not np.all(np.isfinite(v)) else math.nan
        w = 1 / e ** 2; m = (w * v).sum() / w.sum(); c = (w * (v - m) ** 2).sum()
        return c, chi2.sf(c, len(v) - 1)
    P("\n## Универсальность: одно ли число у H, Ps, Mu (χ² равенства, p; p<0.01 = различие)")
    for Q in ('Q1', 'Q2', 'Q3', 'Q5'):
        P(f"\n### {Q}")
        for spec in specs:
            for V, S in sets.items():
                vals = [tuple(D[(spec, s, v)][Q]) for s, v in S]
                c, p = chi_eq(vals)
                txt = '  '.join(f"{s[:3]}={x[0]:.4g}±{x[1]:.2g}" for (s, v), x in zip(S, vals))
                fin = [x[0] for x in vals if np.isfinite(x[0])]
                spread = max(fin) / min(fin) if len(fin) == 3 and min(fin) > 0 else math.inf
                P(f"{str(spec):22s} {V}: {txt}   max/min={spread:.3g}  χ²={c:.3g} p={p:.2g}")
    P("\n## Зависимость от размера (один тип расположения, разные N): χ² равенства по N")
    groups = {}
    for spec in specs:
        if spec[0] in ('D2', 'D3'):
            groups.setdefault(spec[0], []).append(spec)
    for g, sp in groups.items():
        for s, v in [('proton', 'V1'), ('positron', 'V1'), ('positron', 'V2'), ('muon', 'V1'), ('muon', 'V2')]:
            for Q in ('Q1', 'Q2'):
                vals = [tuple(D[(x, s, v)][Q]) for x in sp]
                c, p = chi_eq(vals)
                P(f"{g} {s:8s} {v} {Q}: " + '  '.join(f"N={x[1]}:{y[0]:.4g}±{y[1]:.2g}" for x, y in zip(sp, vals)) + f"   p={p:.2g}")
    P("\n## Сравнение с 274.071998 и 137.036 (z = (Q-цель)/σ), только конечные значения")
    best = []
    for r in T:
        for Q in ('Q1', 'Q2', 'Q3', 'Q5'):
            v, e = r[Q]
            if np.isfinite(v) and np.isfinite(e) and e > 0:
                for tgt in (P_TARGET, ALPHA_INV):
                    best.append((abs(v - tgt) / tgt, (v - tgt) / e, Q, tuple(r['spec']), r['sys'], r['V'], v, e, tgt))
    best.sort()
    for b in best[:15]:
        P(f"{b[2]} {str(b[3]):22s} {b[4]:8s} {b[5]}: {b[6]:.4g} ± {b[7]:.2g}  vs {b[8]:.3f}: отклонение {100*(b[6]/b[8]-1):+.1f}%  z={b[1]:+.2f}")
    with open('junction_count_summary.txt', 'w') as fh:
        fh.write('\n'.join(out) + '\n')

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'summary':
    summary(); sys.exit(0)

if __name__ == '__main__':
    import multiprocessing as mp
    out = []
    def P(*a):
        s = ' '.join(str(x) for x in a); print(s); out.append(s); sys.stdout.flush()

    P(f"# junction_count.py  {time.strftime('%Y-%m-%d %H:%M')}")
    P(f"P_target = 2/alpha = {P_TARGET:.6f};  1/alpha = {ALPHA_INV}")
    P("\n## Разрешённые общие акты")
    for s in HEAVY:
        for V in ('V1', 'V2'):
            P(f"{s:9s} {V}: {allowed_pairs(s, V)}")

    P("\n## Q4 (только расписания)")
    q4res = {}
    for s in HEAVY:
        for V in ('V1', 'V2'):
            a, b, T = q4(s, V); q4res[(s, V)] = (a, b)
            P(f"{s:9s} {V}: Q4a = {a:.6g}   Q4b = {b:.6g}   (joint period {T} steps)")

    P("\n## Проверка: точный перебор против оценки Кнута (одно размещение, все фазы)")
    for spec, s, V, tr in [(('square', 8), 'proton', 'V1', 'cons'), (('tri', 8), 'positron', 'V2', 'cons'),
                           (('square', 8), 'positron', 'V2', 'all')]:
        ex, kn = exact_vs_knuth(spec, s, V, tr)
        P(f"{spec} {s} {V} {tr}: exact(valid,shared,other) = {ex}   knuth = {np.round(kn, 1)}")

    NS = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    specs = [('D2', 144, 11), ('D2', 400, 12), ('D2', 1600, 13),
             ('D3', 216, 21), ('D3', 512, 22), ('D3', 1000, 23),
             ('square', 20), ('tri', 20), ('cubic', 8), ('fcc', 8)]
    jobs = []
    k = 0
    for spec in specs:
        for s in HEAVY:
            for V in ('V1', 'V2'):
                if s == 'proton' and V == 'V2':
                    continue  # V1 = V2 for the proton (no pure (−) point)
                for tr in ('cons', 'all'):
                    k += 1
                    jobs.append((spec, s, V, tr, NS, 1000 + k))
    t0 = time.time()
    with mp.Pool(2) as pool:
        res = pool.map(job, jobs, chunksize=1)
    P(f"\n## Оценки Кнута: {NS} выборок на прогон, {len(jobs)} прогонов, {time.time()-t0:.0f} c")
    R = {}
    for spec, s, V, tr, S, md, dt in res:
        R[(spec, s, V, tr)] = (np.array(S), md)
    table = []
    for spec in specs:
        md = R[(spec, 'proton', 'V1', 'cons')][1]
        P(f"\n### {spec}  mean degree = {md:.3f}")
        for s in HEAVY:
            for V in ('V1', 'V2'):
                if s == 'proton' and V == 'V2':
                    continue
                Sc = R[(spec, s, V, 'cons')][0]; Sa = R[(spec, s, V, 'all')][0]
                q1 = ratio_err(Sc[:, 0], Sc[:, 1])
                q5 = ratio_err(Sc[:, 0], Sc[:, 2])
                q2 = ratio_err(Sa[:, 0], Sc[:, 0])
                q3 = (q1[0] * q2[0], q1[0] * q2[0] * math.hypot(q1[1] / q1[0], q2[1] / q2[0]) if q1[0] < math.inf and q2[0] < math.inf else math.nan)
                fa = Sa[:, 1].sum() / Sa[:, 0].sum() if Sa[:, 0].sum() else math.nan
                P(f"{s:9s} {V}: Q1 = {q1[0]:9.4f} ± {q1[1]:.4f}   Q2 = {q2[0]:9.4f} ± {q2[1]:.4f}   "
                  f"Q3 = {q3[0]:9.4f} ± {q3[1]:.4f}   Q5 = {q5[0]:10.4f} ± {q5[1]:.4f}   "
                  f"(доля Ω_all с любым перекрёстным актом = {fa:.4f})")
                table.append(dict(spec=list(spec), sys=s, V=V, deg=md, Q1=q1, Q2=q2, Q3=q3, Q5=q5,
                                  Q4a=q4res[(s, V)][0], Q4b=q4res[(s, V)][1]))
    with open('junction_count_results.json', 'w') as fh:
        json.dump(table, fh, indent=1, default=float)

    P("\n## Проверка локального правила среды полным паросочетанием (блоссом) на листьях Ω")
    for spec in [('D2', 144, 1), ('D3', 216, 2)]:
        for s in HEAVY:
            ok, bad = medium_check(spec, s, 'V2', nleaf=30)
            P(f"{spec} {s:9s}: листьев {ok+bad}, полное паросочетание есть: {ok}, нет: {bad}")

    with open('junction_count_output.txt', 'w') as fh:
        fh.write('\n'.join(out) + '\n')
