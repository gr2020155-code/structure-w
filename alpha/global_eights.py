"""Global eights on small systems - formalisation F1-F8 fixed BEFORE running in global_eights_rules_2026-11-18.md.

F1 Poisson cloud (2D square / 3D cube), Delaunay neighbours, O = place nearest the centre.
F2 every tick a perfect matching (every point exchanges, A-43).   F3 exchange = swap along the edge (A-03).
F4 period T = 2n, every point back at its start.                  F5 every point at O exactly twice per period.
F6 the two loops of a point have opposite orientation (2D: signed areas about O of opposite sign; 3D: area vectors
   with negative dot product).                                     F7 (variant b) each loop visits a hull vertex.
F8 no edge works two ticks running (cyclically): no return in the medium (A-40, A-41).
Q1 existence, Q2 number of solutions up to time shift and reversal, Q3 min / max repeats without F8 vs random matchings.
Every solution found is re-checked by an independent simulation (verify()).
"""
import numpy as np, sys, time, itertools, json
import networkx as nx
from scipy.spatial import Delaunay
from ortools.sat.python import cp_model

GRID3 = 200            # 3D: coordinates on a grid 1/GRID3 so that areas are exact integers at scale GRID3**2

def cloud(n, seed, dim):
    rng = np.random.default_rng(seed)
    while True:
        X = rng.random((n, dim))
        if dim == 3: X = np.round(X * GRID3) / GRID3
        tri = Delaunay(X)
        E = sorted({tuple(sorted((int(s[a]), int(s[b])))) for s in tri.simplices for a, b in itertools.combinations(range(dim + 1), 2)})
        G = nx.Graph(E); O = int(np.argmin(((X - 0.5) ** 2).sum(1)))
        hull = sorted({int(v) for f in tri.convex_hull for v in f})
        if O not in hull and len(nx.max_weight_matching(G, maxcardinality=True)) == n // 2:
            return X, E, O, hull

def area(X, O, u, v):
    a, b = X[u] - X[O], X[v] - X[O]                   # 2D: scalar (doubled area), 3D: vector
    return a[0] * b[1] - a[1] * b[0] if len(a) == 2 else np.cross(a, b)

def build(X, E, O, hull, f7, f8, scale, objective=None, blocks=(), at_least=None, f10=False):
    n, dim = len(X), X.shape[1]; T = 2 * n; m = cp_model.CpModel()
    x = {(t, e): m.NewBoolVar('') for t in range(T) for e in E}
    inc = {v: [e for e in E if v in e] for v in range(n)}
    for t in range(T):
        for v in range(n): m.AddExactlyOne(x[t, e] for e in inc[v])
    p = {(i, t, v): m.NewBoolVar('') for i in range(n) for t in range(T + 1) for v in range(n)}
    for i in range(n):
        for v in range(n):
            m.Add(p[i, 0, v] == int(v == i)); m.Add(p[i, T, v] == int(v == i))
        for t in range(1, T): m.AddExactlyOne(p[i, t, v] for v in range(n))
    for t in range(1, T):
        for v in range(n): m.AddExactlyOne(p[i, t, v] for i in range(n))
    D = [(u, v) for (a, b) in E for (u, v) in ((a, b), (b, a))]
    z = {}
    for i in range(n):
        for t in range(T):
            for (u, v) in D:
                e = (min(u, v), max(u, v)); zz = m.NewBoolVar(''); z[i, t, u, v] = zz
                m.AddBoolOr([p[i, t, u].Not(), x[t, e].Not(), p[i, t + 1, v]])
                m.AddImplication(zz, p[i, t, u]); m.AddImplication(zz, x[t, e])
                m.AddBoolOr([p[i, t, u].Not(), x[t, e].Not(), zz])
    rep = {}
    for t in range(T):
        for e in E:
            r = m.NewBoolVar(''); rep[t, e] = r
            m.AddImplication(r, x[t, e]); m.AddImplication(r, x[(t + 1) % T, e])
            m.AddBoolOr([x[t, e].Not(), x[(t + 1) % T, e].Not(), r])
    if f8: m.Add(sum(rep.values()) == 0)
    if at_least is not None: m.Add(sum(rep.values()) >= at_least)
    A = {(u, v): np.round(scale * area(X, O, u, v)).astype(int) for (u, v) in D}
    big = T * int(np.abs(np.array(list(A.values()))).max() + 1)
    for i in range(n):
        m.Add(sum(p[i, t, O] for t in range(T)) == 2)
        if f10:
            for t in range(T): m.AddImplication(p[i, t, O], p[i, (t + n) % T, O])
        c = [m.NewIntVar(0, 2, '') for t in range(T)]; b = [m.NewBoolVar('') for t in range(T)]
        for t in range(T):
            m.Add(c[t] == sum(p[i, s, O] for s in range(t + 1)))
            m.Add(c[t] == 1).OnlyEnforceIf(b[t]); m.Add(c[t] != 1).OnlyEnforceIf(b[t].Not())
        comps = range(1 if dim == 2 else 3)
        S1, S = [], []
        for k in comps:
            s1 = []
            for t in range(T):
                a = sum((int(A[u, v]) if dim == 2 else int(A[u, v][k])) * z[i, t, u, v] for (u, v) in D)
                w = m.NewIntVar(-big, big, ''); m.Add(w == a).OnlyEnforceIf(b[t]); m.Add(w == 0).OnlyEnforceIf(b[t].Not())
                s1.append(w)
            tot = sum((int(A[u, v]) if dim == 2 else int(A[u, v][k])) * z[i, t, u, v] for t in range(T) for (u, v) in D)
            v1 = m.NewIntVar(-big, big, ''); m.Add(v1 == sum(s1)); v2 = m.NewIntVar(-big, big, ''); m.Add(v2 == tot - v1)
            S1.append(v1); S.append(v2)
        if dim == 2:
            o = m.NewBoolVar('')
            m.Add(S1[0] >= 1).OnlyEnforceIf(o); m.Add(S[0] <= -1).OnlyEnforceIf(o)
            m.Add(S1[0] <= -1).OnlyEnforceIf(o.Not()); m.Add(S[0] >= 1).OnlyEnforceIf(o.Not())
        else:
            pr = []
            for k in range(3):
                q = m.NewIntVar(-big * big, big * big, ''); m.AddMultiplicationEquality(q, [S1[k], S[k]]); pr.append(q)
            m.Add(sum(pr) <= -1)
        if f7:
            for want in (True, False):
                g = []
                for t in range(T):
                    gg = m.NewBoolVar('')
                    m.AddImplication(gg, b[t] if want else b[t].Not())
                    m.AddBoolOr([gg.Not()] + [p[i, t, h] for h in hull]); g.append(gg)
                m.AddBoolOr(g)
    for sol in blocks:
        m.AddBoolOr([x[t, e].Not() for t in range(T) for e in sol[t]])
    m.p_vars = p                      # exposed for extra conditions (ge_symmetric.py)
    if objective == 'max_rep': m.Maximize(sum(rep.values()))
    if objective == 'min_rep': m.Minimize(sum(rep.values()))
    return m, x, rep, T

def extract(solver, x, E, T):
    return [frozenset(e for e in E if solver.Value(x[t, e])) for t in range(T)]

def transforms(M):
    T = len(M); out = []
    for s in range(T):
        out.append([M[(s + t) % T] for t in range(T)])
        out.append([M[(s - t - 1) % T] for t in range(T)])
    return out

def canon(M):
    return min(tuple(tuple(sorted(Mt)) for Mt in Mx) for Mx in transforms(M))

def verify(X, E, O, hull, M, f7, f8):
    n, dim = len(X), X.shape[1]; T = len(M); Es = set(E)
    pos = list(range(n)); traj = [[i] for i in range(n)]; where = {v: v for v in range(n)}
    for t in range(T):
        assert all(e in Es for e in M[t]) and sorted(v for e in M[t] for v in e) == list(range(n)), 'not a perfect matching'
        new = pos[:]
        for (a, b) in M[t]:
            ia, ib = pos.index(a), pos.index(b); new[ia], new[ib] = b, a
        pos = new
        for i in range(n): traj[i].append(pos[i])
    assert pos == list(range(n)), 'not periodic'
    rep = sum(len(M[t] & M[(t + 1) % T]) for t in range(T))
    if f8: assert rep == 0, 'repeat'
    for i in range(n):
        tr = traj[i]; vis = [t for t in range(T) if tr[t] == O]; assert len(vis) == 2, 'visits'
        t1, t2 = vis
        L1 = list(range(t1, t2)); L2 = [t for t in range(T) if t not in L1]
        a1 = sum(area(X, O, tr[t], tr[t + 1]) for t in L1); a2 = sum(area(X, O, tr[t], tr[t + 1]) for t in L2)
        assert (a1 * a2 < 0) if dim == 2 else (np.dot(a1, a2) < 0), 'orientation'
        if f7:
            assert any(tr[t] in hull for t in L1[1:]) and any(tr[t] in hull for t in L2 if t != t2), 'hull'
    return rep

def random_repeat_rate(E, n, ticks, seed):
    rng = np.random.default_rng(seed); G = nx.Graph(E); prev = None; hits = 0
    for k in range(ticks):
        for (a, b) in G.edges: G[a][b]['w'] = rng.random()
        Mk = {tuple(sorted(e)) for e in nx.max_weight_matching(G, maxcardinality=True, weight='w')}
        if prev is not None: hits += len(Mk & prev)
        prev = Mk
    return 2 * hits / (n * (ticks - 1))

def solve(m, limit):
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = limit; s.parameters.num_workers = 2
    st = s.Solve(m); return s, st

if __name__ == '__main__':
    dim = int(sys.argv[1]); ns = [int(a) for a in sys.argv[2].split(',')]; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    limit = float(sys.argv[4]) if len(sys.argv) > 4 else 300; cap = int(sys.argv[5]) if len(sys.argv) > 5 else 12
    scale = 1000 if dim == 2 else 40000
    for n in ns:
        X, E, O, hull = cloud(n, seed, dim); T = 2 * n
        degs = np.bincount(np.array(E).ravel(), minlength=n)
        print(f'== dim {dim}, n {n}, seed {seed}: edges {len(E)}, mean degree {degs.mean():.2f}, O degree {degs[O]}, hull {len(hull)}, T {T}')
        rr = random_repeat_rate(E, n, 400, seed); print(f'   random perfect matchings: repeat (same partner next tick) per point-tick {rr:.3f}')
        for f7 in (False, True):
            tag = 'F1-F8' + ('+F7' if f7 else '')
            t0 = time.time(); orbits = []; blocks = []; status = None
            while len(orbits) < cap:
                m, x, rep, T = build(X, E, O, hull, f7, True, scale, blocks=blocks)
                s, st = solve(m, limit)
                if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): status = s.StatusName(st); break
                M = extract(s, x, E, T); verify(X, E, O, hull, M, f7, True)
                orbits.append(canon(M)); blocks += transforms(M)
            print(f'   {tag}: orbits found {len(orbits)}' + (f' (stopped: {status})' if status else ' (cap)') + f'  [{time.time()-t0:.0f}s]')
            sys.stdout.flush()
            # Q3: min repeats = 0 if F1-F8 has a solution; max repeats by feasibility search R = 1, 2, 4, ... then bisection
            lo, hi, proved = 0, None, False
            R = 1
            t0 = time.time()
            while True:
                m, x, rep, T = build(X, E, O, hull, f7, False, scale, at_least=R); s, st = solve(m, limit)
                if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                    M = extract(s, x, E, T); r = verify(X, E, O, hull, M, f7, False); lo = max(lo, r); R = max(2 * R, r + 1)
                    if R > n * T // 2: break
                else:
                    hi = R; proved = (st == cp_model.INFEASIBLE); break
            while hi is not None and hi - lo > 1:
                R = (lo + hi) // 2
                m, x, rep, T = build(X, E, O, hull, f7, False, scale, at_least=R); s, st = solve(m, limit)
                if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                    M = extract(s, x, E, T); lo = max(lo, verify(X, E, O, hull, M, f7, False))
                else:
                    hi = R; proved = proved and (st == cp_model.INFEASIBLE)
            print(f'   {tag} without F8: max repeats >= {lo} = {2*lo/(n*T):.3f} per point-tick; '
                  f'{"none with " + str(hi) + " (proved)" if (hi is not None and proved) else "upper bound not proved"}  [{time.time()-t0:.0f}s]')
            sys.stdout.flush()
