"""Symmetric global eights - FS1, FS2 fixed BEFORE running in global_eights_rules_2026-11-18.md (addendum FS).

FS1: 2D places symmetric about the line x = 0.5 through O = (0.5, 0.5): O and one more place on the line, the rest mirror
pairs; neighbours = Delaunay edges whose mirror image is also a Delaunay edge.
FS2: after n ticks every point stands at the mirror image of its place: pos(t + n) = sigma(pos(t)) (second loop = mirror
copy of the first; equal lengths and opposite orientation follow). Other conditions F2-F5, F8 as in global_eights.py.
Reports per arrangement: number of solutions up to time shift and reversal, number of distinct paths and reversed pairs.
"""
import sys, time, itertools, numpy as np, networkx as nx
from scipy.spatial import Delaunay
from ortools.sat.python import cp_model
from global_eights import build, extract, transforms, canon, verify, solve

def sym_cloud(n, seed, a=2):
    rng = np.random.default_rng(seed); k = (n - a) // 2
    while True:
        L = np.column_stack([rng.uniform(0.02, 0.48, k), rng.uniform(0, 1, k)])
        X = np.vstack([[0.5, 0.5]] + [[0.5, y] for y in rng.uniform(0, 1, a - 1)] + [L, np.column_stack([1 - L[:, 0], L[:, 1]])])
        sig = list(range(a)) + [a + k + j for j in range(k)] + [a + j for j in range(k)]
        tri = Delaunay(X)
        D = {tuple(sorted((int(s[a]), int(s[b])))) for s in tri.simplices for a, b in itertools.combinations(range(3), 2)}
        E = sorted(e for e in D if tuple(sorted((sig[e[0]], sig[e[1]]))) in D)
        hull = sorted({int(v) for f in tri.convex_hull for v in f})
        G = nx.Graph(E); G.add_nodes_from(range(n))
        if 0 not in hull and nx.is_connected(G) and len(nx.max_weight_matching(G, maxcardinality=True)) == n // 2:
            return X, E, 0, hull, sig

def paths_info(M, n):
    T = len(M); where = list(range(n)); tr = [[i] for i in range(n)]
    for t in range(T):
        owner = {where[i]: i for i in range(n)}
        for (a, b) in M[t]: where[owner[a]], where[owner[b]] = b, a
        for i in range(n): tr[i].append(where[i])
    P = [p[:-1] for p in tr]; c = lambda p: min(tuple(p[s:] + p[:s]) for s in range(len(p)))
    g = {}
    for i, p in enumerate(P): g.setdefault(c(p), []).append(i)
    rev = sum(1 for a in g for b in g if a < b and c(list(a)[::-1]) == b)
    return len(g), sorted(len(v) for v in g.values()), rev, P

n = int(sys.argv[1]); seeds = [int(s) for s in sys.argv[2].split(',')]; limit = float(sys.argv[3]); cap = int(sys.argv[4])
a = int(sys.argv[5]) if len(sys.argv) > 5 else 2
fs2 = 'nofs2' not in sys.argv[6:]
for seed in seeds:
    X, E, O, hull, sig = sym_cloud(n, seed, a); T = 2 * n
    print(f'== {"FS1+FS2" if fs2 else "FS1 only (symmetric places, motion free)"}, n {n}, places on the mirror line {a}, mirror pairs {(n - a) // 2}, symmetric cloud {seed}: edges {len(E)}, O degree {sum(O in e for e in E)}'); sys.stdout.flush()
    t0 = time.time(); orbits, blocks, status, info = [], [], None, []
    while len(orbits) < cap:
        m, x, rep, T = build(X, E, O, hull, False, True, 1000, blocks=blocks); p = m.p_vars
        if fs2:
            for i in range(n): m.Add(p[i, n, sig[i]] == 1)                                   # FS2: config(n) = sigma(config(0))
            for t in range(n):
                for e in E: m.Add(x[t, e] == x[t + n, tuple(sorted((sig[e[0]], sig[e[1]])))])  # FS2: M(t+n) = sigma(M(t))
        s, st = solve(m, limit)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): status = s.StatusName(st); break
        M = extract(s, x, E, T); verify(X, E, O, hull, M, False, True)
        npaths, sizes, rev, P = paths_info([sorted(Mt) for Mt in M], n)
        if fs2:
            for i in range(n): assert all(P[i][(t + n) % T] == sig[P[i][t]] for t in range(T)), 'FS2 violated'
        orbits.append(canon(M)); blocks += transforms(M); info.append((npaths, sizes, rev))
    print(f'   solutions (up to time shift and reversal): {len(orbits)}' + (f', search ended: {status}' if status else ' (cap)')
          + f'  [{time.time()-t0:.0f}s]; paths per solution {info}'); sys.stdout.flush()
