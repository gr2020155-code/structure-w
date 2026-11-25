"""Stage 1: sides of exchange fixed by zero contradictions over ALL neighbours - rules fixed BEFORE running (17.11.2026).

Author's positions: A-03 (exchange = swap of neighbours along a circle; side - property of the act), A-43 (every point
exchanges every tick), A-63 (the side is the one with zero contradictions), A-42 (one consistent trajectory),
rule of loop directions (neighbouring loops of one direction - contradiction; V-175, eight = loops of opposite direction A-07).

Formalisation (my choices, stated before running):
  F1 places: Poisson cloud of n points in a unit cube, neighbours = Delaunay (P16, as in all previous runs).
  F2 one tick: every point in exactly one exchange -> a perfect matching of the Delaunay graph (A-43).
  F3 side of the exchange on edge p-q: the circle (A-03) passes through one sector of the ring of points around the edge
     (vertices of Delaunay tetrahedra containing p-q), in one of two orientations. The loop's direction is the oriented normal
     n = +-(u x w), u - bond direction, w - direction to the ring point perpendicular to u. Options: 2k (k - ring length).
  F4 contradiction: two exchanges of the same tick on neighbouring edges (an end of one is a Delaunay neighbour of an end of
     the other) whose loops have one direction: n_e . n_f > 0. Zero contradictions = every neighbouring pair has n_e . n_f < 0.
Questions (descriptive, no comparison with 1/alpha in this stage):
  Q1 does a zero-contradiction tick exist?   Q2 how many (enumerate up to a cap)?   Q3 given all other loops of a solution,
  how many side options of one exchange remain consistent (is the side determined locally, A-63)?
"""
import numpy as np, sys, itertools, time
from scipy.spatial import Delaunay
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

def build(n, seed):
    rng = np.random.default_rng(seed); X = rng.random((n, 3))
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
    adj = [set(nb[ind[i]:ind[i + 1]]) for i in range(n)]
    edges = sorted({(min(i, j), max(i, j)) for i in range(n) for j in adj[i]})
    vt = [[] for _ in range(n)]
    for t, s in enumerate(tri.simplices):
        for v in s: vt[v].append(t)
    opts = {}
    for (p, q) in edges:
        u = X[q] - X[p]; u /= np.linalg.norm(u)
        ringv = {v for t in vt[p] if q in tri.simplices[t] for v in tri.simplices[t]} - {p, q}
        N = []
        for r in ringv:
            w = X[r] - X[p]; w -= (w @ u) * u; w /= np.linalg.norm(w)
            nv = np.cross(u, w); N += [nv, -nv]
        opts[(p, q)] = np.array(N)
    return X, adj, edges, opts

def solve(n, seed, cap=50):
    X, adj, edges, opts = build(n, seed)
    pool = IDPool(); xv = {e: pool.id(('x', e)) for e in edges}
    yv = {(e, o): pool.id(('y', e, o)) for e in edges for o in range(len(opts[e]))}
    cl = []
    inc = [[] for _ in range(n)]
    for e in edges: inc[e[0]].append(xv[e]); inc[e[1]].append(xv[e])
    for i in range(n):
        cl += CardEnc.equals(lits=inc[i], bound=1, vpool=pool, encoding=EncType.seqcounter).clauses
    for e in edges:
        ys = [yv[(e, o)] for o in range(len(opts[e]))]
        cl.append([-xv[e]] + ys)
        for y in ys: cl.append([xv[e], -y])
        cl += CardEnc.atmost(lits=ys, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses
    npairs = 0
    for a, b in itertools.combinations(edges, 2):
        if len({*a, *b}) < 4: continue
        if not any(v in adj[w] for v in a for w in b): continue
        npairs += 1
        D = opts[a] @ opts[b].T
        for o, o2 in zip(*np.nonzero(D > 0)):
            cl.append([-yv[(a, o)], -yv[(b, o2)]])
    s = Cadical153(bootstrap_with=cl)
    sols = []
    while len(sols) < cap and s.solve():
        m = set(l for l in s.get_model() if l > 0)
        sol = {e: o for e in edges if xv[e] in m for o in range(len(opts[e])) if yv[(e, o)] in m}
        sols.append(sol)
        s.add_clause([-yv[(e, o)] for e, o in sol.items()] + [-xv[e] for e in sol])
    exhausted = len(sols) < cap
    # Q3: local freedom of the side, other loops fixed (first solution)
    free = []
    if sols:
        sol = sols[0]
        for e, o in sol.items():
            nbr = [f for f in sol if f != e and len({*e, *f}) == 4 and any(v in adj[w] for v in e for w in f)]
            ok = sum(all(opts[e][k] @ opts[f][sol[f]] < 0 for f in nbr) for k in range(len(opts[e])))
            free.append(ok)
    return dict(n=n, seed=seed, edges=len(edges), pairs=npairs, clauses=len(cl), solutions=len(sols),
                exhausted=exhausted, free=np.array(free))

if __name__ == '__main__':
    for n, seed in ((20, 1), (20, 2), (30, 1), (30, 2), (40, 1), (60, 1)):
        t = time.time(); r = solve(n, seed)
        fr = r['free']
        print(f"n={n} seed={seed}: edges {r['edges']}, neighbouring pairs {r['pairs']}, clauses {r['clauses']}; "
              f"zero-contradiction ticks found {r['solutions']}{' (all)' if r['exhausted'] else '+ (cap)'}; "
              + (f"side options consistent with the rest: mean {fr.mean():.2f}, =1 for {np.mean(fr == 1):.2f} of exchanges" if len(fr) else "no solution")
              + f"  [{time.time() - t:.0f}s]")
        sys.stdout.flush()
