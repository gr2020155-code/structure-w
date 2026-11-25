"""Uniqueness check of the vacuum trajectory (W-Structure), 15.11.2026.

Constraints, fixed BEFORE running (author's positions):
  C1  every point exchanges in every step            -> each step is a perfect matching of the graph
  C2  K1: a point takes part in one exchange per step -> (same: a matching)
  C3  K2: an exchange does not open a closure        -> vacuous here: the region holds no particle
  C4  global self-consistency: the whole trajectory closes on itself after T steps
      (every point is back at its place after T steps)
  O12 indistinguishable trajectories count as one   -> solutions are counted up to graph automorphisms
Two readings of the medium (open question in electron v1.1, section 9):
  A  "there and back": M1 M1 M2 M2 (T=4), M1, M2 disjoint            (medium without a lattice)
  B  no repeated pairing on consecutive steps (cyclically), T=4        (medium on the square lattice)
Number 1/137 is not used anywhere.
"""
import itertools, sys, random
import networkx as nx
from pysat.solvers import Minisat22
from pysat.card import CardEnc, EncType

CAP = 200_000

def torus(dims):
    G = nx.grid_graph(dim=list(dims), periodic=True)
    return nx.convert_node_labels_to_integers(G, label_attribute='pos')

def cloud(n, seed):
    import numpy as np
    from scipy.spatial import Delaunay
    rng = np.random.default_rng(seed); P = rng.random((n, 3))
    G = nx.Graph()
    for s in Delaunay(P).simplices:
        for a, b in itertools.combinations(s, 2): G.add_edge(int(a), int(b))
    return G

class Enc:
    def __init__(s): s.n = 0; s.cl = []
    def var(s): s.n += 1; return s.n
    def exactly1(s, lits):
        r = CardEnc.equals(lits=lits, bound=1, top_id=s.n, encoding=EncType.seqcounter)
        s.n = max(s.n, r.nv); s.cl += r.clauses

def perfect_matching_vars(e, G):
    x = {}
    for u, v in G.edges(): x[frozenset((u, v))] = e.var()
    for u in G.nodes(): e.exactly1([x[frozenset((u, w))] for w in G[u]])
    return x

def model_A(G):
    e = Enc(); m1 = perfect_matching_vars(e, G); m2 = perfect_matching_vars(e, G)
    for k in m1: e.cl.append([-m1[k], -m2[k]])          # disjoint: the second partner differs
    return e, [m1, m2]

def model_B(G, T=4):
    e = Enc(); M = [perfect_matching_vars(e, G) for _ in range(T)]
    for t in range(T):                                   # no repeated pairing on consecutive steps
        for k in M[t]: e.cl.append([-M[t][k], -M[(t + 1) % T][k]])
    # position of each point after each step: p[t][point][site]; start p[0][i][i]
    nodes = list(G.nodes()); N = len(nodes)
    # track only "where is point i at time t" among sites reachable in t steps
    pos = [{(i, i): None} for i in nodes]
    P = {}
    for i in nodes:
        reach = {i}; P[(0, i, i)] = True
        for t in range(T):
            nxt = set()
            for s in reach: nxt |= set(G[s])
            for s2 in nxt: P[(t + 1, i, s2)] = e.var()
            # moves: if i at s at time t and edge (s,s2) in M[t] -> i at s2 at t+1
            for s in reach:
                at = P[(t, i, s)]
                for s2 in G[s]:
                    lit = M[t][frozenset((s, s2))]
                    if at is True: e.cl.append([-lit, P[(t + 1, i, s2)]])
                    else: e.cl.append([-at, -lit, P[(t + 1, i, s2)]])
            e.exactly1([P[(t + 1, i, s2)] for s2 in nxt])
            reach = nxt
        if (T, i, i) not in P: return None, None
        e.cl.append([P[(T, i, i)]])                       # C4: back at its place
    return e, M

def enumerate_solutions(e, Ms):
    sols = []
    with Minisat22(bootstrap_with=e.cl) as s:
        allv = [v for M in Ms for v in M.values()]
        while s.solve():
            m = set(l for l in s.get_model() if l > 0)
            sol = tuple(frozenset(k for k, v in M.items() if v in m) for M in Ms)
            sols.append(sol)
            if len(sols) >= CAP: break
            s.add_clause([-v if v in m else v for v in allv])
    return sols

def automorphisms(G, limit=100_000):
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G)
    out = []
    for a in GM.isomorphisms_iter():
        out.append(a)
        if len(out) >= limit: break
    return out

def classes(sols, auts, cyc):
    """orbits under graph automorphisms (O12); for model B also under a cyclic shift of time"""
    seen = set(); n = 0
    def img(sol, a, sh):
        mapped = [frozenset(frozenset(a[x] for x in k) for k in M) for M in sol]
        mapped = mapped[sh:] + mapped[:sh]
        return tuple(mapped)
    for sol in sols:
        key = tuple(sol)
        if key in seen: continue
        n += 1
        for a in auts:
            for sh in (range(len(sol)) if cyc else [0]):
                seen.add(img(sol, a, sh))
    return n

if __name__ == '__main__':
    cases = [('square torus 4x4', torus((4, 4))), ('square torus 6x6', torus((6, 6))),
             ('cubic torus 4x4x4', torus((4, 4, 4))), ('cloud 3D, 24 points, no lattice', cloud(24, 1)),
             ('cloud 3D, 40 points, no lattice', cloud(40, 2))]
    for name, G in cases:
        auts = automorphisms(G)
        print(f'\n== {name}: {G.number_of_nodes()} points, {G.number_of_edges()} links, automorphisms {len(auts)}')
        for lab, mk, cyc in (('A  M1 M1 M2 M2', model_A, False), ('B  no repetition, T=4', model_B, True)):
            e, Ms = mk(G)
            if e is None: print(f'   {lab}: no closed trajectory possible'); continue
            sols = enumerate_solutions(e, Ms)
            capped = len(sols) >= CAP
            cl = classes(sols, auts, cyc) if len(sols) * len(auts) <= 3_000_000 else None
            print(f'   {lab}: trajectories {">=" if capped else ""}{len(sols)}'
                  + (f'; distinguishable (O12) {">=" if capped else ""}{cl}' if cl is not None else '; classes not counted (too many)'))
        sys.stdout.flush()
