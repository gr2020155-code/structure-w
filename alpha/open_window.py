"""Is a local region open when everything moves? (author 15.11: new points are already there; motion brings them)
Fixed before running: big cloud N=1000 (Delaunay, no lattice), same rules as forms.py
(D4; new pair = change; repeat only where unavoidable). Window = 60 places nearest the centre.
Compared with the closed 60-point cloud under the same rules.
Motion here is only what the exchanges themselves produce (no global eight imposed) - the weakest case."""
import random, sys, collections
import numpy as np, networkx as nx
from uniq import cloud
from forms import step, run as run_closed

def run_big(N=1000, seed=1, T=150, W=60):
    rng0 = np.random.default_rng(seed); X = rng0.random((N, 3))
    from scipy.spatial import Delaunay
    import itertools
    G = nx.Graph()
    for s in Delaunay(X).simplices:
        for a, b in itertools.combinations(s, 2): G.add_edge(int(a), int(b))
    d = np.linalg.norm(X - 0.5, axis=1); window = set(np.argsort(d)[:W].tolist())
    rng = random.Random(seed); at = {v: v for v in G}; used = set()
    seen = set(at[v] for v in window); rep_w = []
    for t in range(T):
        M = step(G, at, used, rng); r = k = 0; nat = dict(at)
        for u, v in M:
            p = frozenset((at[u], at[v])); inw = u in window or v in window
            if inw: k += 1; r += p in used
            used.add(p); nat[u], nat[v] = at[v], at[u]
        at = nat; rep_w.append(r / k); seen |= set(at[v] for v in window)
        if t % 25 == 24:
            print(f'   step {t+1}: repeats in window (last 25 steps) {np.mean(rep_w[-25:]):.3f}; distinct points that have been in the window {len(seen)}'); sys.stdout.flush()
    return rep_w

print('closed 60-point cloud, same rules:')
G, partner, rep, used = run_closed(60, 4, 150)
for a in range(0, 150, 25): print(f'   steps {a+1}-{a+25}: repeats {np.mean(rep[a:a+25]):.3f}')
print('window of 60 places inside 1000 points:')
run_big()
