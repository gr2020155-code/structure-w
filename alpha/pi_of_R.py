"""pi of closure as a function of distance R (in steps) around the axis of a core - step 1 of
'find what fits, then why, then test elsewhere'.
Definition: shell S_R = points at R steps (hops) from O; poles = the shell point closest in direction to O->A and
to O->B; L_R = minimal number of shell points separating the poles within the shell (links = Delaunay links
between shell points); pi(R) = L_R / (2R) (diameter through O = 2R steps)."""
import numpy as np, sys
import networkx as nx
from scipy.spatial import Delaunay
rng = np.random.default_rng(7); X = rng.random((200000, 3))
tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
N = lambda v: nb[ind[v]:ind[v + 1]]
inner = np.where(np.all((X > 0.3) & (X < 0.7), axis=1))[0]
picks = rng.choice(inner, 400, replace=False)
res = {R: [] for R in range(1, 8)}
for o in picks:
    nbrs = N(o); d = X[nbrs] - X[o]; d /= np.linalg.norm(d, axis=1)[:, None]
    C = d @ d.T; np.fill_diagonal(C, 2); i, j = np.unravel_index(np.argmin(C), C.shape)
    ua, ub = d[i], d[j]
    dist = {int(o): 0}; front = [int(o)]
    for R in range(1, 8):
        nxt = []
        for u in front:
            for v in N(u):
                v = int(v)
                if v not in dist: dist[v] = R; nxt.append(v)
        front = nxt
        S = set(front)
        if len(S) < 4: continue
        G = nx.Graph(); G.add_nodes_from(S)
        for u in S:
            for v in N(u):
                if int(v) in S: G.add_edge(u, int(v))
        P = np.array(list(S)); dirs = X[P] - X[o]; dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        pa, pb = int(P[np.argmax(dirs @ ua)]), int(P[np.argmax(dirs @ ub)])
        if pa == pb or G.has_edge(pa, pb) or not nx.has_path(G, pa, pb): continue
        res[R].append(nx.node_connectivity(G, pa, pb) / (2 * R))
for R in range(1, 8):
    v = np.array(res[R])
    if len(v): print(f'R = {R} steps: cores {len(v)}, pi(R) = {v.mean():.4f} ± {v.std(ddof=1)/np.sqrt(len(v)):.4f}')
    sys.stdout.flush()
