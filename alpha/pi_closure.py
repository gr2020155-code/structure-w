"""pi of closure at the distance of one step, and alpha = s/pi^2 - definitions written BEFORE running (15.11.2026).

Positions used: length is a number of steps (A-01); pi is the closure operator, a rational number depending on the
closed path (A-51); acceptance of a 1D core O with axis A-O-B (B-148, B-164).
Definition (per core, no free parameter):
  closure path = the shortest closed path of steps around the axis at the distance of one step from O:
                 the shortest cycle among the neighbours of O (links = Delaunay links between neighbours of O)
                 that separates A from B (every path from A to B over the neighbours of O must cross it);
                 its length L is the minimal number of neighbours of O separating A from B.
  diameter    = A-O-B = 2 steps.
  pi_i        = L_i / 2   (a rational number for every core)
  share s_i   = (Omega_A + Omega_B) / (4 pi_euclid) / 2  (as in phase.py: directions through the axis faces, phase 1/2)
  alpha_W     = mean over cores of s_i / pi_i^2
Reported whatever the result. 1/137 is not used in the computation.
"""
import numpy as np, sys
import networkx as nx
from scipy.spatial import Delaunay, Voronoi
from core1d import solid_angle_poly

def run(n, seed):
    rng = np.random.default_rng(seed); X = rng.random((n, 3))
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices; vor = Voronoi(X)
    ridge = {}
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        ridge[(p, q)] = rv; ridge[(q, p)] = rv
    nbset = lambda v: set(nb[ind[v]:ind[v + 1]].tolist())
    inner = np.where(np.all((X > 0.2) & (X < 0.8), axis=1))[0]
    out = []
    for o in inner:
        nbrs = nb[ind[o]:ind[o + 1]]
        dirs = X[nbrs] - X[o]; dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        C = dirs @ dirs.T; np.fill_diagonal(C, 2)
        i, j = np.unravel_index(np.argmin(C), C.shape); a, b = int(nbrs[i]), int(nbrs[j])
        om = 0.0; ok = True
        for q in (a, b):
            rv = ridge.get((o, q))
            if rv is None or -1 in rv: ok = False; break
            om += solid_angle_poly(vor.vertices[rv], X[o])
        if not ok: continue
        Nset = nbset(o)
        L = nx.Graph(); L.add_nodes_from(Nset)
        for u in Nset:
            for v in nbset(u) & Nset:
                L.add_edge(u, v)
        if L.has_edge(a, b): continue                       # axis ends adjacent: no separating ring
        k = nx.node_connectivity(L, a, b)                    # min number of neighbours separating A from B
        s = om / (4 * np.pi) / 2
        out.append((k, s, len(Nset)))
    return out

if __name__ == '__main__':
    allr = []
    for k, seed in enumerate((21, 22, 23, 24)):
        allr += run(120000, seed)
        Ls = np.array([r[0] for r in allr]); S = np.array([r[1] for r in allr])
        pis = Ls / 2; aw = S / pis ** 2
        m, se = aw.mean(), aw.std(ddof=1) / np.sqrt(len(aw))
        print(f'clouds {k+1}: cores {len(allr)}; ring length L: mean {Ls.mean():.3f} (min {Ls.min()}, max {Ls.max()}); '
              f'pi = L/2 mean {pis.mean():.4f}; s mean {S.mean():.5f}; alpha_W = {m:.6f} ± {se:.6f}, 1/alpha_W = {1/m:.2f} ± {se/m**2:.2f}')
        sys.stdout.flush()
    vals, cnt = np.unique(Ls, return_counts=True)
    print('distribution of L:', dict(zip(vals.tolist(), (cnt / cnt.sum()).round(4).tolist())))
    print('data: 1/alpha = 137.036')
