"""(b) placement: rule R-52 (A-70) on Matern II hard-core clouds (rules: narrow_band_rules_2026-11-25.md)."""
import numpy as np, sys
from math import pi
from scipy.spatial import Delaunay, cKDTree

def matern2(N0, keep, seed):
    rng = np.random.default_rng(seed)
    lam = N0; lv = {0.5: 1.5936, 0.25: 3.9207}[keep]          # lambda*V solving (1-exp(-x))/x = keep
    d = (lv / lam / (4 / 3 * pi)) ** (1 / 3)
    X = rng.random((N0, 3)); mk = rng.random(N0); tree = cKDTree(X)
    keepm = np.ones(N0, bool)
    for i, nb in enumerate(tree.query_ball_point(X, d)):
        for j in nb:
            if j != i and mk[j] < mk[i]: keepm[i] = False; break
    return X[keepm], rng, d

def run(X, rng, T=60, ends_n=3000):
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
    vt = [[] for _ in range(len(X))]
    for t, s in enumerate(tri.simplices):
        for v in s: vt[v].append(t)
    def options(o, q):
        u = X[q] - X[o]; u /= np.linalg.norm(u)
        ring = {v for t in vt[o] if q in tri.simplices[t] for v in tri.simplices[t]} - {o, q}
        N = []
        for r in ring:
            w = X[r] - X[o]; w -= (w @ u) * u; w /= np.linalg.norm(w); n = np.cross(u, w); N += [n, -n]
        return np.array(N)
    inner = np.where(np.all((X > 0.25) & (X < 0.75), axis=1))[0]
    ends = rng.choice(inner, ends_n, replace=False); res = []
    for o in ends:
        nbrs = nb[ind[o]:ind[o + 1]]; prev_q = None; prev_n = None
        for t in range(T):
            q = rng.choice([x for x in nbrs if x != prev_q]); opt = options(o, q)
            n = opt[rng.integers(len(opt))] if prev_n is None else opt[np.argmin(opt @ prev_n)]
            if prev_n is not None: res.append((1 + n @ prev_n) / 2)
            prev_n = n; prev_q = q
    res = np.array(res); return res.mean(), res.std(ddof=1) / np.sqrt(len(res))

for keep, seed in ((0.5, 71), (0.5, 72), (0.25, 81)):
    X, rng, d = matern2(240000 if keep == 0.5 else 480000, keep, seed)
    m, se = run(X, rng)
    print(f'Matern II keep {keep}: points {len(X)}, min distance d = {d:.5f} (in units of mean spacing {d * len(X) ** (1/3):.3f}): q0 = {m:.5f} +- {se:.5f}, 1/q0^2 = {1/m**2:.2f}, 1/alpha_W(r_Z) = {1/m**2 + 63.249:.2f} (data 129.72)')
    sys.stdout.flush()
