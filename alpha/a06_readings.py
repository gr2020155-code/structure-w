"""Two readings of A-06 (rules: a06_readings_rules_2026-11-25.md). Rule R-52 (A-70) on Poisson vs equalised clouds."""
import numpy as np, sys
from scipy.spatial import Delaunay, cKDTree

def equalise(X, iters):
    s = len(X) ** (-1 / 3)                     # mean spacing in the unit cube
    for _ in range(iters):
        tree = cKDTree(X); pairs = tree.query_pairs(s, output_type='ndarray')
        d = X[pairs[:, 1]] - X[pairs[:, 0]]; r = np.linalg.norm(d, axis=1); r[r == 0] = 1e-12
        push = ((s - r) / 2)[:, None] * d / r[:, None] * 0.5
        disp = np.zeros_like(X); np.add.at(disp, pairs[:, 1], push); np.add.at(disp, pairs[:, 0], -push)
        X = np.clip(X + disp, 0, 1)
    return X

def nn_cv(X):
    dd, _ = cKDTree(X).query(X, 2); d = dd[:, 1]
    inner = np.all((X > 0.25) & (X < 0.75), axis=1); d = d[inner]
    return d.std() / d.mean()

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

for label, iters, seed in (('P1 Poisson', 0, 91), ('P2 equalised 10', 10, 91), ('P2 equalised 40', 40, 91)):
    rng = np.random.default_rng(seed); X = rng.random((120000, 3))
    if iters: X = equalise(X, iters)
    cv = nn_cv(X); m, se = run(X, rng)
    print(f'{label:16s}: NN-distance CV = {cv:.3f}; q0 = {m:.5f} +- {se:.5f}; 1/q0^2 = {1/m**2:.2f}; '
          f'1/alpha_W(r_Z) = {1/m**2 + 63.249:.2f} (data 129.72, diff {1/m**2 + 63.249 - 129.72:+.2f})')
    sys.stdout.flush()
