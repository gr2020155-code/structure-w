"""Poisson placement derived from positions (rules: placement_rules_2026-11-25.md): precise q0 over 8 new clouds."""
import numpy as np, sys
from scipy.spatial import Delaunay
def run(seed, T=60, ends_n=3000):
    rng = np.random.default_rng(seed); X = rng.random((120000, 3)); tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
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
    return np.mean(res)
vals = []
for seed in range(101, 109):
    q = run(seed); vals.append(q)
    print(f'seed {seed}: q0 = {q:.5f}'); sys.stdout.flush()
v = np.array(vals); m = v.mean(); sd = v.std(ddof=1); se = sd / np.sqrt(len(v))
inv_rZ = 1 / m**2 + 63.249; d_inv = 2 / m**3 * se
print(f'\nmean q0 = {m:.5f}, spread between clouds sd = {sd:.5f} ({sd/m*100:.2f} %), error of mean = {se:.5f}')
print(f'1/alpha_W(r_Z) = {inv_rZ:.2f} +- {d_inv:.2f}   (data 129.72 +- 0.05; diff {inv_rZ-129.72:+.2f} = {(inv_rZ-129.72)/np.hypot(d_inv,0.05):+.1f} sigma)')
print(f'1/alpha far = {inv_rZ + 7.32:.2f} +- {d_inv:.2f}   (data 137.036)')
print(f'q0 needed for exact agreement: 0.12266')
