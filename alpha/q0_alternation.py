"""Robustness of R-52 (alternation, A-70) - rules: q0_alternation_rules_2026-11-25.md. Same rule as alternation.py,
only the clouds change. Reports the mean pair remainder (1 + n.n')/2 = candidate q0."""
import numpy as np, sys
from scipy.spatial import Delaunay

def cloud(kind, seed, N=120000):
    rng = np.random.default_rng(seed)
    if kind == 'poisson':
        return rng.random((N, 3)), rng
    m = round(N ** (1 / 3)); g = np.stack(np.meshgrid(*[np.arange(m)] * 3, indexing='ij'), -1).reshape(-1, 3).astype(float)
    g += rng.normal(0, 0.3, g.shape)
    return g / m, rng

def run(kind, seed, T=60, ends_n=3000):
    X, rng = cloud(kind, seed); tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
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
    lo, hi = X.min(0) + 0.25 * np.ptp(X, 0), X.max(0) - 0.25 * np.ptp(X, 0)
    inner = np.where(np.all((X > lo) & (X < hi), axis=1))[0]
    ends = rng.choice(inner, ends_n, replace=False); res = []
    for o in ends:
        nbrs = nb[ind[o]:ind[o + 1]]; prev_q = None; prev_n = None
        for t in range(T):
            q = rng.choice([x for x in nbrs if x != prev_q]); opt = options(o, q)
            n = opt[rng.integers(len(opt))] if prev_n is None else opt[np.argmin(opt @ prev_n)]
            if prev_n is not None: res.append((1 + n @ prev_n) / 2)
            prev_n = n; prev_q = q
    res = np.array(res); return res.mean(), res.std(ddof=1) / np.sqrt(len(res)), len(res)

cases = [('poisson', 11, 60), ('poisson', 21, 60), ('poisson', 31, 60), ('poisson', 41, 60),
         ('lattice s=0.3', 51, 60), ('poisson', 61, 200)]
vals = []
for kind, seed, T in cases:
    m, se, n = run('poisson' if kind == 'poisson' else 'lattice', seed, T, 3000 if T == 60 else 1000)
    vals.append(m)
    print(f'{kind:14s} seed {seed} T {T:3d}: pairs {n:6d}  q0 = {m:.5f} +- {se:.5f}  1/q0 = {1/m:.3f}  1/alpha_step = 1/q0^2 = {1/m**2:.2f}')
    sys.stdout.flush()
v = np.array(vals); print(f'\nspread over clouds: min {v.min():.5f}, max {v.max():.5f}, relative spread {(v.max()-v.min())/v.mean()*100:.2f} %')
for F, lab in ((71.71, 'current quark masses'), (69.94, 'light quarks 300 MeV')):
    print(f'1/alpha far = 1/q0^2 + {F} ({lab}): {1/v.max()**2 + F:.2f} .. {1/v.min()**2 + F:.2f}   (data 137.036)')
