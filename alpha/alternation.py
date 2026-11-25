"""Side of the free end's outward exchanges must alternate - rule fixed BEFORE running (18.11.2026).

Author (18.11): the side of each outward exchange of the free end must alternate with the previous one, as the layers
alternate for gravity (each layer passes to the next, A-64; one layer per tick, V-117).
Other positions: A-03 (exchange = circle; side - property of the act), A-63 (the side taken is the contradiction-free one),
A-40 (no return: consecutive partners differ), A-65 (charge = uncompensated part), A-39 (combination by squares).

Model:
  * places: Poisson cloud 120 000 (no lattice), neighbours - Delaunay;
  * a free end A makes one outward exchange per tick with a neighbour different from the previous partner (random);
  * side of an exchange on bond A-X: oriented normal of its circle through one sector of the ring around the bond
    (2n options, n - ring length), as in R-45/R-51;
  * alternation: the side taken is the available option most opposite to the previous side (min n_t . n_{t-1});
  * uncompensated part of a pair (by squares, A-39): |(n_{t-1} + n_t)/2|^2 = (1 + n_{t-1}.n_t)/2.
Primary result: alpha_W = mean of (1 + n_{t-1}.n_t)/2 over all pairs; compare 1/alpha_W with 137.036.
Also reported (descriptive, not a test): static remainder |mean of n_t| over the whole sequence of each end.
"""
import numpy as np, sys
from scipy.spatial import Delaunay
rng = np.random.default_rng(11)
X = rng.random((120000, 3)); tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
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
ends = rng.choice(inner, 3000, replace=False); T = 60
res, stat = [], []
for o in ends:
    nbrs = nb[ind[o]:ind[o + 1]]; prev_q = None; prev_n = None; seq = []
    for t in range(T):
        q = rng.choice([x for x in nbrs if x != prev_q]); opt = options(o, q)
        n = opt[rng.integers(len(opt))] if prev_n is None else opt[np.argmin(opt @ prev_n)]
        if prev_n is not None: res.append((1 + n @ prev_n) / 2)
        seq.append(n); prev_n = n; prev_q = q
    stat.append(np.linalg.norm(np.mean(seq, axis=0)))
res = np.array(res); m = res.mean(); se = res.std(ddof=1) / np.sqrt(len(res))
print(f'pairs {len(res)}: alpha_W = mean (1 + n.n\')/2 = {m:.5f} +- {se:.5f}  ->  1/alpha_W = {1/m:.2f} +- {se/m**2:.2f}   (data 137.036)')
print(f'distribution of the pair remainder: median {np.median(res):.4f}, 10% {np.quantile(res,0.1):.4f}, 90% {np.quantile(res,0.9):.4f}')
print(f'static remainder |mean n| over {T} ticks: mean {np.mean(stat):.4f} (random walk of residuals would give ~{np.sqrt(np.mean(res))/np.sqrt(T):.4f})')
