"""Charge sign in the q0 rule (rules: charge_sign_rules_2026-11-25.md). Same rule as R-52/R-71 (A-70)."""
import numpy as np
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
ends = rng.choice(inner, 1000, replace=False); R = []
for o in ends:
    nbrs = nb[ind[o]:ind[o + 1]]; prev_q = None; prev_n = None
    for t in range(60):
        q = rng.choice([x for x in nbrs if x != prev_q]); opt = options(o, q)
        n = opt[rng.integers(len(opt))] if prev_n is None else opt[np.argmin(opt @ prev_n)]
        if prev_n is not None: R.append((n + prev_n) / 2)
        prev_n = n; prev_q = q
R = np.array(R)
z1 = (R * R).sum(1).mean()
U = rng.normal(size=(20, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
z3 = np.mean([(R @ u) ** 2 for u in U])
for lab, q0 in (('Z1 magnitude, sign from loop (A-15, A-17)', z1), ('Z3 projection on direction to observer', z3)):
    print(f'{lab:45s}: q0 = {q0:.5f}; 1/alpha_W(r_Z) = {1/q0**2 + 63.249:.1f} (data 129.72)')
print(f'Z3/Z1 = {z3/z1:.3f} (isotropic residuals give 1/3)')
print(f'mean residual vector over all pairs |<r>| = {np.linalg.norm(R.mean(0)):.4f} (random-walk direction, B-248 p.2)')
