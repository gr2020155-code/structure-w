"""q from the side of the outward exchange - rule fixed BEFORE running (18.11.2026).
Positions: A-65 (charge = uncompensated outward part), A-63 (side = the contradiction-free one), V-231 (charge = loop
orientation), A-31 (indistinguishable options - equal weight), A-43 (every point exchanges every tick).
Reading: the free end of a charged form makes one outward exchange per tick (electron: one free end per tick).
Its side (ring sector around the outward bond x 2 orientations: 2n options, n = ring length) is fixed by zero
contradictions; with no preferred option it coincides with the charge loop's orientation with probability 1/(2n);
the other outcomes cancel each other (A-65). Net outward share per tick q = E[1/(2n)] over outward bonds
(random place, random neighbour). With unit coupling (V-235) alpha = q^2. Report whatever comes out. 1/137 not used.
"""
import numpy as np
from scipy.spatial import Delaunay
rng = np.random.default_rng(7)
X = rng.random((120000, 3)); tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
vt = [[] for _ in range(len(X))]
for t, s in enumerate(tri.simplices):
    for v in s: vt[v].append(t)
inner = np.where(np.all((X > 0.2) & (X < 0.8), axis=1))[0]
pick = inner
vals = []
for o in pick:
    nbrs = nb[ind[o]:ind[o + 1]]; q = rng.choice(nbrs)
    n = sum(1 for t in vt[o] if q in tri.simplices[t])
    vals.append(1 / (2 * n))
v = np.array(vals); q = v.mean(); se = v.std(ddof=1) / np.sqrt(len(v))
print(f'outward bonds {len(v)}: q = E[1/(2n)] = {q:.5f} +- {se:.5f} (1/q = {1/q:.3f});  alpha = q^2 -> 1/alpha = {1/q**2:.1f} +- {2*se/q**3:.1f}')
print('needed for 1/alpha = 137.036: q = 0.08542 (1/q = 11.706)')
