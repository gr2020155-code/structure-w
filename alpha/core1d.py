"""Share of directions a 1D core (A-O-B) accepts among 3D directions, places without a lattice (P16).
Fixed before running:
  places  - Poisson cloud in 3D (no lattice), neighbours = Delaunay (as in electron v1.1, section 9)
  core    - O with exactly two links A, B on one axis: the pair of neighbours of O closest to opposite
  share 1 - links: 2 / (number of neighbours of O)  (which of O's possible exchange partners the form uses)
  share 2 - directions: solid angle of the two Voronoi faces O|A and O|B seen from O, divided by 4*pi
            (the set of directions in which an arriving attempt meets A or B rather than another neighbour)
Only interior points are used (no boundary effects). 1/137 is not used.
"""
import numpy as np
from scipy.spatial import Delaunay, Voronoi, ConvexHull

def solid_angle_poly(P, O):
    """solid angle of a planar convex polygon P (k x 3) seen from O (fan of triangles, Van Oosterom-Strackee)"""
    c = P.mean(axis=0); tot = 0.0
    k = len(P)
    # order vertices around centroid
    n = np.cross(P[1] - P[0], P[2] - P[0]); n /= np.linalg.norm(n)
    u = P[0] - c; u /= np.linalg.norm(u); v = np.cross(n, u)
    ang = np.arctan2((P - c) @ v, (P - c) @ u); P = P[np.argsort(ang)]
    for i in range(k):
        a, b, cc = c - O, P[i] - O, P[(i + 1) % k] - O
        la, lb, lc = map(np.linalg.norm, (a, b, cc))
        num = abs(np.dot(a, np.cross(b, cc)))
        den = la * lb * lc + np.dot(a, b) * lc + np.dot(a, cc) * lb + np.dot(b, cc) * la
        tot += 2 * np.arctan2(num, den)
    return tot

def run(n=60000, seed=0, samples=3000):
    rng = np.random.default_rng(seed); X = rng.random((n, 3))
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
    vor = Voronoi(X)
    ridge = {}
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        ridge[(p, q)] = rv; ridge[(q, p)] = rv
    inner = np.where(np.all((X > 0.2) & (X < 0.8), axis=1))[0]
    pick = rng.choice(inner, size=min(samples, len(inner)), replace=False)
    s1, s2, deg, cosab = [], [], [], []
    for o in pick:
        nbrs = nb[ind[o]:ind[o + 1]]; d = len(nbrs)
        dirs = X[nbrs] - X[o]; dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        C = dirs @ dirs.T; np.fill_diagonal(C, 2)
        i, j = np.unravel_index(np.argmin(C), C.shape)          # most opposite pair = the axis
        a, b = nbrs[i], nbrs[j]
        om = 0.0; ok = True
        for q in (a, b):
            rv = ridge.get((o, q))
            if rv is None or -1 in rv: ok = False; break
            om += solid_angle_poly(vor.vertices[rv], X[o])
        if not ok: continue
        deg.append(d); s1.append(2 / d); s2.append(om / (4 * np.pi)); cosab.append(C[i, j])
    return np.array(deg), np.array(s1), np.array(s2), np.array(cosab)

if __name__ == '__main__':
    for n, seed in ((30000, 1), (60000, 2), (120000, 3)):
        deg, s1, s2, cab = run(n, seed)
        sem = lambda x: x.std(ddof=1) / np.sqrt(len(x))
        print(f'N={n:6d}  cores={len(deg)}  neighbours {deg.mean():.2f}  '
              f'share(links) {s1.mean():.4f}±{sem(s1):.4f}  share(directions) {s2.mean():.4f}±{sem(s2):.4f}  '
              f'angle A-O-B {np.degrees(np.arccos(-cab)).mean():.1f} deg off straight  '
              f'min share(dir) {s2.min():.4f}')
    print('1/share(directions) for comparison with 1/alpha:', round(1 / s2.mean(), 2))
