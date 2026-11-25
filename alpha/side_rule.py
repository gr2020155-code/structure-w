"""Side of exchange as a factor of the junction share - rule fixed BEFORE running (17.11.2026).

Author's positions used:
  A-03  exchange = swap of neighbouring points; the circle is the exchange path; the side is a property of the act.
  A-63  an act "knows" its side because only one side gives zero contradictions.
  A-31  indistinguishable options cannot have different weight.
  V-219 alpha = share of an act of the light part that realizes the link (scale cancels).
Already computed (R-38, same definitions reused from core1d.py): directions through faces O|A, O|B with the eight's phase,
  s = 0.07101 +- 0.00015, 1/s = 14.083.

Rule (primary, the only one):
  * The exchange O<->Q (Q = A or B) goes along a circle with the bond O-Q as diameter (A-03).
    Without a lattice the directions around the bond are the ring of points that form Delaunay tetrahedra with the edge O-Q
    (ring length n = number of tetrahedra around the edge). The circle passes through one sector of this ring, in one of two
    orientations: 2n options for the side.
  * Exactly one option is contradiction-free (A-63); a priori no option is preferred (A-31), so over the ensemble the
    contradiction-free side coincides with the side the junction needs ((-) loop to (+) end) with probability 1/(2n).
  * Share per core: alpha_W = [ Omega_A/(2 n_A) + Omega_B/(2 n_B) ] / 2 / (4 pi)   (phase 1/2 as in R-33).
Report alpha_W and 1/alpha_W whatever they are. Known before running (stated honestly): the target 1/alpha = 137.036,
  i.e. a missing factor 9.73; rough expectation before running: mean ring length about 5-6, so 1/alpha_W of order 130-170.
Not used in the code: 1/137.
"""
import numpy as np, sys
from scipy.spatial import Delaunay, Voronoi
from core1d import solid_angle_poly

def run(n, seed):
    rng = np.random.default_rng(seed); X = rng.random((n, 3))
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
    vor = Voronoi(X)
    ridge = {}
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        ridge[(p, q)] = rv; ridge[(q, p)] = rv
    # tetrahedra containing each vertex, to get the ring around an edge
    vt = [[] for _ in range(n)]
    for t, simp in enumerate(tri.simplices):
        for v in simp: vt[v].append(t)
    def ring(o, q):
        ts = [t for t in vt[o] if q in tri.simplices[t]]
        return len(ts)                      # number of tetrahedra around the edge = ring length
    inner = np.where(np.all((X > 0.2) & (X < 0.8), axis=1))[0]
    out = []
    for o in inner:
        nbrs = nb[ind[o]:ind[o + 1]]
        dirs = X[nbrs] - X[o]; dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        C = dirs @ dirs.T; np.fill_diagonal(C, 2)
        i, j = np.unravel_index(np.argmin(C), C.shape)
        a, b = nbrs[i], nbrs[j]
        ok = True; om = []; rn = []
        for q in (a, b):
            rv = ridge.get((o, q))
            if rv is None or -1 in rv: ok = False; break
            om.append(solid_angle_poly(vor.vertices[rv], X[o])); rn.append(ring(o, q))
        if not ok: continue
        s = (om[0] + om[1]) / 2 / (4 * np.pi)
        aw = (om[0] / (2 * rn[0]) + om[1] / (2 * rn[1])) / 2 / (4 * np.pi)
        out.append((s, aw, rn[0], rn[1]))
    return np.array(out)

if __name__ == '__main__':
    allv = []
    for k, seed in enumerate((11, 12, 13)):
        r = run(120000, seed); allv.append(r); A = np.concatenate(allv)
        s, aw, ring_n = A[:, 0], A[:, 1], A[:, 2:].ravel()
        se = aw.std(ddof=1) / np.sqrt(len(aw))
        print(f'clouds {k+1}: cores {len(A)}  s = {s.mean():.5f} (1/s {1/s.mean():.3f})  '
              f'ring length mean {ring_n.mean():.3f} (min {ring_n.min():.0f}, max {ring_n.max():.0f})  '
              f'alpha_W = {aw.mean():.6f} +- {se:.6f}  1/alpha_W = {1/aw.mean():.2f} +- {se/aw.mean()**2:.2f}')
        sys.stdout.flush()
    vals, cnt = np.unique(ring_n, return_counts=True)
    print('ring length distribution:', {int(v): round(c / len(ring_n), 4) for v, c in zip(vals, cnt)})
