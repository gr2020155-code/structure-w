"""alpha from streams through two 1D cores - written BEFORE running (15.11.2026).

Chain (each link from a position):
  streams   - points go along their global eights from all directions (A-47, A-48), one step per tick (A-01);
              locally straight and isotropic (A-27); density: one point per cell, cell volume = a^3 (A-06: a point
              has size one step, no gaps)  -> rho * a^3 = 1
  acceptance- a 1D core O (axis A-O-B) interacts only with streams crossing its faces O|A, O|B (Voronoi faces, C1, B-148);
              by K1 face O|A is open only in the 2 ticks of 4 when O exchanges with A (phase, B-148) -> factor 1/2 per core
  shared act- a stream point that leaves the light core through an open face and enters the heavy core through an open
              face is an act of both parts (A-45). Equal charges |q_p| = |q_e|: both cores have the same 1D acceptance.
  momentum  - one act carries action hbar (A-04, E3); over one step a: momentum hbar/a per act
  force     - F = (shared acts per tick) * (hbar/a) * (c/a);  Coulomb F = alpha*hbar*c/r^2
              =>  alpha = (shared acts per tick) * (r/a)^2
  counting  - isotropic intensity I = rho*a/(4 pi) per area per sr per tick; straight lines from face f to face g:
              N = I * S_f * S_g * cos_f * cos_g / r^2  (cos_f>0 outward from e toward p, cos_g>0 inward into p)
              =>  alpha = (rho a^3) * <sum_fg S_f S_g cos_f cos_g> / (4 pi a^4) * (1/2)^2
Number 1/137 is not used. Result is reported whatever it is.
"""
import numpy as np
from scipy.spatial import Delaunay, Voronoi

def faces_of_cores(n, seed, samples):
    rng = np.random.default_rng(seed); X = rng.random((n, 3))
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices; vor = Voronoi(X)
    ridge = {}
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        ridge[(p, q)] = rv; ridge[(q, p)] = rv
    inner = np.where(np.all((X > 0.2) & (X < 0.8), axis=1))[0]
    pick = rng.choice(inner, size=min(samples, len(inner)), replace=False)
    out = []
    for o in pick:
        nbrs = nb[ind[o]:ind[o + 1]]
        dirs = X[nbrs] - X[o]; dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        C = dirs @ dirs.T; np.fill_diagonal(C, 2)
        i, j = np.unravel_index(np.argmin(C), C.shape)
        fs = []
        for q in (nbrs[i], nbrs[j]):
            rv = ridge.get((o, q))
            if rv is None or -1 in rv: fs = None; break
            P = vor.vertices[rv]; c = P.mean(0)
            nrm = X[q] - X[o]; nrm /= np.linalg.norm(nrm)            # face normal, outward from O
            # polygon area (fan around centroid, vertices ordered by angle in the face plane)
            u = P[0] - c; u /= np.linalg.norm(u); v = np.cross(nrm, u)
            ang = np.arctan2((P - c) @ v, (P - c) @ u); P = P[np.argsort(ang)]
            area = 0.5 * np.linalg.norm(sum(np.cross(P[k] - c, P[(k + 1) % len(P)] - c) for k in range(len(P))))
            fs.append((area, nrm))
        if fs: out.append(fs)
    return out, n          # unit cube: rho = n

if __name__ == '__main__':
    rng = np.random.default_rng(99)
    for n, seed in ((30000, 1), (60000, 2), (120000, 3)):
        cores, rho = faces_of_cores(n, seed, 3000)
        a = rho ** (-1 / 3)                                     # rho a^3 = 1
        vals = []
        for _ in range(200000):
            e = cores[rng.integers(len(cores))]; p = cores[rng.integers(len(cores))]
            d = rng.normal(size=3); d /= np.linalg.norm(d)      # direction from e to p (isotropic)
            s = 0.0
            for Sf, nf in e:
                cf = nf @ d
                if cf <= 0: continue
                for Sg, ng in p:
                    cg = -(ng @ d)
                    if cg <= 0: continue
                    s += Sf * Sg * cf * cg
            vals.append(s / a ** 4)
        vals = np.array(vals)
        alpha = vals.mean() / (4 * np.pi) * 0.25
        sem = vals.std(ddof=1) / np.sqrt(len(vals)) / (4 * np.pi) * 0.25
        mean_face = np.mean([S for c in cores for S, _ in c]) / a ** 2
        print(f'N={n:6d}: mean area of an axis face {mean_face:.4f} a^2;  alpha_W = {alpha:.6f} ± {sem:.6f}  (1/alpha_W = {1/alpha:.1f})')
    print('data: alpha = 0.0072973526 (1/137.036)')
