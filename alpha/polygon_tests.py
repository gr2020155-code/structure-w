"""Independent checks of the candidate V-223: alpha = s / pi_d^2, pi_d = d*sin(pi/d) - rules fixed BEFORE running (18.11.2026).

Candidate (found after looking at mean values on 17.11, R-47: 1/alpha_W = 136.49 +- 0.36 per core).
Checks declared now, each can refute the candidate:
  T1 universality across forms. alpha is the same for hydrogen and muonic hydrogen (data), so the light part's form must not
     change it. A cross (muon, A-07): two axes through O; O makes one act per tick (A-43), so one face is active per tick:
     s_cross = (Om_A + Om_B + Om_C + Om_D) / 4 / (4 pi). Axis (electron): s_axis = (Om_A + Om_B) / 2 / (4 pi) (as R-38).
     Second axis C-D: among the remaining neighbours the pair minimising cos(C,D) + |cos(axis CD, axis AB)|.
     Pass: 1/alpha with s_cross equals 1/alpha with s_axis within 1 %.
  T2 precision. Poisson clouds, all interior cores, per core alpha_i = s_i / pi_{d_i}^2, mean over cores.
     Pass: 1/alpha_W within 2 sigma of 137.036 at the final precision.
  T3 arrangement. The same formula on other lattice-free-like arrangements (jittered cubic, hard-core Matern II).
     Not pass/fail for W: shows whether the match is a property of the Poisson arrangement only.
Not used in the code: 1/137.
"""
import numpy as np, sys, time
from scipy.spatial import Delaunay, Voronoi, cKDTree
from core1d import solid_angle_poly

def cores(X, inner):
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
    vor = Voronoi(X)
    ridge = {}
    for (p, q), rv in zip(vor.ridge_points, vor.ridge_vertices):
        ridge[(p, q)] = rv; ridge[(q, p)] = rv
    out = []
    for o in np.where(inner)[0]:
        nbrs = nb[ind[o]:ind[o + 1]]; d = len(nbrs)
        dirs = X[nbrs] - X[o]; dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        C = dirs @ dirs.T; np.fill_diagonal(C, 2)
        i, j = np.unravel_index(np.argmin(C), C.shape)
        uab = dirs[i] - dirs[j]; uab /= np.linalg.norm(uab)
        R = np.array([k for k in range(d) if k != i and k != j])
        V = dirs[R][:, None, :] - dirs[R][None, :, :]
        nv = np.linalg.norm(V, axis=2); np.fill_diagonal(nv, 1)
        sc = C[np.ix_(R, R)] + np.abs(V @ uab) / nv
        np.fill_diagonal(sc, 9)
        a, b = np.unravel_index(np.argmin(sc), sc.shape); k, l = R[a], R[b]
        om = []
        for q in (nbrs[i], nbrs[j], nbrs[k], nbrs[l]):
            rv = ridge.get((o, q))
            if rv is None or -1 in rv: om = None; break
            om.append(solid_angle_poly(vor.vertices[rv], X[o]))
        if om is None: continue
        s_ax = (om[0] + om[1]) / 2 / (4 * np.pi); s_cr = sum(om) / 4 / (4 * np.pi)
        pid2 = (d * np.sin(np.pi / d)) ** 2
        out.append((d, s_ax, s_cr, s_ax / pid2, s_cr / pid2))
    return np.array(out)

def report(tag, A):
    d, sa, scr, aa, ac = A.T
    se = lambda x: x.std(ddof=1) / np.sqrt(len(x))
    print(f'{tag}: cores {len(A)}  mean degree {d.mean():.3f}  1/s_axis {1/sa.mean():.3f}  1/s_cross {1/scr.mean():.3f}  '
          f'1/alpha(axis) = {1/aa.mean():.3f} +- {se(aa)/aa.mean()**2:.3f}  1/alpha(cross) = {1/ac.mean():.3f} +- {se(ac)/ac.mean()**2:.3f}  '
          f'ratio cross/axis = {aa.mean()/ac.mean():.4f}')
    sys.stdout.flush()

def poisson(n, seed):
    X = np.random.default_rng(seed).random((n, 3)); return X, np.all((X > 0.2) & (X < 0.8), axis=1)

def jitter(m, sig, seed):
    g = (np.indices((m, m, m)).reshape(3, -1).T + 0.5) / m
    X = g + np.random.default_rng(seed).normal(0, sig / m, g.shape)
    return X, np.all((X > 0.2) & (X < 0.8), axis=1)

def matern2(n, r, seed):
    rng = np.random.default_rng(seed); X = rng.random((n, 3)); mark = rng.random(n)
    keep = np.ones(n, bool)
    for p, q in cKDTree(X).query_pairs(r):
        if mark[p] < mark[q]: keep[q] = False
        else: keep[p] = False
    X = X[keep]; return X, np.all((X > 0.2) & (X < 0.8), axis=1)

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'poisson':
        seeds = range(int(sys.argv[2]), int(sys.argv[3])); allA = []
        for s in seeds:
            t = time.time(); X, inn = poisson(120000, s); A = cores(X, inn); allA.append(A)
            np.save(f'poly_poisson_{s}.npy', A)
            report(f'Poisson clouds {len(allA)} (seed {s}, {time.time()-t:.0f}s)', np.concatenate(allA))
    elif mode == 'arrangements':
        for tag, (X, inn) in (('jittered cubic sigma=0.15', jitter(50, 0.15, 1)), ('jittered cubic sigma=0.30', jitter(50, 0.30, 2)),
                              ('jittered cubic sigma=0.60', jitter(50, 0.60, 3)), ('hard-core Matern II', matern2(400000, 0.012, 4))):
            t = time.time(); A = cores(X, inn); report(f'{tag} ({len(X)} points, {time.time()-t:.0f}s)', A)
