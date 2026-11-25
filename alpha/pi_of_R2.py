"""pi(R) by winding cycles (corrected definition; the shell-cut definition in pi_of_R.py does not measure an equator:
the minimal separator inside a one-hop shell stays ~6 for every R, so pi(R) fell as 1/R - a broken measure).
Definition: points at >= R steps from O; the shortest closed path of steps that winds once around the axis line A-O-B;
L_R = its number of steps; pi(R) = L_R / (2R)."""
import numpy as np, sys, collections
from scipy.spatial import Delaunay
rng = np.random.default_rng(7); X = rng.random((60000, 3))
tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
N = lambda v: nb[ind[v]:ind[v + 1]]
edge = np.mean([np.linalg.norm(X[v] - X[0]) for v in N(0)])
inner = np.where(np.all((X > 0.35) & (X < 0.65), axis=1))[0]
picks = rng.choice(inner, 60, replace=False)
RMAX = 8
res = collections.defaultdict(list)
for o in picks:
    nbrs = N(o); d = X[nbrs] - X[o]; d /= np.linalg.norm(d, axis=1)[:, None]
    C = d @ d.T; np.fill_diagonal(C, 2); i, j = np.unravel_index(np.argmin(C), C.shape)
    u = d[i] - d[j]; u /= np.linalg.norm(u)                     # axis direction
    e1 = np.cross(u, [1, 0, 0]); e1 = e1 if np.linalg.norm(e1) > .1 else np.cross(u, [0, 1, 0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    dist = {int(o): 0}; front = [int(o)]
    for R in range(1, RMAX + 6):
        nxt = []
        for a in front:
            for b in N(a):
                b = int(b)
                if b not in dist: dist[b] = R; nxt.append(b)
        front = nxt
    def theta(v):
        y = X[v] - X[o]; return np.arctan2(y @ e2, y @ e1), np.linalg.norm(y - (y @ u) * u)
    for R in range(1, RMAX + 1):
        allowed = {v for v, k in dist.items() if R <= k <= R + 5}
        th = {}
        for v in allowed:
            t, rp = theta(v)
            if rp > 0.5 * edge: th[v] = t
        allowed = set(th)
        def step(a, b):
            ta, tb = th[a], th[b]
            if ta > np.pi / 2 and tb < -np.pi / 2: return 1
            if ta < -np.pi / 2 and tb > np.pi / 2: return -1
            return 0
        starts = [v for v in allowed if dist[v] == R]
        starts = sorted(starts, key=lambda v: abs((X[v] - X[o]) @ u))[:6]
        best = None
        for s in starts:
            seen = {(s, 0): 0}; q = collections.deque([(s, 0)])
            while q:
                a, lv = q.popleft(); da = seen[(a, lv)]
                if best is not None and da + 1 >= best: break
                for b in N(a):
                    b = int(b)
                    if b not in allowed: continue
                    nl = lv + step(a, b)
                    if abs(nl) > 1: continue
                    if (b, nl) == (s, 1) or (b, nl) == (s, -1):
                        best = da + 1 if best is None else min(best, da + 1); q.clear(); break
                    if (b, nl) not in seen:
                        seen[(b, nl)] = da + 1; q.append((b, nl))
        if best: res[R].append(best / (2 * R))
    print('core done', len(res[1]), {R: round(float(np.mean(v)),3) for R, v in res.items()}, flush=True)
for R in sorted(res):
    v = np.array(res[R]); print(f'R = {R}: cores {len(v)}, pi(R) = {v.mean():.4f} ± {v.std(ddof=1)/np.sqrt(len(v)):.4f}')
