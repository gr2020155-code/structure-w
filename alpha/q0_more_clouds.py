"""q0 over more Poisson clouds - rules: q0_more_clouds_rules_2026-11-26.md (R-113).
Same rule as placement_poisson.py (R-76); new seeds 109-140; main result over 44 clouds (101-108, 11, 21, 31, 41, 109-140)."""
import numpy as np, sys
from multiprocessing import Pool
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
    res = np.array(res)
    return seed, res.mean(), res.std(ddof=1) / np.sqrt(len(res))

if __name__ == '__main__':
    old = {101: 0.12245, 102: 0.12310, 103: 0.12266, 104: 0.12302, 105: 0.12227, 106: 0.12271, 107: 0.12323, 108: 0.12316,
           11: 0.12336, 21: 0.12326, 31: 0.12315, 41: 0.12262}            # R-76, R-73 (printed values)
    new, ses = {}, []
    with Pool(6) as p:
        for seed, m, se in p.imap_unordered(run, range(109, 141)):
            new[seed] = m; ses.append(se); print(f'seed {seed}: q0 = {m:.5f} +- {se:.5f}'); sys.stdout.flush()
    def stats(v):
        v = np.array(v); return v.mean(), v.std(ddof=1), v.std(ddof=1) / np.sqrt(len(v))
    mn, sdn, sen = stats(list(new.values()))
    m76, sd76, se76 = stats([old[s] for s in range(101, 109)])
    m, sd, se = stats(list(old.values()) + list(new.values()))
    print(f'\n32 new clouds: mean {mn:.5f}, sd {sdn:.5f}, error of mean {sen:.5f}')
    print(f'R-76 (8 clouds): mean {m76:.5f} +- {se76:.5f};  new - R-76 = {mn - m76:+.5f} = {(mn - m76) / np.hypot(sen, se76):+.1f} sigma')
    print(f'spread between clouds {sdn:.5f} vs error inside one cloud {np.mean(ses):.5f} (ratio {sdn / np.mean(ses):.2f})')
    print(f'MAIN (44 clouds): q0 = {m:.5f} +- {se:.5f}, sd {sd:.5f} ({sd / m * 100:.2f} %)')
    inv = 137.354 + (1 / m**2 - 66.291); dq = 2 / m**3 * se; err = np.hypot(0.121, dq)
    print(f'1/alpha far = {inv:.3f} +- {err:.3f}  (q0 part +- {dq:.3f})  vs 137.036: {inv - 137.036:+.3f} = {(inv - 137.036) / err:+.1f} sigma')
    print(f'seed 61 (T = 200, 1000 ends, not in main): 0.12299')
