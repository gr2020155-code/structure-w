"""Driver for global_eights.py (same formalisation F1-F8): orbit counts over several clouds, optional repeat search.
Usage: python3 global_eights_runs.py dim n seeds limit cap [q3]
"""
import sys, time
import numpy as np
from ortools.sat.python import cp_model
import os, global_eights
if os.environ.get('GRID3'): global_eights.GRID3 = int(os.environ['GRID3'])
from global_eights import cloud, build, extract, transforms, canon, verify, solve, random_repeat_rate

dim, n = int(sys.argv[1]), int(sys.argv[2]); seeds = [int(s) for s in sys.argv[3].split(',')]
limit, cap = float(sys.argv[4]), int(sys.argv[5]); q3 = 'q3' in sys.argv[6:]; f10 = 'f10' in sys.argv[6:]
scale = 1000 if dim == 2 else global_eights.GRID3 ** 2

def visit_gaps(M, n, O):
    pos = list(range(n)); vis = {i: [] for i in range(n)}
    for t in range(len(M)):
        for i in range(n):
            if pos[i] == O: vis[i].append(t)
        new = pos[:]
        for (a, b) in M[t]:
            ia, ib = pos.index(a), pos.index(b); new[ia], new[ib] = b, a
        pos = new
    return sorted({v[1] - v[0] for v in vis.values()})
for seed in seeds:
    X, E, O, hull = cloud(n, seed, dim); T = 2 * n
    degs = np.bincount(np.array(E).ravel(), minlength=n)
    print(f'== dim {dim}, n {n}, seed {seed}: edges {len(E)}, mean degree {degs.mean():.2f}, O degree {degs[O]}, T {T}'); sys.stdout.flush()
    t0 = time.time(); orbits, blocks, status = [], [], None
    while len(orbits) < cap:
        m, x, rep, T = build(X, E, O, hull, False, True, scale, blocks=blocks, f10=f10); s, st = solve(m, limit)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): status = s.StatusName(st); break
        M = extract(s, x, E, T); verify(X, E, O, hull, M, False, True)
        if f10: assert visit_gaps(M, n, O) == [n], 'F10 violated'
        orbits.append(canon(M)); blocks += transforms(M)
    print(f'   F1-F8{"+F10" if f10 else ""}: distinct solutions (up to time shift and reversal) {len(orbits)}'
          + (f', search ended: {status}' if status else ' (cap reached)') + f'  [{time.time()-t0:.0f}s]'); sys.stdout.flush()
    if q3 and orbits:
        rr = random_repeat_rate(E, n, 400, seed)
        lo, hi, proved, R, t0 = 0, None, False, 1, time.time()
        while True:
            m, x, rep, T = build(X, E, O, hull, False, False, scale, at_least=R); s, st = solve(m, limit)
            if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                M = extract(s, x, E, T); lo = max(lo, verify(X, E, O, hull, M, False, False)); R = max(2 * R, lo + 1)
                if R > n * T // 2: break
            else:
                hi, proved = R, (st == cp_model.INFEASIBLE); break
        while hi is not None and hi - lo > 1:
            R = (lo + hi) // 2
            m, x, rep, T = build(X, E, O, hull, False, False, scale, at_least=R); s, st = solve(m, limit)
            if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                M = extract(s, x, E, T); lo = max(lo, verify(X, E, O, hull, M, False, False))
            else:
                hi = R; proved = proved and (st == cp_model.INFEASIBLE)
        print(f'   Q3 (F1-F7 without F8): repeats from 0 up to >= {lo} = {2*lo/(n*T):.3f} per point-tick'
              + (f', {hi} impossible (proved)' if hi is not None and proved else ', upper bound not proved')
              + f'; random matchings {rr:.3f}  [{time.time()-t0:.0f}s]'); sys.stdout.flush()
