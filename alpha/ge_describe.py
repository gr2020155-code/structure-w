"""Descriptive statistics of all F1-F8 solutions of one cloud (no new rule; for the report)."""
import sys, json, numpy as np
from ortools.sat.python import cp_model
from global_eights import cloud, build, extract, transforms, canon, verify, solve, area
from global_eights_f9 import f9_violations
dim, n, seed, limit = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
X, E, O, hull = cloud(n, seed, dim); scale = 1000 if dim == 2 else 40000
orb, blocks, sols = [], [], []
while True:
    m, x, rep, T = build(X, E, O, hull, False, True, scale, blocks=blocks); s, st = solve(m, limit)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): print('search ended:', s.StatusName(st)); break
    M = extract(s, x, E, T); verify(X, E, O, hull, M, False, True); sols.append(M); blocks += transforms(M)
for k, M in enumerate(sols):
    pos = list(range(n)); traj = [[i] for i in range(n)]
    for t in range(T):
        new = pos[:]
        for (a, b) in M[t]:
            ia, ib = pos.index(a), pos.index(b); new[ia], new[ib] = b, a
        pos = new
        for i in range(n): traj[i].append(pos[i])
    lens = []
    for i in range(n):
        v = [t for t in range(T) if traj[i][t] == O]; lens.append(tuple(sorted((v[1] - v[0], T - (v[1] - v[0])))))
    occ = [next(i for i in range(n) if traj[i][t] == O) for t in range(T)]
    bad, tot = f9_violations(X, O, M) if dim == 2 else (None, None)
    print(f'solution {k+1}: loop lengths per point {sorted(lens)}; points at O by tick {occ}; '
          f'acts with partners on loops of opposite orientation {bad} of {tot}')
json.dump([[sorted(map(list, Mt)) for Mt in M] for M in sols], open(f'ge_solutions_{dim}d_n{n}_s{seed}.json', 'w'))
