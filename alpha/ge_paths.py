"""How many distinct global eights (closed paths up to time shift) the points follow, and whether shared paths are
counter-streams (one path = the other traversed backwards). Saved 2D solutions, n = 8; descriptive, no new rule."""
import json
from global_eights import cloud
def paths(M, n):
    T = len(M); where = list(range(n)); tr = [[i] for i in range(n)]
    for t in range(T):
        owner = {where[i]: i for i in range(n)}
        for (a, b) in map(tuple, M[t]): where[owner[a]], where[owner[b]] = b, a
        for i in range(n): tr[i].append(where[i])
    return [p[:-1] for p in tr]
def canon(p): return min(tuple(p[s:] + p[:s]) for s in range(len(p)))
for seed in (1, 3, 5, 6):
    X, E, O, hull = cloud(8, seed, 2)
    for k, M in enumerate(json.load(open(f'ge_solutions_2d_n8_s{seed}.json'))):
        P = paths(M, 8); groups = {}
        for i, p in enumerate(P): groups.setdefault(canon(p), []).append(i)
        keys = list(groups)
        rev = sum(1 for a in keys for b in keys if a < b and canon(list(a)[::-1]) == b)
        lag = []
        for g, members in groups.items():
            starts = sorted(next(t for t in range(len(P[i])) if tuple(P[i][t:] + P[i][:t]) == g) for i in members)
            lag.append([starts[j + 1] - starts[j] for j in range(len(starts) - 1)])
        print(f'cloud {seed} sol {k+1}: distinct paths {len(groups)}, points per path {[len(v) for v in groups.values()]}, '
              f'pairs of paths that are each other reversed {rev}, time lags between points on one path {lag}')
