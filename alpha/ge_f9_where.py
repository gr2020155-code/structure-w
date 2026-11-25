"""Where do same-orientation acts happen? (descriptive, saved solutions)"""
import sys, json, numpy as np
from global_eights import cloud, area
for seed in (1, 3, 5, 6):
    X, E, O, hull = cloud(8, seed, 2); n = 8
    sols = json.load(open(f'ge_solutions_2d_n8_s{seed}.json'))
    for k, Ms in enumerate(sols):
        M = [set(tuple(e) for e in Mt) for Mt in Ms]; T = len(M)
        pos = list(range(n)); traj = [[i] for i in range(n)]
        for t in range(T):
            new = pos[:]
            for (a, b) in M[t]:
                ia, ib = pos.index(a), pos.index(b); new[ia], new[ib] = b, a
            pos = new
            for i in range(n): traj[i].append(pos[i])
        sgn = np.zeros((n, T)); star = True
        for i in range(n):
            tr = traj[i]; t1, t2 = [t for t in range(T) if tr[t] == O]
            s1 = np.sign(sum(area(X, O, tr[t], tr[t + 1]) for t in range(t1, t2)))
            for t in range(T):
                s = s1 if t1 <= t < t2 else -s1; sgn[i, t] = s
                a = area(X, O, tr[t], tr[t + 1])
                if a != 0 and np.sign(a) != s: star = False
        same_at_O = same_off = 0
        for t in range(T):
            at = {traj[i][t]: i for i in range(n)}
            for (a, b) in M[t]:
                if sgn[at[a], t] == sgn[at[b], t]:
                    if O in (a, b): same_at_O += 1
                    else: same_off += 1
        print(f'seed {seed} sol {k+1}: same-orientation acts at O {same_at_O}, away from O {same_off}; '
              f'every step of every loop turns the same way about O (star-shaped loops): {star}')
