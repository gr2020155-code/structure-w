"""Medium under A-40/A-41: points never repeat a pair (every medium exchange is a change).
Question: does such a medium exist for many steps without a lattice? Depth-first search with random order.
Fixed before running: Delaunay cloud, every place in one exchange per step (D4), no pair of POINTS twice."""
import random, sys
import networkx as nx
from uniq import cloud
sys.setrecursionlimit(10000)

def random_matchings(G, at, used, rng, tries=200):
    nodes = list(G.nodes())
    for _ in range(tries):
        free = set(nodes); M = []
        order = nodes[:]; rng.shuffle(order)
        ok = True
        while free:
            u = min(free, key=lambda x: sum(1 for v in G[x] if v in free))   # most constrained first
            cand = [v for v in G[u] if v in free and frozenset((at[u], at[v])) not in used]
            if not cand: ok = False; break
            v = rng.choice(cand); M.append((u, v)); free -= {u, v}
        if ok: yield M

def run(n, seed, T):
    G = cloud(n, seed); rng = random.Random(seed)
    at = {v: v for v in G}; used = set(); steps = 0
    def rec(t, at, used):
        if t == T: return t
        best = t
        for M in random_matchings(G, at, used, rng, tries=30):
            nat = dict(at); nu = set(used)
            for u, v in M: nu.add(frozenset((at[u], at[v]))); nat[u], nat[v] = at[v], at[u]
            r = rec(t + 1, nat, nu)
            best = max(best, r)
            if r == T: return T
            break   # greedy with limited backtracking
        return best
    return rec(0, at, used), n * (n - 1) // 2

for n, seed in ((40, 1), (80, 2), (160, 3)):
    reached, pairs = run(n, seed, 200)
    print(f'cloud {n} points: no-repeat medium reached {reached} steps (pairs of points possible {pairs}, used {reached * n // 2})'); sys.stdout.flush()
