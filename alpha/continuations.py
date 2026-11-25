"""How many continuations does the medium have at each step? (A-40, A-41: medium never repeats a pair of points)
Fixed before running:
  places: Delaunay cloud without lattice, n points, closed region
  step t: count C_t = number of perfect matchings of places whose point pairs are all new
  trajectory: continue with a random one of the C_t continuations, until C_t = 0 (saturation)
Hypothesis to test: C_t falls with the history and reaches exactly 1 (a unique trajectory)."""
import random, sys
from functools import lru_cache
from uniq import cloud

def count_and_pick(G, at, used, rng, cap=3_000_000):
    nodes = sorted(G.nodes())
    adj = {u: [v for v in G[u] if frozenset((at[u], at[v])) not in used] for u in nodes}
    idx = {u: i for i, u in enumerate(nodes)}
    @lru_cache(maxsize=None)
    def cnt(mask):
        if mask == 0: return 1
        # pick free place with fewest free options
        best = None
        for u in nodes:
            if mask >> idx[u] & 1:
                k = sum(1 for v in adj[u] if mask >> idx[v] & 1)
                if best is None or k < best[0]: best = (k, u)
                if k == 0: return 0
        u = best[1]; tot = 0
        for v in adj[u]:
            if mask >> idx[v] & 1:
                tot += cnt(mask & ~(1 << idx[u]) & ~(1 << idx[v]))
        return tot
    full = (1 << len(nodes)) - 1
    C = cnt(full)
    if C == 0: return 0, None
    # uniform random continuation by sequential sampling
    mask = full; M = []
    while mask:
        u = min((u for u in nodes if mask >> idx[u] & 1), key=lambda x: sum(1 for v in adj[x] if mask >> idx[v] & 1))
        opts = [v for v in adj[u] if mask >> idx[v] & 1]
        ws = [cnt(mask & ~(1 << idx[u]) & ~(1 << idx[v])) for v in opts]
        r = rng.random() * sum(ws); acc = 0
        for v, w in zip(opts, ws):
            acc += w
            if r <= acc: break
        M.append((u, v)); mask &= ~(1 << idx[u]) & ~(1 << idx[v])
    cnt.cache_clear()
    return C, M

if __name__ == '__main__':
    for n, seed in ((16, 7), (20, 3), (24, 5)):
        G = cloud(n, seed)
        for run in range(3):
            rng = random.Random(100 * seed + run)
            at = {v: v for v in G}; used = set(); seq = []
            while True:
                C, M = count_and_pick(G, at, used, rng)
                seq.append(C)
                if C == 0: break
                nat = dict(at)
                for u, v in M: used.add(frozenset((at[u], at[v]))); nat[u], nat[v] = at[v], at[u]
                at = nat
            print(f'cloud {n} pts, run {run}: C_t = {seq}'); sys.stdout.flush()
