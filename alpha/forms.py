"""What forms appear when a closed region runs out of new pairs?
Rules fixed BEFORE running (A-40, A-41, A-43, D4):
  every place in exactly one exchange per step (perfect matching of the Delaunay places graph);
  a pair of points that has never exchanged = change; an already used pair = repeat (act along an existing link);
  at each step the number of repeats is the smallest possible (repeat only where unavoidable);
  ties broken at random (seeded).
Recorded: fraction of repeats per step; per point the sequence of partners; patterns:
  point loop = same pair at t and t+1 ("there and back");
  eight      = point loop with A then point loop with B (B != A) for the same point: X<->A, X<->A, X<->B, X<->B;
  persistent = the same point keeps repeating pairs for many steps.
"""
import random, sys, collections
import networkx as nx
from uniq import cloud

def step(G, at, used, rng):
    H = nx.Graph()
    for u, v in G.edges():
        new = frozenset((at[u], at[v])) not in used
        H.add_edge(u, v, weight=(2.0 if new else 1.0) + rng.random() * 1e-3)
    M = nx.max_weight_matching(H, maxcardinality=True)
    assert 2 * len(M) == G.number_of_nodes(), 'no perfect matching of places'
    return M

def run(n, seed, T):
    G = cloud(n, seed)
    rng = random.Random(seed); at = {v: v for v in G}; used = set()
    partner = collections.defaultdict(list)      # point -> partner point per step
    rep = []
    for t in range(T):
        M = step(G, at, used, rng); r = 0
        nat = dict(at)
        for u, v in M:
            a, b = at[u], at[v]; p = frozenset((a, b))
            if p in used: r += 1
            used.add(p); partner[a].append(b); partner[b].append(a)
            nat[u], nat[v] = b, a
        at = nat; rep.append(r / (n / 2))
    return G, partner, rep, used

def patterns(partner, T, t0):
    loops = eights = 0; eight_pts = set(); loop_pts = set()
    for x, seq in partner.items():
        for t in range(t0, T - 1):
            if seq[t] == seq[t + 1]: loops += 1; loop_pts.add(x)
        for t in range(t0, T - 3):
            if seq[t] == seq[t + 1] and seq[t + 2] == seq[t + 3] and seq[t] != seq[t + 2]:
                eights += 1; eight_pts.add(x)
    return loops, eights, loop_pts, eight_pts

if __name__ == '__main__':
    for n, seed in ((24, 5), (40, 1), (80, 2)):
        T = 3 * n
        G, partner, rep, used = run(n, seed, T)
        first = next((t for t, r in enumerate(rep) if r > 0), None)
        q = lambda a, b: sum(rep[a:b]) / max(1, b - a)
        print(f'\n== cloud {n} points, {T} steps; pairs of points possible {n*(n-1)//2}, realized {len(used)}')
        print(f'   first repeat at step {first} (bound n-1 = {n-1}); share of repeats: '
              f'steps 0-{n//2}: {q(0,n//2):.2f}, {n//2}-{n}: {q(n//2,n):.2f}, {n}-{2*n}: {q(n,2*n):.2f}, {2*n}-{3*n}: {q(2*n,3*n):.2f}')
        L, E, lp, ep = patterns(partner, T, first or 0)
        steps_after = T - (first or 0)
        print(f'   after the first repeat: point loops {L}, eights {E}; points that made an eight {len(ep)} of {n}')
        # persistence: longest run of consecutive repeated partners for a point
        longest = 0
        for x, seq in partner.items():
            run_ = 0
            for t in range(1, T):
                run_ = run_ + 1 if seq[t] == seq[t - 1] else 0
                longest = max(longest, run_)
        print(f'   longest run of the same pair on consecutive steps: {longest + 1}')
        sys.stdout.flush()
