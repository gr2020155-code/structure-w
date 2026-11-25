"""Arrow of time (author, 15.11.2026): no return; every exchange changes the global system.
Reading tested: a realized link is a pair of POINTS (not places). Points move: an exchange swaps their places.
Constraints fixed before running:
  C1/K1  every step is a perfect matching of places
  ARROW  no pair of points exchanges twice (every exchange realizes a new link; saturation grows)
  no return condition (C4 removed)
Counts all trajectories of T steps from one start; O12 classes are not taken (start fixes labels).
"""
import itertools, sys
import networkx as nx
from uniq import torus, cloud

def matchings(G):
    nodes = sorted(G.nodes()); out = []
    def rec(free, cur):
        if not free: out.append(tuple(cur)); return
        u = min(free)
        for v in G[u]:
            if v in free:
                rec(free - {u, v}, cur + [(u, v)])
    rec(frozenset(nodes), []); return out

def count(G, T, cap=10**7):
    Ms = matchings(G)
    n = 0
    def rec(t, at, used):   # at: place -> point
        nonlocal n
        if n >= cap: return
        if t == T: n += 1; return
        for M in Ms:
            pairs = [frozenset((at[u], at[v])) for u, v in M]
            if any(p in used for p in pairs): continue
            nat = dict(at)
            for u, v in M: nat[u], nat[v] = at[v], at[u]
            rec(t + 1, nat, used | set(pairs))
    rec(0, {v: v for v in G.nodes()}, frozenset()); return n, len(Ms)

if __name__ == '__main__':
    for name, G in (('square torus 4x4', torus((4, 4))), ('cloud 3D, 16 points', cloud(16, 7))):
        for T in (1, 2, 3, 4):
            n, m = count(G, T)
            print(f'{name}: matchings {m}; trajectories of {T} steps without repeated pairs: {n}'); sys.stdout.flush()
