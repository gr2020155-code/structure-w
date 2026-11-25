"""Pairs +/- with the medium point holding the side.
Rules: pair_hold_rules_2026-11-25.md (written before this run).
Units: step per tact = 1 (c = 1); circle radius R = N/(2*pi) steps, D = 2R."""
from math import pi, cos, sin, dist
import random

N = 2000
R = N / (2 * pi); D = 2 * R

def path_G1():
    # both circuits on the same circle; circuit 1 (+) steps 0..N-1, circuit 2 (-) steps N..2N-1
    p = [(R * cos(2 * pi * k / N), R * sin(2 * pi * k / N), 0.0) for k in range(N)]
    return p, p
def path_G2():
    # figure eight: two circles mirrored in the plane x = 0, junction at O = origin
    c1 = [(R - R * cos(2 * pi * k / N), R * sin(2 * pi * k / N), 0.0) for k in range(N)]   # starts and ends at O
    c2 = [(-x, y, z) for x, y, z in c1]
    return c1, c2

random.seed(25112026)
def rand_X():
    while True:
        v = [random.uniform(-3 * D, 3 * D) for _ in range(3)]
        if sum(t * t for t in v) <= (3 * D) ** 2:
            return v

for name, mk in (('G1 same circle', path_G1), ('G2 figure eight', path_G2)):
    plus, minus = mk()
    order_viol = 0; spans = []; gaps = []
    for _ in range(20000):
        X = rand_X()
        ap = [k + dist(X, plus[k]) for k in range(N)]            # arrival of + born at tact k
        am = [N + k + dist(X, minus[k]) for k in range(N)]        # arrival of - born at tact N+k
        arr = ap + am
        order_viol += sum(1 for i in range(2 * N - 1) if arr[i + 1] < arr[i] - 1e-9)
        first_minus = min(range(N), key=lambda k: am[k])
        held = max((k for k in range(N) if ap[k] <= am[first_minus]), key=lambda k: ap[k])
        spans.append(dist(plus[held], minus[first_minus]))
        gaps.append((N + first_minus) - held)                     # tacts between the two births
    spans.sort()
    print(f'{name}: arrival-order violations {order_viol}; '
          f'held + born at tact {held} of 0..{N-1}, first - born at tact {N + first_minus} (last X)')
    print(f'   span S in steps: min {spans[0]:.3f}, median {spans[len(spans)//2]:.3f}, max {spans[-1]:.3f}; '
          f'S/D max {spans[-1]/D:.2e};  tacts between births: min {min(gaps)}, max {max(gaps)}')
print(f'\nFor comparison (data, R-61): half of spans below 0.15*D = {0.15*D:.0f} steps at this N.')
