"""One charge loop without a second particle - count fixed BEFORE running (18.11.2026).
Positions: A-03 (exchange = swap of neighbours), A-43 (every point exchanges every tick), A-07 + V-147 (eight: O<->A, O<->A,
O<->B, O<->B), A-41 (a repeated pair inside a closed form is not a change - goes by the existing link), A-58 (points change).
Count per cycle of the eight (4 ticks):
  inner  - acts between points of the form (closure, 'magnetic side', circulation);
  outer  - acts of form points with medium points (link outward, 'electric side');
  return - outer acts that are undone by the next act of the same place with the same partner (the point goes out and
           comes back: self-link). With d neighbours of the end place and no preferred partner (A-31): probability 1/d.
Question: is any of these shares of order alpha = 1/137? d from the Poisson clouds already computed (poly_poisson_*.npy).
"""
import numpy as np, glob
sched = [('O', 'A'), ('O', 'A'), ('O', 'B'), ('O', 'B')]
inner = outer = 0; outer_by = {'A': 0, 'B': 0, 'O': 0}
for t, (p, q) in enumerate(sched):
    for pt in 'OAB':
        if pt in (p, q): inner += 1
        else: outer += 1; outer_by[pt] += 1
tot = inner + outer
print(f'per eight cycle: point-acts {tot}, inner {inner} ({inner/tot:.3f}), outer {outer} ({outer/tot:.3f}); outer by point {outer_by}')
print(f'end point: outer share {outer_by["A"]/4:.3f} of its ticks; centre O: {outer_by["O"]/4:.3f}')
d = np.concatenate([np.load(f)[:, 0] for f in sorted(glob.glob('poly_poisson_*.npy'))])
print(f'self-link (out and back with the same partner) per outer act: E[1/d] = {np.mean(1/d):.5f} = 1/{1/np.mean(1/d):.2f} (Poisson, {len(d)} places)')
print(f'self-link per eight cycle (2 ends, 1 out-and-back pair each): {2*np.mean(1/d)/tot:.5f} of point-acts = 1/{tot/(2*np.mean(1/d)):.1f}')
print('alpha = 1/137.036 = 0.0072974 for comparison')
