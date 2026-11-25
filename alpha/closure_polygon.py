"""alpha = s / pi_O^2 with pi_O the closure operator of the core's own exchange closure - rule fixed BEFORE this run (17.11.2026).

Positions: A-03 (exchange = swap of neighbours; the circle is the exchange path), A-51 (pi is not a constant: the closure
operator of the path being closed, a rational number), A-43 (every point exchanges every tick), V-219 (alpha = share of an act
that realizes the link; scale cancels), R-38 (s = directions through the active face with the eight's phase, 1/s = 14.083).

Reading (stated before this run; found by looking at mean values, so it is a fit to be tested, not a derivation):
  The closure of O is made of its exchanges with its neighbours: d neighbours = d equal chords of one closed path.
  The closure operator of a closed path of d equal chords on its circle is pi_d = d * sin(pi/d) (perimeter / diameter of the
  regular d-gon). Then alpha_i = s_i / pi_{d_i}^2 for each core with its own d_i, and alpha_W = mean of alpha_i.
Already known before this run (honest): with the MEAN degree 15.54 and the MEAN s it gives 1/alpha = 137.11 +- 0.30.
The test here: per core (own d_i, own s_i - correlated), more clouds. Pass: 1/alpha_W within 2 sigma of 137.036.
Not used in the code: 1/137.
"""
import numpy as np, sys
from core1d import run
A = []
for k, seed in enumerate((21, 22, 23, 24)):
    deg, s1, s2, cab = run(120000, seed, samples=10**9)
    s = s2 / 2
    pid = deg * np.sin(np.pi / deg)
    A.append(np.column_stack([s, deg, s / pid**2])); X = np.concatenate(A)
    aw = X[:, 2]; m = aw.mean(); se = aw.std(ddof=1) / np.sqrt(len(aw))
    print(f'clouds {k+1}: cores {len(X)}  1/s = {1/X[:,0].mean():.3f}  mean degree {X[:,1].mean():.3f}  '
          f'alpha_W = {m:.6f} +- {se:.6f}  1/alpha_W = {1/m:.3f} +- {se/m**2:.3f}'); sys.stdout.flush()
