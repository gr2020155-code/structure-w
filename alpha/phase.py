"""Phase rule for the 1D core (fixed before running).
Eight of O: O<->A, O<->A, O<->B, O<->B (electron v1.1, definition 'Eight').
K1: at a tick O exchanges with one partner only. An attempt that reaches O through the face O|A
can become a change only at the two ticks when O exchanges with A; through O|B only at the other two.
Share per tick, averaged over the 4 ticks = (Omega_A + Omega_B) / 2 / (4 pi).
Also recorded (input for the next question, not used in the share): the number of points around the
links O-A and O-B (common Delaunay neighbours of the two ends) = candidates for the 'side' of an exchange (P3).
"""
import numpy as np
from scipy.spatial import Delaunay
from core1d import run

for n, seed in ((30000, 1), (60000, 2), (120000, 3)):
    deg, s1, s2, cab = run(n, seed)
    ph = s2 / 2
    sem = ph.std(ddof=1) / np.sqrt(len(ph))
    print(f'N={n:6d}  share with phase {ph.mean():.4f} ± {sem:.4f}   1/share {1/ph.mean():.1f}')

# ring size of a link (common neighbours of its two ends) in the same kind of cloud
rng = np.random.default_rng(5); X = rng.random((60000, 3)); tri = Delaunay(X)
ind, nb = tri.vertex_neighbor_vertices
inner = np.where(np.all((X > 0.2) & (X < 0.8), axis=1))[0][:5000]
rs = []
for o in inner:
    no = set(nb[ind[o]:ind[o+1]])
    for q in list(no)[:3]:
        nq = set(nb[ind[q]:ind[q+1]]); rs.append(len(no & nq))
rs = np.array(rs)
print(f'points around a link: mean {rs.mean():.2f}, min {rs.min()}, max {rs.max()}')
