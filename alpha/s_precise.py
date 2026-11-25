"""Precise share s (directions accepted by a 1D core with the eight's phase) - same definitions as core1d.py/phase.py,
all interior points of each cloud instead of 3000 samples. Fixed before running. 1/alpha is not used in the code."""
import numpy as np, sys
from core1d import run
vals = []
for k, seed in enumerate((11, 12, 13, 14, 15, 16)):
    deg, s1, s2, cab = run(120000, seed, samples=10**9)
    ph = s2 / 2; vals.append(ph)
    allv = np.concatenate(vals)
    m = allv.mean(); se = allv.std(ddof=1) / np.sqrt(len(allv))
    print(f'clouds {k+1}: cores {len(allv)}, s = {m:.5f} ± {se:.5f}, 1/s = {1/m:.3f} ± {se/m**2:.3f}'); sys.stdout.flush()
