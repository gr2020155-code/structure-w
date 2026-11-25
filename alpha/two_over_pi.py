"""Factor 2/pi of the slope (rules: two_over_pi_rules_2026-11-25.md).
Loss per e-fold = alpha * H / (2 pi D^2), H = c * p * <nA*nB>; compared with alpha/(3 pi).
Part 1 checks the geometry factor 1/(2 pi) per e-fold by Monte Carlo; part 2 counts the branches."""
from math import pi, log, sqrt, e
import random

# 1. geometry: H points uniform on a sphere of diameter D through the charge;
#    share per cell summed over the layers of one e-fold = sum over points in the shell of 1/(4 pi d^2)
random.seed(7); D = 2000.0; R = D / 2; H = 2 * 10**6
for x1 in (0.003, 0.01, 0.1):
    r1, r2 = x1 * D, x1 * D * e
    s = 0.0; k = 0
    for _ in range(H):
        v = [random.gauss(0, 1) for _ in range(3)]; nv = sqrt(sum(t * t for t in v))
        d = sqrt((R * v[0] / nv) ** 2 + (R * v[1] / nv) ** 2 + (R * v[2] / nv - R) ** 2)
        if r1 <= d < r2: s += 1 / (4 * pi * d * d); k += 1
    print(f'geometry, e-fold from {x1} D: measured {s:.4f} ({k} points), predicted H/(2 pi D^2) = {H / (2 * pi * D * D):.4f}')

# 2. branches of the count of held points
n = 2000
mean_pairs = sum(a * (n - a) for a in range(1, n)) / (n - 1)
print(f'\n<nA*nB>/D^2 = {mean_pairs / n**2:.6f} (1/6 = {1/6:.6f})')
target = 1 / (3 * pi)
for c in (1, 2):
    for p in (1, 2):
        loss = c * p * mean_pairs / (2 * pi * n**2)
        print(f'   c = {c}, p = {p}: loss per e-fold / alpha = {loss:.6f}  (data 1/(3 pi) = {target:.6f})  '
              f'ratio {loss / target:.4f}  {"HIT" if abs(loss / target - 1) < 0.01 else "miss"}')
