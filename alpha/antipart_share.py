"""Anti-part share in a W nucleon: taking ends of transfers (R-153).
Rules: antipart_share_rules_2026-11-28.md (committed before the run, 6e79b28; pre-run arithmetic note ab49263).
usage: python -P -E -X utf8 antipart_share.py"""
from math import sqrt

K = {'u': 2, 'd': 1}                       # axes of a part (V-304)
NUCLEONS = {'p': 'uud', 'n': 'udd'}         # A-23, V-305
FD, FE = 0.161, 0.017                       # needed by data (R-152)


def counts(parts):
    axes = sum(K[q] for q in parts)
    missing = sum(3 - K[q] for q in parts)  # one transfer per missing axis (V-343 p.3, V-355 p.2)
    return len(parts), axes, missing


def share(branch, phi=1.0):
    """isoscalar average; phi = fraction of the full complement."""
    num = den = 0.0
    for parts in NUCLEONS.values():
        npart, axes, miss = counts(parts)
        t = phi * miss                      # transfers present
        if branch == 'F1':                  # measure: axis; each end half an axis
            anti, total = t / 2, axes + t
        elif branch == 'F2':                # measure: object; end = half a transfer
            anti, total = t / 2, npart + t
        else:                               # F3: each end a full object
            anti, total = t, npart + 2 * t
        num += anti / 2; den += total / 2
    return num / den


for name, parts in NUCLEONS.items():
    n, a, mi = counts(parts)
    print(f'{name} ({parts}): parts {n}, axes in parts {a}, missing axes (transfers) {mi}, total {a + mi}')
print(f'needed by data f = {FD} +- {FE}')
for br in ('F1', 'F2', 'F3'):
    f = share(br)
    R = (1 + 2 * f) / (3 - 2 * f)
    pull = (f - FD) / FE
    lo, hi = 0.0, 1.0                        # phi needed by data (bisection; share grows with phi)
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if share(br, mid) < FD else (lo, mid)
    print(f'{br}: full complement f = {f:.4f} ({pull:+.1f} sigma), nubar/nu = {R:.3f} (data 0.493 +- 0.016) -> '
          f'{"PASS" if abs(pull) <= 2 else "FAIL"}; complement share needed by data phi = {lo:.3f}')
