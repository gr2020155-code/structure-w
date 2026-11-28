"""Ratio 'transfer acting on a transfer / transfer acting on a part' in W vs C_A/C_F (R-150).
Rules: colour_ratio_rules_2026-11-28.md (committed before the run, 0fcf6d7).
usage: python -P -E -X utf8 colour_ratio.py"""
from itertools import permutations

AX = (1, 2, 3)
TR = list(permutations(AX, 2))          # 6 oriented transfers a -> b, no diagonal (V-355)
assert len(TR) == 6


def acts_on_set(t, held):
    a, b = t
    return a in held and b not in held   # takes an axis the object has, puts it where it has none


def c_part(held):
    return sum(acts_on_set(t, held) for t in TR)


cp = {frozenset(h): c_part(set(h)) for h in [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3)]}
print('C_part by held axes:', {tuple(sorted(k)): v for k, v in cp.items()})
assert len(set(cp.values())) == 1
CP = next(iter(cp.values()))


def compose(T, S):
    """S acting on T = (i -> j). Returns list of results ('nothing' or new transfer)."""
    i, j = T; a, b = S; out = []
    if a == j:                           # S moves the delivered axis
        out.append('nothing' if b == i else (i, b))
    if b == i:                           # S fills the place T took the axis from
        out.append('nothing' if a == j else (a, j))
    return out


res = {}
for T in TR:
    i, j = T
    other = ({1, 2, 3} - {i, j}).pop()
    b1 = c_part({j}) + c_part({j, other})            # B1: two ends as two places of action
    acting = {S: compose(T, S) for S in TR if compose(T, S)}
    b3 = len(acting)                                  # B3: each acting variant once
    b2 = sum(1 for S, r in acting.items() if any(x != 'nothing' for x in r))   # B2: drop the one that annuls
    res[T] = (b1, b2, b3)
assert len(set(res.values())) == 1
B1, B2, B3 = next(iter(res.values()))
print(f'C_part = {CP}; C_transfer: B1 = {B1}, B2 = {B2}, B3 = {B3}  (same for all 6 transfers)')

DATA = [('ALEPH 4-jet', 2.20, (0.09**2 + 0.13**2) ** 0.5), ('LEP combined', 2.89 / 1.30, 2.89 / 1.30 * ((0.21 / 2.89)**2 + (0.09 / 1.30)**2) ** 0.5)]
for name, ct in (('B1 by ends', B1), ('B2 variants, keeping transfer', B2), ('B3 variants incl. annulling', B3)):
    r = ct / CP
    pulls = [(r - v) / e for _, v, e in DATA]
    ok = all(abs(p) <= 2 for p in pulls)
    print(f'{name:32s}: ratio {r:.3f}; ' + ', '.join(f'{n} {v:.2f}+-{e:.2f}: {p:+.1f} sigma' for (n, v, e), p in zip(DATA, pulls))
          + f' -> {"PASS" if ok else "FAIL"}')
print('for comparison: SU(3) Yang-Mills 9/4 = 2.25; abelian model 0')
