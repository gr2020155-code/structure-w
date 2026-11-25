"""Combinatorial route to alpha - checks fixed BEFORE running (18.11.2026).
C1 Small integers of the junction itself (W positions): 2 loops of the eight, 3 points of the electron, 4 ticks of the eight /
   points of the proton, 5 loops of the atom, 7 points of the atom. How many expressions of up to 4 of them (+ - * / ^)
   land on 1/alpha = 137.036 within the precision W can claim (2e-5, terms of order alpha^2)? Many hits = this route
   cannot give evidence, whatever expression one "finds".
C2 Strict constraint (V-198 + A-39): if beta and iota are ratios of whole counts with beta^2 + iota^2 = 1 exactly,
   (beta, iota, 1) is a Pythagorean triple: beta = 2mn/(m^2+n^2) or (m^2-n^2)/(m^2+n^2). Simplest families vs 137.036,
   and how many primitive triples with hypotenuse <= C hit within 2e-5 (density = no selection without a rule for m, n).
"""
import itertools, math
target = 137.035999084; tol = 2e-5
S = [2, 3, 4, 5, 7]
def combine(a, b):
    out = {a + b, a - b, a * b}
    if b != 0: out.add(a / b)
    if 0 < abs(a) < 50 and abs(b) <= 8 and float(b).is_integer(): out.add(a ** b)
    return out
vals = {1: {x: 1 for x in S}}
exprs = {k: {} for k in range(1, 5)}
exprs[1] = {float(x) for x in S}
for k in range(2, 5):
    acc = set()
    for i in range(1, k):
        for a in exprs[i]:
            for b in exprs[k - i]:
                for v in combine(a, b):
                    if abs(v) < 1e7: acc.add(round(v, 9))
    exprs[k] = acc
for k in range(1, 5):
    hits = sorted(v for v in exprs[k] if abs(v / target - 1) < tol)
    near = sum(1 for v in exprs[k] if abs(v / target - 1) < 1e-3)
    print(f'C1 up to {k} numbers: distinct values {len(exprs[k])}, within 1e-3: {near}, within 2e-5: {len(hits)} {hits[:5]}')
print('C1 example: 2^7 + 3^2 + 1/(4*7) =', 2**7 + 3**2 + 1/(4*7), ' relative error', (2**7 + 3**2 + 1/28) / target - 1)
# C2
fam1 = [(m, (m*m + 1) / (2*m)) for m in (273, 274, 275)]
print('C2 family n=1 (beta = 2m/(m^2+1)):', [(m, round(x, 4)) for m, x in fam1])
fam2 = [(n, (2*n*n + 2*n + 1) / (2*n + 1)) for n in (136, 137)]
print('C2 family m=n+1 (beta = (2n+1)/(2n^2+2n+1)):', [(n, round(x, 4)) for n, x in fam2])
for C in (10**4, 10**5, 10**6, 10**7):
    cnt = 0; best = None
    mmax = int(math.isqrt(C)) + 1
    for m in range(2, mmax):
        for n in range(1, m):
            if (m - n) % 2 == 0 or math.gcd(m, n) != 1: continue
            c = m*m + n*n
            if c > C: break
            for a in (2*m*n, m*m - n*n):
                inv = c / a
                if abs(inv / target - 1) < tol:
                    cnt += 1
                    if best is None or abs(inv - target) < abs(best[0] - target): best = (inv, m, n, a, c)
    print(f'C2 primitive triples with c <= {C:>8}: hits within 2e-5: {cnt}; best {best}')
