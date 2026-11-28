"""p-bar + n at rest (in deuterium): prediction by the rule of R-148 (meeting region D = 2 hbar/(m_pi c)) (R-149).
Rules: pbarn_rules_2026-11-28.md (committed before the run, acb0475); counting as in ppbar.py.
usage: python -P -E -X utf8 pbarn.py"""
from math import factorial, pi
import numpy as np

MP, ME, MPC, MP0 = 0.93827208943, 0.00051099895069, 0.13957039, 0.1349768   # GeV
MN, BD = 0.93956542052, 0.00222457
E = MP + MN - BD
NEV = 200_000
rng = np.random.default_rng(147)
DATA = {'ch': (3.0, 0.2), 'neu': (2.0, 0.2)}
PRONG = {0: 4.1, 2: 43.2, 4: 48.6, 6: 4.1}


def pstar(M, m1, m2):
    a = (M**2 - (m1 + m2)**2) * (M**2 - (m1 - m2)**2)
    return np.sqrt(np.clip(a, 0, None)) / (2 * M)


def rand_dir(n):
    c = rng.uniform(-1, 1, n); phi = rng.uniform(0, 2 * pi, n); s = np.sqrt(1 - c**2)
    return np.stack([s * np.cos(phi), s * np.sin(phi), c], axis=1)


def boost(P, beta):
    """P: (N, k, 4) four-vectors; beta: (N, 3)."""
    b2 = np.sum(beta**2, axis=1)[:, None]
    g = 1 / np.sqrt(1 - b2)
    bp = np.einsum('nkj,nj->nk', P[:, :, 1:], beta)
    En = P[:, :, 0]
    coef = np.where(b2 > 0, (g - 1) * bp / np.where(b2 > 0, b2, 1), 0) + g * En
    out = np.empty_like(P)
    out[:, :, 0] = g * (En + bp)
    out[:, :, 1:] = P[:, :, 1:] + coef[:, :, None] * beta[:, None, :]
    return out


def integral(masses):
    """I_n = integral of prod d3p_i delta(E - sum E_i) delta3(sum p_i), in GeV^(3n-4)."""
    m = np.array(masses); n = len(m); T = E - m.sum()
    if T <= 0:
        return 0.0
    if n == 2:
        p = pstar(E, m[0], m[1]); e1, e2 = np.sqrt(p**2 + m[0]**2), np.sqrt(p**2 + m[1]**2)
        return float(4 * pi * p * e1 * e2 / E)
    r = np.sort(rng.uniform(0, 1, (NEV, n - 2)), axis=1)
    cm = np.cumsum(m)
    M = np.empty((NEV, n)); M[:, 0] = m[0]; M[:, -1] = E
    M[:, 1:n - 1] = cm[1:n - 1] + r * T
    p = pstar(M[:, 1:], M[:, :-1], m[1:])               # (NEV, n-1): p_i for i = 2..n
    w = pi**(n - 1) * 2**(n - 2) / E * T**(n - 2) / factorial(n - 2) * np.prod(p, axis=1)
    # build momenta: start with particles 1, 2 in the rest frame of M_2
    P = np.zeros((NEV, n, 4))
    u = rand_dir(NEV); q = p[:, 0][:, None] * u
    P[:, 0, 1:] = q; P[:, 0, 0] = np.sqrt(p[:, 0]**2 + m[0]**2)
    P[:, 1, 1:] = -q; P[:, 1, 0] = np.sqrt(p[:, 0]**2 + m[1]**2)
    for i in range(2, n):                                # add particle i+1 (0-based i)
        u = rand_dir(NEV); pi_ = p[:, i - 1]
        esub = np.sqrt(pi_**2 + M[:, i - 1]**2)
        beta = (pi_ / esub)[:, None] * u
        P[:, :i] = boost(P[:, :i], beta)
        P[:, i, 1:] = -pi_[:, None] * u; P[:, i, 0] = np.sqrt(pi_**2 + m[i]**2)
    e = P[:, :, 0]
    assert np.allclose(e.sum(axis=1), E, atol=1e-9) and np.allclose(P[:, :, 1:].sum(axis=1), 0, atol=1e-9)
    return float(np.mean(w * np.prod(2 * e, axis=1)))


SETS = []
for k in range(0, 7):
    for j in range(0, 14):
        if (2 * k + 1) * MPC + j * MP0 < E and 2 * k + 1 + j >= 2:
            SETS.append((k, j, integral([MPC] * (2 * k + 1) + [MP0] * j)))
print(f'E = {E * 1000:.1f} MeV, sets counted: {len(SETS)}, events per set: {NEV}')
D = 2 / MPC
V = 4 / 3 * pi * D**3
S = np.array([(V / (2 * pi)**3)**(2 * k + j) * I / (factorial(k) * factorial(k + 1) * factorial(j)) for k, j, I in SETS])
S /= S.sum()
npl = sum(s * k for s, (k, j, _) in zip(S, SETS)); nmi = sum(s * (k + 1) for s, (k, j, _) in zip(S, SETS))
ne = sum(s * j for s, (k, j, _) in zip(S, SETS))
print(f'D = {D * 0.1973269804:.3f} fm (as R-148)')
print(f'PREDICTION: <n+> = {npl:.2f}, <n-> = {nmi:.2f}, <n+-> = {npl + nmi:.2f}, <n0> = {ne:.2f}, total {npl + nmi + ne:.2f}')
pr = {c: 100 * sum(s for s, (k, j, _) in zip(S, SETS) if 2 * k + 1 == c) for c in (1, 3, 5, 7)}
print('  charged prongs, %: ' + ', '.join(f'{c}: {v:.1f}' for c, v in pr.items()))
nt = {}
for s, (k, j, _) in zip(S, SETS):
    nt[2 * k + 1 + j] = nt.get(2 * k + 1 + j, 0) + s
print('  total pions, %: ' + ', '.join(f'{n}: {100 * v:.1f}' for n, v in sorted(nt.items()) if v > 0.001))
