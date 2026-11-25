"""End-to-end check: from the positions of Structure W to 1/alpha, in one program (25.11.2026).

Every number is either DERIVED here from a position (A-xx) / matrix record (V-xxx), or taken from DATA (D-xx).
Run:  python -P -E -X utf8 end_to_end.py          (fast: 2 Poisson clouds for q0)
      python -P -E -X utf8 end_to_end.py --full   (8 clouds, as R-76)
Needs numpy, scipy.
"""
import sys
from math import pi, sqrt, log, exp, comb
import numpy as np
from scipy.spatial import Delaunay

FULL = '--full' in sys.argv
def say(tag, text): print(f'[{tag}] {text}')

# ---------------------------------------------------------------- DATA (only these numbers come from measurement)
hbar = 1.054571817e-34; c = 299792458.0; G = 6.67430e-11            # CODATA 2022
MeV = 1.602176634e-13 / c**2
M_Z = 91187.6; m_t = 172570.0                                         # MeV (PDG)
inv_alpha_far_data = 137.035999177                                    # D-48
inv_alpha_MZ_data = 128.9                                             # D-48, "effective" (by energy), +-0.05
GAMMA = 0.5772156649015329                                            # Euler constant (mathematics, not data)
say('DATA', f'G, hbar, c (CODATA); M_Z = {M_Z} MeV, m_t = {m_t} MeV (PDG); 1/alpha far = {inv_alpha_far_data}; 1/alpha(M_Z) = {inv_alpha_MZ_data} (D-48)')

# ---------------------------------------------------------------- 1. step of space (A-80 -> c_g = 1, V-271)
l_P = sqrt(hbar * G / c**3); a = 2 * pi * l_P
say('1 step', f'c_g = 1 (A-80, V-271) -> a = 2*pi*l_P = {a:.4e} m')

# ---------------------------------------------------------------- 2. slope per unit-charge species (V-288, V-289, A-31, A-88)
n = 2000
mean_pairs_n = sum(k * (n - k) for k in range(1, n)) / (n - 1) / n**2         # <nA nB>/D^2 at n = 2000: (n+1)/(6n)
mean_pairs = 1 / 6                                                             # limit for a real form: D/a ~ 1e22 steps
c_roles, p_points = 2, 2                                                       # A-31 partners equal; A-88 act holds both points
geom = 1 / (2 * pi)                                                            # ring 2*pi*r on the sphere of the form (V-283, R-68)
kappa = c_roles * p_points * mean_pairs * geom * 2                             # *2: 1/alpha = q^-2 (V-287)
say('2 slope', f'<nA nB>/D^2 = {mean_pairs_n:.6f} at n = {n} (finite size), limit 1/6 for D/a ~ 1e22;'
               f' kappa = 2*2*(1/6)*(1/2pi)*2 = {kappa:.6f}  (2/(3pi) = {2/(3*pi):.6f})')

# ---------------------------------------------------------------- 3. edge of screening (V-286, A-31, A-39, A-84)
def lam(nA, nB):                                               # decay length: 1/sqrt(<k^2>) over the time in each loop
    rA, rB = nA / 2, nB / 2; xA, xB = nA / n, nB / n
    return 1 / sqrt(xA / rA**2 + xB / rB**2)
items = [(k * (n - k), lam(k, n - k)) for k in range(1, n)]
W = sum(w for w, _ in items)
E = sum(w * log(l / n) for w, l in items) / W - GAMMA                          # ln(r_eff/D); exponential gives -gamma
say('3 edge', f'ln(r_eff/D) = {E:.4f}  (real-distance target -(5/6+gamma+ln2) = {-(5/6 + GAMMA + log(2)):.4f})')

# ---------------------------------------------------------------- 4. charges and colours from three axes (V-304, V-305)
axes = {k: (comb(3, k), k / 3) for k in range(4)}            # k axes: variants, |Q|
NcQ2 = {k: v * q * q for k, (v, q) in axes.items()}
say('4 axes', '; '.join(f'k={k}: variants {v}, |Q| {q:.3f}, Nc*Q^2 {NcQ2[k]:.4f}' for k, (v, q) in axes.items()))
gen = {'lepton': NcQ2[3], 'up': NcQ2[2], 'down': NcQ2[1]}
sum_active_at_rZ = 3 * gen['lepton'] + 2 * gen['up'] + 3 * gen['down']       # e,mu,tau; u,c; d,s,b (t is above M_Z)
say('4 axes', f'sum Nc*Q^2 of species active at the Z distance = {sum_active_at_rZ:.4f} (20/3 = {20/3:.4f}); top: {gen["up"]:.4f}')

# ---------------------------------------------------------------- 5. placement and q0 (V-307, V-295, V-303, A-70, A-39)
def q0_cloud(seed, T=60, ends_n=3000):
    rng = np.random.default_rng(seed); X = rng.random((120000, 3))              # Poisson placement (A-89, A-02, A-43, A-31)
    tri = Delaunay(X); ind, nb = tri.vertex_neighbor_vertices
    vt = [[] for _ in range(len(X))]
    for t, s in enumerate(tri.simplices):
        for v in s: vt[v].append(t)
    def options(o, q):
        u = X[q] - X[o]; u /= np.linalg.norm(u)
        ring = {v for t in vt[o] if q in tri.simplices[t] for v in tri.simplices[t]} - {o, q}
        N = []
        for r in ring:
            w = X[r] - X[o]; w -= (w @ u) * u; w /= np.linalg.norm(w); nn = np.cross(u, w); N += [nn, -nn]
        return np.array(N)
    inner = np.where(np.all((X > 0.25) & (X < 0.75), axis=1))[0]
    ends = rng.choice(inner, ends_n, replace=False); res = []
    for o in ends:
        nbrs = nb[ind[o]:ind[o + 1]]; prev_q = None; prev_n = None
        for t in range(T):
            q = rng.choice([x for x in nbrs if x != prev_q]); opt = options(o, q)
            nv = opt[rng.integers(len(opt))] if prev_n is None else opt[np.argmin(opt @ prev_n)]   # A-70: most opposite side
            if prev_n is not None: res.append((1 + nv @ prev_n) / 2)                               # A-39: by squares; sign from loop (V-303)
            prev_n = nv; prev_q = q
    return float(np.mean(res))
seeds = list(range(101, 109)) if FULL else [101, 102]
qs = []
for s in seeds:
    qs.append(q0_cloud(s)); say('5 q0', f'cloud {s}: q0 = {qs[-1]:.5f}'); sys.stdout.flush()
q0 = float(np.mean(qs)); q0_se = float(np.std(qs, ddof=1) / sqrt(len(qs))) if len(qs) > 1 else float('nan')
say('5 q0', f'q0 = {q0:.5f} +- {q0_se:.5f} over {len(qs)} clouds  ->  1/alpha at the step = 1/q0^2 = {1/q0**2:.2f}')

# ---------------------------------------------------------------- 6. running from the step to the Z distance (W boson does not screen, V-294)
r_Z = hbar / (M_Z * MeV * c); D_t = 2 * hbar / (m_t * MeV * c)
run_W = kappa * (sum_active_at_rZ * log(r_Z / a) + gen['up'] * (log(D_t / a) + E))
inv_W_rZ = 1 / q0**2 + run_W
d_inv_W = 2 / q0**3 * q0_se if len(qs) > 1 else float('nan')
say('6 run', f'r_Z = {r_Z:.4e} m; running step -> r_Z without q0 = {run_W:.3f}; 1/alpha_W(r_Z) = {inv_W_rZ:.2f} +- {d_inv_W:.2f}')

# ---------------------------------------------------------------- 7. compare with data at the Z distance and far away
inv_data_rZ = inv_alpha_MZ_data + kappa * GAMMA * sum_active_at_rZ                 # energy -> real distance (log regime)
drop_data = inv_alpha_far_data - inv_data_rZ
say('7 data', f'1/alpha(r_Z) from D-48 in real distance = {inv_data_rZ:.2f} +- 0.05; measured drop r_Z -> far = {drop_data:.2f}')
diff = inv_W_rZ - inv_data_rZ
print()
print(f'RESULT: 1/alpha at the Z distance: W {inv_W_rZ:.2f} +- {d_inv_W:.2f}  vs  data {inv_data_rZ:.2f} +- 0.05  (diff {diff:+.2f})')
print(f'        1/alpha far away:          W {inv_W_rZ + drop_data:.2f} +- {d_inv_W:.2f}  vs  data {inv_alpha_far_data:.3f}  ({(inv_W_rZ + drop_data)/inv_alpha_far_data*100 - 100:+.2f} %)')
print('Derived from positions: a (c_g=1), kappa, edge E, Nc*Q^2, placement, q0.  From data: G, hbar, c, M_Z, m_t, 1/alpha(M_Z), 1/alpha far.')
print('Assumed, not checked: W boson does not screen (V-294, V-299); "pair remainder = q0" (V-295); energy -> distance conversion (gamma).')
