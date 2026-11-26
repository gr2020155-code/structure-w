"""Screening boundary for the return of the self-coupling (R-118). Rules: screen_boundary_rules_2026-11-26.md.
w(rho) = int_0^1 6x(1-x) exp(-rho/sqrt(x(1-x))) dx, rho = r/R (soft profile V-286, f = 1); sharp profile: w = 1 for rho < 2."""
from math import pi, log, exp, sqrt
alpha = 1 / 137.035999177; kappa = 2 / (3 * pi)
delta_data = 0.001513
hbar, c, G = 1.054571817e-34, 299792458.0, 6.67430e-11
a = 2 * pi * sqrt(hbar * G / c**3)
R = hbar / (0.51099895 * 1.602176634e-13 / c)            # hbar/(m_e c), m
eps = a / R

NX = 400
xs = [(k + 0.5) / NX for k in range(NX)]
lam = [sqrt(x * (1 - x)) for x in xs]; wt = [6 * x * (1 - x) / NX for x in xs]
def w_soft(rho): return sum(q * exp(-rho / l) for q, l in zip(wt, lam))
def w_sharp(rho): return 1.0 if rho < 2 else 0.0

def integ(f, lo, hi, n=4000):                            # in ln(rho), midpoint
    if hi <= lo: return 0.0
    L = log(hi / lo); h = L / n
    return sum(f(lo * exp((i + 0.5) * h)) for i in range(n)) * h   # f already includes the 1/rho*rho = 1 of d ln rho

rows = []
for pname, w in (('P1 sharp 2R (excluded, V-283 p.5)', w_sharp), ('P2 soft V-286', w_soft)):
    # part of the path inside the orbit for H2: from the step eps to 1
    inner = integ(lambda r: w(r), eps, 1.0)
    # beyond the orbit, averaged over the reach: M1 weight P(r>rho) = 1/rho, M2 weight 2/rho - 1 on [1,2]
    outer = {'M1': integ(lambda r: w(r) / r, 1.0, 40.0), 'M2': integ(lambda r: w(r) * (2 / r - 1), 1.0, 2.0)}
    print(f'{pname}: int w dln(rho) from step to orbit = {inner:.4f}; beyond orbit <...> M1 = {outer["M1"]:.5f}, M2 = {outer["M2"]:.5f}')
    for start in ('H1 from orbit', 'H2 from step'):
        for legs in (1, 2):
            for m in ('M1', 'M2'):
                L = outer[m] + (inner if start.startswith('H2') else 0.0)
                rows.append((f'{pname[:2]} {start}, legs {legs}, {m}', legs * kappa * alpha * L, pname.startswith('P2')))
print(f'check: soft profile, int w dln(rho) from step to far - ln(D/a) = {integ(w_soft, eps, 40.0) - log(2 / eps):+.4f} (edge E = -2.1037)')
print(f'step/R = {eps:.3e}; ln(D/a) = {log(2/eps):.2f}\n')
for name, d, allowed in rows:
    dev = d / delta_data - 1
    ok = abs(dev) <= 0.05
    verdict = ('PASS' if allowed else 'pass, but profile excluded') if ok else 'FAIL'
    print(f'{name:38s}: delta = {d:.6f}  ({dev*100:+9.1f} %)  -> {verdict}')
