"""Edge of screening for the electron form in the medium.
Rules: electron_edge_rules_2026-11-25.md (written before this run).
w(r) = P(S >= r); smooth fall-off equals a sharp edge at r_eff = D*exp(E[ln(S/D)]).
Data side: Uehling potential (standard vacuum polarization), checked numerically here."""
from math import pi, log, sqrt, exp, sin, sinh, cosh, acosh
import random, bisect

GAMMA = 0.5772156649015329
kappa = 2 / (3 * pi)
TARGET_POS = -(5 / 6 + GAMMA + log(2))   # real distance between charges
TARGET_MOM = -(5 / 6 + log(2))           # convention r = hbar/(Qc), as in V-276

# ---------- data side: Uehling ----------
def uehling(rho, n=40000):
    """I(rho) and w(rho) = -dI/dln(rho); rho = m c r / hbar.  t = cosh u."""
    U = acosh(max(2.0, 80.0 / (2 * rho)))
    h = U / n
    sI = sW = 0.0
    for k in range(n + 1):
        u = k * h
        t = cosh(u); sh = sinh(u)
        if t == 1.0:
            fI = fW = 0.0
        else:
            e = exp(-2 * rho * t) * (1 + 1 / (2 * t * t)) * sh * sh
            fI = e / (t * t)          # integrand of I times dt/du = sinh u
            fW = e / t
        wgt = 1 if k in (0, n) else (4 if k % 2 else 2)
        sI += wgt * fI; sW += wgt * fW
    return sI * h / 3, 2 * rho * sW * h / 3

print('Data side (Uehling), check of the target number')
for rho in (1e-3, 1e-4, 1e-5):
    I, w = uehling(rho)
    print(f'   rho = {rho:.0e}: I + ln(rho) = {I + log(rho):.5f}   (limit -(gamma+5/6) = {-(GAMMA + 5/6):.5f}),  w = {w:.5f}')
I, _ = uehling(1e-5)
print(f'   ln(r_eff/D) from Uehling = {I + log(1e-5) - log(2):.4f};  target in rules {TARGET_POS:.4f}')
print(f'   (convention r = hbar/(Qc), as in V-276: {TARGET_MOM:.4f})')

xs = [10 ** (k / 40) for k in range(-120, 21)]          # x = r/D from 1e-3 to 3.16
wU = [uehling(2 * x, 8000)[1] for x in xs]
# x where Uehling w = 1/2
xh_U = next(xs[i - 1] * (xs[i] / xs[i - 1]) ** ((wU[i - 1] - 0.5) / (wU[i - 1] - wU[i]))
            for i in range(1, len(xs)) if wU[i] < 0.5 <= wU[i - 1])
print(f'   Uehling: w = 1/2 at r/D = {xh_U:.4f};  w at r = D: {wU[xs.index(1.0)]:.4f}')

# ---------- readings ----------
random.seed(20261125)
N = 600000
def iso():
    while True:
        x, y, z = (random.uniform(-1, 1) for _ in range(3))
        q = x * x + y * y + z * z
        if 1e-12 < q <= 1:
            q = sqrt(q); return x / q, y / q, z / q
def in_ball():
    while True:
        x, y, z = (random.uniform(-1, 1) for _ in range(3))
        if x * x + y * y + z * z <= 1:
            return x / 2, y / 2, z / 2      # ball of diameter 1
def chord_circle():
    a, b = random.uniform(0, 2 * pi), random.uniform(0, 2 * pi)
    return abs(sin((a - b) / 2))
def cosabs():
    return abs(random.uniform(-1, 1))     # |cos| of isotropic direction to the line from the charge

readings = {
    'Ch1 chord of circle':                 lambda: chord_circle(),
    'Ch2 segment A-O-B, projected':        lambda: cosabs(),
    'Ch3 whole circle, extent along line': lambda: sqrt(1 - random.random() ** 2),
    'Ch4 chord of circle, projected':      lambda: chord_circle() * cosabs(),
    'Ch5 chord of sphere':                 lambda: sqrt(random.random()),
    'Ch6 chord of sphere, projected':      lambda: sqrt(random.random()) * cosabs(),
    'Ch7 two points in ball':              lambda: (lambda p, q: sqrt(sum((p[i] - q[i]) ** 2 for i in range(3))))(in_ball(), in_ball()),
    'Ch8 two points in ball, projected':   lambda: (lambda p, q: abs(p[2] - q[2]))(in_ball(), in_ball()),
}
exact = {'Ch1 chord of circle': -log(2), 'Ch2 segment A-O-B, projected': -1.0,
         'Ch3 whole circle, extent along line': log(2) - 1, 'Ch4 chord of circle, projected': -log(2) - 1,
         'Ch5 chord of sphere': -0.5, 'Ch6 chord of sphere, projected': -1.5}

r_Z_over = {'e': 0.51099895, 'mu': 105.6583755, 'tau': 1776.86}   # MeV; r_Z = hbar/(M_Z c)
MZ = 91187.6
# drop of 1/alpha from far away to r_Z for leptons with a sharp edge at D*exp(E): kappa*sum ln(D e^E / r_Z), D/r_Z = 2 M_Z/m
def lepton_drop(E):
    return sum(kappa * (log(2 * MZ / m) + E) for m in r_Z_over.values())
drop_pos = lepton_drop(TARGET_POS)

print(f'\nReadings: E[ln(S/D)], edge r_eff/D, shape (x1/2, max|w - w_U|), lepton drop of 1/alpha to r_Z')
print(f'   target (real distance): E = {TARGET_POS:.4f}, r_eff/D = {exp(TARGET_POS):.4f}, x1/2 = {xh_U:.4f}, drop {drop_pos:.2f}')
hits = []; Es = {}
for name, f in readings.items():
    S = sorted(max(f(), 1e-300) for _ in range(N))
    E = sum(log(s) for s in S) / N
    med = S[N // 2]
    dev = max(abs((1 - bisect.bisect_left(S, x) / N) - wu) for x, wu in zip(xs, wU))
    ex = f' (exact {exact[name]:.4f})' if name in exact else ''
    Es[name] = E
    hit = abs(E - TARGET_POS) < 0.05
    if hit: hits.append(name)
    print(f'   {name:37s} E = {E:7.4f}{ex:17s} r_eff/D = {exp(E):.4f}  x1/2 = {med:.4f}  max|dw| = {dev:.3f}  '
          f'drop {lepton_drop(E):.2f} ({lepton_drop(E) - drop_pos:+.2f})  {"HIT" if hit else "miss"}')

print(f'\nResult: hits within 0.05 of {TARGET_POS:.4f}: {", ".join(hits) if hits else "none"}')
for lab, T in (('real distance', TARGET_POS), ('V-276 convention r = hbar/(Qc)', TARGET_MOM)):
    best = min(Es, key=lambda k: abs(Es[k] - T))
    print(f'Closest to {lab} ({T:.4f}): {best}, E = {Es[best]:.4f}, difference {Es[best] - T:+.4f}')
