"""(normalisation at x = 0.003..0.01; below 0.003 g assumed flat when the log test passes)
Soft edge of screening: variants K1-K5.  Rules: soft_edge_rules_2026-11-25.md (written before this run).
f(r) = share of points per cell with the form's side; g = f*r normalised to 1 at small r.
ln(r_eff/D) = integral (g - 1{x<1}) dx/x, x = r/D.  Data (R-61): -2.1037, half slope at x = 0.150."""
import random, math
random.seed(20261125)
TARGET = -2.1037; HALF_DATA = 0.150
N = 3_000_000
lo, hi, nb = -2.5, 0.6, 31                       # log10 bins of x
edges = [10 ** (lo + (hi - lo) * k / nb) for k in range(nb + 1)]

def iso():
    while True:
        v = [random.gauss(0, 1) for _ in range(3)]
        n = math.sqrt(sum(t * t for t in v))
        if n > 0: return [t / n for t in v]
def rot_from(u):
    """random rotation matrix (columns orthonormal) from three gaussian vectors"""
    a = iso(); b = iso()
    b = [b[i] - sum(a[j] * b[j] for j in range(3)) * a[i] for i in range(3)]
    nb_ = math.sqrt(sum(t * t for t in b)); b = [t / nb_ for t in b]
    c = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    return a, b, c

def g_profile(sample):
    cnt = [0] * nb
    for _ in range(N):
        d = sample()
        if d <= 0: continue
        k = int((math.log10(d) - lo) / (hi - lo) * nb)
        if 0 <= k < nb: cnt[k] += 1
    xs, gs = [], []
    for k in range(nb):
        x1, x2 = edges[k], edges[k + 1]
        p = cnt[k] / N / (x2 - x1)               # density of distance
        xs.append(math.sqrt(x1 * x2)); gs.append(p / math.sqrt(x1 * x2))   # g ~ p(x)/x
    return xs, gs

def summarize(name, xs, gs):
    small = [g for x, g in zip(xs, gs) if 3e-3 < x < 1e-2]
    if not small or sum(small) == 0:
        print(f'{name}: no form points near the charge'); return None
    g0 = sum(small) / len(small)
    g_at = lambda X: next(g for x, g in zip(xs, gs) if x >= X) / g0
    flat = g_at(3e-2)
    dl = (hi - lo) / nb * math.log(10)
    E = sum(((g / g0) - (1 if x < 1 else 0)) * dl for x, g in zip(xs, gs) if x > 3e-3)
    half = next((x for i, (x, g) in enumerate(zip(xs, gs)) if x > 3e-3 and g / g0 < 0.5 and i + 1 < len(gs) and gs[i + 1] / g0 < 0.5), None)
    beyond = sum(g / g0 * dl for x, g in zip(xs, gs) if x > 1) / max(1e-12, sum(g / g0 * dl for x, g in zip(xs, gs) if x > 3e-3))
    logok = abs(flat - 1) < 0.15
    hit = logok and abs(E - TARGET) < 0.05 and half is not None and abs(half - HALF_DATA) < 0.03
    print(f'{name}: log inside D: {"yes" if logok else "NO"} (g(0.03)/g(0.003..0.01) = {flat:.2f});  '
          f'ln(r_eff/D) = {E:+.3f} (target {TARGET});  half at x = {half if half is None else round(half, 3)} (data 0.150);  '
          f'{"HIT" if hit else "miss"}')
    return E

R = 0.5   # D = 1
# K1: uniform sphere of diameter 1, charge on it
def k1():
    u = iso(); return math.dist((0, 0, R), [R * t for t in u])
xs, gs = g_profile(k1); E1 = summarize('K1 sharp sphere through charge', xs, gs)
print(f'K2 competition for the held "-": g = g_K1 * 1/(1 + 2x^2): ln(r_eff/D) = {-0.5 * math.log(3):+.3f} (exact, -ln3/2);  '
      f'half at x = {math.sqrt(0.5):.3f} (g = 1/2 at 2x^2 = 1)  miss')
# K3: charge uniform inside the ball of the form
def k3():
    while True:
        c = [random.uniform(-R, R) for _ in range(3)]
        if sum(t * t for t in c) <= R * R: break
    u = iso(); return math.dist(c, [R * t for t in u])
xs, gs = g_profile(k3); summarize('K3 charge anywhere inside the ball', xs, gs)
# K4: circle rotating about one diameter axis; charge = a random point of the same swept set
def k4pt():
    phi = random.uniform(0, 2 * math.pi); psi = random.uniform(0, 2 * math.pi)
    return (R * math.sin(phi) * math.cos(psi), R * math.sin(phi) * math.sin(psi), R * math.cos(phi))
def k4():
    return math.dist(k4pt(), k4pt())
xs, gs = g_profile(k4); summarize('K4 rotation about one axis', xs, gs)
# K5: figure eight (two circles of diameter 1 through O), all orientations about O; charge at distance 1 from O
def k5():
    a, b, c = rot_from(None)
    phi = random.uniform(0, 2 * math.pi); side = 1 if random.random() < 0.5 else -1
    px, py = side * (R - R * math.cos(phi)), R * math.sin(phi)          # circles centred at (+-R, 0), through O
    P = [px * a[i] + py * b[i] for i in range(3)]
    return math.dist((0, 0, 2 * R), P)
xs, gs = g_profile(k5); summarize('K5 figure eight swept about O, charge on its outer sphere', xs, gs)
