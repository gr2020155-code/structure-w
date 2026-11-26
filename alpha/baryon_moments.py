"""Baryon magnetic moments. Rules: baryon_moments_rules_2026-11-26.md (pushed as 62e41e4 before running).

Part i contributes w_i * Q_i * (m_p / M_i) mu_N  (V-14, R-78); weights w_i from V-312.
"""
from fractions import Fraction as F

m_p = 938.272
Q = {'u': F(2, 3), 'd': F(-1, 3), 's': F(-1, 3), 'c': F(2, 3)}
M = {'p': 938.272, 'L': 1115.683, 'S+': 1189.37, 'S-': 1197.449, 'X0': 1314.86, 'X-': 1321.71,
     'O-': 1672.45, 'D++': 1232.0, 'Lc+': 2286.46, 'Sc++': 2453.97, 'Sc0': 2453.75, 'Xc+': 2467.71,
     'Xc0': 2470.44, 'Oc0': 2695.2, 'Xcc++': 3621.6, 'Xcc+': 3621.6}
own = {'u': m_p / 3, 'd': m_p / 3}
own['s'] = M['L'] - 2 * m_p / 3
own['c'] = M['Lc+'] - 2 * m_p / 3

# name: (parts, spin 3/2?)
B = {'L': ('uds', 0), 'S+': ('uus', 0), 'S-': ('dds', 0), 'X0': ('ssu', 0), 'X-': ('ssd', 0), 'O-': ('sss', 1),
     'D++': ('uuu', 1),
     'Lc+': ('udc', 0), 'Sc++': ('uuc', 0), 'Sc0': ('ddc', 0), 'Xc+': ('usc', 0), 'Xc0': ('dsc', 0),
     'Oc0': ('ssc', 0), 'Xcc++': ('ccu', 0), 'Xcc+': ('ccd', 0), 'Occ+': ('ccs', 0), 'Occc++': ('ccc', 1)}
data = {'L': (-0.613, 0.004), 'S+': (2.458, 0.010), 'S-': (-1.160, 0.025), 'X0': (-1.250, 0.014),
        'X-': (-0.6507, 0.0025), 'O-': (-2.02, 0.05)}
HYP = ['L', 'S+', 'S-', 'X0', 'X-', 'O-']
CHARM = ['Lc+', 'Sc++', 'Sc0', 'Xc+', 'Xc0', 'Oc0', 'Xcc++', 'Xcc+', 'Occ+', 'Occc++']
P, PAIR, ANTI, ZERO = F(2, 3), F(2, 3), F(1, 3), F(0)

def weights(parts, s32, L):
    """list of weights in the order of `parts`."""
    if s32:
        return [F(1)] * 3
    for i in range(3):                       # identical pair (V-312, A-31)
        rest = [j for j in range(3) if j != i]
        if parts[rest[0]] == parts[rest[1]]:
            w = [None] * 3; w[rest[0]] = w[rest[1]] = PAIR; w[i] = F(-1, 3)
            return w
    order = sorted(range(3), key=lambda j: own[parts[j]])   # all three different
    light, heavy = order[:2], order[2]
    w = [None] * 3
    if L == 'L1':
        w[light[0]] = w[light[1]] = ZERO; w[heavy] = F(1)
    elif L == 'L2':
        w[light[0]] = w[light[1]] = PAIR; w[heavy] = F(-1, 3)
    elif L == 'L3':
        for i in range(3):
            rest = [j for j in range(3) if j != i]
            if (Q[parts[rest[0]]] > 0) == (Q[parts[rest[1]]] > 0):
                w[rest[0]] = w[rest[1]] = PAIR; w[i] = F(-1, 3)
    elif L == 'L4':
        w[light[0]] = w[light[1]] = ANTI; w[heavy] = ANTI
    return w

def part_mass(parts, name, R):
    if R == 'Rb':
        return [own[x] for x in parts]
    if name not in M:
        return None
    if R == 'Ra':
        return [M[name] / 3] * 3
    tot = sum(own[x] for x in parts)          # Rc
    return [M[name] * own[x] / tot for x in parts]

def moment(name, R, L):
    parts, s32 = B[name]
    mm = part_mass(parts, name, R)
    if mm is None:
        return None
    w = weights(parts, s32, L)
    return sum(float(wi * Q[x]) * m_p / mi for wi, x, mi in zip(w, parts, mm))

def score(f):
    errs = [abs(f(n) / data[n][0] - 1) for n in HYP]
    return sum(errs) / len(errs), max(errs), errs

print(f'own part masses (Rb): u = d = {own["u"]:.1f}, s = {own["s"]:.1f}, c = {own["c"]:.1f} MeV')
print('check: proton, neutron-like nucleon under Ra/Rb/Rc =',
      [round(float(2 * PAIR * Q['u'] - F(1, 3) * Q['d']) * 3, 3)], '(V-308: 3.000)')
print()
hdr = ''.join(f'{n:>8}' for n in HYP)
print(f'{"reading":10}{hdr}  {"mean":>6} {"max":>6}  pass')
print(f'{"data":10}' + ''.join(f'{data[n][0]:8.3f}' for n in HYP))
passed = []
for R in ('Ra', 'Rb', 'Rc'):
    for L in ('L1', 'L2', 'L3', 'L4'):
        f = lambda n: moment(n, R, L)
        mean, mx, errs = score(f)
        ok = mean <= 0.10 and mx <= 0.20
        if ok:
            passed.append((R, L))
        print(f'{R + "*" + L:10}' + ''.join(f'{f(n):8.3f}' for n in HYP) + f'  {mean:6.1%} {mx:6.1%}  {"YES" if ok else "no"}')
        print(f'{"  error":10}' + ''.join(f'{e:8.0%}' for e in errs))
# reference: PDG simple quark model with 3 fitted moments
mq = {'u': 1.852, 'd': -0.972, 's': -0.613}
ref = {'L': mq['s'], 'S+': (4 * mq['u'] - mq['s']) / 3, 'S-': (4 * mq['d'] - mq['s']) / 3,
       'X0': (4 * mq['s'] - mq['u']) / 3, 'X-': (4 * mq['s'] - mq['d']) / 3, 'O-': 3 * mq['s']}
mean, mx, errs = score(lambda n: ref[n])
print(f'{"ref PDG":10}' + ''.join(f'{ref[n]:8.3f}' for n in HYP) + f'  {mean:6.1%} {mx:6.1%}  {"YES" if mean <= .1 and mx <= .2 else "no"}  (3 moments fitted, for scale)')
print(f'{"  error":10}' + ''.join(f'{e:8.0%}' for e in errs))
print()
print('Delta++ (data ~6.1 +- 0.5, model-dependent):', {R: round(moment('D++', R, 'L1'), 2) for R in ('Ra', 'Rb', 'Rc')})
print()
print('passed:', passed if passed else 'none')
print()
print('charmed baryons (consequences of each reading; predictions of W only for passed readings)')
print(f'{"reading":10}' + ''.join(f'{n:>8}' for n in CHARM))
for R in ('Ra', 'Rb', 'Rc'):
    for L in ('L1', 'L2', 'L3', 'L4'):
        row = []
        for n in CHARM:
            v = moment(n, R, L)
            row.append(f'{v:8.3f}' if v is not None else f'{"-":>8}')
        print(f'{R + "*" + L:10}' + ''.join(row))
