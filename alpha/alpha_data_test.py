"""Check of the alpha chain against data only (26.11.2026). Rules: alpha_data_test_rules_2026-11-26.md (pushed before this run).
Run: python -X utf8 alpha_data_test.py     (needs numpy, scipy)
"""
from math import pi, sqrt, log
from scipy.integrate import quad
from scipy.special import exp1

GAMMA = 0.5772156649015329
hbar = 1.054571817e-34; c = 299792458.0; G = 6.67430e-11
MeV = 1.602176634e-13 / c**2
a = 2 * pi * sqrt(hbar * G / c**3)
kappa = 2 / (3 * pi)
inv_q0sq, s_W = 1 / 0.12282**2, 2 * 0.00013 / 0.12282**3          # 66.29 +- 0.14

# data (PDG, MeV)
m_e, m_mu, m_tau = 0.51099895, 105.6583755, 1776.86
m_p, M_L, M_Lc, M_Lb, m_t = 938.27209, 1115.683, 2286.46, 5619.60, 172570.0
target_far = {'Rb': 137.035999206, 'Cs': 137.035999046}
opal, opal_err = 440e-5, sqrt(58**2 + 43**2 + 30**2) * 1e-5
t1, t2 = 6.07e6, 1.81e6                                             # |t| in MeV^2

species = {  # name: (mass, Nc*Q^2)
    'e': (m_e, 1), 'mu': (m_mu, 1), 'tau': (m_tau, 1),
    'u': (m_p / 3, 4 / 3), 'd': (m_p / 3, 1 / 3), 's': (M_L - 2 * m_p / 3, 1 / 3),
    'c': (M_Lc - 2 * m_p / 3, 4 / 3), 'b': (M_Lb - 2 * m_p / 3, 1 / 3), 't': (m_t, 4 / 3)}
print('fragment masses (Rb), MeV:', ', '.join(f'{k} {v[0]:.1f}' for k, v in species.items()))

# edge readings: weight over the split x (normalised), symmetric S1 is a point x = 1/2
edges = {
    'S1': None,
    'S2': lambda x: 1.0,
    'S3': lambda x: 6 * x * (1 - x)}


def avg(edge, f):
    if edge == 'S1':
        return f(0.5)
    w = edges[edge]
    return quad(lambda x: w(x) * f(x), 0, 1, limit=400, points=[0.5])[0]


def E_of(edge):  # ln(r_eff/D), D = 2*hbar/(m c), lambda_x = (hbar/(m c))*sqrt(x(1-x))
    return avg(edge, lambda x: log(sqrt(x * (1 - x)) / 2)) - GAMMA


def inv_far(edge):
    E = E_of(edge)
    tot = 0.0
    for m, w in species.values():
        D = 2 * hbar / (m * MeV * c)
        tot += kappa * w * (log(D / a) + E)
    return inv_q0sq + tot, E


def delta_inv(edge, reading, Q2):  # 1/alpha(Q) - 1/alpha(far), Q2 in MeV^2
    Q = sqrt(Q2); tot = 0.0
    for m, w in species.values():
        if reading == 'M1':     # point: r = hbar/(Q c); r/lambda_x = m/(Q sqrt(x(1-x)))
            f = lambda x: exp1(m / (Q * sqrt(x * (1 - x))))
        else:                   # amplitude: (1/2) ln(1 + Q^2 lambda_x^2)
            f = lambda x: 0.5 * log(1 + Q2 * x * (1 - x) / m**2)
        tot += kappa * w * avg(edge, f)
    return -tot


print('\n== C1: 1/alpha far vs atom recoil (137.036), pass if |dev| <= 2*sigma_W =', f'{2*s_W:.2f}')
alpha0 = 1 / target_far['Rb']
res = {}
for ed in edges:
    inv, E = inv_far(ed)
    dev = inv - target_far['Rb']
    ok1 = abs(dev) <= 2 * s_W
    res[ed] = ok1
    print(f'{ed}: E = {E:.4f}; 1/alpha_W(0) = {inv:.3f} +- {s_W:.2f}; dev {dev:+.3f} ({dev/target_far["Rb"]*100:+.2f} %, {dev/s_W:+.1f} sigma) -> {"PASS" if ok1 else "FAIL"}')

print('\n== C2: OPAL  Dalpha(-6.07) - Dalpha(-1.81 GeV^2) = 440 +- 78 (1e-5), pass if |dev| <= 156')
for ed in edges:
    for rd in ('M1', 'M2'):
        d1, d2 = delta_inv(ed, rd, t1), delta_inv(ed, rd, t2)
        diff = alpha0 * (d2 - d1)            # alpha0*(1/alpha(t2) - 1/alpha(t1))
        dev = diff - opal
        ok2 = abs(dev) <= 2 * opal_err
        both = res[ed] and ok2
        print(f'{ed}+{rd}: {diff*1e5:.0f} (1e-5), dev {dev*1e5:+.0f} ({dev/opal_err:+.1f} sigma) -> {"PASS" if ok2 else "FAIL"}; both targets: {"yes" if both else "no"}')

# split of the far value, S3, for the record
print('\n== breakdown S3: contribution of each species, kappa*NcQ^2*(ln(D/a)+E)')
E3 = E_of('S3')
for k, (m, w) in species.items():
    D = 2 * hbar / (m * MeV * c)
    print(f'  {k}: {kappa*w*(log(D/a)+E3):.3f}')
need = target_far['Rb'] - inv_far('S3')[0]
print(f'missing to 137.036 at S3: {need:+.3f}; that equals kappa*sum(NcQ2) * {need/(kappa*8):+.3f} of edge, or 1/q0^2 = {inv_q0sq+need:.2f} (q0 = {1/sqrt(inv_q0sq+need):.5f})')
