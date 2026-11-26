"""Screening part of the 2S-2P Lamb shift in muonic (and ordinary) hydrogen for decay-length readings L1-L4 (26.11.2026).
Rules: muonic_h_rules_2026-11-26.md (pushed before run). Bridge: non-relativistic hydrogen 2S, 2P functions (QM) - model, not yardstick.
W profile S3 (P-64) = sum of exponentials; in closed form per species: d alpha/alpha = alpha*(2/3pi)*Nc Q^2 * int_1^inf e^{-2rt/(f lam)} g(t) dt.
Run: python -X utf8 muonic_h.py   (needs scipy)
"""
from math import pi, sqrt
from scipy.integrate import quad

hbarc = 197.3269804                      # MeV fm
alpha = 1 / 137.035999206
m_e, m_mu, m_p = 0.51099895, 105.6583755, 938.27209
species = [(m_e, 1), (m_mu, 1), (1776.86, 1), (m_p / 3, 4 / 3), (m_p / 3, 1 / 3), (1115.683 - 2 * m_p / 3, 1 / 3),
           (2286.46 - 2 * m_p / 3, 4 / 3), (5619.60 - 2 * m_p / 3, 1 / 3)]          # t negligible
g = lambda t: (1 + 1 / (2 * t * t)) * sqrt(t * t - 1) / (t * t)


def shift(m_orb, f):                      # E(2P) - E(2S) from screening, MeV
    m_r = m_orb * m_p / (m_orb + m_p); a = hbarc / (alpha * m_r)
    def kern(t, lam):
        c = 1 / a + 2 * t / (f * lam)
        I2S = (1 / (2 * a**3)) * (1 / c**2 - 2 / (a * c**3) + 1.5 / (a * a * c**4))
        I2P = 1 / (4 * a**5 * c**4)
        return g(t) * (I2P - I2S)
    tot = 0.0
    for m, w in species:
        lam = hbarc / m
        tot += w * quad(lambda t: kern(t, lam), 1, float('inf'), limit=400)[0]
    return -alpha * hbarc * alpha * (2 / (3 * pi)) * tot


data, lo, hi = 202.3706, 202.37 + 3.35 - 3, 202.37 + 4.23 + 3
print(f'muonic H: measured 2P-2S {data} meV; pass band for screening part {lo:.1f} - {hi:.1f} meV (finite size 3.35-4.23, rest +-3)')
for name, f in (('L1 f=1', 1.0), ('L2 f=sqrt3', sqrt(3)), ('L3 f=sqrt(3/2)', sqrt(1.5)), ('L4 f=sqrt2', sqrt(2))):
    mu = shift(m_mu, f) * 1e9; h = shift(m_e, f) * 1e9 * 241798.93      # meV ; MHz
    print(f'{name}: muonic H screening {mu:.2f} meV -> {"PASS" if lo <= mu <= hi else "FAIL"} (needed {data+3.35:.1f}-{data+4.23:.1f}, dev {mu-(data+3.79):+.1f});'
          f'  ordinary H screening part of 2P-2S: {h:+.1f} MHz (info only)')
