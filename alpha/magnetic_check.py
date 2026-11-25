"""Magnetic moments of composite forms: rules M0, M1 fixed in magnetic_rules_2026-11-18.md BEFORE comparison.

Data (CODATA 2018): mu_p = 2.79284734 mu_N, mu_n = -1.91304273 mu_N, mu_d = 0.85743823 mu_N,
a_e = 0.00115965218, a_mu = 0.00116592061; reduced Compton wavelength of the proton 0.21031 fm.
Also printed: the data translated into W terms - radius of a unit charge circulating at c that gives each moment
(mu = e c r / 2, so r / R_p = mu / mu_N). This is a translation of data, not a derivation.
"""
import math
mu_p, mu_n, mu_d = 2.79284734, -1.91304273, 0.85743823
a_e, a_mu, alpha = 0.00115965218, 0.00116592061, 1 / 137.035999084
Rp_fm = 0.21030891        # hbar/(m_p c), fm
kp, kn = mu_p - 1, mu_n     # anomalous parts (Dirac part: charge 1 and 0)
print('M0 (one circulation of the total charge):   mu_p = 1.000, mu_n = 0.000, ratio 0.000;'
      f'  data {mu_p:.4f}, {mu_n:.4f}, ratio {mu_n/mu_p:.5f}')
print(f'M1 leptons: anomaly of order alpha only.  data a_e = {a_e:.8f}, a_mu = {a_mu:.8f}, alpha/2pi = {alpha/2/math.pi:.8f}')
print(f'M1 proton: anomaly of order 1.             data kappa_p = {kp:.4f}')
print(f'M1 neutron: mu_n = kappa_p = {kp:+.4f} (sign +).  data {mu_n:+.4f} (sign -)  ->  sign contradiction')
print(f'what the data need from the eight inside the neutron under M1: mu_8 = mu_n - kappa_p = {mu_n - kp:+.4f} mu_N')
print(f'isoscalar kappa_p + kappa_n = {kp + kn:+.4f}, isovector kappa_p - kappa_n = {kp - kn:+.4f}')
print('radius of a unit charge at speed c giving the moment (fm):')
for name, m in (('mu_N (Dirac circulation of the proton)', 1.0), ('mu_p', mu_p), ('|mu_n|', -mu_n), ('kappa_p', kp),
                ('kappa_p - kappa_n (eight inside the neutron, M1)', kp - kn)):
    print(f'   {name:52s} {m * Rp_fm:.3f}')
print('measured sizes (fm): proton charge radius 0.841, magnetic radii p 0.851, n 0.864; W core (4 points) < 3e-5 fm')
print(f'deuteron: mu_p + mu_n = {mu_p + mu_n:.4f} vs data {mu_d:.4f} ({(mu_d/(mu_p + mu_n) - 1)*100:+.1f} %)')
