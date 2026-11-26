"""Which data really test spin hbar/2? (26.11.2026)
Rules: spin_data_rules_2026-11-26.md (pushed before run, commit 258f7e3). Record: Р-103.
Run: python -I -X utf8 spin_data.py
"""
from math import pi

hbar = 1.054571817e-34; e = 1.602176634e-19; m = 9.1093837015e-31; c = 299792458.0
alpha = 1 / 137.035999
R = hbar / (m * c)
muB = e * hbar / (2 * m)
B = 5.0  # tesla, any value; cancels

print('1. Penning trap (D-01)')
mu_orbit = e * c * R / 2
for L_name, L in (('orbit L = hbar', hbar), ('L = hbar/2', hbar / 2)):
    w_s = 2 * mu_orbit * B / hbar          # A-04: flip frequency = energy difference / hbar
    w_c = e * B / m
    w_gyro = mu_orbit * B / L              # classical gyroscope precession (not what A-04 gives)
    print(f'   {L_name:15s}: mu/muB = {mu_orbit/muB:.6f}; w_s/w_c (A-04) = {w_s/w_c:.6f}; classical precession / w_c = {w_gyro/w_c:.3f}')
print(f'   data g/2 = 1.00115965; residual 1.15965e-3; alpha/(2 pi) = {alpha/(2*pi):.5e} (noticed before rules, not an argument)')
print('   -> w_s/w_c does not contain L: the trap measures mu_z, not the mechanical moment')

print('2. Einstein - de Haas (D-05)')
dmu = 2 * mu_orbit
for L_name, dL in (('orbit reverses, dL = 2 hbar', 2 * hbar), ('dL = hbar', hbar)):
    gp = (dmu / dL) / (e / (2 * m))
    print(f'   {L_name:28s}: g\' = {gp:.3f}')
for name, gp, s in (('Fe 1959', 1.919, 0.002), ('Co 1956', 1.850, 0.004)):
    print(f'   {name}: g\' = {gp} +- {s}; W orbit (1.000) off by {(gp-1)/s:.0f} sigma; per flip data need dL = {2/gp:.3f} hbar')

print('3. Can the medium carry hbar/2? (B-245 p.3: share alpha/2 of energy outside the orbit)')
E_out = alpha / 2 * m * c**2
r_need = (hbar / 2) * c / E_out
print(f'   L_field <= E_out * r / c; for hbar/2 need r >= {r_need/R:.1f} R  (R = {R:.3e} m)')
