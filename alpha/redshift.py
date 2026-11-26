"""Redshift without expansion. Rules: redshift_rules_2026-11-26.md (pushed as c02425e before running)."""
from math import sqrt

# data
b_dat, b_sig = 1.003, sqrt(0.005**2 + 0.010**2)          # D-63 DES SN Ia time dilation
z_m, T_m, T_sig = 0.89, 5.08, 0.10                        # D-64 Muller 2013
beta_dat, beta_sig = 0.006, 0.013                         # D-64 Hurier 2014
T0 = 2.7255
tU, tU_sig = 13.5, sqrt(0.15**2 + 0.23**2)      # D-65 Valcin 2021, Gyr (stat avg, sys)
H_loc = {'SH0ES 2022': (73.04, 1.04), 'CCHP 2025': (70.39, sqrt(1.22**2 + 1.33**2 + 0.70**2))}
KMSMPC_GYR = 977.79                                       # 1/Gyr in km/s/Mpc

print('Fork 1: what the factor acts on')
for name, b, beta in (('M1 photon only', 0.0, 1.0), ('M2 one factor (A-33)', 1.0, 0.0)):
    Tz = T0 * (1 + z_m) ** (1 - beta)
    nb, nbeta, nT = abs(b - b_dat) / b_sig, abs(beta - beta_dat) / beta_sig, abs(Tz - T_m) / T_sig
    ok = nb <= 3 and nbeta <= 3
    print(f'  {name:22} b = {b:.0f} ({nb:5.1f} sigma)   beta = {beta:.0f} ({nbeta:5.1f} sigma)   '
          f'T(0.89) = {Tz:.2f} K vs {T_m} ({nT:4.1f} sigma)   {"PASS" if ok else "FAIL"}')

print()
print('Fork 2, Z1: H0 = 1/t_U')
H = KMSMPC_GYR / tU
H_sig = H * tU_sig / tU
print(f'  t_U = {tU} +- {tU_sig:.2f} Gyr  ->  H0 = {H:.2f} +- {H_sig:.2f} km/s/Mpc')
allok = True
for k, (h, s) in H_loc.items():
    n = abs(H - h) / sqrt(s**2 + H_sig**2)
    allok &= n <= 2
    print(f'  vs {k:11} {h:.2f} +- {s:.2f}: {n:.1f} sigma  {"pass" if n <= 2 else "fail"}')
print('  Z1:', 'PASS' if allok else 'FAIL')
print(f'  (for scale, not a criterion: Planck LCDM-derived 67.4 -> {abs(H - 67.4) / sqrt(0.5**2 + H_sig**2):.1f} sigma)')

print()
print('Consequence of Z1 (not a criterion): time at T = 0.1 MeV')
kB = 8.617333e-11                                         # MeV/K
onepz = 0.1 / kB / T0
t_s = tU * 1e9 * 3.15576e7 / onepz
print(f'  1+z = {onepz:.3g},  t = t_U/(1+z) = {t_s:.3g} s = {t_s / 3.15576e7:.0f} yr;  neutron lifetime 879 s, '
      f'ratio {t_s / 879:.2g}')
