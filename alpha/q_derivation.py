"""Budget of screening in the medium (rules: q_derivation_rules_2026-11-25.md).
Each charged form s screens inside its diameter D_s = 2*hbar/(m_s c), slope k_s = kappa*Nc*Q^2 per e-fold."""
from math import pi, log, sqrt
hbar = 1.054571817e-34; c = 299792458; G = 6.67430e-11; MeV = 1.602176634e-13 / c**2  # kg per MeV
lP = sqrt(hbar * G / c**3); a = 2 * pi * lP
kappa = 2 / (3 * pi)
# name, mass MeV, Nc, Q^2  (PDG; light-quark masses MSbar at 2 GeV, top pole)
forms = [('e', 0.51099895, 1, 1), ('mu', 105.6583755, 1, 1), ('tau', 1776.86, 1, 1),
         ('u', 2.16, 3, 4/9), ('d', 4.70, 3, 1/9), ('s', 93.5, 3, 1/9),
         ('c', 1273.0, 3, 4/9), ('b', 4183.0, 3, 1/9), ('t', 172570.0, 3, 4/9)]
D = lambda m: 2 * hbar / (m * MeV * c)
inv_far = 137.035999177

print(f'step a = {a:.4e} m, kappa = {kappa:.5f}')
print('\nQ1: 1/alpha at r = hbar/(M_Z c) from 1/alpha = 137.036 far away, edge at diameter')
rZ = hbar / (91187.6 * MeV * c)
dz = sum(kappa * Nc * Q2 * log(D(m) / rZ) for _, m, Nc, Q2 in forms if D(m) > rZ)
print(f'   r_Z = {rZ:.3e} m; drop {dz:.2f}; 1/alpha = {inv_far - dz:.2f}  (data: 127.95 MSbar, 128.9 effective)')
# same with QED-like edge at hbar/(m c) (reduced Compton) and at hbar/(2mc) for comparison
for lab, f in (('radius hbar/mc', 0.5), ('hbar/2mc', 0.25)):
    dz2 = sum(kappa * Nc * Q2 * log(f * D(m) / rZ) for _, m, Nc, Q2 in forms if f * D(m) > rZ)
    print(f'   (for comparison, edge at {lab}: 1/alpha = {inv_far - dz2:.2f})')

print('\nQ2: budget from the step a to each diameter')
tot = 0
for name, m, Nc, Q2 in forms:
    k = kappa * Nc * Q2; L = log(D(m) / a); tot += k * L
    print(f'   {name:4s} D = {D(m):.3e} m  ln(D/a) = {L:6.2f}  k = {k:.4f}  contribution {k*L:6.2f}')
print(f'   sum of fermions = {tot:.2f};  needed 136.04;  remainder R = {136.036 - tot:.2f}')
# W boson in standard counting: d(1/alpha)/dlnmu contribution -7/(2pi)... (one-loop: b_W = -7 in units 1/(2pi)? report as standard)
kW = -7 / (2 * pi)  # W loop in the standard one-loop count: b_W = -7 vs (4/3)NcQ^2 for fermions -> antiscreening
LW = log(D(80369.2) / a)
print(f'   W boson (standard sign, antiscreening): contribution {kW*LW:.2f};  with W: sum {tot + kW*LW:.2f}, R = {136.036 - tot - kW*LW:.2f}')

print('\nQ1 for leptons only (precise data, no hadron uncertainty)')
al = 1 / inv_far
dlep_data = sum(al / (3 * pi) * (2 * log(91187.6 / m) - 5 / 3) for _, m, _, _ in forms[:3])  # standard one-loop QED, matches measured lepton part
print(f'   standard leptonic Delta alpha(M_Z) = {dlep_data:.5f} (published 0.031498) -> drop in 1/alpha {dlep_data*inv_far:.2f}')
for lab, f in (('diameter 2hbar/mc', 1.0), ('radius hbar/mc', 0.5), ('effective edge from data hbar/mc*e^(-5/6)', 0.5 * 2.718281828 ** (-5 / 6))):
    d = sum(kappa * log(f * D(m) / rZ) for _, m, _, _ in forms[:3])
    print(f'   edge at {lab}: drop {d:.2f}')
