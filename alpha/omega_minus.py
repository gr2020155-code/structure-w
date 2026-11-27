"""Structure of Omega- as a whole (R-146). Rules: omega_minus_rules_2026-11-27.md (committed before the run).
usage: python -P -E -X utf8 omega_minus.py"""
ME, IA, M_OMEGA = 0.51099895069, 137.035999177, 1672.43
mass = lambda x: (1 + 1.5 * IA * x) * ME
xdata = (M_OMEGA / ME - 1) / (1.5 * IA)
k = 3
B = [('Omega1 one plane (as proton)', k**2), ('Omega2 two planes at once', k**4), ('Omega3 two planes in turn', 2 * k**2),
     ('Omega4 three planes in turn', 3 * k**2), ('Omega5 3D as tau', 2 * k**4)]
print(f'data: m = {M_OMEGA} MeV, X = {xdata:.4f}')
for n, x in B:
    d = mass(x) / M_OMEGA - 1
    print(f'  {n:30s} X = {x:4d}: m = {mass(x):8.1f} MeV ({100 * d:+.1f} %) -> {"PASS" if abs(d) <= 0.01 else "FAIL"}')
print(f'noticed before the rules (not a branch): X = 16 -> m = {mass(16):.1f} MeV ({100 * (mass(16) / M_OMEGA - 1):+.2f} %)')
