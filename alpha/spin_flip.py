"""What does not take part in the orbit at a spin flip? All branches (26.11.2026).
Rules: spin_flip_rules_2026-11-26.md (pushed before run). Record: Р-104.
Run: python -I -X utf8 spin_flip.py
"""
from fractions import Fraction as F

gp_min = 1.919 - 2 * 0.002          # Einstein - de Haas, Fe (D-05); orbital part only lowers g'
branches = [
    ('a1 all inertia in the core',                     F(1),    True,  'from positions'),
    ('a2 core + one global eight, equal (A-60, B-212)', F(1, 2), True,  'from positions; "equal" by A-31 is weak'),
    ('a3 core + three global eights (A-73)',            F(1, 4), True,  'from positions'),
    ('a4 all inertia in the global eight',              F(0),    True,  'from positions'),
    ('b1 two eights = two systems, one flips (B-42)',   F(1, 2), False, 'conflicts with B-198: the form moves as a whole'),
    ('b2 orientation moment hbar*j, j = 1/2 (B-27)',    F(1, 2), True,  'premise outside the positions'),
]
passed = []
for name, f, ok_pos, why in branches:
    dL = 2 * f
    c1 = dL == 1                                     # one photon (A-04, D-07, A-28)
    gp = float('inf') if f == 0 else 1 / float(f)
    c2 = gp >= gp_min
    ok = c1 and c2 and ok_pos
    if ok:
        passed.append(name)
    print(f'{name:50s} f = {str(f):4s} dL = {str(dL):3s} hbar  g\' = {gp:6.3f}  '
          f'photon {"ok" if c1 else "NO"}, EdH {"ok" if c2 else "NO"}, positions {"ok" if ok_pos else "NO"}  -> {"PASS" if ok else "fail"}  [{why}]')
print(f'\npass: {len(passed)} -> ' + ('; '.join(passed)))
