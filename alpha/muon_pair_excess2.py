"""Excess of mu+mu- pairs (B-355) vs CMD-3, SND, BaBar (26.11.2026).
Rules: muon_pair_excess2_rules_2026-11-26.md (pushed before run). Record: Р-111.
Run: python -I -X utf8 muon_pair_excess2.py
"""
from math import pi, sqrt
a = 1 / 137.036
d = a - 3 * a / (4 * pi)
sets = {'CMD-3 mumu/ee (Bhabha)': (1.0017, (0.0016, 0.0025, 0.0056), 'ee'),
        'SND mumu (gamma gamma)': (1.005, (0.007, 0.018), 'gg'),
        'BaBar mumu gamma (Bhabha)': (1.0040, (0.0020, 0.0055, 0.0094), 'ee')}
for name, (v, errs, norm) in sets.items():
    s = sqrt(sum(e * e for e in errs))
    out = []
    for br in ('B1', 'B2'):
        pred = 1.0 if (norm == 'ee' and br == 'B1') else 1 + d
        out.append(f'{br}: W {pred:.4f} ({(pred - v) / s:+.2f} s)')
    print(f'{name:26s} data {v:.4f} +- {s:.4f}: QED 1.0000 ({(1 - v) / s:+.2f} s); ' + '; '.join(out))
print(f'W - QED = {d*100:.3f} %; needed 2-sigma precision {d/2*100:.2f} %')
