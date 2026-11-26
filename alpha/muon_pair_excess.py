"""Excess of mu+mu- pair production predicted by B-355 vs BaBar (26.11.2026).
Rules: muon_pair_excess_rules_2026-11-26.md (pushed before run). Record: Р-110.
Run: python -I -X utf8 muon_pair_excess.py
"""
from math import pi, sqrt
alpha = 1 / 137.036
eps_W, eps_QED = alpha, 3 * alpha / (4 * pi)
diff = eps_W - eps_QED
data, errs = 1 + 40e-4, (20e-4, 55e-4, 94e-4)
s = sqrt(sum(e * e for e in errs))
for name, pred in (('W (K1, eps = alpha)', 1 + diff), ('QED (NLO)', 1.0)):
    print(f'{name:20s} ratio to QED {pred:.4f}; data {data:.4f} +- {s:.4f} -> {(pred - data) / s:+.2f} sigma')
print(f'eps_W = {eps_W*100:.3f} %, eps_QED = {eps_QED*100:.3f} %, W - QED = {diff*100:.3f} %')
print(f'needed total precision to separate at 2 sigma: {diff/2*100:.2f} %  (now {s*100:.2f} %, luminosity alone {errs[2]*100:.2f} %)')
