"""W and Z decay fractions in the two-part reading (26.11.2026).
Rules: wz_decays_rules_2026-11-26.md (pushed before run). Record: Р-106.
Run: python -I -X utf8 wz_decays.py
"""
from math import log, sqrt

c, dc, E0 = 0.251, 0.017, 0.3128
def eps(M):
    L = log(M / E0); return c / L, dc / L

print('W (M = 80.3625 GeV)')
data_W = {'lnu': (10.86, 0.09), 'had': (67.41, 0.27)}
flav = {'e nu': (10.71, 0.16), 'mu nu': (10.63, 0.15), 'tau nu': (11.38, 0.21)}
e, de = eps(80.3625)
for name, ee, dee in (('W0 no excess', 0.0, 0.0), ('W1 with excess', e, de)):
    q = 6 * (1 + ee); T = 3 + q
    lnu, had = 100 / T, 100 * q / T
    dl = 100 * 6 * dee / T**2; dh = 100 * 3 * 6 * dee / T**2
    p1 = (lnu - data_W['lnu'][0]) / sqrt(data_W['lnu'][1]**2 + dl**2)
    p2 = (had - data_W['had'][0]) / sqrt(data_W['had'][1]**2 + dh**2)
    ok = abs(p1) <= 2 and abs(p2) <= 2
    print(f'  {name:15s} eps = {ee:.4f} +- {dee:.4f}: l nu = {lnu:.3f} % ({p1:+.1f} sigma), hadrons = {had:.2f} % ({p2:+.1f} sigma) -> {"PASS" if ok else "fail"}')
    print('     by flavour (info): ' + ', '.join(f'{k} {(lnu-v)/s:+.1f} sigma' for k, (v, s) in flav.items()))
    ex = ((100 / data_W['lnu'][0]) - 3) / 6 - 1
    if ee: print(f'     excess from data: eps = {ex:.3f} +- {(100/data_W["lnu"][0]**2*data_W["lnu"][1])/6:.3f}')

print('\nZ (M = 91.19 GeV)')
data_Z = {'inv': (20.000, 0.055), 'l': (3.3658, 0.0023), 'had': (69.911, 0.056)}
e, de = eps(91.19)
weights = {
    'Z1 equal per variant':           lambda k: 1.0,
    'Z2 (k/3)^2 (photon-like)':       lambda k: (k / 3) ** 2,
    'Z3 (k/3)^2 + ((3-k)/3)^2':       lambda k: (k / 3) ** 2 + ((3 - k) / 3) ** 2,
    'Z4 (k/3)(1-k/3)':                lambda k: (k / 3) * (1 - k / 3),
}
for name, w in weights.items():
    inv = 3 * w(0); lep1 = w(3); up = 2 * 3 * w(2) * (1 + e); dn = 3 * 3 * w(1) * (1 + e)
    T = inv + 3 * lep1 + up + dn
    if inv == 0 or lep1 == 0:
        print(f'  {name:30s} invisible {100*inv/T:.2f} %, l {100*lep1/T:.3f} % -> fail (a seen channel is zero)')
        continue
    r = {'inv': 100 * inv / T, 'l': 100 * lep1 / T, 'had': 100 * (up + dn) / T}
    pulls = {k: (r[k] - data_Z[k][0]) / data_Z[k][1] for k in r}
    ok = all(abs(p) <= 2 for p in pulls.values())
    print(f'  {name:30s} ' + ', '.join(f'{k} {r[k]:.3f} % ({pulls[k]:+.0f} sigma)' for k in r) + f' -> {"PASS" if ok else "fail"}')
