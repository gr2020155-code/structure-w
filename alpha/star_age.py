"""Age of the oldest stars vs today's rate: the whole cannot be younger than its stars (R-159).
Rules: star_age_rules_2026-11-28.md (committed before the run, b9e7b51). Data: D-96 (Valcin et al. III, IV).
usage: python -P -E -X utf8 star_age.py"""
from math import sqrt

GYR = 977.79222 / 100                     # 1/(100 km/s/Mpc) in Gyr
RATES = [('Cepheids', 0.7354, 0.011), ('red giants', 0.710, 0.011), ('far rate', 0.685, 0.010), ('clocks', 0.687, 0.040)]
STARS = [('III (main)', 13.39, sqrt(0.10**2 + 0.23**2)), ('IV', 13.61, sqrt(0.25**2 + 0.23**2))]
SHAPES = [('ruler, Om = 0.2973 (main)', 0.9663), ('supernovae, Om ~ 0.33', 0.938)]

for sname, tg in SHAPES:
    for gname, t, st in STARS:
        main = sname.endswith('(main)') and gname.endswith('(main)')
        hmax = tg * GYR / t; shmax = hmax * st / t
        print(f'{"MAIN  " if main else ""}shape {sname}; clusters {gname}: t = {t} +- {st:.2f} Gyr -> h_max = {hmax:.4f} +- {shmax:.4f}')
        for lab, h, sh in RATES:
            age = tg * GYR / h; sage = age * sh / h
            s = (t - age) / sqrt(st**2 + sage**2)
            verdict = 'EXCLUDED' if s > 2 else 'not excluded'
            print(f'   {lab:10s} h = {h}: age {age:.2f} +- {sage:.2f} Gyr; younger than clusters by {s:+.2f} sigma -> {verdict}')
        print()
