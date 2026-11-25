"""Spin addition of form parts. Rules: spin_addition_rules_2026-11-26.md (written before running).

w_i = (s_i . J) / |J|^2 : projection of part i on the whole, in units where a fully aligned spin-1/2 part has w = 1.
"""
from fractions import Fraction as F

m_p, m_n = 938.272, 939.565
mu_p_data, mu_n_data = 2.79284734, -1.91304273
u, d = F(2, 3), F(-1, 3)
q = F(1, 4)                                   # (1/2)^2: projection squared on one axis

readings = {
    'M1 (along axis only)':           dict(s2=q,     dot=q,     J2=q),
    'M2*P1 (3 axes, independent)':    dict(s2=3 * q, dot=q,     J2=3 * q),
    'M2*P2 (3 axes, rigid pair)':     dict(s2=3 * q, dot=3 * q, J2=3 * q),
}

def weights(s2, dot, J2):
    S2 = 2 * s2 + 2 * dot                     # |s1 + s2|^2
    SJ = (J2 + S2 - s2) / 2
    s3J = (J2 + s2 - S2) / 2
    return SJ / J2 / 2, s3J / J2              # each pair member, third part

def moment(w_pair, w3, pair_q, third_q, M):
    unit = 3 * m_p / M                        # R2 (V-308): each of 3 parts has M/3
    return float(2 * w_pair * pair_q + w3 * third_q) * unit

print(f'{"reading":32} {"w pair":>7} {"w third":>8} {"mu_p":>7} {"mu_n":>7} {"ratio":>7}')
for name, r in readings.items():
    wp, w3 = weights(**r)
    mp = moment(wp, w3, u, d, m_p)            # proton: pair u,u; third d
    mn = moment(wp, w3, d, u, m_n)            # neutron (V-310): pair d,d; third u
    print(f'{name:32} {str(wp):>7} {str(w3):>8} {mp:7.3f} {mn:7.3f} {mn / mp:7.3f}')
print(f'{"data":32} {"":7} {"":8} {mu_p_data:7.3f} {mu_n_data:7.3f} {mu_n_data / mu_p_data:7.3f}')
print()
print('observation (not proof): |J|^2 of N aligned spin-1/2 parts vs standard N/2*(N/2+1)')
for name, r in readings.items():
    row = []
    for N in (1, 2, 3, 4):
        J2 = N * r['s2'] + N * (N - 1) * r['dot']
        std = F(N, 2) * (F(N, 2) + 1)
        row.append(f'N={N}: {J2} vs {std}{"" if J2 == std else " x"}')
    print(f'{name:32} ' + ';  '.join(row))
