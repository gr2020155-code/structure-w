import sys
src = open(sys.argv[2], encoding='utf-8').read().split("print(f'SN reference")[0]
sys.argv = [sys.argv[2], sys.argv[1]]; exec(compile(src, 'tv', 'exec'))
for v0 in (2.675, 2.230, 0.0):
    S, dr = rays(v0)
    j = B.size - 1
    th = [rs / np.interp(Dr, sgrid, S[:, j]) for rs, Dr in relic]
    print(f'v0 {v0}: meridian (beta=pi) relic theta {min(th):.4g}..{max(th):.4g}')
    for deg in (30, 60, 90, 120, 150, 180):
        jj = int(round(deg / 2.5)); Sz = [np.interp(D(zz), sgrid, S[:, jj]) for zz in (0.1, 0.5, 1.0)]
        print(f'   beta {deg:3d}: S at z = 0.1, 0.5, 1: ' + ', '.join(f'{x:.3g}' for x in Sz) + f'  (flat {D(0.1):.3g}, {D(0.5):.3g}, {D(1.0):.3g})')
