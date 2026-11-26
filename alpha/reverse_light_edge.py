import sys
src = open(sys.argv[2], encoding='utf-8').read()
src = src.replace("np.linspace(0, pi - 0.02, 40)", "[0.0]").replace("np.arange(2.55, 2.905, 0.005)", "np.arange(2.665, 2.711, 0.005)")
src = src.replace("i = int(np.argmin", "[print(f'   v0 {a:.3f}: theta {t[0]:.4g}..{t[1]:.4g} SN dchi2 {x:+.1f}') for a,t,x in fine]\ni = int(np.argmin")
sys.argv = [sys.argv[2], sys.argv[1]]; exec(compile(src, 'rl', 'exec'))
