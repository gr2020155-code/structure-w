"""Pre-registration check (17.11.2026): do the rules R1-R6 fix alpha?
R1 cycle = lag of one exchange behind the global eight (A-60, V-212); R2 junction (-)loop to (+)end, transition around,
eight of two loops (A-24, V-170); R3 whole-form motion, each act a whole step (beta, iota) (V-198);
R4 energy = lag at a place, clock = lag along the particle (V-213); R5 one shared cycle per transition K(1-iota)=1 (A-45, V-158);
R6 unique trajectory (A-42). Units: a = c = 1, h = m*N (V-17).
If every relation holds for any beta, the rules do not fix alpha and running a model would only return the extra rule chosen."""
import math
for beta in (1/137.035999084, 0.01, 0.05, 0.2):
    iota = math.sqrt(1 - beta**2); K = 1/(1 - iota); N = 1.0; m = 1.0
    T = K*N/iota            # ticks per transition
    L = beta*T              # transfer steps per transition
    p = beta*m/iota         # gamma*beta*m (V-213)
    phase = p*L/(m*N)       # phase of the lag field around the transition, in cycles
    print(f"beta={beta:.6f}  K={K:10.1f}  L/N={L/N:9.2f}  2/beta={2/beta:9.2f}  phase around transition = {phase:.5f} cycles")
