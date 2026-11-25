"""Global eights with F9 (side of the act fixed without contradiction) - F9 and Q4 fixed BEFORE running
in global_eights_rules_2026-11-18.md (addendum).

F1-F8 as in global_eights.py (2D only here). F9: in every exchange the two partners move on loops of the same
orientation (same sign of the signed area about O). Q2: number of solutions F1-F9 up to time shift and reversal.
Q4: if F9 is infeasible, the least share of acts with contradictory orientations.
"""
import numpy as np, sys, time
from ortools.sat.python import cp_model
from global_eights import cloud, area, extract, transforms, canon, verify, solve

def build_f9(X, E, O, hull, scale, f9=True, blocks=(), max_viol=None):
    n = len(X); T = 2 * n; m = cp_model.CpModel()
    x = {(t, e): m.NewBoolVar('') for t in range(T) for e in E}
    inc = {v: [e for e in E if v in e] for v in range(n)}
    for t in range(T):
        for v in range(n): m.AddExactlyOne(x[t, e] for e in inc[v])
    p = {(i, t, v): m.NewBoolVar('') for i in range(n) for t in range(T + 1) for v in range(n)}
    for i in range(n):
        for v in range(n):
            m.Add(p[i, 0, v] == int(v == i)); m.Add(p[i, T, v] == int(v == i))
        for t in range(1, T): m.AddExactlyOne(p[i, t, v] for v in range(n))
    for t in range(1, T):
        for v in range(n): m.AddExactlyOne(p[i, t, v] for i in range(n))
    D = [(u, v) for (a, b) in E for (u, v) in ((a, b), (b, a))]
    z = {}
    for i in range(n):
        for t in range(T):
            for (u, v) in D:
                e = (min(u, v), max(u, v)); zz = m.NewBoolVar(''); z[i, t, u, v] = zz
                m.AddBoolOr([p[i, t, u].Not(), x[t, e].Not(), p[i, t + 1, v]])
                m.AddImplication(zz, p[i, t, u]); m.AddImplication(zz, x[t, e])
                m.AddBoolOr([p[i, t, u].Not(), x[t, e].Not(), zz])
    for t in range(T):
        for e in E: m.AddBoolOr([x[t, e].Not(), x[(t + 1) % T, e].Not()])          # F8
    A = {(u, v): int(round(scale * area(X, O, u, v))) for (u, v) in D}
    big = T * (max(abs(a) for a in A.values()) + 1)
    sig = {}
    for i in range(n):
        m.Add(sum(p[i, t, O] for t in range(T)) == 2)
        c = [m.NewIntVar(0, 2, '') for t in range(T)]; b = [m.NewBoolVar('') for t in range(T)]
        for t in range(T):
            m.Add(c[t] == sum(p[i, s, O] for s in range(t + 1)))
            m.Add(c[t] == 1).OnlyEnforceIf(b[t]); m.Add(c[t] != 1).OnlyEnforceIf(b[t].Not())
        s1 = []
        for t in range(T):
            a = sum(A[u, v] * z[i, t, u, v] for (u, v) in D)
            w = m.NewIntVar(-big, big, ''); m.Add(w == a).OnlyEnforceIf(b[t]); m.Add(w == 0).OnlyEnforceIf(b[t].Not()); s1.append(w)
        tot = sum(A[u, v] * z[i, t, u, v] for t in range(T) for (u, v) in D)
        v1 = m.NewIntVar(-big, big, ''); m.Add(v1 == sum(s1)); v2 = m.NewIntVar(-big, big, ''); m.Add(v2 == tot - v1)
        o = m.NewBoolVar('')                                   # o: loop 1 positive
        m.Add(v1 >= 1).OnlyEnforceIf(o); m.Add(v2 <= -1).OnlyEnforceIf(o)
        m.Add(v1 <= -1).OnlyEnforceIf(o.Not()); m.Add(v2 >= 1).OnlyEnforceIf(o.Not())
        for t in range(T):                                     # sig = orientation of the current loop is positive
            s = m.NewBoolVar(''); sig[i, t] = s
            m.AddBoolOr([o.Not(), b[t].Not(), s]); m.AddBoolOr([o, b[t], s])
            m.AddBoolOr([o.Not(), b[t], s.Not()]); m.AddBoolOr([o, b[t].Not(), s.Not()])
    rho = {(t, v): m.NewBoolVar('') for t in range(T) for v in range(n)}
    for i in range(n):
        for t in range(T):
            for v in range(n):
                m.AddBoolOr([p[i, t, v].Not(), rho[t, v].Not(), sig[i, t]])
                m.AddBoolOr([p[i, t, v].Not(), rho[t, v], sig[i, t].Not()])
    viol = []
    for t in range(T):
        for (u, v) in E:
            if f9:
                m.AddBoolOr([x[t, (u, v)].Not(), rho[t, u].Not(), rho[t, v]])
                m.AddBoolOr([x[t, (u, v)].Not(), rho[t, u], rho[t, v].Not()])
            else:
                q = m.NewBoolVar(''); viol.append(q)
                m.AddBoolOr([x[t, (u, v)].Not(), rho[t, u].Not(), rho[t, v], q])
                m.AddBoolOr([x[t, (u, v)].Not(), rho[t, u], rho[t, v].Not(), q])
    if max_viol is not None: m.Add(sum(viol) <= max_viol)
    for sol in blocks: m.AddBoolOr([x[t, e].Not() for t in range(T) for e in sol[t]])
    return m, x, T

def f9_violations(X, O, M):
    """independent count: acts whose partners move on loops of opposite orientation"""
    n = len(X); T = len(M); pos = list(range(n)); traj = [[i] for i in range(n)]
    for t in range(T):
        new = pos[:]
        for (a, b) in M[t]:
            ia, ib = pos.index(a), pos.index(b); new[ia], new[ib] = b, a
        pos = new
        for i in range(n): traj[i].append(pos[i])
    sgn = np.zeros((n, T))
    for i in range(n):
        tr = traj[i]; t1, t2 = [t for t in range(T) if tr[t] == O]
        a1 = sum(area(X, O, tr[t], tr[t + 1]) for t in range(t1, t2))
        for t in range(T): sgn[i, t] = np.sign(a1) if t1 <= t < t2 else -np.sign(a1)
    bad = 0; where = None
    for t in range(T):
        at = {traj[i][t]: i for i in range(n)}
        for (a, b) in M[t]:
            if sgn[at[a], t] != sgn[at[b], t]: bad += 1
    return bad, n * T // 2

if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1].split(',')]; seed = int(sys.argv[2]); limit = float(sys.argv[3]); cap = int(sys.argv[4])
    for n in ns:
        X, E, O, hull = cloud(n, seed, 2); T = 2 * n
        print(f'== n {n}, seed {seed}: edges {len(E)}, T {T}, acts per period {n*T//2}')
        t0 = time.time(); orbits = []; blocks = []; status = None
        while len(orbits) < cap:
            m, x, T = build_f9(X, E, O, hull, 1000, f9=True, blocks=blocks); s, st = solve(m, limit)
            if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): status = s.StatusName(st); break
            M = extract(s, x, E, T); verify(X, E, O, hull, M, False, True)
            bad, tot = f9_violations(X, O, M); assert bad == 0, 'F9 not satisfied'
            orbits.append(canon(M)); blocks += transforms(M)
        print(f'   F1-F9: orbits {len(orbits)}' + (f' (stopped: {status})' if status else ' (cap)') + f'  [{time.time()-t0:.0f}s]'); sys.stdout.flush()
        if not orbits:          # Q4: least share of contradictory acts
            t0 = time.time(); best = None; K = n * T // 2
            lo, hi = 0, K
            while lo < hi:
                mid = (lo + hi) // 2
                m, x, T = build_f9(X, E, O, hull, 1000, f9=False, max_viol=mid); s, st = solve(m, limit)
                if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                    M = extract(s, x, E, T); verify(X, E, O, hull, M, False, True); bad, tot = f9_violations(X, O, M)
                    best = bad; hi = min(mid, bad)
                else:
                    lo = mid + 1 if st == cp_model.INFEASIBLE else hi
                    if st != cp_model.INFEASIBLE: break
            print(f'   Q4: least contradictory acts {best} of {K} = {best / K if best is not None else float("nan"):.4f}  [{time.time()-t0:.0f}s]')
            sys.stdout.flush()
