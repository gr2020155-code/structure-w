# -*- coding: utf-8 -*-
"""Проверка матрицы «утверждение — доказательство — следствие».

Читает matrix.json (по умолчанию рядом с программой) и проверяет:
  1. ID уникальны, префиксы и статусы допустимы;
  2. все ссылки существуют (proof, cond, cons, replaced_by, противоречия);
  3. снятые положения (С-) нигде не служат основанием;
  4. граф доказательств (proof + cond) без циклов;
  5. следствия (cons) — точное обращение proof + cond: ни одно утверждение
     не опирается на своё следствие;
  6. ни один [вывод]/[расчёт] не опирается (proof) на [кандидат]/[открыто]/[противоречие];
  7. условия (cond), полученные по цепочке proof, объявлены явно — скрытых условий нет;
  8. у [вывод] есть основание, у каждого утверждения — источник;
  9. каждое [противоречие] записано в разделе противоречий;
 10. калибровки по измеренным величинам не имеют статуса [вывод]/[расчёт];
 11. справка: неиспользованные данные/расчёты; что зависит от посылок не автора.
Запуск: python3 check_matrix.py [matrix.json]
"""
import json, sys, os

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "matrix.json")
M = json.load(open(path, encoding="utf-8"))

STAT = set(M["meta"]["statuses"])
WEAK = {"кандидат", "открыто", "противоречие"}
STRONG = {"вывод", "расчёт"}

fails = []
def check(name, problems):
    if problems:
        fails.append(name)
        print(f"FAIL {name}")
        for p in problems[:40]:
            print("     ", p)
        if len(problems) > 40:
            print(f"      … ещё {len(problems) - 40}")
    else:
        print(f"PASS {name}")

ax = {a["id"]: a for a in M["axioms"]}
da = {d["id"]: d for d in M["data"]}
ca = {c["id"]: c for c in M["calcs"]}
st = {s["id"]: s for s in M["statements"]}
wd = {w["id"]: w for w in M["withdrawn"]}
xs = {x["id"]: x for x in M["contradictions"]}
allid = {}
for grp in (ax, da, ca, st, wd, xs):
    allid.update(grp)

print(f"Файл: {path}")
print(f"Положений {len(ax)}, данных {len(da)}, расчётов {len(ca)}, утверждений {len(st)}, снято {len(wd)}, противоречий {len(xs)}")
print("Статусы утверждений:", ", ".join(f"{k} {v}" for k, v in sorted(
    {s: sum(1 for x in st.values() if x["status"] == s) for s in STAT}.items(), key=lambda t: -t[1]) if v))
print()

# 1
prob = []
n_all = sum(len(M[k]) for k in ("axioms", "data", "calcs", "statements", "withdrawn", "contradictions"))
if n_all != len(allid):
    prob.append(f"повторяющиеся ID: {n_all - len(allid)}")
for grp, pre in ((ax, "А-"), (da, "Д-"), (ca, "Р-"), (st, "В-"), (wd, "С-"), (xs, "Х-")):
    for i in grp:
        if not i.startswith(pre):
            prob.append(f"{i}: префикс не {pre}")
for s in st.values():
    if s["status"] not in STAT:
        prob.append(f"{s['id']}: статус «{s['status']}» не допустим")
check("1. ID уникальны, префиксы и статусы допустимы", prob)

# 2
prob = []
for s in st.values():
    for f in ("proof", "cond", "cons"):
        for r in s[f]:
            if r not in allid:
                prob.append(f"{s['id']}.{f}: нет {r}")
for grp in (ax, da, ca):
    for a in grp.values():
        for r in a.get("cons", []):
            if r not in st:
                prob.append(f"{a['id']}.cons: нет {r}")
for w in wd.values():
    for r in w["replaced_by"]:
        if r not in allid:
            prob.append(f"{w['id']}.replaced_by: нет {r}")
for x in xs.values():
    for r in x["ids"]:
        if r not in allid:
            prob.append(f"{x['id']}: нет {r}")
check("2. Все ссылки существуют", prob)

# 3
prob = [f"{s['id']} опирается на снятое {r}" for s in st.values() for r in s["proof"] + s["cond"] if r.startswith("С-")]
check("3. Снятые положения не служат основанием", prob)

# 4
def deps(i):
    s = st.get(i)
    return (s["proof"] + s["cond"]) if s else []
color, stack, cycles = {}, [], []
def dfs(u):
    color[u] = 1; stack.append(u)
    for w in deps(u):
        if color.get(w, 0) == 0:
            dfs(w)
        elif color.get(w) == 1:
            cycles.append(" → ".join(stack[stack.index(w):] + [w]))
    stack.pop(); color[u] = 2
sys.setrecursionlimit(10000)
for i in st:
    if color.get(i, 0) == 0:
        dfs(i)
check("4. Граф доказательств без циклов", cycles)

# 5
inv = {i: set() for i in allid}
for s in st.values():
    for r in s["proof"] + s["cond"]:
        if r in inv:
            inv[r].add(s["id"])
prob = []
for grp in (ax, da, ca, st):
    for a in grp.values():
        if set(a.get("cons", [])) != inv[a["id"]]:
            prob.append(f"{a['id']}: cons {sorted(a.get('cons', []))} ≠ обращение {sorted(inv[a['id']])}")
anc_memo = {}
def ancestors(i):
    if i in anc_memo:
        return anc_memo[i]
    out = set()
    for r in deps(i):
        out.add(r); out |= ancestors(r)
    anc_memo[i] = out
    return out
desc_memo = {}
def descendants(i):
    if i in desc_memo:
        return desc_memo[i]
    out = set()
    for r in inv.get(i, ()):
        out.add(r); out |= descendants(r)
    desc_memo[i] = out
    return out
if not cycles:
    for i in st:
        bad = ancestors(i) & descendants(i)
        if bad:
            prob.append(f"{i} опирается на своё следствие: {sorted(bad)}")
check("5. Следствия — обращение оснований; нет опоры на своё следствие", prob)

# 6
prob = []
for s in st.values():
    if s["status"] in STRONG:
        for r in s["proof"]:
            if r in st and st[r]["status"] in WEAK:
                prob.append(f"{s['id']} [{s['status']}] опирается на {r} [{st[r]['status']}]")
check("6. [вывод]/[расчёт] не опираются на [кандидат]/[открыто]/[противоречие]", prob)

# 7
inh_memo = {}
def inherited(i):
    """условия, накопленные по цепочке proof (cond — явные посылки)."""
    if i in inh_memo:
        return inh_memo[i]
    s = st[i]
    out = set(s["cond"])
    for r in s["proof"]:
        if r in st:
            out |= inherited(r)
    inh_memo[i] = out
    return out
prob = []
for s in st.values():
    if s["status"] in STRONG:
        hidden = inherited(s["id"]) - set(s["cond"])
        # условие, которое уже стоит в proof явно, не скрыто
        hidden -= set(s["proof"])
        if hidden:
            prob.append(f"{s['id']}: условия по цепочке не объявлены: {sorted(hidden)}")
check("7. Скрытых условий нет (всё, что принято условно по цепочке, объявлено в cond)", prob)

# 8
prob = [f"{s['id']}: [вывод] без основания" for s in st.values() if s["status"] == "вывод" and not s["proof"]]
prob += [f"{s['id']}: нет источника" for s in st.values() if not s["src"].strip()]
prob += [f"{a['id']}: нет источника" for g in (ax, da, ca) for a in g.values() if not a["src"].strip()]
check("8. У [вывод] есть основание, у всех пунктов — источник", prob)

# 9
inx = {r for x in xs.values() for r in x["ids"]}
prob = [f"{s['id']} [противоречие] не записано в разделе противоречий" for s in st.values() if s["status"] == "противоречие" and s["id"] not in inx]
check("9. Каждое [противоречие] записано в разделе противоречий", prob)

# 10
calib = [s["id"] for s in st.values() if s.get("calibration")]
prob = [f"{i}: калибровка со статусом [{st[i]['status']}]" for i in calib if st[i]["status"] in STRONG]
for s in st.values():
    if s["status"] in STRONG:
        hit = (set(s["proof"]) | {a for r in s["proof"] if r in st for a in ancestors(r)}) & set(calib)
        if hit:
            prob.append(f"{s['id']} [{s['status']}] опирается на калибровку {sorted(hit)}")
check("10. Калибровки по измеренным величинам не выданы за вывод", prob)
print("      калибровки:", ", ".join(f"{i} [{st[i]['status']}]" for i in calib))

# 11 справка
print()
print("Справка:")
used = {r for s in st.values() for r in s["proof"] + s["cond"]}
for grp, name in ((da, "данные"), (ca, "расчёты"), (ax, "положения")):
    un = sorted(i for i in grp if i not in used)
    print(f"  не использованы ({name}):", ", ".join(un) if un else "нет")
nonauth = [a["id"] for a in ax.values() if "не положение автора" in a["label"]]
for a in nonauth:
    dep = sorted((i for i in st if a in ancestors(i)), key=lambda t: int(t[2:]))
    print(f"  от посылки не автора {a} зависят ({len(dep)}):", ", ".join(dep))
conds = sorted({c for s in st.values() for c in s["cond"]}, key=lambda t: (t[0], int(t[2:])))
for c in conds:
    users = sorted((s["id"] for s in st.values() if c in s["cond"]), key=lambda t: int(t[2:]))
    print(f"  условно принято {c} [{st[c]['status'] if c in st else 'положение'}]: {', '.join(users)}")
down = [s["id"] for s in st.values() if s.get("was")]
print(f"  понижено против источника ({len(down)}):", ", ".join(down))

print()
print("ИТОГ:", "все проверки пройдены" if not fails else "не пройдено: " + "; ".join(fails))
sys.exit(1 if fails else 0)
