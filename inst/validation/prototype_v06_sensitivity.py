"""Prototype of cm_sensitivity() (v0.6 item 3; design note section 17).
Builds on reference_v041.py and reference_v05_input_check.py; does not use
the package.

For each input, the model is rerun with that input at a low and a high value
and everything else at its central value (pooling of association rows frozen
at the central fit, as in the draws):
- entered inputs with a distribution: -1 and +1 SD on the normal score
  u = Phi^-1(F(x)) (as in the diagnosis);
- entered inputs without one: x (1 - change) and x (1 + change), change 0.2
  by default (probabilities kept inside (0, 1));
- unknowns (not entered):
  pair:d1:d2         a pair that no association row mentions: OR 0.5 and 3
  interaction:d1;d2  a pair without a combination term: an impact row for
                     the combination adjusted for all, at -/+ half the
                     smaller absolute adjusted impact of the two (event
                     impacts: HR 0.67 and 1.5)
  or_ratio:d1:d2|d3  only on request: OR_ratio 0.5 and 2.
Reported: the aggregate (attributable risk for event impacts) at low and
high, the swing, each disease's adjusted value and the largest change in a
disease's rank by contribution.

Draws (all inputs with a distribution at once, n draws of the normal
scores): each input's share of the variance of the aggregate, as the squared
standardised regression coefficient of the aggregate on the normal scores
(with R^2), and the probability that each disease ranks 1st, 2nd, ... by
contribution.

Run from this directory: python3 prototype_v06_sensitivity.py [case]
"""
import itertools
import math
import sys

import numpy as np
from scipy import stats

import reference_v041 as R
from reference_v05_input_check import Case, fit_frozen, event_adjust, read, u_of
from reference_v06_event_parity import shapley_event

CHANGE = 0.2


def run(case, changes, extra_assoc=(), extra_impacts=()):
    """Case.run with extra rows (unknowns); also returns contributions."""
    P = []
    for r in case.drows:
        v = changes.get(f"prob:{r['disease']}", float(r["value"]))
        P.append(1 - math.exp(-v) if r["measure"] == "incidence_rate" else v)
    pop = R.Pop(case.ids, P)
    arows = []
    for r in case.arows:
        rr = dict(r)
        rr["value"] = repr(changes.get(f"associations:{r['row']}", float(r["value"])))
        arows.append(rr)
    n0 = len(arows)
    for i, r in enumerate(extra_assoc):
        arows.append(dict(r, row=n0 + i + 1))
    specs = [R.spec(pop, r) for r in arows]
    p, groups, info = fit_frozen(pop, specs, case.frozen)
    rows = []
    for r in case.irows:
        rr = dict(r)
        rr["value"] = repr(changes.get(f"impacts:{r['row']}", float(r["value"])))
        rows.append(rr)
    n0 = len(rows)
    for i, r in enumerate(extra_impacts):
        rows.append(dict(r, row=n0 + i + 1))
    if case.event:
        risk = changes.get("risk", case.risk)
        hr, ar, alone = event_adjust(pop, p, rows, risk)
        terms = list(hr)
        theta = [math.log(hr[t]) for t in terms]
        h0 = -math.log(1 - (risk - ar))
        contrib = shapley_event(pop, p, terms, theta, h0)
        return dict(adj=hr, agg=ar, contrib=contrib)
    b, contrib, total, *_ = R.adjust_additive(pop, p, rows)
    return dict(adj=b, agg=total, contrib=contrib)


def ranks(contrib, ids):
    order = sorted(ids, key=lambda d: -contrib[d])
    return {d: order.index(d) + 1 for d in ids}


def entered_ranges(case, change=CHANGE):
    out = []
    for key, lab, v, qf, kind, meas in case.inputs:
        if qf is not None:
            u0 = u_of(qf, v)
            lo, hi = float(qf(stats.norm.cdf(u0 - 1))), float(qf(stats.norm.cdf(u0 + 1)))
            how = "1 SD"
        else:
            lo, hi = v * (1 - change), v * (1 + change)
            if kind == "prob" and meas != "incidence_rate":
                lo, hi = max(lo, 1e-6), min(hi, 1 - 1e-6)
            how = f"{change:.0%}"
        out.append(dict(key=key, label=lab, kind="entered", how=how, value=v, low=lo, high=hi))
    if case.event and not any(i["key"] == "risk" for i in out):
        v = case.risk
        out.append(dict(key="risk", label="overall risk", kind="entered", how=f"{change:.0%}", value=v,
                        low=v * (1 - change), high=v * (1 + change)))
    return out


def unknown_ranges(case, central, three_way=()):
    ids = case.ids
    out = []
    # Pairs that no association row mentions.
    mentioned = set()
    for r in case.arows:
        ds = set(R.split(r["disease"]) + R.split(r.get("associated_with", "")))
        ds |= {x.replace("no ", "").strip() for x in R.split(r.get("among", "")) if x.strip() and x.strip() != "no others"}
        for a, b in itertools.combinations(sorted(ds & set(ids)), 2):
            mentioned.add((a, b))
    for a, b in itertools.combinations(ids, 2):
        if tuple(sorted((a, b))) not in mentioned:
            out.append(dict(key=f"pair:{a}:{b}", label=f"{a}-{b} OR (unknown)", kind="unknown", how="OR 0.5-3",
                            value=None, low=0.5, high=3.0,
                            row=lambda v, a=a, b=b: ("assoc", dict(disease=a, associated_with=b, among="",
                                                                   adjusted_for="", measure="OR", value=repr(float(v))))))
    # Interactions for pairs without a combination term.
    terms = set(central["adj"])
    for a, b in itertools.combinations(ids, 2):
        if f"{a};{b}" in terms:
            continue
        if case.event:
            lo, hi, how = 0.67, 1.5, "HR 0.67-1.5"
            mk = lambda v, a=a, b=b: ("impact", dict(disease=f"{a};{b}", among="", adjusted_for="all",
                                                     measure="HR", value=repr(float(v))))
        else:
            half = 0.5 * min(abs(central["adj"][a]), abs(central["adj"][b]))
            lo, hi, how = -half, half, "+/-50% of smaller"
            mk = lambda v, a=a, b=b: ("impact", dict(disease=f"{a};{b}", among="", adjusted_for="all",
                                                     value=repr(float(v))))
        out.append(dict(key=f"interaction:{a};{b}", label=f"{a};{b} interaction (unknown)", kind="unknown",
                        how=how, value=None, low=lo, high=hi, row=mk))
    for a, b, c in three_way:
        out.append(dict(key=f"or_ratio:{a}:{b}|{c}", label=f"{a}-{b} OR ratio among {c} (unknown)",
                        kind="unknown", how="0.5-2", value=None, low=0.5, high=2.0,
                        row=lambda v, a=a, b=b, c=c: ("assoc", dict(disease=a, associated_with=b, among=c,
                                                                    adjusted_for="", measure="OR_ratio",
                                                                    value=repr(float(v))))))
    return out


def at(case, inp, x):
    if inp["kind"] == "entered":
        return run(case, {inp["key"]: x})
    tab, row = inp["row"](x)
    return run(case, {}, extra_assoc=[row] if tab == "assoc" else [], extra_impacts=[row] if tab == "impact" else [])


def one_at_a_time(case, three_way=(), change=CHANGE):
    central = run(case, {})
    r0 = ranks(central["contrib"], case.ids)
    rows = []
    for inp in entered_ranges(case, change) + unknown_ranges(case, central, three_way):
        res = []
        for x in (inp["low"], inp["high"]):
            try:
                res.append(at(case, inp, x))
            except Exception as e:  # infeasible inputs
                res.append(None)
        lo, hi = res
        aggs = [r["agg"] if r else float("nan") for r in res]
        shift = max((max(abs(ranks(r["contrib"], case.ids)[d] - r0[d]) for d in case.ids) for r in res if r),
                    default=float("nan"))
        rows.append(dict(inp, agg_low=aggs[0], agg_high=aggs[1], swing=abs(aggs[1] - aggs[0]),
                         adj_low={d: lo["adj"][d] for d in case.ids} if lo else None,
                         adj_high={d: hi["adj"][d] for d in case.ids} if hi else None, rank_shift=shift,
                         failed=sum(r is None for r in res)))
    rows.sort(key=lambda r: -(r["swing"] if np.isfinite(r["swing"]) else np.inf))
    return central, rows


def draws(case, n=400, seed=1):
    unc = [(key, qf) for key, lab, v, qf, kind, meas in case.inputs if qf is not None]
    rng = np.random.default_rng(seed)
    U = rng.standard_normal((n, len(unc)))
    agg, rk, ok = [], [], []
    for i in range(n):
        ch = {key: float(qf(stats.norm.cdf(U[i, j]))) for j, (key, qf) in enumerate(unc)}
        try:
            r = run(case, ch)
        except Exception:
            continue
        if not np.isfinite(r["agg"]):
            continue
        ok.append(i); agg.append(r["agg"]); rk.append(ranks(r["contrib"], case.ids))
    U, agg = U[ok], np.array(agg)
    X = np.column_stack([np.ones(len(agg)), U])
    beta = np.linalg.lstsq(X, agg, rcond=None)[0]
    fitted = X @ beta
    r2 = 1 - np.sum((agg - fitted) ** 2) / np.sum((agg - agg.mean()) ** 2)
    src = beta[1:] * U.std(axis=0) / agg.std()
    share = dict(zip([k for k, _ in unc], src ** 2))
    k = len(case.ids)
    rank_p = {d: np.bincount([r[d] for r in rk], minlength=k + 1)[1:] / len(rk) for d in case.ids}
    mean_rank = {d: float(np.mean([r[d] for r in rk])) for d in case.ids}
    return dict(n=n, accepted=len(ok), r2=float(r2), share=share, rank_p=rank_p, mean_rank=mean_rank,
                agg_mean=float(agg.mean()), agg_sd=float(agg.std()))


def report(name, case, three_way=(), n_draws=400, top=15):
    central, rows = one_at_a_time(case, three_way)
    lab = "attributable risk" if case.event else "aggregate"
    print(f"\n== {name}\ncentral {lab} {central['agg']:.10f}; inputs varied: {len(rows)} "
          f"({sum(r['kind'] == 'unknown' for r in rows)} unknown)")
    dr = draws(case, n_draws) if n_draws else None
    print(f"{'input':24s} {'kind':8s} {'how':18s} {'low':>10s} {'high':>10s} {lab + ' low':>14s} {'high':>12s}"
          f" {'swing':>10s} {'rank':>4s} {'var share':>9s}")
    for r in rows[:top]:
        sh = dr["share"].get(r["key"]) if dr else None
        print(f"{r['key']:24s} {r['kind']:8s} {r['how']:18s} {r['low']:10.5g} {r['high']:10.5g} "
              f"{r['agg_low']:14.8f} {r['agg_high']:12.8f} {r['swing']:10.6f} {r['rank_shift']:4.0f} "
              f"{'' if sh is None else f'{sh:9.1%}'}")
    if dr:
        print(f"draws: {dr['accepted']}/{dr['n']} accepted; R^2 of the linear fit {dr['r2']:.3f}; "
              f"{lab} mean {dr['agg_mean']:.5f}, sd {dr['agg_sd']:.5f}")
        top_share = sorted(dr["share"].items(), key=lambda kv: -kv[1])[:8]
        print("  variance shares:", ", ".join(f"{k} {v:.1%}" for k, v in top_share))
        print("  rank probabilities (P rank 1, P rank 2, ...; mean rank):")
        for d in sorted(case.ids, key=lambda d: dr["mean_rank"][d]):
            print(f"    {d:5s} " + " ".join(f"{x:5.2f}" for x in dr["rank_p"][d]) + f"   mean {dr['mean_rank'][d]:.2f}")
    return central, rows, dr


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    np.set_printoptions(linewidth=160)
    if which in ("all", "A"):
        report("A. five diseases, yield", Case("five_diseases", "impacts_yield.csv"))
    if which in ("all", "B"):
        arows = read("five_diseases", "associations.csv")
        arows = [r for r in arows if not (r["disease"] == "MAS" and r["associated_with"] == "RP")]
        for i, r in enumerate(arows):
            r["row"] = i + 1
        report("B. five diseases, yield, MAS-RP unknown, three-way MAS:RP|MET on request",
               Case("five_diseases", "impacts_yield.csv", assoc_rows=arows), three_way=[("MAS", "RP", "MET")])
    if which in ("all", "C"):
        report("C. five diseases, culling (overall risk 0.25)", Case("five_diseases", "culling.csv", event=True, risk=0.25),
               n_draws=200)
    if which in ("all", "D"):
        report("D. global dairy 2024, yield", Case("global_dairy_2024", "impacts_yield.csv"), n_draws=300, top=20)
