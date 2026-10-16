"""Independent reference values for cm_threshold() targets (v0.6 item 12;
design note section 19). Builds on prototype_v06_sensitivity.py (which uses
reference_v041.py and reference_v05_input_check.py), not on the package.

The loss of a set of diseases is the sum of their contributions (Shapley
values; for event impacts, of the attributable risk). A target is a value,
or a fraction of the set's adjusted value (the model as given) or of its
raw value (the same model with the diseases unrelated: the independent
joint distribution, the same impact rows). One input is varied, everything
else at its central value with the pooling of association rows frozen at
the central fit (the package refits it; the two agree here, since no pooled
row is varied); the threshold is the input value where the set's loss
equals the level (scan on a grid, then Brent's method).

Printed: the reference value (adjusted or raw) of the set, the level, and
each threshold.

Run from this directory: python3 reference_v06_thresholds.py
"""
import math

import numpy as np
from scipy.optimize import brentq

import prototype_v06_sensitivity as P
import reference_v041 as R
from reference_v05_input_check import Case, event_adjust
from reference_v06_event_parity import shapley_event


def set_loss(out, ds):
    return out["agg"] if ds is None else sum(out["contrib"][d] for d in ds)


def raw_contrib(case):
    """Each disease's contribution with the diseases unrelated."""
    probs = [1 - math.exp(-float(r["value"])) if r["measure"] == "incidence_rate" else float(r["value"])
             for r in case.drows]
    pop = R.Pop(case.ids, probs)
    p = np.prod(np.where(pop.cells == 1, pop.P, 1 - pop.P), axis=1)
    rows = [dict(r) for r in case.irows]
    if case.event:
        hr, ar, alone = event_adjust(pop, p, rows, case.risk)
        terms = list(hr)
        theta = [math.log(hr[t]) for t in terms]
        h0 = -math.log(1 - (case.risk - ar))
        return shapley_event(pop, p, terms, theta, h0)
    b, contrib, total, *_ = R.adjust_additive(pop, p, rows)
    return contrib


def crossings(case, key, lo, hi, ds, level, n=41, log=False):
    xs = np.exp(np.linspace(math.log(lo), math.log(hi), n)) if log else np.linspace(lo, hi, n)
    f = lambda x: set_loss(P.run(case, {key: float(x)}), ds) - level
    fx = [f(x) for x in xs]
    out = []
    for a, b, fa, fb in zip(xs[:-1], xs[1:], fx[:-1], fx[1:]):
        if fa == 0:
            out.append(a)
        elif fa * fb < 0:
            out.append(brentq(f, a, b, xtol=1e-12, rtol=1e-12))
    return out


def show(name, case, key, lo, hi, ds, target, relative_to=None, log=False):
    central = P.run(case, {})
    if relative_to == "adjusted":
        ref = set_loss(central, ds)
    elif relative_to == "raw":
        rc = raw_contrib(case)
        ref = sum(rc[d] for d in (ds or case.ids))
    else:
        ref = None
    level = target * ref if ref is not None else target
    th = crossings(case, key, lo, hi, ds, level, log=log)
    print(f"\n== {name}: {key} in [{lo}, {hi}], set {ds or 'all'}, target {target}"
          + (f" x {relative_to}" if relative_to else ""))
    print(f"central set loss {set_loss(central, ds):.10f}" + (f"; reference {ref:.10f}" if ref is not None else "")
          + f"; level {level:.10f}")
    print("thresholds", [round(float(t), 10) for t in th])


if __name__ == "__main__":
    yield_case = Case("five_diseases", "impacts_yield.csv")
    welfare = Case("five_diseases", "impacts_welfare.csv")
    cull = Case("five_diseases", "culling.csv", event=True, risk=0.25)
    # A. Yield, LAM-SCK OR (associations:3), loss of LAM, MAS and SCK.
    show("A1", yield_case, "associations:3", 1.0, 4.0, ["LAM", "MAS", "SCK"], 1.95)
    show("A2", yield_case, "associations:3", 1.0, 4.0, ["LAM", "MAS", "SCK"], 0.97, "adjusted")
    show("A3", yield_case, "prob:SCK", 0.2, 0.6, ["LAM", "MAS", "SCK"], 0.75, "raw")
    show("A4", yield_case, "prob:SCK", 0.2, 0.6, None, 0.75, "raw")
    # B. Welfare (interactions): raw contributions from the refit on the
    # independent joint distribution.
    print("\n== B0: raw contributions, welfare", {d: round(v, 10) for d, v in raw_contrib(welfare).items()})
    show("B1", welfare, "impacts:6", 0.0, 6.0, ["LAM", "MAS"], 0.8, "raw")
    show("B1b", welfare, "impacts:1", 5.0, 15.0, ["LAM", "MAS"], 0.85, "raw")
    show("B2", welfare, "prob:LAM", 0.1, 0.4, ["LAM", "MAS"], 3.0)
    # C. Culling (event impacts), attributable risk of LAM and MAS.
    print("\n== C0: raw contributions, culling", {d: round(v, 10) for d, v in raw_contrib(cull).items()})
    show("C1", cull, "impacts:2", 0.8, 2.4, ["LAM", "MAS"], 0.03)
    show("C2", cull, "impacts:6", 2.0, 5.0, ["LAM", "MAS"], 0.9, "adjusted", log=True)
    show("C3", cull, "risk", 0.15, 0.35, ["LAM", "MAS"], 0.6, "raw")
