"""Independent reference values for cm_threshold_prob() (v0.6 item 12;
design note section 19). Builds on reference_v041.py, not on the package.

A case with a closed-form answer: the supplement model (three diseases,
crude pair rows, one crude impact row per disease) with one uncertain input,
the impact row of d1 ~ normal(2.5, 0.5), and the varied input the d2-d3
odds ratio x. For fixed x the adjusted impacts are linear in the raw
values, so the aggregate (and the loss of any set) is
    L(x) = mu(x) + s(x) * (I1 - 2.5),
with mu(x) the loss at the point values and s(x) its derivative with
respect to the d1 row. Then
    P(L < t) = Phi((t - mu(x)) / (|s(x)| * 0.5))      (either sign of s),
and the 50% crossing is where mu(x) = t (the deterministic threshold).

Printed, for the aggregate (t = 2.15) and for the loss of d1 and d2
(t = 0.75): mu, s and P at x = 1, 2, 3, 5, 8, and the x where P = 0.5, 0.8
and 0.95.

Run from this directory: python3 reference_v06_threshold_prob.py
"""
import math

import numpy as np
from scipy import stats
from scipy.optimize import brentq

import reference_v041 as R

IDS = ["d1", "d2", "d3"]
PREV = [0.10, 0.15, 0.20]
RAW = [2.5, 5.0, 7.5]
SD1 = 0.5


def loss(x, raw1, ds):
    pop = R.Pop(IDS, PREV)
    arows = [dict(row=1, disease="d1", associated_with="d2", among="", adjusted_for="", measure="OR", value="2"),
             dict(row=2, disease="d1", associated_with="d3", among="", adjusted_for="", measure="OR", value="1"),
             dict(row=3, disease="d2", associated_with="d3", among="", adjusted_for="", measure="OR",
                  value=repr(float(x)))]
    p, groups, info = R.fit(pop, [R.spec(pop, r) for r in arows])
    rows = [dict(row=i + 1, disease=d, among="", adjusted_for="", value=repr(float(v)))
            for i, (d, v) in enumerate(zip(IDS, [raw1] + RAW[1:]))]
    b, contrib, total, *_ = R.adjust_additive(pop, p, rows)
    return total if ds is None else sum(contrib[d] for d in ds)


def mu_s(x, ds):
    mu = loss(x, RAW[0], ds)
    h = 1e-4
    s = (loss(x, RAW[0] + h, ds) - loss(x, RAW[0] - h, ds)) / (2 * h)
    return mu, s


def prob(x, ds, t):
    mu, s = mu_s(x, ds)
    return stats.norm.cdf((t - mu) / (abs(s) * SD1))


def show(ds, t):
    print(f"\n== set {ds or 'all'}, target {t}")
    for x in (1, 2, 3, 5, 8):
        mu, s = mu_s(x, ds)
        print(f"x {x}: mu {mu:.10f} s {s:.10f} P {prob(x, ds, t):.10f}")
    for lv in (0.5, 0.8, 0.95):
        f = lambda x: prob(x, ds, t) - lv
        try:
            print(f"level {lv}: x {brentq(f, 1.0, 10.0, xtol=1e-12):.10f}")
        except ValueError:
            print(f"level {lv}: not crossed in [1, 10]")


if __name__ == "__main__":
    show(None, 2.15)
    show(["d1", "d2"], 0.75)
