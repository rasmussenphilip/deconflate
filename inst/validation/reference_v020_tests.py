"""Reference values for the 0.2.0 tests of hazard ratios, comparisons,
reports, sensitivity and simulation (impacts in percent: 2.5, 5, 7.5).
Run from this directory: python3 reference_v020_tests.py
Requires numpy and scipy. Values printed here are hard-coded in
tests/testthat/test-hazard.R, test-compare.R, test-report.R,
test-shapley-sensitivity.R and test-simulation.R.
"""
import itertools
import numpy as np
from scipy.optimize import fsolve
from python_reference import or_to_joint, excess_matrix, published, simultaneous, ipf, crude

P = np.array([.10, .15, .20])
OR = np.ones((3, 3)); OR[0, 1] = OR[1, 0] = 2; OR[1, 2] = OR[2, 1] = 3
m = np.array([2.5, 5, 7.5])
A = np.eye(3) + excess_matrix(P, OR).T
CELLS = np.array(list(itertools.product([0, 1], repeat=3)))


def fl(x):
    return [float(v) for v in np.atleast_1d(x)]


def fit_three_way(ratio, sweeps=5000):
    """Maximum-entropy joint with a three-way term: IPF from a start with
    the term ratio ** (D1 D2 D3), as fit_joint() does."""
    p = np.prod(np.where(CELLS == 1, P, 1 - P), axis=1) * np.where(CELLS.sum(1) == 3, ratio, 1.0)
    p /= p.sum()
    tg = []
    for i, k in itertools.combinations(range(3), 2):
        p11 = or_to_joint(P[i], P[k], OR[i, k])
        tg.append((i, k, {(1, 1): p11, (1, 0): P[i] - p11, (0, 1): P[k] - p11,
                          (0, 0): 1 - P[i] - P[k] + p11}))
    for _ in range(sweeps):
        for i, k, t in tg:
            for (a, b), v in t.items():
                msk = (CELLS[:, i] == a) & (CELLS[:, k] == b)
                p[msk] *= v / p[msk].sum()
    return p


if __name__ == "__main__":
    # Methods, aggregates and productivity gaps (observed 10,000, percent).
    for name, adj in [("published", published(P, OR, m)), ("simultaneous", simultaneous(P, OR, m))]:
        L = adj @ P
        xh = 10000 / (1 - L / 100)
        print(name, "adjusted", fl(adj), "aggregate", float(L))
        print("  decrease: disease-free", float(xh), "gap", float(xh - 10000),
              "attribution", fl(xh * adj * P / 100))
        print("  increase: disease-free", float(10000 / (1 + L / 100)))
    print("sign change, raw (0.1, 5, 7.5): simultaneous", fl(simultaneous(P, OR, np.array([.1, 5, 7.5]))),
          "published", fl(published(P, OR, np.array([.1, 5, 7.5]))))
    print("published denominator offset of d2:", float((A[1] - np.eye(3)[1]) @ m))

    # Association screen (simultaneous aggregate with one OR scaled).
    for (i, k) in [(0, 1), (0, 2), (1, 2)]:
        for f in (0.5, 2):
            O2 = OR.copy(); O2[i, k] = O2[k, i] = OR[i, k] * f
            print("screen d%d:d%d x %g:" % (i + 1, k + 1, f), float(simultaneous(P, O2, m) @ P))

    # Global method with an interaction d1:d2 = 1, with three-way ratios;
    # and the snapshot hazard-ratio model with a three-way ratio.
    HR = np.array([1.5, 2.0, 1.3])
    for r in (1, 0.5, 2):
        p = fit_three_way(r)
        J = (CELLS * p[:, None]).T @ CELLS
        S = J - np.outer(P, P); np.fill_diagonal(S, P * (1 - P))
        Ar = S / np.diag(S)[:, None]
        g = CELLS[:, 0] * CELLS[:, 1] * 1.0
        off = (CELLS.T @ (p * g) - P * (p @ g)) / np.diag(S)
        b = np.linalg.solve(Ar, m - off)
        print("three-way ratio %g: global total with interaction" % r, float(P @ b + J[0, 1]))

        def F(beta):
            w = p * np.exp(CELLS @ beta)
            return [np.log(w[CELLS[:, i] == 1].sum() / p[CELLS[:, i] == 1].sum())
                    - np.log(w[CELLS[:, i] == 0].sum() / p[CELLS[:, i] == 0].sum()) - np.log(HR[i])
                    for i in range(3)]
        print("  snapshot hazard ratios", fl(np.exp(fsolve(F, np.log(HR), xtol=1e-13))))

    # Snapshot hazard ratios with d1 adjusted for d2 (Mantel-Haenszel-type
    # weighted log ratio within strata of d2).
    cells, p = ipf(P, OR, iters=4000)

    def F_strat(beta):
        w = p * np.exp(cells @ beta)
        x, s = cells[:, 0], cells[:, 1]
        num = den = 0.0
        for k in (0, 1):
            msk = s == k
            d1 = p[msk & (x == 1)].sum(); d0 = p[msk & (x == 0)].sum()
            r = np.log(w[msk & (x == 1)].sum() / d1) - np.log(w[msk & (x == 0)].sum() / d0)
            wt = d1 * d0 / (d1 + d0); num += wt * r; den += wt
        out = [num / den - np.log(HR[0])]
        for i in (1, 2):
            x = cells[:, i]
            out.append(np.log(w[x == 1].sum() / p[x == 1].sum())
                       - np.log(w[x == 0].sum() / p[x == 0].sum()) - np.log(HR[i]))
        return out
    print("snapshot, d1 adjusted for d2:", fl(np.exp(fsolve(F_strat, np.log(HR), xtol=1e-13))))

    # Simulated raw impacts (true impacts 2, 4, 6; interactions 1 and 1.5).
    truth = np.array([2., 4, 6])
    y = cells @ truth
    print("simulated crude raw", fl(crude(cells, p, y)))
    yi = y + cells[:, 0] * cells[:, 1] + 1.5 * cells[:, 1] * cells[:, 2]
    print("simulated crude raw with interactions", fl(crude(cells, p, yi)))

    # Legacy conversions (Rasmussen et al. 2022, Table 6; 2024 overall odds).
    c = 0.27
    for prev, hr in [(.03, 3.83), (.30, 3.40)]:
        p11 = or_to_joint(c, prev, hr)
        print("legacy excess HR %g, P %g:" % (hr, prev), float(p11 / prev - (c - p11) / (1 - prev)))
    print("overall-odds excess HR 2.75:", 2.75 * .27 / (2.75 * .27 + .73) - .27)
