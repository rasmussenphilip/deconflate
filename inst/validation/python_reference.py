"""Independent Python reference calculations for the deconflate test suite.

Run with: python3 python_reference.py   (requires numpy, scipy)
Values printed here are hard-coded in tests/testthat/.
"""
import itertools
import numpy as np
from scipy.optimize import brentq


def or_to_joint(p1, p2, odds):
    """P(d1 and d2) for a 2x2 table with marginals p1, p2 and odds ratio."""
    if abs(odds - 1) < 1e-12:
        return p1 * p2
    a = (odds - 1) * p2
    b = -((odds - 1) * (p1 + p2) + 1)
    c = odds * p1
    g = (-b - np.sqrt(b * b - 4 * a * c)) / (2 * a)   # P(d1 | d2), eq. 11
    return g * p2


def excess_matrix(P, OR):
    n = len(P)
    E = np.zeros((n, n))  # E[k, i] = P(k|i) - P(k|not i)
    for i in range(n):
        for k in range(n):
            if i != k:
                p11 = or_to_joint(P[i], P[k], OR[i, k])
                E[k, i] = p11 / P[i] - (P[k] - p11) / (1 - P[i])
    return E


def published(P, OR, m):
    E = excess_matrix(P, OR)
    return np.array([m[i] ** 2 / (m[i] + sum(E[k, i] * m[k] for k in range(len(P)) if k != i))
                     for i in range(len(P))])


def simultaneous(P, OR, m):
    return np.linalg.solve(np.eye(len(P)) + excess_matrix(P, OR).T, m)


def ipf(P, OR, iters=3000):
    n = len(P)
    cells = np.array(list(itertools.product([0, 1], repeat=n)))
    p = np.prod(np.where(cells == 1, P, 1 - P), axis=1)
    targets = []
    for i, k in itertools.combinations(range(n), 2):
        p11 = or_to_joint(P[i], P[k], OR[i, k])
        targets.append((i, k, {(1, 1): p11, (1, 0): P[i] - p11, (0, 1): P[k] - p11,
                               (0, 0): 1 - P[i] - P[k] + p11}))
    for _ in range(iters):
        for i, k, t in targets:
            for (a, b), v in t.items():
                msk = (cells[:, i] == a) & (cells[:, k] == b)
                p[msk] *= v / p[msk].sum()
    return cells, p


def crude(cells, p, f):
    return np.array([(p[cells[:, i] == 1] @ f[cells[:, i] == 1]) / p[cells[:, i] == 1].sum()
                     - (p[cells[:, i] == 0] @ f[cells[:, i] == 0]) / p[cells[:, i] == 0].sum()
                     for i in range(cells.shape[1])])


if __name__ == "__main__":
    np.set_printoptions(precision=8, suppress=True)
    # Supplementary File example
    P = np.array([.10, .15, .20])
    OR = np.ones((3, 3)); OR[0, 1] = OR[1, 0] = 2; OR[1, 2] = OR[2, 1] = 3
    m = np.array([.025, .05, .075])
    print("E[k,i]\n", excess_matrix(P, OR))
    for name, adj in [("published", published(P, OR, m)), ("simultaneous", simultaneous(P, OR, m))]:
        L = adj @ P; xh = 10000 / (1 - L)
        print(name, adj, "xh", repr(xh), "gaps", adj * P / L * (xh - 10000))
    # Simulated population, true impacts 2/4/6%, with and without interactions
    cells, p = ipf(P, OR)
    true = np.array([.02, .04, .06])
    add = cells @ true
    inter = add + .01 * cells[:, 0] * cells[:, 1] + .015 * cells[:, 1] * cells[:, 2]
    print("raw additive", crude(cells, p, add))
    print("raw interactions", crude(cells, p, inter))
    # Hazard ratios (2022 Table 6), culling rate 0.27
    c = 0.27
    for d, prev, hr in [("DA", .03, 3.83), ("LAM", .30, 3.40)]:
        p11 = or_to_joint(c, prev, hr)
        p0 = brentq(lambda x: (1 - prev) * x + prev * (1 - (1 - x) ** hr) - c, 1e-12, 1 - 1e-12, xtol=1e-15)
        print(d, "or_approx", repr(p11 / prev), repr((c - p11) / (1 - prev)),
              "PH", repr(1 - (1 - p0) ** hr), repr(p0))
