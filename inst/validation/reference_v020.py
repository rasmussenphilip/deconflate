"""Reference values for deconflate 0.2.0 (estimands, interactions,
feasibility, three-way terms). Run from this directory:
    python3 reference_v020.py
Requires numpy. Values printed here are hard-coded in
tests/testthat/test-estimands.R, test-interactions.R and test-feasibility.R.
"""
import itertools
import numpy as np
from python_reference import or_to_joint, excess_matrix, ipf

np.set_printoptions(precision=12, suppress=True, linewidth=150)

# Supplementary File population: P = (0.10, 0.15, 0.20), OR d1:d2 = 2,
# d2:d3 = 3, d1:d3 independent. True additive impacts 0.02, 0.04, 0.06.
P = np.array([.10, .15, .20])
OR = np.ones((3, 3)); OR[0, 1] = OR[1, 0] = 2; OR[1, 2] = OR[2, 1] = 3
b = np.array([.02, .04, .06])
n = 3

# Pairwise covariance matrix of the disease indicators.
J = np.diag(P).astype(float)
for i, k in itertools.combinations(range(n), 2):
    J[i, k] = J[k, i] = or_to_joint(P[i], P[k], OR[i, k])
Sigma = J - np.outer(P, P)
np.fill_diagonal(Sigma, P * (1 - P))
E = excess_matrix(P, OR)          # E[k, i]
A_crude = np.eye(n) + E.T         # row i: crude estimate of disease i


def projection_row(i, S):
    X = [i] + list(S)
    coef = np.linalg.solve(Sigma[np.ix_(X, X)], Sigma[np.ix_(X, range(n))])
    row = coef[0].copy()
    row[list(S)] = 0
    return row


# 1. Bob's counterexample: d1 adjusted for d2 (additive regression), d2 and
#    d3 crude.
A = A_crude.copy()
A[0] = projection_row(0, [1])
raw = A @ b
print("projection raw (d1 adj. for d2; d2, d3 crude):", repr(raw))
print("projection recovers:", np.linalg.solve(A, raw))
# The shortcut: crude row of d1 with the d2 entry set to zero.
A_short = A_crude.copy(); A_short[0, 1] = 0
print("shortcut (zero A[1,2]) gives:", np.linalg.solve(A_short, raw))

# Independent check of the raw value by weighted least squares on the
# maximum-entropy joint distribution (as simulate_raw_impacts() does).
cells, p = ipf(P, OR)
Y = cells @ b
X = np.column_stack([np.ones(len(p)), cells[:, 0], cells[:, 1]])
W = np.diag(p)
beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ Y)
print("WLS coefficient of d1 adjusted for d2:", repr(beta[1]))

# 2. d1 adjusted for all (d2, d3) and d3 adjusted for d1: raw values.
A2 = A_crude.copy()
A2[0] = projection_row(0, [1, 2])
A2[2] = projection_row(2, [0])
print("mixed raw (d1 adj. all, d2 crude, d3 adj. d1):", repr(A2 @ b))

# 3. Interactions: delta(d1, d2) = 0.01, delta(d2, d3) = 0.015 (crude raw).
inter = Y + .01 * cells[:, 0] * cells[:, 1] + .015 * cells[:, 1] * cells[:, 2]


def crude(f):
    return np.array([(p[cells[:, i] == 1] @ f[cells[:, i] == 1]) / p[cells[:, i] == 1].sum()
                     - (p[cells[:, i] == 0] @ f[cells[:, i] == 0]) / p[cells[:, i] == 0].sum()
                     for i in range(n)])


print("raw crude with interactions:", repr(crude(inter)))
pjk = [J[0, 1], J[1, 2]]
print("aggregate with interactions:", repr(P @ b + .01 * pjk[0] + .015 * pjk[1]))
print("contributions (main + half of each interaction):",
      repr(P * b + np.array([.005 * pjk[0], .005 * pjk[0] + .0075 * pjk[1], .0075 * pjk[1]])))

# 4. Bob's infeasible triple: p = 0.5, ORs 20, 20 and 0.05.
P3 = np.array([.5, .5, .5])
OR3 = np.ones((3, 3))
OR3[0, 1] = OR3[1, 0] = 20; OR3[0, 2] = OR3[2, 0] = 20; OR3[1, 2] = OR3[2, 1] = .05
A3 = np.eye(3) + excess_matrix(P3, OR3).T
print("infeasible triple: cond(A) =", repr(np.linalg.cond(A3)))
p11 = {(i, k): or_to_joint(.5, .5, OR3[i, k]) for i, k in itertools.combinations(range(3), 2)}
# Necessary condition: p11(0,1) + p11(0,2) - p11(1,2) <= p0.
print("triple screen: p12 + p13 - p23 =", p11[(0, 1)] + p11[(0, 2)] - p11[(1, 2)], "> p1 = 0.5")

# 5. Scaling: results scale with the units of the impacts.
print("simultaneous, impacts in percent (2.5, 5, 7.5):", repr(np.linalg.solve(A_crude, [2.5, 5, 7.5])))
print("simultaneous, impacts in kg (250, 500, 750):", repr(np.linalg.solve(A_crude, [250, 500, 750])))
