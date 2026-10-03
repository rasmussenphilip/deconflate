"""Reference values for deconflate 0.3.0: threshold searches, importance-
sampling support and three-way terms with unknown pairs.
Run from this directory: python3 reference_v030.py  (numpy, scipy)
Values printed here are hard-coded in tests/testthat/test-threshold.R and
test-review-v03.R.
"""
import itertools
import numpy as np
from scipy.optimize import brentq
from python_reference import excess_matrix

np.set_printoptions(precision=12, suppress=True, linewidth=150)
P = np.array([.10, .15, .20])
m0 = np.array([2.5, 5, 7.5])


def model(or12=2.0, or13=1.0, or23=3.0):
    OR = np.ones((3, 3))
    OR[0, 1] = OR[1, 0] = or12
    OR[0, 2] = OR[2, 0] = or13
    OR[1, 2] = OR[2, 1] = or23
    return np.eye(3) + excess_matrix(P, OR).T     # A[i, k] = E[k, i]


def simultaneous(A, m):
    return np.linalg.solve(A, m)


def published(A, m):
    conf = (A - np.eye(3)) @ m
    return m ** 2 / (m + conf)


A0 = model()
b0 = simultaneous(A0, m0)
total0 = P @ b0
print("baseline adjusted:", repr(b0), "total:", repr(total0))

# 1. Change: association strengths at which the total departs from the
#    baseline by 10%.
def change_roots(key, lo, hi, t):
    f = lambda o: P @ simultaneous(model(**{key: o}), m0) - (1 + t) * total0
    g = np.exp(np.linspace(np.log(lo), np.log(hi), 101))
    v = [f(x) for x in g]
    return [brentq(f, g[i], g[i + 1], xtol=1e-14) for i in range(100) if np.sign(v[i]) != np.sign(v[i + 1])]
print("change -10% via OR(d1,d3) in [1, 50]:", change_roots("or13", 1, 50, -0.1))
print("change +10% via OR(d2,d3) in [0.2, 20]:", change_roots("or23", 0.2, 20, 0.1))

# 2. Sign: raw impact of d1 at which its adjusted impact is zero (linear).
Ainv = np.linalg.inv(A0)
m1_star = -(Ainv[0, 1] * 5 + Ainv[0, 2] * 7.5) / Ainv[0, 0]
print("sign of d1 (simultaneous), impact:d1 root:", repr(m1_star))
# Published: pole at m1 = -c, c = A[0,1]*5 + A[0,2]*7.5.
c = A0[0, 1] * 5 + A0[0, 2] * 7.5
print("published pole of d1 at impact:d1 =", repr(-c))

# 3. Total: raw impact of d2 at which the total equals 3.
w = P @ Ainv
m2_star = (3 - w[0] * 2.5 - w[2] * 7.5) / w[1]
print("total = 3 via impact:d2:", repr(m2_star))

# 4. Rank by contribution as OR(d2, d3) varies (pairs d1-d2, d1-d3, d2-d3).
def contrib(o):
    return P * simultaneous(model(or23=o), m0)
grid3 = np.exp(np.linspace(np.log(0.1), np.log(100), 201))
C = np.array([contrib(o) for o in grid3])
for a, b in [(0, 1), (0, 2), (1, 2)]:
    d = C[:, a] - C[:, b]
    r = [brentq(lambda o: contrib(o)[a] - contrib(o)[b], grid3[i], grid3[i + 1], xtol=1e-14)
         for i in range(200) if np.sign(d[i]) != np.sign(d[i + 1])]
    print(f"rank d{a+1} vs d{b+1} via OR(d2,d3) in [0.1, 100]:", r)

# 5. Importance sampling with a gap in the proposal's support:
# target U(0,1), proposal 0.5 U(0,0.25) + 0.5 U(0.5,1): the self-normalised
# estimate converges to the mean over the covered part.
print("IS with a gap converges to:", (0.25 ** 2 / 2 + (1 - 0.25) / 2) / 0.75)

# 6. Three-way term with all pairs unknown (Bob): p = 0.5, crude impacts 1,2,3.
cells = np.array(list(itertools.product([0, 1], repeat=3)))
for r in [1, 4]:
    p = np.exp(np.log(r) * cells.prod(1))
    for _ in range(3000):
        for i in range(3):
            cur = p[cells[:, i] == 1].sum() / p.sum()
            p = np.where(cells[:, i] == 1, p * .5 / cur, p * .5 / (1 - cur))
            p /= p.sum()
    J = (cells * p[:, None]).T @ cells
    Pm = np.diag(J)
    S = J - np.outer(Pm, Pm)
    np.fill_diagonal(S, Pm * (1 - Pm))
    A = S / np.diag(S)[:, None]
    b = np.linalg.solve(A, [1, 2, 3])
    print(f"three-way ratio {r}: P(1 and 2) = {J[0, 1]!r}, aggregate = {Pm @ b!r}, adjusted = {b!r}")
