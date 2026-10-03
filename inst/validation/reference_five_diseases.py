"""Checks and reference values for inst/extdata/five_diseases.
Run from this directory: python3 reference_five_diseases.py  (numpy, scipy)
* the pairwise tables of every association measure, joint feasibility (triple
  screen and LP) and the IPF fit with the three-way term;
* adjusted impacts for the yield analysis (simultaneous and published), the
  calving-interval analysis (simultaneous, with adjusted_linear estimands) and
  the welfare analysis (global, with interactions);
* the share of Monte Carlo draws that give valid, jointly feasible inputs.
"""
import itertools
import numpy as np
from scipy.optimize import linprog
from python_reference import or_to_joint

ids = ["LAM", "MAS", "MET", "SCK", "RP"]
n = 5
ix = {d: i for i, d in enumerate(ids)}
value = {"LAM": 0.25, "MAS": 0.30, "MET": 0.10, "SCK": 0.35, "RP": 0.06}
typ = {"LAM": "prevalence", "MAS": "incidence_rate", "MET": "prevalence",
       "SCK": "probability", "RP": "prevalence"}


def prob(d, v):
    return 1 - np.exp(-v) if typ[d] == "incidence_rate" else v


assoc = [  # disease1, disease2, measure, value
    ("LAM", "MAS", "OR", 1.8), ("LAM", "SCK", "OR", 2.0), ("MET", "RP", "OR", 6.0),
    ("MET", "SCK", "RR", 1.5), ("MAS", "SCK", "cond_prob", 0.35), ("MAS", "MET", "RD", 0.08),
    ("LAM", "MET", "phi", 0.05), ("RP", "SCK", "table", 30 * 610 / (40 * 320)),
    ("LAM", "RP", "independent", 1.0), ("MAS", "RP", "OR", 1.5)]


def to_joint(meas, v, p1, p2):
    if meas == "independent":
        return p1 * p2
    if meas in ("OR", "table"):
        return or_to_joint(p1, p2, v)
    if meas == "RR":
        return v * p1 / (1 - p2 + v * p2) * p2
    if meas == "RD":
        return (v * (1 - p2) + p1) * p2
    if meas == "cond_prob":
        return v * p2
    if meas == "phi":
        return p1 * p2 + v * np.sqrt(p1 * (1 - p1) * p2 * (1 - p2))
    raise ValueError(meas)


def joint_matrix(P, avals):
    J = np.diag(P).astype(float)
    for (a, b, meas, _), v in zip(assoc, avals):
        i, k = ix[a], ix[b]
        p11 = to_joint(meas, v, P[i], P[k])
        if not (max(0, P[i] + P[k] - 1) - 1e-12 <= p11 <= min(P[i], P[k]) + 1e-12):
            return None
        J[i, k] = J[k, i] = p11
    return J


def triple_ok(J, P):
    for a, b, c in itertools.combinations(range(n), 3):
        lo = max(0, J[a, b] + J[a, c] - P[a], J[a, b] + J[b, c] - P[b], J[a, c] + J[b, c] - P[c])
        hi = min(J[a, b], J[a, c], J[b, c], 1 - P[a] - P[b] - P[c] + J[a, b] + J[a, c] + J[b, c])
        if lo > hi + 1e-9:
            return False
    return True


cells = np.array(list(itertools.product([0, 1], repeat=n)), float)


def lp_feasible(J, P):
    rows, rhs = [np.ones(len(cells))], [1.0]
    for i in range(n):
        rows.append(cells[:, i]); rhs.append(P[i])
    for i, k in itertools.combinations(range(n), 2):
        rows.append(cells[:, i] * cells[:, k]); rhs.append(J[i, k])
    r = linprog(np.zeros(len(cells)), A_eq=np.array(rows), b_eq=rhs, bounds=(0, None), method="highs")
    return r.status == 0


def ipf(P, J, three=(("LAM", "MAS", "SCK"), 1.5), iters=3000):
    p = np.prod(np.where(cells == 1, P, 1 - P), axis=1)
    (a, b, c), r = three
    p = p * np.exp(np.log(r) * cells[:, ix[a]] * cells[:, ix[b]] * cells[:, ix[c]])
    p /= p.sum()
    for _ in range(iters):
        for i, k in itertools.combinations(range(n), 2):
            t = {(1, 1): J[i, k], (1, 0): P[i] - J[i, k], (0, 1): P[k] - J[i, k],
                 (0, 0): 1 - P[i] - P[k] + J[i, k]}
            for (u, v), tv in t.items():
                m = (cells[:, i] == u) & (cells[:, k] == v)
                p[m] *= tv / p[m].sum()
    return p


def cov(J, P):
    S = J - np.outer(P, P)
    np.fill_diagonal(S, P * (1 - P))
    return S


def conflation(S, adj):
    A = np.eye(n)
    for i in range(n):
        Sset = adj.get(i, [])
        X = [i] + Sset
        others = [k for k in range(n) if k not in X]
        g = np.linalg.solve(S[np.ix_(X, X)], S[np.ix_(X, others)])
        A[i, others] = g[0]
        A[i, Sset] = 0
    return A


P = np.array([prob(d, value[d]) for d in ids])
J = joint_matrix(P, [a[3] for a in assoc])
print("probabilities:", dict(zip(ids, np.round(P, 6))))
print("pairwise tables valid:", J is not None, "| triple screen:", triple_ok(J, P), "| LP feasible:", lp_feasible(J, P))
S = cov(J, P)
A = conflation(S, {})
yield_raw = np.array([4.8, 3.3, 5.6, 2.5, 4.2])
b = np.linalg.solve(A, yield_raw)
pub = yield_raw ** 2 / (yield_raw + (A - np.eye(n)) @ yield_raw)
print("yield simultaneous:", repr(b), "total", repr(P @ b))
print("yield published:   ", repr(pub), "total", repr(P @ pub), "cond(A)", np.linalg.cond(A))
ci_raw = np.array([12, 6, 18, 4, 10])
Aci = conflation(S, {ix["MAS"]: [ix["SCK"]], ix["SCK"]: [i for i in range(n) if i != ix["SCK"]]})
bci = np.linalg.solve(Aci, ci_raw)
print("calving interval simultaneous (MAS adj. SCK, SCK adj. all):", repr(bci), "total", repr(P @ bci))
# Welfare: interactions LAM:MAS 2, MET:RP 1.5, LAM:SCK -0.5; global with three-way 1.5
p = ipf(P, J)
Jg = (cells * p[:, None]).T @ cells
print("IPF residual:", np.abs(Jg - J).max())
Sg = cov(Jg, P)
Ag = conflation(Sg, {})
w_raw = np.array([10, 6, 5, 3, 4.0])
g = 2.0 * cells[:, 0] * cells[:, 1] + 1.5 * cells[:, 2] * cells[:, 4] - 0.5 * cells[:, 0] * cells[:, 3]
covg = cells.T @ (p * g) - P * (p @ g)
off = covg / np.diag(Sg)
bw = np.linalg.solve(Ag, w_raw - off)
Dm = np.zeros((n, n)); Dm[0, 1] = Dm[1, 0] = 2; Dm[2, 4] = Dm[4, 2] = 1.5; Dm[0, 3] = Dm[3, 0] = -0.5
contrib = P * bw + 0.5 * (Dm * Jg).sum(1)
print("welfare global:", repr(bw), "total", repr(contrib.sum()))
# Thresholds used by run_all_features.R and the tests.
from scipy.optimize import brentq
def yield_adj(i, v):
    r = yield_raw.copy(); r[i] = v
    return np.linalg.solve(A, r)
print("sign threshold, raw SCK yield impact:", repr(brentq(lambda v: yield_adj(3, v)[3], 0, 2.5, xtol=1e-14)))
print("total = 3 threshold, raw MET yield impact:", repr(brentq(lambda v: P @ yield_adj(2, v) - 3, 0, 15, xtol=1e-14)))
def welfare_total(delta):
    g = delta * cells[:, 0] * cells[:, 1] + 1.5 * cells[:, 2] * cells[:, 4] - 0.5 * cells[:, 0] * cells[:, 3]
    cg = cells.T @ (p * g) - P * (p @ g)
    bb = np.linalg.solve(Ag, w_raw - cg / np.diag(Sg))
    D2 = Dm.copy(); D2[0, 1] = D2[1, 0] = delta
    return (P * bb + 0.5 * (D2 * Jg).sum(1)).sum()
print("change = -10% threshold, LAM:MAS welfare interaction:",
      repr(brentq(lambda d: welfare_total(d) / welfare_total(2.0) - 0.9, -5, 20, xtol=1e-14)))
def yield_contrib_psck(psck):
    P2 = P.copy(); P2[ix["SCK"]] = psck
    A2 = conflation(cov(joint_matrix(P2, [a[3] for a in assoc]), P2), {})
    return P2 * np.linalg.solve(A2, yield_raw)
print("rank threshold MET vs SCK (yield contributions), prob:SCK:",
      repr(brentq(lambda v: np.diff(yield_contrib_psck(v)[[2, 3]])[0], 0.1, 0.35, xtol=1e-14)))
print("published pole of SCK (raw yield impact):", repr(-((A - np.eye(n)) @ yield_raw)[ix["SCK"]]))
# Monte Carlo: share of valid and feasible draws.
rng = np.random.default_rng(1)
def pert(mn, mode, mx, size):
    a = 1 + 4 * (mode - mn) / (mx - mn); b_ = 1 + 4 * (mx - mode) / (mx - mn)
    return mn + (mx - mn) * rng.beta(a, b_, size)
N = 20000
draw_p = {"LAM": rng.beta(78.29, 227.42, N), "MAS": pert(0.20, 0.30, 0.45, N),
          "MET": rng.uniform(0.07, 0.13, N), "SCK": pert(0.20, 0.35, 0.50, N), "RP": np.full(N, 0.06)}
z = 1.959963984540054
from scipy.stats import truncnorm
draw_a = [np.exp(rng.normal(np.log(1.8), (np.log(2.5) - np.log(1.3)) / (2 * z), N)),
          truncnorm.rvs((0 - 2.0) / 0.3, np.inf, loc=2.0, scale=0.3, size=N, random_state=2),
          np.exp(rng.normal(np.log(6), 0.3, N)),
          pert(1.1, 1.5, 2.2, N),
          rng.beta(35, 65, N),
          rng.normal(0.08, 0.02, N),
          rng.uniform(0.02, 0.08, N),
          np.full(N, 30 * 610 / (40 * 320)), np.ones(N), np.full(N, 1.5)]
ok = 0
for r in range(N):
    Pr = np.array([prob(d, draw_p[d][r]) for d in ids])
    Jr = joint_matrix(Pr, [x[r] for x in draw_a])
    if Jr is not None and triple_ok(Jr, Pr):
        ok += 1
print("Monte Carlo: share of draws with valid, triple-feasible inputs:", ok / N)
