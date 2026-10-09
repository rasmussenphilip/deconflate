"""Prototype for v0.4 items 3-4, part 1 (five-disease example, no three-way term
unless stated). Run from this directory: python3 prototype_v04_part1.py
1a. Several impact estimates per disease: weighted least squares (GLS),
    standardised residuals and the heterogeneity statistic.
1b. Several association estimates per pair: pooling on the log-OR scale.
1c. An adjusted odds ratio as a conditional association in the global method
    (logistic projection), compared with using it as a marginal OR.
"""
import itertools
import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, norm
from python_reference import or_to_joint

ids = ["LAM", "MAS", "MET", "SCK", "RP"]
n = 5
ix = {d: i for i, d in enumerate(ids)}
P = np.array([0.25, 1 - np.exp(-0.30), 0.10, 0.35, 0.06])
assoc = [("LAM", "MAS", "OR", 1.8), ("LAM", "SCK", "OR", 2.0), ("MET", "RP", "OR", 6.0),
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


def joint_matrix(P, assoc):
    J = np.diag(P).astype(float)
    for a, b, meas, v in assoc:
        i, k = ix[a], ix[b]
        J[i, k] = J[k, i] = to_joint(meas, v, P[i], P[k])
    return J


def cov(J, P):
    S = J - np.outer(P, P)
    np.fill_diagonal(S, P * (1 - P))
    return S


def conf_row(S, i, adj):
    """Row of the conflation matrix for an estimate of disease i adjusted for adj."""
    X = [i] + adj
    others = [k for k in range(n) if k not in X]
    row = np.zeros(n)
    row[i] = 1
    if others:
        row[others] = np.linalg.solve(S[np.ix_(X, X)], S[np.ix_(X, others)])[0]
    return row


def pert_sd(a, m, b):
    mu = (a + 4 * m + b) / 6
    return np.sqrt((mu - a) * (b - mu) / 7)


def lognormal_sd(ml, sl):
    return np.sqrt(np.exp(sl ** 2) - 1) * np.exp(ml + sl ** 2 / 2)


J = joint_matrix(P, assoc)
S = cov(J, P)

# ---------------------------------------------------------------- 1a
print("=== 1a. Calving interval with two estimates for MAS ===")
rows = [  # disease, value, sd, adjusted_for (None = crude), label
    ("LAM", 12.0, 3.0, None, "LAM crude"),
    ("MAS", 9.0, 3.0, None, "MAS crude (new)"),
    ("MAS", 6.0, 2.0, ["SCK"], "MAS adj. SCK"),
    ("MET", 18.0, pert_sd(8, 18, 30), None, "MET crude"),
    ("SCK", 4.0, 1.0, "all", "SCK adj. all"),
    ("RP", 10.0, lognormal_sd(np.log(10), 0.25), None, "RP crude")]
A = np.array([conf_row(S, ix[d], ([k for k in range(n) if k != ix[d]] if adj == "all"
                                   else [ix[a] for a in adj] if adj else []))
              for d, _, _, adj, _ in rows])
y = np.array([r[1] for r in rows])
s = np.array([r[2] for r in rows])
W = np.diag(1 / s ** 2)
V = np.linalg.inv(A.T @ W @ A)            # covariance of b-hat
b = V @ A.T @ W @ y
res = y - A @ b
H = A @ V @ A.T @ W                       # hat matrix (weighted)
z = res / (s * np.sqrt(1 - np.diag(H)))   # standardised residuals
Q = np.sum((res / s) ** 2)
df = len(y) - n
print("sd used as weights:", np.round(s, 4))
print("GLS adjusted impacts:", dict(zip(ids, np.round(b, 4))), "se", np.round(np.sqrt(np.diag(V)), 4))
print("aggregate:", round(P @ b, 4), "se", round(np.sqrt(P @ V @ P), 4))
for r, zz, rr in zip(rows, z, res):
    print(f"  {r[4]:16s} raw {r[1]:6.2f}  fitted {r[1] - rr:7.3f}  std. residual {zz:6.3f}")
print(f"Q = {Q:.4f} on {df} df, p = {chi2.sf(Q, df):.4f}")
# Each MAS estimate alone (exact solve) for comparison
for keep in (1, 2):
    sel = [0, keep, 3, 4, 5]
    bb = np.linalg.solve(A[sel], y[sel])
    print(f"only '{rows[keep][4]}':", dict(zip(ids, np.round(bb, 4))), "aggregate", round(P @ bb, 4))
# What the MAS crude estimate implies: the MAS-SCK OR at which both MAS estimates agree
def gap_at(o):
    a2 = [x if (x[0], x[1]) != ("MAS", "SCK") else ("MAS", "SCK", "OR", o) for x in assoc]
    S2 = cov(joint_matrix(P, a2), P)
    A2 = np.array([conf_row(S2, ix[d], ([k for k in range(n) if k != ix[d]] if adj == "all"
                                        else [ix[a] for a in adj] if adj else []))
                   for d, _, _, adj, _ in rows])
    sel = [0, 2, 3, 4, 5]
    b2 = np.linalg.solve(A2[sel], y[sel])
    return (A2[1] @ b2) - 9.0
p_ms = J[ix["MAS"], ix["SCK"]]
or_ms = p_ms * (1 - P[1] - P[3] + p_ms) / ((P[1] - p_ms) * (P[3] - p_ms))
print(f"MAS-SCK odds ratio implied by cond_prob 0.35: {or_ms:.4f}")
print(f"MAS-SCK odds ratio at which both MAS estimates agree exactly: {brentq(gap_at, 1.0, 50):.4f}")

# ---------------------------------------------------------------- 1b
print("\n=== 1b. Pooling two LAM-MAS odds ratios ===")
z975 = norm.ppf(0.975)
est = [(np.log(1.8), (np.log(2.5) - np.log(1.3)) / (2 * z975), "OR 1.8 (95% CI 1.3-2.5)"),
       (np.log(60 * 650 / (140 * 150)), np.sqrt(1 / 60 + 1 / 140 + 1 / 150 + 1 / 650),
        "table 60/140/150/650")]
yy = np.array([e[0] for e in est]); ss = np.array([e[1] for e in est])
w = 1 / ss ** 2
fe = np.sum(w * yy) / np.sum(w)
Qa = np.sum(w * (yy - fe) ** 2)
tau2 = max(0, (Qa - (len(yy) - 1)) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
wr = 1 / (ss ** 2 + tau2)
re = np.sum(wr * yy) / np.sum(wr)
I2 = max(0, (Qa - 1) / Qa) if Qa > 0 else 0
for e in est:
    print(f"  {e[2]:26s} log OR {e[0]:.4f} se {e[1]:.4f}")
print(f"fixed effect OR {np.exp(fe):.4f} (se log {np.sqrt(1 / w.sum()):.4f}); Q {Qa:.4f}; tau2 {tau2:.5f}; I2 {I2:.3f}")
print(f"random effects OR {np.exp(re):.4f} (se log {np.sqrt(1 / wr.sum()):.4f})")

# ---------------------------------------------------------------- 1c
print("\n=== 1c. LAM-SCK odds ratio 2.0 entered as adjusted (conditional) ===")
cells = np.array(list(itertools.product([0, 1], repeat=n)), float)


def ipf(P, J, pairs, fixed=None, iters=4000, tol=1e-13):
    """Max-entropy joint matching P and the listed pairs. fixed = (i, k, lam):
    a pairwise log-linear term held fixed (its pair is not fitted)."""
    p = np.prod(np.where(cells == 1, P, 1 - P), axis=1)
    if fixed:
        i, k, lam = fixed
        p = p * np.exp(lam * cells[:, i] * cells[:, k]); p /= p.sum()
    for _ in range(iters):
        for i in range(n):
            m1 = cells[:, i] == 1
            s1 = p[m1].sum()
            p[m1] *= P[i] / s1; p[~m1] *= (1 - P[i]) / (1 - s1)
        for i, k in pairs:
            t = {(1, 1): J[i, k], (1, 0): P[i] - J[i, k], (0, 1): P[k] - J[i, k],
                 (0, 0): 1 - P[i] - P[k] + J[i, k]}
            for (u, v), tv in t.items():
                m = (cells[:, i] == u) & (cells[:, k] == v)
                p[m] *= tv / p[m].sum()
        Jg = (cells * p[:, None]).T @ cells
        err = max(abs(np.diag(Jg) - P).max(), max(abs(Jg[i, k] - J[i, k]) for i, k in pairs))
        if err < tol:
            break
    return p


def logistic_coef(p, j, X):
    """Population logistic projection of D_j on [1, D_X]: coefficients (IRLS)."""
    Z = np.column_stack([np.ones(len(cells)), cells[:, X]])
    yv = cells[:, j]
    beta = np.zeros(Z.shape[1])
    for _ in range(100):
        mu = 1 / (1 + np.exp(-Z @ beta))
        g = Z.T @ (p * (yv - mu))
        Hm = (Z * (p * mu * (1 - mu))[:, None]).T @ Z
        step = np.linalg.solve(Hm, g)
        beta += step
        if np.abs(step).max() < 1e-13:
            break
    return beta


iL, iS = ix["LAM"], ix["SCK"]
pairs_all = [(ix[a], ix[b]) for a, b, _, _ in assoc]
pairs_wo = [pr for pr in pairs_all if set(pr) != {iL, iS}]


def fit_adjusted(target_logor, adj):
    X = [iL] + adj
    def f(lam):
        p = ipf(P, J, pairs_wo, fixed=(iL, iS, lam))
        return logistic_coef(p, iS, X)[1] - target_logor
    lam = brentq(f, -3, 5, xtol=1e-12)
    return lam, ipf(P, J, pairs_wo, fixed=(iL, iS, lam))


def global_yield(p):
    Jg = (cells * p[:, None]).T @ cells
    Sg = cov(Jg, P)
    Ag = np.array([conf_row(Sg, i, []) for i in range(n)])
    return Jg, np.linalg.solve(Ag, np.array([4.8, 3.3, 5.6, 2.5, 4.2]))


def marg_or(Jg, i, k):
    p11 = Jg[i, k]
    return p11 * (1 - P[i] - P[k] + p11) / ((P[i] - p11) * (P[k] - p11))


p_marg = ipf(P, J, pairs_all)
Jm, bm = global_yield(p_marg)
print("as marginal OR 2.0 : yield", dict(zip(ids, np.round(bm, 4))), "aggregate", round(P @ bm, 4))
for lab, adj in (("adjusted for all", [k for k in range(n) if k not in (iL, iS)]),
                 ("adjusted for MAS", [ix["MAS"]])):
    lam, p_a = fit_adjusted(np.log(2.0), adj)
    Ja, ba = global_yield(p_a)
    print(f"{lab:17s}: pairwise term exp(lambda) {np.exp(lam):.4f}; implied marginal OR {marg_or(Ja, iL, iS):.4f}")
    print("                    yield", dict(zip(ids, np.round(ba, 4))), "aggregate", round(P @ ba, 4))

# ---------------------------------------------------------------- 1d
print("\n=== 1d. The 1a example under the three settings of the toggle ===")
from scipy.optimize import minimize
from scipy.stats import truncnorm
rs = np.random.default_rng(11)
M = 400000
def pert_draw(a, m, b, size):
    al = 1 + 4 * (m - a) / (b - a); be = 1 + 4 * (b - m) / (b - a)
    return a + (b - a) * rs.beta(al, be, size)
# input SD of each association, on the scale it is entered (log scale for OR-type rows)
sd_assoc = {("LAM", "MAS"): (np.log(2.5) - np.log(1.3)) / (2 * z975),            # lognormal_ci, log scale
            ("LAM", "SCK"): np.log(truncnorm.rvs(-2.0 / 0.3, np.inf, loc=2.0, scale=0.3, size=M, random_state=3)).std(),
            ("MET", "RP"): 0.3,
            ("MET", "SCK"): pert_draw(1.1, 1.5, 2.2, M).std(),                     # RR, natural scale
            ("MAS", "SCK"): rs.beta(35, 65, M).std(),                                # cond_prob
            ("MAS", "MET"): 0.02,                                                    # RD
            ("LAM", "MET"): (0.08 - 0.02) / np.sqrt(12),                             # phi
            ("RP", "SCK"): np.sqrt(1 / 30 + 1 / 40 + 1 / 320 + 1 / 610)}              # table: Woolf, log scale
logscale = {("LAM", "MAS"), ("LAM", "SCK"), ("MET", "RP"), ("RP", "SCK")}
movable = [q for q, a in enumerate(assoc) if (a[0], a[1]) in sd_assoc]
def assoc_at(u):
    out = list(assoc)
    for q, uq in zip(movable, u):
        a, b_, meas, v = assoc[q]
        sd = sd_assoc[(a, b_)]
        out[q] = (a, b_, meas, np.exp(np.log(v) + sd * uq) if (a, b_) in logscale else v + sd * uq)
    return out
def A_at(u):
    S2 = cov(joint_matrix(P, assoc_at(u)), P)
    return np.array([conf_row(S2, ix[d], ([k for k in range(n) if k != ix[d]] if adj == "all"
                                          else [ix[a] for a in adj] if adj else []))
                     for d, _, _, adj, _ in rows])
nu = len(movable)
def report(lab, u, b):
    Au = A_at(u)
    zr = (y - Au @ b) / s
    shifts = ', '.join(f"{assoc[q][0]}-{assoc[q][1]} {uq:+.2f}" for q, uq in zip(movable, u) if abs(uq) >= 0.05)
    print(f"  {lab:12s} b = {np.round(b, 3)}  aggregate {P @ b:.4f}")
    print(f"  {'':12s} impact residuals (input SDs): {np.round(zr, 3)}; association shifts (input SDs): {shifts or 'none'}")
# associations: ORs fixed, GLS (as 1a)
report("associations", np.zeros(nu), b)
# balanced: minimise impact misfit + association shifts, both in input SDs
def obj(u):
    Au = A_at(u)
    Vu = np.linalg.inv(Au.T @ W @ Au)
    bu = Vu @ Au.T @ W @ y
    return np.sum(((y - Au @ bu) / s) ** 2) + u @ u, bu
r_bal = minimize(lambda u: obj(u)[0], np.zeros(nu), method="BFGS", options={"gtol": 1e-10})
report("balanced", r_bal.x, obj(r_bal.x)[1])
# impacts: the limit of 'balanced' as the impact SDs shrink (here by a factor 1000):
# the smallest association shift (in input SDs) that fits every impact row.
K = 1e6
def obj_imp(u):
    Au = A_at(u)
    Vu = np.linalg.inv(Au.T @ W @ Au)
    bu = Vu @ Au.T @ W @ y
    return K * np.sum(((y - Au @ bu) / s) ** 2) + u @ u, bu
r_imp = minimize(lambda u: obj_imp(u)[0], r_bal.x, method="BFGS", options={"gtol": 1e-8})
for _ in range(3):
    r_imp = minimize(lambda u: obj_imp(u)[0], r_imp.x, method="Nelder-Mead" if not r_imp.success else "BFGS",
                     options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 20000} if not r_imp.success else {"gtol": 1e-8})
report("impacts", r_imp.x, obj_imp(r_imp.x)[1])
print("  (impacts: converged", r_imp.success, ")")
