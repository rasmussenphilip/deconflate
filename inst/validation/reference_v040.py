"""Reference values for deconflate 0.4 (item 3: unknown pairs by default,
event impacts through deconflate()). Run from this directory:
python3 reference_v040.py   (numpy, scipy)

1. 2024 global dairy (inputs = "analysis"): pairs not in Table 3 unknown, so
   the global method: adjusted yield and fertility impacts and aggregates.
2. 2024 culling hazard ratios as event impacts (snapshot_crude), unknown
   pairs, overall risk 0.2366: adjusted hazard ratios and attributable risk.
3. UK 2022 (example_uk_dairy_2022()): yield with unknown pairs (global).
4. Five-disease example: culling.csv (event impacts with mixed measures).
5. Supplement example (all pairs given): risk-based rows in the event model
   reproduce the hazard-ratio solution.
"""
import itertools
import numpy as np
from scipy.optimize import brentq, fsolve
from python_reference import or_to_joint

np.set_printoptions(precision=10, linewidth=160)


def cells_of(n):
    return np.array(list(itertools.product([0, 1], repeat=n)), float)[:, ::-1]


def ipf(P, pairs, three=(), tol=1e-13, iters=20000):
    """Max-ent joint matching P and pairs [(i, k, p11)]; three: [(i, j, k, ratio)]."""
    n = len(P)
    cells = cells_of(n)
    logp = cells @ np.log(P) + (1 - cells) @ np.log(1 - P)
    for i, j, k, r in three:
        logp = logp + np.log(r) * cells[:, i] * cells[:, j] * cells[:, k]
    p = np.exp(logp - logp.max()); p /= p.sum()
    groups = [(2 * cells[:, i] + cells[:, k]).astype(int) for i, k, _ in pairs]
    targets = [np.array([1 - P[i] - P[k] + q, P[k] - q, P[i] - q, q]) for i, k, q in pairs]
    for it in range(iters):
        for g, t in zip(groups, targets):
            cur = np.bincount(g, weights=p, minlength=4)
            p = p * (t / cur)[g]
        for i in range(n):
            m = cells[:, i] == 1
            cur = p[m].sum()
            p = np.where(m, p * P[i] / cur, p * (1 - P[i]) / (1 - cur))
        res = max(max(np.abs(np.bincount(g, weights=p, minlength=4) - t).max() for g, t in zip(groups, targets)) if pairs else 0,
                  np.abs(cells.T @ p - P).max())
        if res < tol:
            break
    return cells, p


def crude_A(cells, p):
    P = cells.T @ p
    J = (cells * p[:, None]).T @ cells
    S = J - np.outer(P, P)
    np.fill_diagonal(S, P * (1 - P))
    return S / np.diag(S)[:, None], J   # A[i, k] = cov(i, k) / var(i)


def or_from(J, i, k):
    p11 = J[i, k]; p1 = J[i, i]; p2 = J[k, k]
    return p11 * (1 - p1 - p2 + p11) / ((p1 - p11) * (p2 - p11))


# ----------------------------------------------------------------- 2024 inputs
D = ['CK', 'CM', 'DA', 'DYS', 'LAM', 'MET', 'MF', 'OC', 'PTB', 'RP', 'SCK', 'SCM']
inc = {'CK': 0.0306287320599544, 'CM': 0.305081708702075, 'DA': 0.0216266798689032,
       'DYS': 0.0611107185570655, 'LAM': 0.254288506511192, 'MET': 0.0958013664473386,
       'MF': 0.0241829408262443, 'OC': 0.112683551256031, 'PTB': 0.103964073833149,
       'RP': 0.12346845541711, 'SCK': 0.479052187708235, 'SCM': 0.409414831921559}
P24 = np.array([inc[d] if d in ('PTB', 'SCM') else 1 - np.exp(-inc[d]) for d in D])
OR24 = {"CK:CM": 2.13, "CK:LAM": 1.65, "CK:MF": 1.60, "CK:OC": 1.97, "CK:RP": 1.55, "CK:SCK": 6.95,
        "CK:SCM": 2.40, "CM:LAM": 2.10, "CM:PTB": 1.89, "CM:RP": 2.70, "CM:SCK": 1.64, "CM:SCM": 3.05,
        "DA:CM": 3.45, "DA:MF": 2.50, "DA:RP": 3.50, "DA:SCK": 3.87, "DA:SCM": 3.60, "DYS:LAM": 2.09,
        "DYS:OC": 0.40, "DYS:RP": 2.74, "LAM:OC": 2.63, "LAM:PTB": 2.70, "LAM:RP": 1.50, "LAM:SCK": 2.01,
        "MET:CK": 2.42, "MET:CM": 2.30, "MET:DA": 3.40, "MET:DYS": 2.95, "MET:LAM": 6.10, "MET:MF": 1.50,
        "MET:OC": 1.94, "MET:RP": 3.53, "MET:SCK": 1.94, "MF:DYS": 9.70, "MF:LAM": 3.60, "MF:RP": 2.40,
        "RP:OC": 2.18, "SCK:RP": 1.52}
YIELD = dict(CK=0.4321944, CM=3.2499, DA=2.83693, DYS=4.919088, LAM=4.8061, MET=5.613085,
             MF=0.5365130859025376, OC=3.747839, PTB=4.3, RP=4.198664, SCK=8.396472, SCM=6.293184)
FERT = dict(CK=1.445122, CM=8.42, DA=1.082641, DYS=2.399177, LAM=3.304898, MET=14.67308,
            MF=2.414949, OC=9.685465, PTB=5.349124, RP=6.760971, SCK=1.122037, SCM=0.2636645)
HR24 = dict(CK=1.5001, CM=2.3, DA=2.851179, DYS=1.258143, LAM=1.744976, MET=1.116444,
            MF=2.999886, OC=1.62, PTB=2.310508, RP=1.599928, SCK=1.92, SCM=1.449996)

pairs24 = []
for key, o in OR24.items():
    a, b = key.split(":")
    i, k = D.index(a), D.index(b)
    pairs24.append((i, k, or_to_joint(P24[i], P24[k], o)))
cells24, p24 = ipf(P24, pairs24)
A24, J24 = crude_A(cells24, p24)
print("=== 1. 2024 inputs, unknown pairs (global) ===")
for lab, imp in (("yield", YIELD), ("fertility", FERT)):
    raw = np.array([imp[d] for d in D])
    b = np.linalg.solve(A24, raw)
    print(f"{lab}: adjusted", dict(zip(D, np.round(b, 6))))
    print(f"{lab}: raw sum {P24 @ raw:.10f}  adjusted total {P24 @ b:.10f}")
listed = {tuple(sorted((D.index(k.split(':')[0]), D.index(k.split(':')[1])))) for k in OR24}
unk = [(i, k) for i, k in itertools.combinations(range(12), 2) if (i, k) not in listed]
fo = [or_from(J24, i, k) for i, k in unk]
print(f"unknown pairs: {len(unk)}; fitted odds ratios from {min(fo):.6f} to {max(fo):.6f}")
print("fitted OR CK:DA", or_from(J24, 0, 2))


# --------------------------------------------------------------- event model
def strata_of(cells, i, adj):
    if not adj:
        return np.zeros(len(cells), int)
    return (cells[:, adj] @ (2 ** np.arange(len(adj)))).astype(int)


def measure_value(cells, p, i, adj, kind, beta, h0):
    st = strata_of(cells, i, adj)
    x = cells[:, i]
    mult = np.exp(cells @ beta)
    val = mult if kind in ("HR", "rate_ratio") else 1 - np.exp(-h0 * mult)
    num, w = [], []
    for s in np.unique(st):
        m = st == s
        d1, d0 = np.sum(p[m] * x[m]), np.sum(p[m] * (1 - x[m]))
        if d1 <= 0 or d0 <= 0:
            continue
        m1 = np.sum(p[m] * x[m] * val[m]) / d1
        m0 = np.sum(p[m] * (1 - x[m]) * val[m]) / d0
        if kind in ("HR", "rate_ratio", "RR"):
            e = np.log(m1 / m0)
        elif kind == "OR":
            e = np.log(m1 / (1 - m1)) - np.log(m0 / (1 - m0))
        else:
            e = m1 - m0
        num.append(e); w.append(d1 * d0 / (d1 + d0))
    return np.average(num, weights=w)


def h0_for(cells, p, beta, r):
    mult = np.exp(cells @ beta)
    return brentq(lambda h: np.sum(p * (1 - np.exp(-h * mult))) - r, 1e-14, 50, xtol=1e-16)


def solve_event(cells, p, rows, r):
    """rows: (disease index, adj list, measure, value)."""
    n = cells.shape[1]
    risk = any(k in ("RR", "OR", "RD") for _, _, k, _ in rows)
    target = np.array([v if k == "RD" else np.log(v) for _, _, k, v in rows])
    def F(beta):
        h0 = h0_for(cells, p, beta, r) if risk else 1.0
        return np.array([measure_value(cells, p, i, a, k, beta, h0) for i, a, k, _ in rows]) - target
    beta = fsolve(F, np.zeros(n), xtol=1e-13)
    return beta, np.abs(F(beta)).max()


def attributable(cells, p, beta, r):
    h0 = h0_for(cells, p, beta, r)
    return r - (1 - np.exp(-h0))


print("\n=== 2. 2024 culling, event model, unknown pairs, overall risk 0.2366 ===")
rows = [(i, [], "HR", HR24[d]) for i, d in enumerate(D)]
beta24, res = solve_event(cells24, p24, rows, 0.2366)
print("adjusted HR", dict(zip(D, np.round(np.exp(beta24), 6))), f"residual {res:.1e}")
ar24 = attributable(cells24, p24, beta24, 0.2366)
print(f"attributable risk {ar24:.10f}; fraction {ar24 / 0.2366:.10f}")

# --------------------------------------------------------------- UK 2022
print("\n=== 3. UK 2022 yield, unknown pairs (global) ===")
U = ["CO", "DA", "DYS", "FAS", "GIN", "LAM", "MAS", "MET", "MF", "NEO", "PTB", "RP", "SCK"]
prevU = dict(CO=0.09, DA=0.03, DYS=0.02, FAS=0.10, GIN=0.21, LAM=0.30, MAS=0.30, MET=0.10, MF=0.08,
             NEO=0.15, PTB=0.07, RP=0.05, SCK=0.22)
ORU = {"RP:MET": 6.20, "DA:SCK": 4.25, "RP:DYS": 4.10, "MET:DA": 3.40, "MET:DYS": 3.20, "LAM:PTB": 2.70,
       "MF:DA": 2.50, "MAS:MET": 2.30, "DA:RP": 2.20, "MAS:DA": 2.10, "SCK:MF": 2.10, "LAM:SCK": 2.01,
       "MAS:MF": 1.90, "MAS:PTB": 1.89, "MAS:CO": 1.65, "MAS:SCK": 1.64, "SCK:CO": 1.60, "MET:SCK": 1.40,
       "SCK:RP": 1.20}
yieldU = dict(CO=0, DA=4.04, DYS=4.05, FAS=7.33, GIN=3.28, LAM=5.54, MAS=4.57, MET=3.95, MF=0.41,
              NEO=4.20, PTB=5.90, RP=7.38, SCK=3.05)
PU = np.array([prevU[d] for d in U])
pairsU = []
for key, o in ORU.items():
    a, b = key.split(":")
    i, k = U.index(a), U.index(b)
    pairsU.append((i, k, or_to_joint(PU[i], PU[k], o)))
cellsU, pU = ipf(PU, pairsU)
AU, JU = crude_A(cellsU, pU)
rawU = np.array([yieldU[d] for d in U])
bU = np.linalg.solve(AU, rawU)
print("adjusted", dict(zip(U, np.round(bU, 6))))
print(f"raw sum {PU @ rawU:.10f}  adjusted total {PU @ bU:.10f}")

# --------------------------------------------------------------- five diseases
print("\n=== 4. Five-disease culling.csv (event impacts), overall risk 0.25 ===")
ids5 = ["LAM", "MAS", "MET", "SCK", "RP"]
P5 = np.array([0.25, 1 - np.exp(-0.30), 0.10, 0.35, 0.06])
def to_joint(meas, v, p1, p2):
    if meas == "OR":
        return or_to_joint(p1, p2, v)
    if meas == "RR":
        return v * p1 / (1 - p2 + v * p2) * p2
    if meas == "RD":
        return (v * (1 - p2) + p1) * p2
    if meas == "cond_prob":
        return v * p2
    if meas == "phi":
        return p1 * p2 + v * np.sqrt(p1 * (1 - p1) * p2 * (1 - p2))
assoc5 = [("LAM", "MAS", "OR", 1.8), ("LAM", "SCK", "OR", 2.0), ("MET", "RP", "OR", 6.0),
          ("MET", "SCK", "RR", 1.5), ("MAS", "SCK", "cond_prob", 0.35), ("MAS", "MET", "RD", 0.08),
          ("LAM", "MET", "phi", 0.05), ("RP", "SCK", "OR", 30 * 610 / (40 * 320)), ("LAM", "RP", "OR", 1.0),
          ("MAS", "RP", "OR", 1.5)]
pairs5 = [(ids5.index(a), ids5.index(b), to_joint(m, v, P5[ids5.index(a)], P5[ids5.index(b)]))
          for a, b, m, v in assoc5]
cells5, p5 = ipf(P5, pairs5, three=[(0, 1, 3, 1.5)])
ix = {d: i for i, d in enumerate(ids5)}
ALL = [None]
rows5 = [(0, [], "HR", 1.74), (1, [0, 3], "HR", 1.6), (2, [], "RR", 1.23),
         (3, [], "OR", 1.45), (4, [0, 1, 2, 3], "RD", 0.10)]
beta5, res5 = solve_event(cells5, p5, rows5, 0.25)
print("adjusted HR", dict(zip(ids5, np.round(np.exp(beta5), 8))), f"residual {res5:.1e}")
ar5 = attributable(cells5, p5, beta5, 0.25)
print(f"attributable risk {ar5:.10f}")
# Shapley allocation of the attributable risk (cell by cell)
def shapley(cells, p, loss):
    n = cells.shape[1]
    phi = np.zeros(n)
    from math import factorial
    for c, pc in zip(cells, p):
        pres = [i for i in range(n) if c[i] == 1]
        k = len(pres)
        for i in pres:
            others = [j for j in pres if j != i]
            for s in range(len(others) + 1):
                for S in itertools.combinations(others, s):
                    x = np.zeros(n); x[list(S)] = 1
                    w = factorial(s) * factorial(k - s - 1) / factorial(k)
                    xi = x.copy(); xi[i] = 1
                    phi[i] += pc * w * (loss(xi) - loss(x))
    return phi
h05 = h0_for(cells5, p5, beta5, 0.25)
r05 = 1 - np.exp(-h05)
sh5 = shapley(cells5, p5, lambda x: 1 - np.exp(-h05 * np.exp(x @ beta5)) - r05)
print("attributable by disease", dict(zip(ids5, np.round(sh5, 10))), "sum", sh5.sum())

# --------------------------------------------------------------- supplement
print("\n=== 5. Supplement example: event model with mixed measures (overall risk 0.25) ===")
Ps = np.array([0.10, 0.15, 0.20])
pairs_s = [(0, 1, or_to_joint(0.10, 0.15, 2.0)), (0, 2, 0.10 * 0.20), (1, 2, or_to_joint(0.15, 0.20, 3.0))]
cells_s, p_s = ipf(Ps, pairs_s)
rows_hr = [(0, [], "HR", 1.5), (1, [], "HR", 2.0), (2, [], "HR", 1.3)]
beta_s, _ = solve_event(cells_s, p_s, rows_hr, 0.25)
print("adjusted HR (HR inputs)", np.round(np.exp(beta_s), 10))
h0s = h0_for(cells_s, p_s, beta_s, 0.25)
rr2 = np.exp(measure_value(cells_s, p_s, 1, [], "RR", beta_s, h0s))
rd3 = measure_value(cells_s, p_s, 2, [], "RD", beta_s, h0s)
or1 = np.exp(measure_value(cells_s, p_s, 0, [], "OR", beta_s, h0s))
print(f"implied crude RR of d2 {rr2:.12f}; RD of d3 {rd3:.12f}; OR of d1 {or1:.12f}")
beta_m, res_m = solve_event(cells_s, p_s, [(0, [], "OR", round(or1, 4)), (1, [], "RR", round(rr2, 4)),
                                           (2, [], "RD", round(rd3, 4))], 0.25)
print(f"rounded inputs OR {round(or1, 4)}, RR {round(rr2, 4)}, RD {round(rd3, 4)}: adjusted HR",
      np.round(np.exp(beta_m), 10), f"residual {res_m:.1e}")
print(f"attributable (HR inputs) {attributable(cells_s, p_s, beta_s, 0.25):.10f}")

# ------------------------------------------------- 6. values for the item-3 tests
# (tests/testthat/test-import.R and test-supplement.R)
print("\n=== 6a. cm_template(): yield.csv (additive) and culling.csv (event, overall risk 0.25) ===")
# diseases.csv: LAM prevalence 0.25, SCK incidence rate 0.48, MET prevalence 0.10;
# associations.csv: LAM:SCK OR 2.01, MET:SCK OR 1.94, MET:LAM OR 6.10 (every pair given).
idsT = ["LAM", "SCK", "MET"]
PT = np.array([0.25, 1 - np.exp(-0.48), 0.10])
pairsT = [(0, 1, or_to_joint(PT[0], PT[1], 2.01)), (2, 1, or_to_joint(PT[2], PT[1], 1.94)),
          (2, 0, or_to_joint(PT[2], PT[0], 6.10))]
JT = np.diag(PT)
for i, k, q in pairsT:
    JT[i, k] = JT[k, i] = q
ST = JT - np.outer(PT, PT)
np.fill_diagonal(ST, PT * (1 - PT))
AT = ST / np.diag(ST)[:, None]
rawT = np.array([4.81, 8.40, 5.61])
bT = np.linalg.solve(AT, rawT)
print("yield simultaneous: adjusted", repr(bT), "raw sum", repr(PT @ rawT), "total", repr(PT @ bT))
cellsT, pTj = ipf(PT, pairsT)
rowsT = [(0, [], "HR", 1.74), (1, [], "HR", 1.92), (2, [], "RR", 1.45)]
betaT, resT = solve_event(cellsT, pTj, rowsT, 0.25)
print("culling: adjusted HR", repr(np.exp(betaT)), f"residual {resT:.1e}")
print(f"culling: attributable risk {attributable(cellsT, pTj, betaT, 0.25):.10f}")

print("\n=== 6b. Supplement population with d1:d3 unknown (d1:d2 OR 2, d2:d3 OR 3), global ===")
pairs_u = [(0, 1, or_to_joint(0.10, 0.15, 2.0)), (1, 2, or_to_joint(0.15, 0.20, 3.0))]
cells_u, p_u = ipf(Ps, pairs_u)
A_u, J_u = crude_A(cells_u, p_u)
raw_u = np.array([2.5, 5, 7.5])
b_u = np.linalg.solve(A_u, raw_u)
print("yield global: adjusted", repr(b_u), "total", repr(Ps @ b_u))
print("fitted OR d1:d3", repr(or_from(J_u, 0, 2)))
