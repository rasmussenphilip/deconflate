"""Prototype for v0.4 item 3, part 3: several kinds of mortality estimate
mapped onto the snapshot hazard model (five-disease example: joint with the
LAM:MAS:SCK three-way ratio 1.5, as in run_all_features.R).
Run from this directory: python3 prototype_v04_part3.py

Model: animal with diseases d has period cumulative hazard h0 * exp(beta'd),
so period risk R(d) = 1 - exp(-h0 exp(beta'd)). The overall risk r fixes h0:
sum_d p(d) R(d) = r. For disease i (within strata of its adjustment set,
combined with weights d1 d0 / (d1 + d0) on the measure's own scale):
  HR, rate_ratio : log( E[exp(beta'D) | i] / E[exp(beta'D) | not i] )
  RR             : log( E[R | i] / E[R | not i] )
  OR             : log( odds(E[R | i]) / odds(E[R | not i]) )
  RD             : E[R | i] - E[R | not i]
"""
import itertools
import numpy as np
from scipy.optimize import brentq, fsolve
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    import reference_five_diseases as f5      # builds the example joint

ids, n, cells = f5.ids, f5.n, f5.cells
p = f5.ipf(f5.P, f5.J)                        # joint with the three-way ratio
ix = f5.ix
ALL = "all"


def strata(i, adj):
    S = [k for k in range(n) if k != i] if adj == ALL else [ix[a] for a in adj]
    if not S:
        return np.zeros(len(cells), int)
    return (cells[:, S] @ (2 ** np.arange(len(S)))).astype(int)


def measure(i, adj, kind, beta, h0):
    st = strata(i, adj)
    x = cells[:, i]
    mult = np.exp(cells @ beta)
    R = 1 - np.exp(-h0 * mult)
    val = mult if kind in ("HR", "rate_ratio") else R
    num, wts = [], []
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
        num.append(e); wts.append(d1 * d0 / (d1 + d0))
    return np.average(num, weights=wts)


def h0_for(beta, r):
    mult = np.exp(cells @ beta)
    return brentq(lambda h: np.sum(p * (1 - np.exp(-h * mult))) - r, 1e-12, 50, xtol=1e-15)


def solve(rows, r=None):
    """rows: (disease, adj, kind, value on natural scale). Returns beta (and h0)."""
    risk_based = any(k in ("RR", "OR", "RD") for _, _, k, _ in rows)
    target = np.array([v if k == "RD" else np.log(v) for _, _, k, v in rows])
    def F(beta):
        h0 = h0_for(beta, r) if risk_based else 1.0
        return np.array([measure(ix[d], a, k, beta, h0) for d, a, k, _ in rows]) - target
    beta = fsolve(F, np.zeros(n), xtol=1e-13)
    return beta, (h0_for(beta, r) if r is not None else None), np.abs(F(beta)).max()


def attributable(beta, r):
    h0 = h0_for(beta, r)
    return r - (1 - np.exp(-h0))


r = 0.25
hr_rows = [("LAM", [], "HR", 1.74), ("MAS", ["LAM", "SCK"], "HR", 1.6), ("MET", [], "HR", 1.3),
           ("SCK", [], "HR", 1.4), ("RP", ALL, "HR", 1.5)]
beta, h0, res = solve(hr_rows, r)
print("=== 3a. Hazard ratios as in hazard_ratios.csv (snapshot) ===")
print("adjusted HR:", {d: round(v, 6) for d, v in zip(ids, np.exp(beta))}, f"max residual {res:.1e}")
print(f"overall risk {r}: h0 {h0:.6f}; disease-free risk {1 - np.exp(-h0):.6f}; attributable {attributable(beta, r):.6f}")
raw_beta = np.log([1.74, 1.6, 1.3, 1.4, 1.5])
print(f"attributable risk using the raw HRs as if adjusted (item 8, raw vs adjusted): {attributable(raw_beta, r):.6f}")

print("\n=== 3b. What the same population implies for other measures (period risk 0.25) ===")
implied = {}
for d in ids:
    i = ix[d]
    row = {k: measure(i, [], k, beta, h0) for k in ("HR", "RR", "OR", "RD")}
    implied[d] = row
    print(f"  {d:4s} crude HR {np.exp(row['HR']):.4f}  RR {np.exp(row['RR']):.4f}  OR {np.exp(row['OR']):.4f}  RD {row['RD']:.4f}")

print("\n=== 3c. Mixed measures recover the same adjusted HRs ===")
mix = [("LAM", [], "HR", 1.74), ("MAS", ["LAM", "SCK"], "HR", 1.6),
       ("MET", [], "RR", np.exp(implied["MET"]["RR"])),
       ("SCK", [], "OR", np.exp(implied["SCK"]["OR"])),
       ("RP", ALL, "RD", measure(ix["RP"], ALL, "RD", beta, h0))]
for d, a, k, v in mix:
    print(f"  {d:4s} {k:3s} {v:.6f} adjusted_for {a if a else '-'}")
b2, h2, res2 = solve(mix, r)
print("adjusted HR:", {d: round(v, 6) for d, v in zip(ids, np.exp(b2))}, f"max residual {res2:.1e}")
print("max |difference| from 3a:", np.abs(b2 - beta).max())

print("\n=== 3d. Entering the MET risk ratio as if it were a hazard ratio ===")
wrong = [x if x[0] != "MET" else ("MET", [], "HR", np.exp(implied["MET"]["RR"])) for x in hr_rows]
b3, _, _ = solve(wrong, r)
print("adjusted HR:", {d: round(v, 4) for d, v in zip(ids, np.exp(b3))})
print(f"attributable risk {attributable(b3, r):.4f} (correct {attributable(beta, r):.4f})")

print("\n=== 3e. The overall risk matters only for risk-based measures ===")
for rr in (0.15, 0.35):
    b4, _, _ = solve(mix, rr)
    print(f"  mixed measures at overall risk {rr}: adjusted HR", {d: round(v, 4) for d, v in zip(ids, np.exp(b4))})
