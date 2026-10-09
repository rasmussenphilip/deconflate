"""Prototype for v0.4 item 4, part 2: the informativeness toggle and
probabilistic constraints on the 2024 yield inputs (analysis inputs as in
reference_2024_analysis.py; simultaneous method; unlisted pairs independent).
Run from this directory: python3 prototype_v04_part2.py [N]

Draw N association sets (ORs) and N impact sets (raw yield impacts) from the
input distributions. Every combination (s, t) gives adjusted impacts
b[s, t] = A(OR_s)^-1 y_t; F[s, t] = 1 when every adjusted impact is >= 0.
  unconstrained : all pairs, equal weight (the current Monte Carlo)
  balanced      : pairs with F = 1, equal weight (joint conditioning: both
                  kinds of input are updated)
  associations  : weight F[s, t] / Z_s, Z_s = mean_t F[s, t] (the OR draws keep
                  their input distribution; only the impacts are updated)
  impacts       : weight F[s, t] / Z_t (the impact draws keep their input
                  distribution; only the ORs are updated)
"""
import sys
import numpy as np
import reference_2024_analysis as r24   # prints its own reference values first

N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
D, n, P = r24.D, r24.n, r24.P
rng = np.random.default_rng(404)
OR = np.ones((N, n, n))
for i, k, spec in r24.pairs:
    v = r24.draw(spec, N, rng, lower=0)
    OR[:, i, k] = OR[:, k, i] = v
E = r24.excess(OR)
A = np.eye(n)[None] + np.transpose(E, (0, 2, 1))      # A[s, i, k]
Ainv = np.linalg.inv(A)
Y = np.column_stack([r24.draw(r24.YIELD[d], N, rng) for d in D])   # percent
logOR = np.array([np.log(OR[:, i, k]) for i, k, _ in r24.pairs]).T  # (N, pairs)

# Feasibility grid, in chunks
def grid(cols):
    F = np.zeros((N, N), bool)
    for s0 in range(0, N, 200):
        B = np.einsum('sik,tk->sti', Ainv[s0:s0 + 200][:, cols], Y)  # (chunk, N, |cols|)
        F[s0:s0 + 200] = (B >= 0).all(axis=2)
    return F


# ------------------------------------------------------------ 2a
F = grid(list(range(n)))
print(f"\n=== 2a. Constraint: every adjusted impact >= 0 ===")
print(f"N = {N} OR sets x {N} impact sets; share of combinations that satisfy it: {F.mean():.5f}")
Bdiag = np.einsum('sik,sk->si', Ainv, Y)        # matched pairs = the usual Monte Carlo
print("unconstrained P(adjusted < 0):", ', '.join(f"{d} {np.mean(Bdiag[:, j] < 0):.3f}" for j, d in enumerate(D)))

# ------------------------------------------------------------ 2b
from scipy.optimize import minimize
print("\n=== 2b. Best-fit diagnostic: smallest shift of the inputs (in input SDs) that satisfies it ===")
y0 = np.array([r24.central(r24.YIELD[d]) for d in D])
# Input SDs from a large separate sample, so that they do not depend on N.
rng_sd = np.random.default_rng(7)
ysd = np.array([r24.draw(r24.YIELD[d], 200000, rng_sd).std() for d in D])
lo0 = np.array([np.log(r24.central(spec)) for _, _, spec in r24.pairs])
losd = np.array([np.log(r24.draw(spec, 200000, rng_sd, lower=0)).std() for _, _, spec in r24.pairs])
ymov = np.nonzero(ysd > 1e-9)[0]
omov = np.nonzero(losd > 1e-9)[0]


def adjusted(uy, uo):
    y = y0.copy(); y[ymov] += ysd[ymov] * uy
    lo = lo0.copy(); lo[omov] += losd[omov] * uo
    O = np.ones((1, n, n))
    for q, (i, k, _) in enumerate(r24.pairs):
        O[0, i, k] = O[0, k, i] = np.exp(lo[q])
    Ax = np.eye(n) + r24.excess(O)[0].T
    return np.linalg.solve(Ax, y)


b0 = adjusted(np.zeros(len(ymov)), np.zeros(len(omov)))
print("central adjusted impacts:", ', '.join(f"{d} {v:.3f}" for d, v in zip(D, b0)))
for mode in ("balanced", "associations", "impacts"):
    my, mo = mode != "impacts", mode != "associations"
    ny, no = len(ymov) * my, len(omov) * mo
    split = lambda u: (u[:ny] if my else np.zeros(len(ymov)), u[ny:] if mo else np.zeros(len(omov)))
    con = [{"type": "ineq", "fun": lambda u: adjusted(*split(u))}]
    r = minimize(lambda u: u @ u, np.zeros(ny + no), jac=lambda u: 2 * u, method="SLSQP",
                 constraints=con, options={"maxiter": 2000, "ftol": 1e-12})
    for _ in range(5):   # restart from the last point until SLSQP reports convergence
        if r.success:
            break
        r = minimize(lambda u: u @ u, r.x, jac=lambda u: 2 * u, method="SLSQP",
                     constraints=con, options={"maxiter": 2000, "ftol": 1e-12})
    uy, uo = split(r.x)
    b = adjusted(uy, uo)
    big = sorted([(abs(v), f"{D[j]} {v:+.2f}") for j, v in zip(ymov, uy) if abs(v) >= 0.5] +
                 [(abs(v), f"{D[r24.pairs[q][0]]}-{D[r24.pairs[q][1]]} {v:+.2f}") for q, v in zip(omov, uo) if abs(v) >= 0.5],
                 reverse=True)
    print(f"\n  {mode} (converged: {r.success}); total shift sqrt(sum z^2) = {np.sqrt(r.fun):.2f} SD")
    print("    shifts >= 0.5 SD:", ', '.join(x[1] for x in big))
    print("    adjusted:", ', '.join(f"{d} {v:.2f}" for d, v in zip(D, b)), f"| aggregate {P @ b:.3f} (central {P @ b0:.3f})")

# ------------------------------------------------------------ 2c
cons = [D.index("CM"), D.index("DA")]
F = grid(cons)
Zs, Zt = F.mean(1), F.mean(0)
print(f"\n=== 2c. Constraint on CM and DA only (adjusted >= 0) ===")
print(f"share of combinations that satisfy it: {F.mean():.4f}; OR sets with none: {np.mean(Zs == 0):.4f}; impact sets with none: {np.mean(Zt == 0):.4f}")
s_idx, t_idx = np.nonzero(F)
modes = {'balanced': np.ones(len(s_idx)), 'associations': 1 / Zs[s_idx], 'impacts': 1 / Zt[t_idx]}
M = 40000
res = {'unconstrained': Bdiag}
ess, sel_pairs = {}, {}
for m, w in modes.items():
    w = w / w.sum()
    ess[m] = 1 / np.sum(w ** 2)
    pick = rng.choice(len(s_idx), size=M, p=w)
    s, t = s_idx[pick], t_idx[pick]
    res[m] = np.einsum('sik,sk->si', Ainv[s], Y[t])
    sel_pairs[m] = (s, t)
print("Adjusted yield impacts (% decrease): mean [2.5%, 97.5%]")
print("      " + "".join(f"{m:>26s}" for m in res))
for j, d in enumerate(D):
    print(f"{d:5s} " + "".join(f"{res[m][:, j].mean():8.3f} [{np.quantile(res[m][:, j], .025):6.2f},{np.quantile(res[m][:, j], .975):6.2f}]"
                             for m in res))
print("aggr  " + "".join(f"{(res[m] @ P).mean():8.3f} [{np.quantile(res[m] @ P, .025):6.2f},{np.quantile(res[m] @ P, .975):6.2f}]"
                         for m in res))
print("effective sample size (pairs):", {m: round(v) for m, v in ess.items()})
print("How far the inputs moved (posterior mean - input mean, in input SDs); |shift| >= 0.1 shown")
for m, (s, t) in sel_pairs.items():
    ysh = (Y[t].mean(0) - Y.mean(0)) / np.where(Y.std(0) > 1e-9, Y.std(0), np.inf)
    osh = (logOR[s].mean(0) - logOR.mean(0)) / np.where(logOR.std(0) > 1e-9, logOR.std(0), np.inf)
    yl = [f"{D[j]} {ysh[j]:+.2f}" for j in np.argsort(-np.abs(ysh)) if abs(ysh[j]) >= 0.1]
    ol = [f"{D[r24.pairs[q][0]]}-{D[r24.pairs[q][1]]} {osh[q]:+.2f}" for q in np.argsort(-np.abs(osh)) if abs(osh[q]) >= 0.1]
    print(f"  {m:12s} impacts: {', '.join(yl) or 'none'}")
    print(f"  {'':12s} ORs:     {', '.join(ol) or 'none'}")
