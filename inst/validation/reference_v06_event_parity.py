"""Independent reference values for the event-model parity work (v0.6;
design note section 15). Builds on reference_v041.py and does not use the
package.

1. Raw results: the same model with the diseases unrelated (same
   prevalences and impact rows, no association rows).
2. Contribution split: each disease's own part is its Shapley value with
   the interaction terms set to zero (baseline hazard kept); the interaction
   share is the rest. For additive impacts this is p_d b_d exactly.
3. Event row equations: the derivative of each row's model value (on its
   measure's scale) with respect to each log-hazard term, at the solution.
4. Risk-scale results: an animal with only disease d vs a disease-free
   animal: risk R_d = 1 - exp(-h0 exp(beta_d)), R0 = 1 - exp(-h0); RR, RD
   and OR.
5. Simulated raw estimates (exact): each row's value given the true terms.
6. baseline_risk: h0 = -log(1 - baseline_risk).

Run from this directory: python3 reference_v06_event_parity.py
"""
import itertools
import math

import numpy as np
from scipy.optimize import brentq, least_squares

import reference_v041 as R
from reference_v05_input_check import fit_frozen, read

FOLDER = "five_diseases"


def population(folder, unrelated=False):
    d = read(folder, "diseases.csv")
    ids = [r["disease"] for r in d]
    P = [1 - math.exp(-float(r["value"])) if r["measure"] == "incidence_rate" else float(r["value"]) for r in d]
    pop = R.Pop(ids, P)
    if unrelated:
        return pop, np.prod(np.where(pop.cells == 1, pop.P, 1 - pop.P), axis=1)
    specs = [R.spec(pop, r) for r in read(folder, "associations.csv")]
    p, groups, info = fit_frozen(pop, specs, None)
    return pop, p


def terms_matrix(pop, terms):
    return np.column_stack([pop.allm(R.split(t)) for t in terms]).astype(float)


def event_solve(pop, p, rows, risk=None, baseline=None):
    """Snapshot model: log-hazard terms, h0 (from the overall risk, or from
    the baseline risk)."""
    specs = [R.impact_spec(pop, r) for r in rows]
    terms = R.terms_of(pop, specs)
    F = terms_matrix(pop, terms)
    sus = [R.event_setup(pop, p, s, terms) for s in specs]
    target = np.array([s["v"] if s["m"] == "RD" else np.log(s["v"]) for s in specs])

    def h0_of(mult):
        if baseline is not None:
            return -math.log(1 - baseline)
        return R.baseline_h0(mult, p, risk)

    def values(theta):
        mult = np.exp(F @ theta)
        h0 = h0_of(mult)
        Rk = 1 - np.exp(-h0 * mult)
        return np.array([R.event_value(su, s["m"], theta, mult, Rk, p) for su, s in zip(sus, specs)])

    sol = least_squares(lambda th: values(th) - target, np.zeros(len(terms)), xtol=1e-15, ftol=1e-15, gtol=1e-15)
    theta = sol.x
    mult = np.exp(F @ theta)
    h0 = h0_of(mult)
    overall = float((p * (1 - np.exp(-h0 * mult))).sum())
    return dict(terms=terms, theta=theta, h0=h0, overall=overall, F=F, values=values, specs=specs,
                resid=float(np.abs(sol.fun).max()))


def shapley_event(pop, p, terms, theta, h0, drop_interactions=False):
    """Cell-wise Shapley allocation of the attributable risk."""
    ids = pop.ids
    th = np.array(theta, float).copy()
    if drop_interactions:
        th[len(ids):] = 0.0
    r0 = 1 - math.exp(-h0)
    tset = [set(R.split(t)) for t in terms]
    phi = {d: 0.0 for d in ids}
    for c, pc in zip(pop.cells, p):
        present = [d for d, x in zip(ids, c) if x == 1]
        if not present or pc == 0:
            continue
        memo = {}

        def loss(sub):
            key = frozenset(sub)
            if key not in memo:
                lp = sum(th[j] for j, s in enumerate(tset) if s <= key)
                memo[key] = 1 - math.exp(-h0 * math.exp(lp)) - r0
            return memo[key]

        k = len(present)
        for d in present:
            rest = [x for x in present if x != d]
            for size in range(k):
                w = math.factorial(size) * math.factorial(k - size - 1) / math.factorial(k)
                for sub in itertools.combinations(rest, size):
                    phi[d] += pc * w * (loss(list(sub) + [d]) - loss(list(sub)))
    return phi


def jacobian(sol, h=1e-6):
    th = sol["theta"]
    base = sol["values"](th)
    J = np.empty((len(base), len(th)))
    for j in range(len(th)):
        e = np.zeros(len(th))
        e[j] = h
        J[:, j] = (sol["values"](th + e) - sol["values"](th - e)) / (2 * h)
    return J


def alone(sol, ids):
    h0 = sol["h0"]
    r0 = 1 - math.exp(-h0)
    out = {}
    for i, d in enumerate(ids):
        rd = 1 - math.exp(-h0 * math.exp(sol["theta"][i]))
        out[d] = dict(risk=rd, RR=rd / r0, RD=rd - r0, OR=(rd / (1 - rd)) / (r0 / (1 - r0)))
    return out, r0


def simulate_exact(pop, p, rows, true_theta, risk=None, baseline=None):
    """What each row would report, given the true log-hazard terms."""
    specs = [R.impact_spec(pop, r) for r in rows]
    terms = R.terms_of(pop, specs)
    F = terms_matrix(pop, terms)
    th = np.array([true_theta.get(t, 0.0) for t in terms])
    mult = np.exp(F @ th)
    h0 = -math.log(1 - baseline) if baseline is not None else R.baseline_h0(mult, p, risk)
    Rk = 1 - np.exp(-h0 * mult)
    out = []
    for s in specs:
        su = R.event_setup(pop, p, s, terms)
        v = R.event_value(su, s["m"], th, mult, Rk, p)
        out.append(v if s["m"] == "RD" else math.exp(v))
    return out


def fmt(d):
    return {k: round(float(v), 10) for k, v in d.items()}


if __name__ == "__main__":
    pop, p = population(FOLDER)
    popu, pu = population(FOLDER, unrelated=True)
    Pm = pop.cells.T @ p

    print("== 1. Raw results (the diseases unrelated)")
    for f in ("impacts_yield.csv", "impacts_calving_interval.csv", "impacts_welfare.csv"):
        rows = read(FOLDER, f)
        b, contrib, total, *_ = R.adjust_additive(pop, p, rows)
        bu, cu, tu, *_ = R.adjust_additive(popu, pu, rows)
        print(f"  {f}: adjusted total {total:.10f}; raw total {tu:.10f}")
    rows = read(FOLDER, "impacts_yield.csv")
    print(f"  yield, sum p * raw: {sum(Pm[i] * float(r['value']) for i, r in enumerate(rows)):.10f}")
    cull = read(FOLDER, "culling.csv")
    sol = event_solve(pop, p, cull, risk=0.25)
    solu = event_solve(popu, pu, cull, risk=0.25)
    ids = pop.ids
    phi = shapley_event(pop, p, sol["terms"], sol["theta"], sol["h0"])
    phiu = shapley_event(popu, pu, solu["terms"], solu["theta"], solu["h0"])
    att = 0.25 - (1 - math.exp(-sol["h0"]))
    attu = 0.25 - (1 - math.exp(-solu["h0"]))
    print(f"  culling (overall risk 0.25): adjusted HRs {fmt(dict(zip(sol['terms'], np.exp(sol['theta']))))}")
    print(f"     raw HRs (unrelated) {fmt(dict(zip(solu['terms'], np.exp(solu['theta']))))}")
    print(f"     attributable: adjusted {att:.10f}, raw {attu:.10f}")
    print(f"     raw Shapley {fmt(phiu)}")

    print("== 2. Contribution split")
    welfare = read(FOLDER, "impacts_welfare.csv")
    b, contrib, total, *_ = R.adjust_additive(pop, p, welfare)
    print(f"  welfare (additive): own p*b {fmt({d: Pm[i] * b[d] for i, d in enumerate(ids)})}")
    print(f"     total {fmt(contrib)}")
    main = shapley_event(pop, p, sol["terms"], sol["theta"], sol["h0"], drop_interactions=True)
    print(f"  culling: total {fmt(phi)} (sum {sum(phi.values()):.10f})")
    print(f"     own {fmt(main)} (sum {sum(main.values()):.10f})")
    print(f"     interaction share {fmt({d: phi[d] - main[d] for d in ids})}")

    print("== 3. Event row equations (at the solution, and at HR = 1)")
    J = jacobian(sol)
    for r, row in zip(cull, J):
        print(f"  row {r['row']} {r['disease']} {r['measure']}: {np.round(row, 10).tolist()}")
    sol0 = dict(sol)
    sol0["theta"] = np.zeros(len(sol["theta"]))
    J0 = jacobian(sol0)
    print("  at zero:", np.round(J0, 6).tolist())
    yrows = [dict(r, measure="", value="1") for r in cull]
    specs = [R.impact_spec(pop, r) for r in yrows]
    A = np.array([R.additive_coefs(pop, p, s, sol["terms"]) for s in specs])
    print("  additive coefficients of the same rows:", np.round(A, 6).tolist())

    print("== 4. Risk scale: an animal with only the disease vs a disease-free animal")
    al, r0 = alone(sol, ids)
    print(f"  disease-free risk {r0:.10f}")
    for d in ids:
        print(f"  {d}: {fmt(al[d])}")

    print("== 5. Simulated raw estimates (exact), then adjusted back")
    true = dict(zip(sol["terms"], sol["theta"]))
    sim = simulate_exact(pop, p, cull, true, risk=0.25)
    print("  simulated:", [round(x, 10) for x in sim], "inputs:", [float(r["value"]) for r in cull])

    print("== 6. baseline_risk")
    solb = event_solve(pop, p, cull, baseline=0.18)
    print(f"  baseline 0.18: HRs {fmt(dict(zip(solb['terms'], np.exp(solb['theta']))))}; overall risk {solb['overall']:.10f};"
          f" attributable {solb['overall'] - 0.18:.10f}")
    solbu = event_solve(popu, pu, cull, baseline=0.18)
    print(f"     raw (unrelated): overall risk {solbu['overall']:.10f}; attributable {solbu['overall'] - 0.18:.10f}")

    print("== 7. Simulated raw estimates (exact): supplement population, HR, RR, OR and RD rows")
    sp = R.Pop(["d1", "d2", "d3"], [0.10, 0.15, 0.20])
    arows = [dict(row=1, disease="d1", associated_with="d2", among="", adjusted_for="", measure="OR", value="2", dist=""),
             dict(row=2, disease="d1", associated_with="d3", among="", adjusted_for="", measure="OR", value="1", dist=""),
             dict(row=3, disease="d2", associated_with="d3", among="", adjusted_for="", measure="OR", value="3", dist="")]
    ps, _, _ = fit_frozen(sp, [R.spec(sp, r) for r in arows], None)
    srows = [dict(row=i + 1, disease=d, among="", adjusted_for="", measure=m, value="1", dist="")
             for i, (d, m) in enumerate((("d1", "HR"), ("d2", "RR"), ("d3", "OR"), ("d1;d2", "RD")))]
    truth = {"d1": math.log(1.2), "d2": math.log(1.5), "d3": math.log(1.1), "d1;d2": math.log(1.4)}
    print("  overall risk 0.2:", [round(x, 10) for x in simulate_exact(sp, ps, srows, truth, risk=0.2)])
    print("  baseline risk 0.15:", [round(x, 10) for x in simulate_exact(sp, ps, srows, truth, baseline=0.15)])
