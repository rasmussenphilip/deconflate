"""Independent reference values for deconflate(diagnose = TRUE) (v0.4,
item 5; design note section 14). Builds on the item 4 reference
implementation (reference_v041.py) and does not use the package.

Definitions (section 14 of the design note):
- A disease term is flagged when c = sign(r) a - 0.1 |r| < 0, with a the
  adjusted value and r the raw reference on the impact's scale (additive:
  units; event: on the rows' measure when they share one, log for HR, rate
  ratio, RR and OR, RD as is; the adjusted HR, or the RR, OR or RD of an
  animal with only the disease against a disease-free animal; with mixed
  measures, the log HR against the ratio rows).
  With several rows, r is their inverse-variance weighted mean on that
  scale (an exact row's value when one is exact; for event impacts, ratio
  rows only when there are any). r stays at the inputs as entered while
  inputs are changed.
- Inputs with a distribution: normal score u = Phi^-1(F(x)) (clamped to
  [-5, 5]); the effect of +1 SD is d(adjusted)/du (natural scale: units or
  HR), and the shift that alone removes the flag is -c / (dc/du).
  Association pooling is frozen at the central fit (as in the draws).
- Inputs without a distribution: the value at which c = 0, on the grid
  value x (0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.6, 2.2, 3) (valid values only;
  RD and phi: value + (-0.5, -0.3, -0.15, 0, 0.15, 0.3, 0.5) within
  (-1, 1)), refined by root finding; the root nearest the value (relative
  distance) is reported. Additive impact rows: solved exactly (linear),
  reported for the flagged disease's own rows, and for other rows when
  within a quarter to three times the value.
- Disagreeing rows: |standardised residual| >= 2; the model is rerun without
  the row (pooling not frozen).

Run from this directory: python3 reference_v05_input_check.py
"""
import csv
import math
import os

import numpy as np
from scipy import stats
from scipy.optimize import brentq, least_squares

import reference_v041 as R

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(HERE, "..", "extdata")
GRID = (0.25, 0.4, 0.6, 0.8, 1, 1.25, 1.6, 2.2, 3)
OFFSETS = (-0.5, -0.3, -0.15, 0, 0.15, 0.3, 0.5)
H = 1e-4


def read(folder, name):
    with open(os.path.join(EXT, folder, name), newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for i, r in enumerate(rows):
        r["row"] = i + 1
    return rows


def u_of(qf, v):
    f = lambda x: qf(stats.norm.cdf(x)) - v
    if f(-5.0) >= 0:
        return -5.0
    if f(5.0) <= 0:
        return 5.0
    return brentq(f, -5.0, 5.0, xtol=1e-14)


def ev_scale(v, m):
    return v if m == "RD" else math.log(v)


# ------------------------------------------------------------------ fitting
def fit_frozen(pop, specs, frozen=None, tol=1e-13):
    """reference_v041.fit, with pooled groups at fixed weights (frozen)."""
    groups = R.groups_of(specs)
    p = np.prod(np.where(pop.cells == 1, pop.P, 1 - pop.P), axis=1)
    masks = [pop.allm(gr["event"][0]) & pop.nonem(gr["event"][1]) for gr in groups]
    margins = [pop.cells[:, i] == 1 for i in range(len(pop.ids))]
    info = None
    for outer in range(5000):
        targets, info = [], []
        for gr in groups:
            zs = [R.row_z(pop, s, s["v"], p, gr["scale"]) for s in gr["specs"]]
            fz = frozen.get(gr["event"]) if frozen else None
            if fz is not None:
                zstar, inf = float(np.dot(fz, zs) / np.sum(fz)), dict(exact=False)
                ses = [1.0] * len(zs)
            else:
                ses = [R.se_on(s["qf"], lambda v, s=s: R.row_z(pop, s, v, p, gr["scale"])) if len(gr["specs"]) > 1 else 1
                       for s in gr["specs"]]
                if len(gr["specs"]) == 1:
                    ses = [0 if gr["specs"][0]["qf"] is None else 1]
                zstar, inf = R.pool(zs, ses)
            inf.update(z=zs, se=ses, estimate=zstar)
            info.append(inf)
            s1 = gr["specs"][int(np.argmax(np.array(ses) == 0))] if inf.get("exact") else gr["specs"][0]
            targets.append(R.z_target(pop, s1, zstar, p, gr["scale"]))
        for sweep in range(10):
            for m, t in zip(masks + margins, targets + list(pop.P)):
                s = p[m].sum()
                p[m] *= t / s
                p[~m] *= (1 - t) / (1 - s)
        err = max(abs(p[m].sum() - t) for m, t in zip(masks + margins, targets + list(pop.P)))
        new = []
        for gr, inf in zip(groups, info):
            s1 = gr["specs"][int(np.argmax(np.array(inf["se"]) == 0))] if inf.get("exact") else gr["specs"][0]
            new.append(R.z_target(pop, s1, inf["estimate"], p, gr["scale"]))
        if err < tol and max(abs(a - b) for a, b in zip(new, targets)) < tol:
            break
    return p, groups, info


def event_adjust(pop, p, rows, risk):
    specs = [R.impact_spec(pop, r) for r in rows]
    terms = R.terms_of(pop, specs)
    F = np.column_stack([pop.allm(R.split(t)) for t in terms]).astype(float)
    sus = [R.event_setup(pop, p, s, terms) for s in specs]
    target = np.array([s["v"] if s["m"] == "RD" else np.log(s["v"]) for s in specs])

    def resid(theta):
        mult = np.exp(F @ theta)
        h0 = R.baseline_h0(mult, p, risk)
        Rk = 1 - np.exp(-h0 * mult)
        return np.array([R.event_value(su, s["m"], theta, mult, Rk, p) for su, s in zip(sus, specs)]) - target

    sol = least_squares(resid, np.zeros(len(terms)), xtol=1e-15, ftol=1e-15, gtol=1e-15)
    theta = sol.x
    mult = np.exp(F @ theta)
    h0 = R.baseline_h0(mult, p, risk)
    r0 = 1 - np.exp(-h0)
    alone = {}
    for i, d in enumerate(pop.ids):
        rd = 1 - np.exp(-h0 * np.exp(theta[i]))
        alone[d] = dict(HR=float(np.exp(theta[i])), RR=rd / r0, RD=rd - r0, OR=(rd / (1 - rd)) / (r0 / (1 - r0)))
    return dict(zip(terms, np.exp(theta))), risk - r0, alone


class Case:
    def __init__(self, folder, impacts, event=False, risk=None, impact_rows=None, assoc_rows=None):
        self.folder = folder
        self.drows = read(folder, "diseases.csv")
        self.arows = assoc_rows if assoc_rows is not None else read(folder, "associations.csv")
        self.irows = impact_rows if impact_rows is not None else read(folder, impacts)
        for i, r in enumerate(self.irows):
            r["row"] = i + 1
        self.ids = [r["disease"] for r in self.drows]
        self.event, self.risk = event, risk
        # inputs: key, label, value, quantile function (None: exact), kind, measure
        self.inputs = []
        for r in self.drows:
            self.inputs.append((f"prob:{r['disease']}", f"{r['disease']} {r['measure']}", float(r["value"]),
                                R.quantile_fn(r), "prob", r["measure"]))
        for r in self.arows:
            lab = f"{r['disease']}-{r['associated_with']} {r['measure']}" + (f" among {r['among']}" if r["among"] else "")
            self.inputs.append((f"associations:{r['row']}", lab, float(r["value"]), R.quantile_fn(r), "assoc", r["measure"]))
        for r in self.irows:
            self.inputs.append((f"impacts:{r['row']}", f"{r['disease']} impact row {r['row']}", float(r["value"]),
                                R.quantile_fn(r), "impact", (r.get("measure") or "").strip()))
        if event:
            self.inputs.append(("risk", "overall risk", risk, None, "risk", "risk"))
        self.central = self.run({}, frozen=None)
        self.frozen = {}
        for gr, inf in zip(self.central["groups"], self.central["info"]):
            if len(gr["specs"]) > 1 and not inf.get("exact"):
                se = np.array(inf["se"])
                self.frozen[gr["event"]] = 1 / (se ** 2 + inf["tau2"])

    def run(self, changes, frozen="central", drop=None):
        """Adjust with some inputs changed (key -> value); `drop`: a key
        whose row is left out."""
        if isinstance(frozen, str):
            frozen = self.frozen
        P = []
        for r in self.drows:
            v = changes.get(f"prob:{r['disease']}", float(r["value"]))
            P.append(1 - math.exp(-v) if r["measure"] == "incidence_rate" else v)
        pop = R.Pop(self.ids, P)
        specs = []
        for r in self.arows:
            k = f"associations:{r['row']}"
            if k == drop:
                continue
            rr = dict(r)
            rr["value"] = repr(changes.get(k, float(r["value"])))
            specs.append(R.spec(pop, rr))
        p, groups, info = fit_frozen(pop, specs, frozen)
        rows = []
        for r in self.irows:
            k = f"impacts:{r['row']}"
            if k == drop:
                continue
            rr = dict(r)
            rr["value"] = repr(changes.get(k, float(r["value"])))
            rows.append(rr)
        if self.event:
            hr, ar, alone = event_adjust(pop, p, rows, changes.get("risk", self.risk))
            return dict(adj=hr, agg=ar, alone=alone, groups=groups, info=info, rows=rows)
        b, contrib, total, resid, Am, y, se = R.adjust_additive(pop, p, rows)
        return dict(adj=b, agg=total, groups=groups, info=info, resid=resid, rows=rows, se=se)

    # ------------------------------------------------------------ flags
    def raw_ref(self, d, rows):
        own = [r for r in rows if r["disease"].strip() == d]
        ms = {(r.get("measure") or "").strip() for r in own}
        if not self.event:
            on = "units"
        elif len(ms) == 1:
            on = "HR" if ms == {"rate_ratio"} else ms.pop()
        else:
            on = "HR"
            own = [r for r in own if (r.get("measure") or "").strip() != "RD"]
        vals, ses, meas = [], [], []
        for r in own:
            m = (r.get("measure") or "").strip()
            scale = (lambda v: v) if not self.event else (lambda v, m=m: ev_scale(v, m))
            vals.append(scale(float(r["value"])))
            ses.append(R.se_on(R.quantile_fn(r), scale))
            meas.append(m)
        if not vals:
            return None, None
        vals, ses = np.array(vals), np.array(ses)
        if (ses == 0).any():
            r0 = vals[np.argmax(ses == 0)]
        else:
            w = 1 / ses ** 2
            r0 = float((w * vals).sum() / w.sum())
        return r0, on

    def compared(self, res, d):
        on = self.raw_ref(d, self.irows)[1]
        return res["adj"][d] if on in ("units", "HR") else res["alone"][d][on]

    def criterion(self, res, d, rows=None):
        # The raw reference stays at the inputs as entered; an event result is
        # compared on its rows' measure (RR, OR and RD: an animal with only
        # the disease against a disease-free animal).
        r, on = self.raw_ref(d, self.irows)
        v = self.compared(res, d)
        a = math.log(v) if self.event and on != "RD" else v
        return np.sign(r) * a - 0.1 * abs(r)

    def flagged(self):
        out = []
        for d in self.ids:
            if self.raw_ref(d, self.irows)[0] is None:
                continue
            if self.criterion(self.central, d) < 0:
                out.append(d)
        return out

    def changed_rows(self, changes):
        rows = []
        for r in self.irows:
            rr = dict(r)
            k = f"impacts:{r['row']}"
            if k in changes:
                rr["value"] = repr(changes[k])
            rows.append(rr)
        return rows

    def c_at(self, key, x, flags):
        ch = {key: x}
        res = self.run(ch)
        rows = self.changed_rows(ch)
        return np.array([self.criterion(res, d, rows) for d in flags]), res

    # ------------------------------------------------------- sensitivities
    def effects(self, flags):
        out = []
        for key, lab, v, qf, kind, meas in self.inputs:
            if qf is None:
                continue
            u0 = u_of(qf, v)
            xp, xm = float(qf(stats.norm.cdf(u0 + H))), float(qf(stats.norm.cdf(u0 - H)))
            cp, rp = self.c_at(key, xp, flags)
            cm, rm = self.c_at(key, xm, flags)
            c0 = np.array([self.criterion(self.central, d) for d in flags])
            for j, d in enumerate(flags):
                eff = (self.compared(rp, d) - self.compared(rm, d)) / (2 * H)
                dc = (cp[j] - cm[j]) / (2 * H)
                need = -c0[j] / dc if abs(dc) > 1e-12 else np.inf
                out.append((key, lab, v, d, eff, need))
        return out

    def thresholds(self, flags):
        out = []
        for key, lab, v, qf, kind, meas in self.inputs:
            if qf is not None:
                continue
            if kind == "impact" and not self.event:
                c0, _ = self.c_at(key, v, flags)
                c1, _ = self.c_at(key, v + 1, flags)
                own = [d for d in flags if any(f"impacts:{r['row']}" == key and r["disease"].strip() == d
                                               for r in self.irows)]
                lo, hi = min(0.25 * v, 3 * v), max(0.25 * v, 3 * v)
                for j, d in enumerate(flags):
                    slope = c1[j] - c0[j]
                    if abs(slope) > 1e-12:
                        x = v - c0[j] / slope
                        if d in own or lo <= x <= hi:
                            out.append((key, lab, v, d, x))
                continue
            if meas in ("RD", "phi"):
                xs = [v + o for o in OFFSETS if -1 < v + o < 1]
            else:
                xs = [v * g for g in GRID]
                if kind == "risk" or meas in ("prevalence", "probability") or (kind == "assoc" and meas == "prevalence"):
                    xs = [x for x in xs if 0 < x < 1]
            cs = []
            for x in xs:
                try:
                    cs.append(self.c_at(key, x, flags)[0])
                except Exception:
                    cs.append(None)
            for j, d in enumerate(flags):
                roots = []
                for k in range(len(xs) - 1):
                    if cs[k] is None or cs[k + 1] is None:
                        continue
                    if np.sign(cs[k][j]) != np.sign(cs[k + 1][j]):
                        f = lambda x: self.c_at(key, x, flags)[0][j]
                        roots.append(brentq(f, xs[k], xs[k + 1], xtol=1e-12))
                if roots:
                    dist = [abs(x - v) / (abs(v) if v != 0 else 1) for x in roots]
                    out.append((key, lab, v, d, roots[int(np.argmin(dist))]))
        return out

    def disagreeing(self):
        out = []
        resid = self.central.get("resid")
        if resid is not None:
            for r, z in zip(self.irows, resid):
                if np.isfinite(z) and abs(z) >= 2:
                    key = f"impacts:{r['row']}"
                    res = self.run({}, frozen=None, drop=key)
                    out.append((key, float(z), res["adj"], res["agg"]))
        for gr, inf in zip(self.central["groups"], self.central["info"]):
            if len(gr["specs"]) < 2 or inf.get("exact"):
                continue
            for s, z, se in zip(gr["specs"], inf["z"], inf["se"]):
                std = (z - inf["estimate"]) / se
                if abs(std) >= 2:
                    key = f"associations:{s['row']}"
                    res = self.run({}, frozen=None, drop=key)
                    out.append((key, float(std), res["adj"], res["agg"]))
        return out


def report(name, case, top=None):
    print(f"\n== {name}")
    flags = case.flagged()
    print("flagged:", [(d, case.raw_ref(d, case.irows)[1], round(case.raw_ref(d, case.irows)[0], 10),
                        round(float(case.compared(case.central, d)), 10),
                        round(float(case.criterion(case.central, d)), 10)) for d in flags])
    print(f"aggregate {case.central['agg']:.10f}")
    if not flags:
        return
    eff = case.effects(flags)
    for d in flags:
        rows = sorted([e for e in eff if e[3] == d], key=lambda e: -abs(e[4]))
        print(f"  {d}: effect per SD (top {top or 6}) and SD needed")
        for key, lab, v, _, e, need in rows[:top or 6]:
            print(f"    {key:18s} {lab:34s} {e:+.8f}  need {need:+.6f}")
    for key, lab, v, d, x in case.thresholds(flags):
        print(f"  threshold {key:18s} {lab:34s} value {v:g} -> {d}: {x:.8f}")


def case_e():
    """E. Five diseases, calving interval, MAS crude row 20 (instead of 9) and
    the second LAM-MAS odds ratio 6.0 (instead of 2.2): disagreeing rows."""
    rows = read("five_diseases", "impacts_calving_interval.csv")
    rows[2]["value"] = "20"
    arows = read("five_diseases", "associations.csv")
    arows[1]["value"] = "6.0"
    case = Case("five_diseases", None, impact_rows=rows, assoc_rows=arows)
    print("\n== E. five diseases, calving interval, MAS crude 20, LAM-MAS OR row 2 = 6.0")
    print("flagged:", case.flagged())
    print("adjusted", {k: round(float(v), 10) for k, v in case.central["adj"].items()}, f"aggregate {case.central['agg']:.10f}")
    for key, z, adj, agg in case.disagreeing():
        print(f"  {key}: std residual {z:+.10f}; without it", {k: round(float(v), 10) for k, v in adj.items()},
              f"aggregate {agg:.10f}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "E":
        case_e()
        sys.exit()
    only_b = len(sys.argv) > 1 and sys.argv[1] == "B"
    # A. Five diseases, yield, SCK raw impact 1.0 (exact): SCK changes sign.
    if not only_b:
        rows = read("five_diseases", "impacts_yield.csv")
        rows[3]["value"] = "1.0"
        report("A. five diseases, yield with SCK raw 1.0", Case("five_diseases", None, impact_rows=rows))

    # B. Five diseases, culling (overall risk 0.25), SCK OR 1.2: SCK on the
    # other side of 1, LAM shrunk by more than 90% on the log scale.
    rows = read("five_diseases", "culling.csv")
    rows[3]["value"] = "1.2"
    report("B. five diseases, culling with SCK OR 1.2", Case("five_diseases", None, event=True, risk=0.25,
                                                               impact_rows=rows))
    if only_b:
        sys.exit()

    # C. 2024 yield (the tables shown on 14 Oct 2026).
    report("C. global dairy 2024, yield", Case("global_dairy_2024", "impacts_yield.csv"))

    # D. 2024 yield with three CM rows (all crude): disagreeing rows.
    rows = [r for r in read("global_dairy_2024", "impacts_yield.csv") if r["disease"] != "CM"]
    cm = [r for r in read("global_dairy_2024", "impacts_yield.csv") if r["disease"] == "CM"][0]
    for v, sd, study in ((6.64, 1.6, "Boujenane 2015"), (4.0, 1.0, "Heikkila 2018"), (0.40, 0.6, "Fukushima 2022")):
        r = dict(cm)
        r.update(value=repr(v), dist="normal", p1=repr(v), p2=repr(sd), p3="", p4="", study=study)
        rows.append(r)
    case = Case("global_dairy_2024", None, impact_rows=rows)
    print("\n== D. 2024 yield with three CM rows")
    print("adjusted", {k: round(float(v), 10) for k, v in case.central["adj"].items()}, f"aggregate {case.central['agg']:.10f}")
    print("std residuals", np.round(case.central["resid"], 10))
    for key, z, adj, agg in case.disagreeing():
        print(f"  {key}: std residual {z:+.10f}; without it CM {adj['CM']:.10f}, aggregate {agg:.10f}")
    case_e()
