"""Independent reference values for cm_sensitivity() (v0.6 item 3; design
note section 17). Uses the prototype (prototype_v06_sensitivity.py), which
builds on reference_v041.py and reference_v05_input_check.py, not on the
package.

One input at a time, low and high value (see the prototype for the ranges),
everything else at its central value with the pooling of association rows
fixed at the central fit. Printed per input: low, high, the aggregate (the
attributable risk for event impacts) at low and high, the largest change in
a disease's rank by contribution, and each disease's adjusted value at low
and high.

Run from this directory: python3 reference_v06_sensitivity.py
"""
import numpy as np

import prototype_v06_sensitivity as P
from reference_v05_input_check import Case, read


def show(name, case, three_way=(), keys=None, adjusted=True):
    central, rows = P.one_at_a_time(case, three_way)
    print(f"\n== {name}\ncentral {central['agg']:.10f}; contributions "
          + ", ".join(f"{d} {central['contrib'][d]:.10f}" for d in case.ids))
    for r in sorted(rows, key=lambda r: r["key"]):
        if keys is not None and r["key"] not in keys:
            continue
        line = (f"{r['key']:24s} low {r['low']:.10g} high {r['high']:.10g} total {r['agg_low']:.10f} "
                f"{r['agg_high']:.10f} rank_shift {r['rank_shift']:.0f}")
        print(line)
        if adjusted and r["adj_low"] is not None:
            print("    low ", {d: round(float(v), 10) for d, v in r["adj_low"].items()})
            print("    high", {d: round(float(v), 10) for d, v in r["adj_high"].items()})


if __name__ == "__main__":
    np.set_printoptions(linewidth=160)
    # A. Five diseases, yield: every input (32, of which 10 interactions).
    show("A. five diseases, yield", Case("five_diseases", "impacts_yield.csv"), adjusted=False)
    show("A (adjusted impacts for three inputs)", Case("five_diseases", "impacts_yield.csv"),
         keys={"impacts:2", "interaction:LAM;MAS", "associations:6"})
    # B. Five diseases, yield, MAS-RP row left out (rows renumbered), and the
    # three-way ratio MAS:RP among MET on request.
    arows = [r for r in read("five_diseases", "associations.csv")
             if not (r["disease"] == "MAS" and r["associated_with"] == "RP")]
    for i, r in enumerate(arows):
        r["row"] = i + 1
    show("B. five diseases, yield, MAS-RP unknown", Case("five_diseases", "impacts_yield.csv", assoc_rows=arows),
         three_way=[("MAS", "RP", "MET")], keys={"pair:MAS:RP", "or_ratio:MAS:RP|MET", "interaction:MAS;RP",
                                                 "associations:5", "prob:RP"})
    # C. Five diseases, culling, overall risk 0.25.
    show("C. five diseases, culling", Case("five_diseases", "culling.csv", event=True, risk=0.25),
         keys={"impacts:6", "impacts:2", "risk", "interaction:MAS;SCK", "associations:1", "prob:SCK"})
