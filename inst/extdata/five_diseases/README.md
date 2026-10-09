# Five-disease example: every feature in one input set

An illustrative dairy population with five diseases and four impact
tables. The values are made up to exercise the package, not taken from
studies. `run_all_features.R` reads the folder and runs every part of the
package on it.

```r
dir <- system.file("extdata", "five_diseases", package = "deconflate")
# One model per impact table; the disease, association and three-way tables
# are shared
yield <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                        associations = file.path(dir, "associations.csv"),
                        impacts = file.path(dir, "impacts_yield.csv"))
deconflate(yield)
culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = file.path(dir, "associations.csv"),
                          impacts = file.path(dir, "culling.csv"),
                          three_way = file.path(dir, "three_way.csv"))
deconflate(culling, event_model = TRUE, overall_risk = 0.25)

source(file.path(dir, "run_all_features.R"), echo = TRUE, max.deparse.length = Inf)

# To keep the printed output and the plots, set a folder first:
out_dir <- "C:/Users/me/Desktop/five_diseases_output"
source(file.path(dir, "run_all_features.R"), echo = TRUE, max.deparse.length = Inf)
```

Diseases: LAM (lameness), MAS (mastitis), MET (metritis), SCK (subclinical
ketosis), RP (retained placenta).

## Files

| File | What it shows |
|---|---|
| `diseases.csv` | All three probability types: `prevalence`, `probability` and `incidence_rate` (MAS, converted with `1 - exp(-rate)`). |
| `associations.csv` | All ten pairs, with every measure: `OR`, `RR`, `RD`, `cond_prob` and `phi`. LAM:RP has an odds ratio of 1 (unrelated); RP:SCK is the odds ratio of a study's 2x2 table. Leaving a row out makes the pair unknown (the script shows this as a variant). |
| `three_way.csv` | A three-way scenario for LAM, MAS and SCK (ratio of conditional odds ratios 1.5). It switches `deconflate()` to the global method; it changes additive results only with interactions or unknown pairs, and it changes event impacts. |
| `impacts_yield.csv` | Crude yield losses in % of yield: the simultaneous method applies, and the published approximation can be shown beside it with `compare_methods()`. |
| `impacts_calving_interval.csv` | Days, with two regression coefficients (`adjusted_linear`): MAS adjusted for SCK, SCK adjusted for `all`. The published approximation does not apply (crude estimates only). |
| `impacts_welfare.csv` + `interactions_welfare.csv` | Score points, with a synergistic (LAM:MAS, MET:RP) and an antagonistic (LAM:SCK) interaction: global method. |
| `culling.csv` | Event impacts (`event_model = TRUE`, with `overall_risk`): hazard ratios (LAM crude, MAS stratified by LAM and SCK), a risk ratio (MET), an odds ratio (SCK) and a risk difference (RP stratified by all other diseases), all on the snapshot hazard model. |

Every table has a free-text `note` column describing its rows.

## Uncertainty in the tables

Uncertain values have a distribution in their own row, in the columns
`dist` and `p1`-`p4`; rows with an empty `dist` are point values. The point
value (`value`, or `ratio` for the three-way term) gives the central
estimate; `deconflate()` draws the distributions (`n_draws`, 1000 by
default) for the intervals. Every distribution type is used at least once:

| `dist` | Used for |
|---|---|
| `beta` | LAM prevalence; MAS:SCK conditional probability |
| `pert` | MAS incidence rate; MET:SCK risk ratio; MAS yield; MET calving interval; MAS welfare; MET culling risk ratio |
| `pert_mean` | SCK probability; MET yield |
| `uniform` | MET prevalence; LAM:MET phi; RP welfare; MET:RP interaction |
| `lognormal_ci` | LAM:MAS odds ratio; the three-way ratio; LAM culling hazard ratio |
| `lognormal` | MET:RP odds ratio; RP calving interval |
| `normal` | LAM:SCK odds ratio (truncated at 0); MAS:MET risk difference; LAM yield; LAM and MAS calving interval; LAM and SCK welfare; LAM:MAS interaction; RP culling risk difference |
| `fixed` | SCK yield (a point mass) |

Point values: RP prevalence, the LAM:RP, RP:SCK and MAS:RP odds ratios, RP
yield, SCK calving interval, MET welfare, the LAM:SCK interaction and the
MAS and SCK culling estimates. The overall risk of culling is an argument
of `deconflate()` and can also be a distribution (e.g. `dist_beta(250, 750)`).

## What the script covers

1. Reading and checking inputs (`cm_check_inputs`, `cm_read_inputs` with one
   impact table per model, `pair_tables`).
2. Feasibility (`check_feasibility`), the exact and sampled joint
   distribution (`fit_joint`), `combination_probs`.
3. `deconflate()` for each additive table: the automatic choice of method
   (and its notes), adjusted estimands, interactions, and the published
   approximation through `compare_methods()` (for comparison only).
4. Event impacts: mixed measures, the attributable risk and its Shapley
   allocation, the `event_model` checks, and `compare_methods()` for event
   impacts.
5. Uncertainty: intervals from draws, rejections, Latin hypercube sampling,
   an uncertain overall risk, and methods compared on the same draws.
6. Reporting (`contribution_table`, `summary`, `attribute_burden`), and
   gaps and values computed from the results in a few lines of base R
   (also over the draws).
7. `sensitivity_oat`, `screen_associations`, `screen_interactions`,
   `screen_three_way` (also for event impacts and without any associations),
   and scenarios built with the `set_*` functions.
8. `cm_threshold` for every conclusion (`rank`, `sign`, `total`, `change`)
   and input type (association, impact, interaction, three-way ratio,
   disease probability, overall risk); poles are reported as discontinuities
   and unusable stretches of a range in `regions`.
9. `simulate_raw_impacts` (crude, adjusted, finite samples) and
   `shapley_by_cell` with additive and multiplicative losses.
10. Variants: an unknown pair (filled in by the global method and listed in
    `$unknown_pairs`), the messages for retired association measures, and a
    covariate-adjusted risk ratio used as marginal.
11. Plots.

Reference values for this set (computed independently in Python) are in
`inst/validation/reference_five_diseases.py` and `reference_v040.py` and
are checked by the tests.
