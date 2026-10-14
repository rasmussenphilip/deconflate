# Five-disease example: every feature in one input set

An illustrative dairy population with five diseases and four impact
tables. The values are made up to exercise the package, not taken from
studies. `run_all_features.R` reads the folder and runs every part of the
package on it.

```r
dir <- system.file("extdata", "five_diseases", package = "deconflate")
# One model per impact table; the disease and association tables are shared
yield <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                        associations = file.path(dir, "associations.csv"),
                        impacts = file.path(dir, "impacts_yield.csv"))
deconflate(yield)
culling <- cm_read_inputs(diseases = file.path(dir, "diseases.csv"),
                          associations = file.path(dir, "associations.csv"),
                          impacts = file.path(dir, "culling.csv"))
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
| `diseases.csv` | All three probability measures: `prevalence`, `probability` and `incidence_rate` (MAS, converted with `1 - exp(-rate)`). |
| `associations.csv` | All ten pairs, with the measures `OR`, `RR`, `RD`, `prevalence` and `phi`. LAM:MAS has two studies (rows 1 and 2), which are pooled (random effects by default). LAM:RP has an odds ratio of 1 (unrelated); RP:SCK is the odds ratio of a study's 2x2 table. Row 12 is an `OR_ratio` row: the LAM-MAS odds ratio is 1.5 times higher among cows with SCK. It changes additive results only with interactions, rows about combinations or subgroups, or unknown pairs, and it changes event impacts. Leaving a row out makes the pair unknown (the script shows this as a variant). |
| `impacts_yield.csv` | Crude yield losses in % of yield, one row per disease (with only the ten crude pair rows, the pairwise shortcut gives the answer). |
| `impacts_calving_interval.csv` | Days: two MAS rows (a coefficient adjusted for SCK and a crude difference), fitted together by weighted least squares; SCK adjusted for `all` (exact); RP among cows with no other disease (`among = "no others"`). |
| `impacts_welfare.csv` | Score points, with interactions: LAM;MAS and MET;RP adjusted for `all` (the interactions themselves, synergistic), and a crude LAM;SCK row (cows with both against cows with neither) that implies an antagonistic LAM:SCK interaction. |
| `culling.csv` | Event impacts (`event_model = TRUE`, with `overall_risk` or `baseline_risk`): hazard ratios (LAM crude, MAS stratified by LAM and SCK), a risk ratio (MET), an odds ratio (SCK), a risk difference (RP stratified by all other diseases) and a crude hazard ratio for cows with both LAM and MAS (an interaction on the log-hazard scale), all on the snapshot hazard model. |

Every table has free-text `study` and `note` columns describing its rows.

## Uncertainty in the tables

Uncertain values have a distribution in their own row, in the columns
`dist` and `p1`-`p4`; rows with an empty `dist` are point values. The point
value (`value`) gives the central estimate; `deconflate()` draws the
distributions (`n_draws`, 1000 by default) for the intervals. Every
distribution type is used at least once:

| `dist` | Used for |
|---|---|
| `beta` | LAM prevalence; MAS prevalence among SCK |
| `pert` | MAS incidence rate; MET:SCK risk ratio; MAS yield; MET calving interval; MAS welfare; MET culling risk ratio |
| `pert_mean` | SCK probability; MET yield |
| `uniform` | MET prevalence; LAM:MET phi; RP welfare; MET;RP interaction |
| `lognormal_ci` | both LAM:MAS odds ratios; the OR_ratio row; LAM culling hazard ratio; LAM;MAS culling hazard ratio |
| `lognormal` | MET:RP odds ratio; RP calving interval |
| `normal` | LAM:SCK odds ratio (truncated at 0); MAS:MET risk difference; LAM yield; LAM and both MAS calving intervals; LAM and SCK welfare; LAM;MAS interaction; RP culling risk difference |
| `fixed` | SCK yield (a point mass) |

Point values: RP prevalence, the LAM:RP, RP:SCK and MAS:RP odds ratios, RP
yield, SCK calving interval, MET welfare, the LAM;SCK welfare row and the
MAS and SCK culling estimates. The overall risk of culling is an argument
of `deconflate()` and can also be a distribution (e.g. `dist_beta(250, 750)`).

## What the script covers

1. Reading and checking inputs (`cm_check_inputs`, `cm_read_inputs` with one
   impact table per model).
2. Feasibility (`check_feasibility`), the exact and sampled joint
   distribution (`fit_joint`), the most common combinations, and the `OR_ratio` row
   in the fitted distribution.
3. `deconflate()` for each additive table: pooled association rows and
   their heterogeneity, the pairwise shortcut, several impact rows (weighted
   least squares), adjusted and subgroup rows, interactions, the diagnosis
   of the inputs (`diagnose = TRUE`).
4. Event impacts: mixed measures, an interaction, the attributable risk and
   its Shapley allocation (own part and interaction share), the raw
   attributable risk, results on the risk scale, `baseline_risk`, the
   diagnosis, and the `event_model` checks.
5. Uncertainty: intervals from draws, rejections, Latin hypercube sampling,
   and an uncertain overall risk.
6. Reporting (`contribution_table`, `summary`, `$contributions`), and
   gaps and values computed from the results in a few lines of base R
   (also over the draws).
7. `sensitivity_oat`, `screen_associations`, `screen_interactions`,
   `screen_three_way` (also for event impacts and without any associations),
   and scenarios built with the `set_*` functions.
8. `cm_threshold` for every conclusion (`rank`, `sign`, `total`, `change`)
   and input type (association, impact, interaction, OR_ratio, any row by
   number, disease probability, overall risk); poles are reported as
   discontinuities and unusable stretches of a range in `regions`.
9. `simulate_raw_impacts` (crude, adjusted, subgroup and combination rows,
   finite samples, event impacts).
10. Variants: an unknown pair (filled in by the fit and listed in
    `$unknown_pairs`), an odds ratio adjusted for other diseases, and the
    messages for files of earlier versions.
11. Plots.

Reference values for this set (computed independently in Python) are in
`inst/validation/reference_v041.py` and are checked by the tests.
