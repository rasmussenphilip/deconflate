# Five-disease example: every feature in one input set

An illustrative dairy population with five diseases and three impact
analyses. The values are made up to exercise the package, not taken from
studies. `run_all_features.R` reads the folder and runs every part of the
package on it.

```r
dir <- system.file("extdata", "five_diseases", package = "deconflate")
inp <- cm_read_inputs(dir = dir)
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
| `associations.csv` | All ten pairs, with every measure: `OR`, `RR`, `RD`, `cond_prob`, `phi`, `table` (counts) and `independent` (the script shows `unknown` as a variant). |
| `three_way.csv` | A three-way scenario for LAM, MAS and SCK (ratio of conditional odds ratios 1.5). It affects the global method (with interactions or unknown pairs), the hazard-ratio snapshot model and `attributable_risk()`. |
| `impacts_yield.csv` | Crude yield losses in % of yield: both methods of `deconflate()` apply, and the published approximation can be shown beside them with `compare_methods()`. |
| `impacts_calving_interval.csv` | Days, with two regression coefficients (`adjusted_linear`): MAS adjusted for SCK, SCK adjusted for `all`. The published approximation does not apply (crude estimates only). |
| `impacts_welfare.csv` + `interactions_welfare.csv` | Score points, with a synergistic (LAM:MAS, MET:RP) and an antagonistic (LAM:SCK) interaction: global method only. |
| `hazard_ratios.csv` | Culling hazard ratios for the snapshot model: `snapshot_crude` and `snapshot_stratified` (MAS stratified by LAM and SCK; RP by all). |

Every table has a free-text `note` column describing its rows.

## Uncertainty in the tables

Uncertain values have a distribution in their own row, in the columns
`dist` and `p1`-`p4`; rows with an empty `dist` are point values. The point
value (`value`, or `ratio` for the three-way term) is what the deterministic
methods use; the distribution is what Monte Carlo runs draw from. Every
distribution type is used at least once:

| `dist` | Used for |
|---|---|
| `beta` | LAM prevalence; MAS:SCK conditional probability |
| `pert` | MAS incidence rate; MET:SCK risk ratio; MAS yield; MET calving interval; MAS welfare |
| `pert_mean` | SCK probability; MET yield |
| `uniform` | MET prevalence; LAM:MET phi; RP welfare; MET:RP interaction |
| `lognormal_ci` | LAM:MAS odds ratio; the three-way ratio |
| `lognormal` | MET:RP odds ratio; RP calving interval |
| `normal` | LAM:SCK odds ratio (truncated at 0); MAS:MET risk difference; LAM yield; LAM and MAS calving interval; LAM and SCK welfare; LAM:MAS interaction |
| `fixed` | SCK yield (a point mass: it cannot be importance-sampled) |

Point values: RP prevalence, the RP:SCK table, the MAS:RP odds ratio, RP
yield, SCK calving interval, MET welfare and the LAM:SCK interaction.
Hazard ratios cannot have distributions.

## What the script covers

1. Reading and checking inputs (`cm_check_inputs`, `cm_read_inputs`, `pair_tables`).
2. Feasibility (`check_feasibility`), the exact and sampled joint
   distribution (`fit_joint`), `combination_probs`.
3. Both methods of `deconflate()` for each analysis, the published
   approximation through `compare_methods()` (for comparison only), and the
   combinations that are not supported (shown as errors or as failed
   methods).
4. Reporting (`contribution_table`, `summary`, `attribute_burden`), and
   gaps and values computed from the results in a few lines of base R.
5. `compare_methods` across analyses.
6. Hazard ratios (`deconflate_hr` with both methods, `compare_methods` with
   the published approach, `attributable_risk`).
7. A batch Monte Carlo run over all analyses on shared population draws.
8. Single-analysis Monte Carlo: several methods, Latin hypercube sampling,
   importance sampling (`cm_suggest_proposal`), `cm_diagnose`,
   `cm_scenario`, `cm_reweight`, gaps and values over the draws in base R,
   and a sampler built in R
   with `cm_dist_table` and a three-way distribution.
9. `sensitivity_oat`, `screen_associations`, `screen_interactions`,
   `screen_three_way`, `compare_scenarios` with the `set_*` functions.
10. `cm_threshold` for every conclusion (`rank`, `sign`, `total`, `change`)
    and input type (association, impact, interaction, three-way ratio,
    disease probability); poles are reported as discontinuities and
    unusable stretches of a range in `regions`.
11. `simulate_raw_impacts` (crude, adjusted, finite samples) and
    `shapley_by_cell` with additive and multiplicative losses.
12. Variants: an unknown pair (global method only), unlisted pairs treated
    as unknown, and a covariate-adjusted risk ratio used as marginal.
13. Plots.

Reference values for this set (computed independently in Python) are in
`inst/validation/reference_five_diseases.py` and are checked by the tests.
