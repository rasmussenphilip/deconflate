# Roadmap

## v0.1 (released as 0.1.0)

- [x] Input constructors with metadata.
- [x] Association measures: OR, RR, RD, conditional probability, phi, contingency table, independent and unknown.
- [x] Three adjustment methods: published, simultaneous (default) and global (IPF/MaxEnt with pairwise interactions).
- [x] Diagnostics: reconstruction residuals, sign changes, condition numbers, IPF convergence.
- [x] Feasibility check: a triple screen, plus an exact LP with `lpSolve`.
- [x] Productivity gaps on the proportion and absolute scales, monetary valuation, and Shapley attribution (closed form, and cell by cell for general losses).
- [x] Culling end to end: hazard ratio → excess risk → adjustment → hazard ratio.
- [x] Reproduction of Rasmussen et al. (2022) Tables 8–10 and (2024) Table 5, with documented differences.
- [x] Input distributions, samplers with outcome correlation, Monte Carlo with economics and rejection reporting, and scenario reweighting.
- [x] Sensitivity: one-at-a-time, association screening, interaction screening and scenario comparison.
- [x] Reporting and base-graphics plots.
- [x] Vignettes and pkgdown configuration.
- [x] Run `devtools::document()`, `test()` and `check()` on 0.0.0.9001 and fix any issues.
- [x] GitHub Actions: `usethis::use_github_action("check-standard")` and `usethis::use_pkgdown_github_pages()`.
- [x] Review the reproduction vignette against the original analysis files (2024: 1st-revision code and inputs; 2022: original files unavailable, corrigendum checked).
- [x] Tag v0.1.0.

## v0.2 (released as 0.2.0)

- [x] Culling hazard ratios: a separate multiplicative adapter (`deconflate_hr()`: snapshot, first-order and published) with stratified adjustment sets.
- [x] `attributable_risk()`: culling attributable to disease without double counting, with Shapley allocation.
- [x] Reading inputs from CSV files or data frames, with all problems reported at once; template files; one impact file per analysis.
- [x] Method comparison tables for models, analyses, hazard-ratio models and Monte Carlo runs.
- [x] Monte Carlo stability diagnostics with suggestions; Latin hypercube (replicate blocks) and importance sampling (support checks).
- [x] External review (Bob): one impact vector per analysis in arbitrary units; explicit estimands with the exact projection mapping for additive-regression coefficients; feasibility screening by default; condition classes; three-way scenarios; complete Shapley allocation; reproduction functions; the reviewer's regression tests.
- [x] Run `devtools::document()`, `test()` and `check()`, and fix any issues.
- [x] Tag v0.2.0.

## v0.3 (released as 0.3.0)

- [x] Second external review: importance-sampling support as intervals, one finiteness check everywhere, `insufficient_info` precision status, three-way documentation, converged joints for simulation, single-analysis template, explicit snapshot hazard-ratio estimands, historical conversions confined to reproduction.
- [x] Threshold searches (`cm_threshold()`): rank, sign, target and relative-change conclusions for one input at a time.
- [x] Sampled backend for the joint distribution (Monte Carlo moment matching, sampling, raking) beyond about 20 diseases.
- [x] Uncertainty given in the input tables (`dist`, `p1`-`p4`), with no separate uncertainty files; three-way ratios can be uncertain.
- [x] Five-disease example that uses every feature (`inst/extdata/five_diseases`, with `run_all_features.R`).
- [x] Run `devtools::document()`, `test()` and `check()` (0 errors, 0 warnings, 0 notes).
- [x] Tag v0.3.0.

## v0.4 (planned)

Items 3-6 follow a separate design document (October 2026); its reference values come from `inst/validation/prototype_v04_part1.py` to `part3.py` and `reference_v040.py`.

In this order:

1. **Test 0.3.0 in use.**
   - [ ] Run `run_all_features.R` end to end and review every printed result as a new user would.
   - [ ] Have someone else work through the CSV workflow with their own data, without help.
   - [ ] Optional: a short third review of the 0.3 additions (in-table uncertainty, `cm_threshold()`, the sampled backend).
   - [x] `run_all_features.R`: optional `out_dir` that saves the printed output and a PDF of all plots.
   - [x] Rewrite `run_all_features.R` for item 3 (one `deconflate()` call per impact table; event impacts; draws inside `deconflate()`; no valuation section; the published method only in `compare_methods()`).
   - [ ] Update `run_all_features.R` again after items 4-6.
2. **Slim the interface to de-conflation itself.**

   Demote the published method (eq. 16 of Rasmussen et al. 2022, and HR - 1 for hazard ratios):
   - [x] Keep it in `compare_methods()`, `cm_monte_carlo()` (to reproduce the 2024 Monte Carlo) and the `reproduce_*()` functions only.
   - [x] Remove it from the method choices of `deconflate()` and `deconflate_hr()`, from `cm_threshold()` and from the screens; simplify the pole handling that exists mainly for it.
   - [x] Documentation: use the simultaneous method when every pair has an association estimate or a defensible independence assumption and impacts are additive; otherwise the global method; the snapshot model for hazards.

   Remove valuation; users compute gaps and values from the adjusted impacts themselves:
   - [x] Remove `productivity_gap()`, `value_losses()`, `cm_mc_gap()` and the exported `uk_dairy_2022_economics()`.
   - [x] Remove every `valuation` argument (`contribution_table()`, `summary()`, `compare_methods()`, `sensitivity_oat()`, the screens, `compare_scenarios()`, `plot_burden()`); sensitivity tools use the adjusted aggregate in the impacts' own units.
   - [x] Remove `unit_value` from `attributable_risk()` (keep the attributable risk itself: it needs the joint distribution).
   - [x] Keep the adjusted aggregate and each disease's contribution, including the Shapley split of interaction effects.
   - [x] Move the gap and valuation code inside `reproduce_rasmussen_2022()` and `reproduce_rasmussen_2024()` as internal helpers, so the published economic tables are still reproduced.
   - [x] A short vignette section showing the gap and value calculation in base R (for percentage impacts, the loss is relative to the disease-free level: observed / (1 - aggregate)).
3. **One function, one impact table.** Everything runs through `deconflate()`.
   - [x] Every input table is a named argument of `cm_read_inputs()` with any file name (`diseases`, `associations`, `impacts`; optionally `interactions` and `three_way`). One impact table per `deconflate()` call; users repeat the call for other outcomes. `cm_analyses()`, multi-analysis reading, `dir` and the fixed file names are removed.
   - [x] Associations: `deconflate()` needs at least one (the sensitivity tools also run without any). Pairs without a row are unknown by default (filled in by the global fit and listed in `$unknown_pairs`); an odds ratio of 1 means unrelated. `missing_associations` and the measures `independent`, `unknown` and `table` are removed, with messages for old files. Associations keep their uncertainty.
   - [x] `method = "auto"` (default): simultaneous when it equals global, otherwise global with a note (interactions, three-way terms, unknown pairs); `"simultaneous"` switches the same way; the sampled backend above 20 diseases.
   - [x] Event impacts through `deconflate(event_model = TRUE, overall_risk = )`: a `measure` column with `HR`, `rate_ratio`, `RR`, `OR` and `RD` (mixable), crude or stratified, all mapped onto the snapshot hazard model; `overall_risk` required (a proportion or a distribution); the attributable risk and its Shapley split in the result. Mismatches between the table and `event_model` are errors. `hazard_ratios.csv` becomes `culling.csv` in the examples.
   - [x] `n_draws` (default 1000; skipped with a note when no input has a distribution): intervals for additive and event impacts, rejections (not replaced) and stability checks, `seed`, Latin hypercube sampling. Uncertainty of event impacts and of the overall risk.
   - [x] Removed from the public interface (kept internally where `reproduce_rasmussen_2024()` needs them): `cm_monte_carlo()`, `cm_sampler()`, `cm_batch_sampler()`, `sampler_global_dairy()`, `deconflate_hr()`, `cm_hr_model()`, `cm_hazard_ratios()`, `attributable_risk()`, `example_global_dairy_hr()`, `cm_analyses()`; removed: `cm_reweight()`, `cm_scenario()`, `compare_scenarios()`, `cm_diagnose()`, importance sampling and `cm_suggest_proposal()`.
   - [x] `compare_methods()` (additive and event impacts, with `n_draws`), `cm_threshold()`, the screens and `sensitivity_oat()` (with `event_model`, point estimates only), the plots, the examples (one model per outcome), the vignettes and `run_all_features.R` updated to the new interface.
4. **Several estimates per disease and per pair**, each with its own estimand, `adjusted_for` and uncertainty.
   - [ ] Impact tables: several rows per disease (optional `study` label), fitted by weighted least squares with weights 1/SD²; a row without uncertainty is exact. Error only when exact rows contradict each other and the anchor does not let the associations reconcile them. For event impacts, several rows per disease on the snapshot model (all measures already supported, one row per disease, in item 3).
   - [ ] Association tables: several rows per pair, mapped to the marginal log odds ratio and pooled (random effects by default, `pool = "fixed"` as an option; Q, I² and tau² reported). Associations adjusted for diseases in the model (ids or `all`) are exact conditional associations (logistic projection, via a fixed-point mapping to the marginal scale); associations adjusted for other covariates are used as marginal, with a note.
5. **The `anchor` argument and constraints.**
   - [ ] `anchor = c("associations", "balanced", "impacts")`: what is held at its input value when the inputs over-determine the model (several impact rows per disease, or constraints). `"associations"` (default, the current behaviour) fits the impacts; `"balanced"` moves both kinds of input in proportion to their SDs; `"impacts"` moves the associations as little as needed. Without over-determination all three give the same results as 0.3.
   - [ ] Diagnostics in the result: a standardised residual for each input row, heterogeneity per disease and pair, and how far each input moved in SD units (shifts of 2 SD or more flagged).
   - [ ] Constraints on adjusted impacts as a `deconflate()` argument only (e.g. `constraints = list(lower = 0)` or `list(lower = c(CM = 0, DA = 0))`; no CSV columns), applied probabilistically with `n_draws`: draw combinations that break them are dropped and the rest weighted according to `anchor`. Reports means and intervals (never exactly at the bound), how often each constraint was broken without it, and the input shifts.
   - [ ] When few draws satisfy the constraints, return the results anyway with a warning, together with the best-fit diagnostic in the result (`$conflict`): the smallest shift of the inputs, in SD units, that satisfies the constraints, and the inputs that have to move most.
6. **Probabilistic machinery** (MCMC or directed importance sampling), for where too few draw combinations satisfy the constraints, constraints on the mortality model, and `anchor = "impacts"` with several impact rows per disease.
   - [ ] Probabilistic thresholds: combine `cm_threshold()` with the Monte Carlo draws, e.g. the input value at which a disease ranks first (or the aggregate exceeds a target) in a given share of draws, with Monte Carlo error.
7. **Revised 2024 inputs**, separate from `reproduce_rasmussen_2024()`, which stays a faithful reproduction.
   - [ ] `example_global_dairy(inputs = "revised")`, using the features of items 4-5.
   - [ ] A vignette section documenting each correction: the SCK conversion from kg to percent, the milk fever and ketosis estimates of Bareille et al. (2003) as direct effects, the adjusted SCK estimates and odds ratios of Raboisson et al. (2014), the clinical mastitis sources by estimand, and the displaced abomasum sources. Unit and measurement-window corrections are made here only, not in the main pipeline.
8. **Simpler, more intuitive input tables.**
   - [ ] Accept missing optional columns everywhere with sensible defaults, and fewer required columns.
   - [ ] Disease labels (full names) next to the short ids, used in printed tables and plots.
   - [ ] Clearer column names and messages; review the template and the example folders accordingly (more important now that tables hold several estimates).
   - [ ] Make uncertainty visible: `deconflate()` notes when the inputs have distributions but `n_draws = 0`.
9. **Clearer output and summary tables.** The package assumes no species, outcome or unit: tables show the user's own labels and units, and the reference population and time horizon from the diseases table when given (otherwise "per average individual in the population").
   - [ ] Consistent column names, order, rounding and layout across `deconflate()`, `compare_methods()`, `summary()` and `contribution_table()`, with and without draws; wide tables do not wrap into several blocks.
   - [ ] Raw vs adjusted attributable risk for mortality files (the same calculation with the raw estimates), so that the double counting removed is visible, as with the raw and adjusted aggregates of impact files.
10. **Better plots.** The same rule as for tables (the user's labels and units only); a consistent look across plot methods, and plots for uncertainty, thresholds and the estimate diagnostics.
    - [ ] Documentation examples from more than dairy cattle (e.g. another livestock system and a human-health example), so that the generality is visible.
11. **Dependencies.** Keep the package on base R (see the Dependencies section of the README); list any new dependency there with its reason (the MCMC in item 6 is the most likely candidate).

Afterwards: a Zenodo DOI for each release, and possibly a CRAN submission and a short software paper.

## To discuss

- Pairwise-only bounds as a complementary sensitivity output for large systems; pseudo-likelihood fitting when individual-level disease records are available.
- Higher-order and non-linear interaction models in the adjustment itself.

Dropped (October 2026): advanced component accounting, and lost productive life as a culling outcome.
