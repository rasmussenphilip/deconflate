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

## v0.5 (released as 0.5.0; planned as v0.4)

Items 3-6 follow a separate design document (October 2026); its reference values come from `inst/validation/prototype_v04_part1.py` to `part3.py` and `reference_v040.py`.

In this order:

1. **Test 0.3.0 in use** (the remaining checks move to v0.6).
   - [x] `run_all_features.R`: optional `out_dir` that saves the printed output and a PDF of all plots.
   - [x] Rewrite `run_all_features.R` for item 3 (one `deconflate()` call per impact table; event impacts; draws inside `deconflate()`; no valuation section; the published method only in `compare_methods()`).
   - [x] Update `run_all_features.R` for item 4.
   - [x] Update `run_all_features.R` for item 5.
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
4. **One model, one input layout, and several estimates per disease, combination and pair.** Prototype: `inst/validation/prototype_v04_item4.py`; reference values: `inst/validation/reference_v041.py` (independent implementation), checked in `test-import.R`, `test-event.R` and `test-item4.R`.
   - [x] One model for the user: the `method` argument goes from `deconflate()` and the sensitivity functions. The package always gives the global model's answer, uses the pairwise shortcut internally when it gives the same answer (recorded in the result), and fits the full joint distribution once for the central estimate (up to about 20 diseases) to confirm the inputs are jointly possible. `compare_methods()` keeps the published approximation beside the package's result.
   - [x] Input layout (the same core columns in every file): `diseases.csv` with `disease`, `measure` (`prevalence`, `probability`, `incidence_rate`), `value`; `associations.csv` with `disease`, `associated_with`, `among`, `adjusted_for`, `measure`, `value`; impact files with `disease`, `among`, `adjusted_for`, `measure` (event impacts only), `value`, `label`, `units`. `disease` and `associated_with` take a disease or a combination (`MAS;SCK`); `among` takes a subgroup (`SCK`, `no SCK`, `SCK;no MET`). Plus `study`, `dist`, `p1`-`p4`, `source`, `note`. A row about two or more diseases belongs in the associations file; the diseases file holds whole-population prevalences only.
   - [x] Association measures: `OR`, `RR`, `RD`, `phi`, `prevalence` (of `disease` in animals with `associated_with`, or of a combination) and `OR_ratio` (the odds ratio among `among` divided by that among its complement). `cond_prob`, `estimand` (decided by `adjusted_for`), the three-way file and the interactions file go; old layouts get messages. Every association row is turned into the probability of one event of disease combinations, and the joint distribution is the maximum-entropy one with all of them (raking on events, with a fixed point for conversions that need other probabilities); for pairs only this is the current fit. Rows that cannot all hold are reported.
   - [x] Impacts of combinations and subgroups: crude or adjusted rows for a combination (`SCK;SCM`) or within a subgroup (`among`). Each disease's own impact and an interaction for each combination listed as a `disease` are estimated; a combination adjusted for `all` is its interaction. Identification check (rows vs unknowns, and which unknowns are not identified). Shapley contributions split each interaction equally among its diseases.
   - [x] Event impacts: combinations as interaction terms on the log-hazard scale, estimated from combination rows (crude, stratified).
   - [x] Sensitivity tools: kept working on the new layout with as few changes as possible (no new features); their redesign is under "To discuss".
   - [x] Several rows per item: impacts fitted by weighted least squares (weights 1/SD²; a row without uncertainty is exact; residuals and heterogeneity reported); several association rows about the same quantity pooled on the log scale (random effects by default, `pool = "fixed"`; Q, I² and tau²); associations adjusted for diseases in the model are exact conditional associations, those adjusted for other factors are used as crude, with a note.
   - [x] Tests, vignettes (`vignette("interactions")` rewritten around combination and subgroup rows), the five-disease example and `run_all_features.R` updated.
5. **Diagnosing the inputs.** `deconflate(diagnose = TRUE)` (off by default) adds `$problems`: which inputs drive problematic adjusted impacts, and rows that disagree with other rows. Design: section 14 of the v0.4 design note; reference values: `inst/validation/reference_v05_input_check.py`.
   - [x] Flagged impacts: an adjusted impact on the other side of zero from its raw impact (for event impacts, of 1), or shrunk by more than 90%.
   - [x] Inputs with a distribution (prevalences, associations, impacts, the overall risk): the change in each flagged adjusted impact when the input moves by one SD, and the shift (in SD) that alone would remove the flag, shown when within 3 SD.
   - [x] Inputs without a distribution: the value at which the flag disappears, searched from a quarter to three times the value.
   - [x] A separate table of rows that disagree with other rows for the same quantity (standardised residual of 2 or more), with the adjusted impacts without each.
   - [x] Printed with the result; the sign-change warning points to it. Additive and event impacts.
6. **Dropped (October 2026): the `anchor` argument, constraints on adjusted impacts, MCMC and importance sampling.** A Python prototype showed that they cannot fix what they were meant to fix: a negative adjusted impact means the inputs contradict each other, a lower bound of zero only parks the impact at zero, and narrowing the distributions makes the negative values more certain. Such conflicts are fixed in the inputs (item 5 shows where). Probabilistic thresholds move to the sensitivity redesign (To discuss).
- [x] Run `devtools::document()`, `test()` and `check()` (R-CMD-check passed on GitHub).
- [x] Tag v0.5.0.

## v0.6 (planned)

0. **Test 0.5.0 in use.**
   - [ ] Run `run_all_features.R` end to end and review every printed result as a new user would.
   - [ ] Have someone else work through the CSV workflow with their own data, without help.
   - [ ] Optional: a short external review of the 0.5 additions (rows about combinations and subgroups, several rows per item, `diagnose = TRUE`).
1. **Event impacts matched to additive impacts.** Design: section 15 of the v0.4 design note; reference values: `inst/validation/reference_v06_event_parity.py`.
   - [x] One raw result for both models: the same model with the diseases unrelated (same prevalences and impact rows, no association rows). `raw_sum` becomes `raw_total` (additive), and event results gain `raw_attributable`.
   - [x] Contributions split into each disease's own part and its share of the interactions (`main`, `interaction`), for additive and event impacts.
   - [x] The local row equations of event impacts (`$conflation$A`, at the adjusted values).
   - [x] Results on the risk scale: risk, risk ratio, risk difference and odds ratio of an animal with only the disease against a disease-free animal; `change` and the sign check on the raw row's own measure.
   - [x] `baseline_risk` (the risk of disease-free animals) as an alternative to `overall_risk`; exactly one is given.
   - [x] Simulated raw estimates of event impacts (`simulate_raw_impacts(event_model = TRUE)`), exact or from sampled animals (rate ratios for HR rows).
   - [x] Interactions of event impacts in `set_interaction()`, `screen_interactions()` and `cm_threshold()`.
   - [x] Run `devtools::document()`, `test()` and `check()`.
2. **One model only.** Remove everything that is not the package's model (the sensitivity tools wait for their redesign).
   - [x] Remove `compare_methods()`, the published approximation (eq. 16, HR - 1), the first-order event method, `reproduce_rasmussen_2022()`/`_2024()`, `example_uk_dairy_2022()`, the reproduction vignette and the reproduction-only code.
   - [x] Remove helpers nothing uses (`attribute_burden()`, `plot_burden()`, `loss_*()`, `combination_probs()`, `excess_matrix()`, `cm_dist_table()`), and stop exporting the internals of the model (`or_to_joint()`, `joint_to_or()`, `association_to_joint()`, `pair_tables()`, `shapley_by_cell()`).
   - [x] Run `devtools::document()`, `test()` and `check()`.
3. **Redesign of the sensitivity tools.** Design: section 17 of the v0.4 design note; reference values: `inst/validation/reference_v06_sensitivity.py`.
   - [x] `cm_sensitivity()`: entered and unknown inputs (pairs, interactions, three-way on request) one at a time, in one table; ranges from each input's SD, a relative change, or defaults for unknowns; variance shares from the draws; print (top 10) and tornado plot.
   - [x] Rank probabilities from the draws (`$ranks`, in `summary()`).
   - [x] `cm_threshold()` takes a result and the same input keys; the diagnosis shares the engine.
   - [x] Remove `sensitivity_oat()` and the three screens; `set_*()` internal.
   - [x] Run `devtools::document()`, `test()` and `check()`.
7. **Revised 2024 inputs**, beside the inputs of the 2024 analysis (`example_global_dairy()`).
   - [ ] `example_global_dairy(inputs = "revised")`, using the features of items 4-5 (several rows per disease with their own estimands, e.g. the clinical mastitis sources).
   - [ ] A vignette section documenting each correction: the SCK conversion from kg to percent, the milk fever and ketosis estimates of Bareille et al. (2003) as direct effects, the adjusted SCK estimates and odds ratios of Raboisson et al. (2014), the clinical mastitis sources by estimand, and the displaced abomasum sources. Unit and measurement-window corrections are made here only, not in the main pipeline.
8. **Simpler, more intuitive input tables.** Design: section 18 of the v0.4 design note.
   - [x] Accept missing optional columns everywhere with sensible defaults, and fewer required columns (only `disease` and `value` are required).
   - [x] Adding a row to a model in R (e.g. an interaction or an `OR_ratio` row, now that the `set_*()` functions are internal) does not need renumbering the `row` column: `cm_model()` numbers rows joined with `rbind()`.
   - [x] Fewer columns: `label` (diseases), `time_horizon`, `reference_population` and `source` are no longer read (a note when present); `note` is the only free-text column. Disease labels were dropped: printed tables show the codes.
   - [x] Clearer messages; the template and the example folders reviewed accordingly; the documentation says that all inputs refer to the same period.
   - [x] Make uncertainty visible: `deconflate()` notes when the inputs have distributions but `n_draws = 0`.
   - [ ] Run `devtools::document()`, `test()` and `check()`.
9. **Clearer output and summary tables.** The package assumes no species, outcome or unit: tables show the user's own disease codes, outcome label and units.
   - [x] Consistent column names, order, rounding (3 significant figures) and layout across `print()`, `summary()` and `contribution_table()`, for both models, with and without draws; tables fit 80 columns and do not wrap into several blocks.
   - [x] The same rules for `cm_sensitivity()`, `cm_threshold()` and the diagnosis.
   - [x] Monte Carlo notes in plain words.
   - [ ] Run `devtools::document()`, `test()` and `check()`.
   - [x] Raw vs adjusted attributable risk for mortality files (the same calculation with the raw estimates), so that the double counting removed is visible, as with the raw and adjusted aggregates of impact files (done in item 1).
10. **Better plots.** The same rule as for tables (the user's labels and units only); a consistent look across plot methods, and plots for uncertainty, thresholds and the estimate diagnostics.
    - [ ] Documentation examples from more than dairy cattle (e.g. another livestock system and a human-health example), so that the generality is visible.
11. **Dependencies.** Keep the package on base R (see the Dependencies section of the README); list any new dependency there with its reason.

Afterwards: a Zenodo DOI for each release, and possibly a CRAN submission and a short software paper.

## To discuss

- Probabilistic thresholds (with the Monte Carlo draws, e.g. the input value at which a disease ranks first in a given share of draws), now that `cm_sensitivity()` and `$ranks` exist (v0.6 item 3).
- Pairwise-only bounds as a complementary sensitivity output for large systems; pseudo-likelihood fitting when individual-level disease records are available.
- Higher-order and non-linear interaction models in the adjustment itself.

Dropped (October 2026): advanced component accounting, and lost productive life as a culling outcome.
