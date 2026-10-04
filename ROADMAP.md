# Roadmap

## v0.1 (released as 0.1.0)

Input constructors with metadata.

Association measures: OR, RR, RD, conditional probability, phi,
contingency table, independent and unknown.

Three adjustment methods: published, simultaneous (default) and global
(IPF/MaxEnt with pairwise interactions).

Diagnostics: reconstruction residuals, sign changes, condition numbers,
IPF convergence.

Feasibility check: a triple screen, plus an exact LP with `lpSolve`.

Productivity gaps on the proportion and absolute scales, monetary
valuation, and Shapley attribution (closed form, and cell by cell for
general losses).

Culling end to end: hazard ratio → excess risk → adjustment → hazard
ratio.

Reproduction of Rasmussen et al. (2022) Tables 8–10 and (2024) Table 5,
with documented differences.

Input distributions, samplers with outcome correlation, Monte Carlo with
economics and rejection reporting, and scenario reweighting.

Sensitivity: one-at-a-time, association screening, interaction screening
and scenario comparison.

Reporting and base-graphics plots.

Vignettes and pkgdown configuration.

Run `devtools::document()`, `test()` and `check()` on 0.0.0.9001 and fix
any issues.

GitHub Actions: `usethis::use_github_action("check-standard")` and
`usethis::use_pkgdown_github_pages()`.

Review the reproduction vignette against the original analysis files
(2024: 1st-revision code and inputs; 2022: original files unavailable,
corrigendum checked).

Tag v0.1.0.

## v0.2 (released as 0.2.0)

Culling hazard ratios: a separate multiplicative adapter
([`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md):
snapshot, first-order and published) with stratified adjustment sets.

[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md):
culling attributable to disease without double counting, with Shapley
allocation.

Reading inputs from CSV files or data frames, with all problems reported
at once; template files; one impact file per analysis.

Method comparison tables for models, analyses, hazard-ratio models and
Monte Carlo runs.

Monte Carlo stability diagnostics with suggestions; Latin hypercube
(replicate blocks) and importance sampling (support checks).

External review (Bob): one impact vector per analysis in arbitrary
units; explicit estimands with the exact projection mapping for
additive-regression coefficients; feasibility screening by default;
condition classes; three-way scenarios; complete Shapley allocation;
reproduction functions; the reviewer’s regression tests.

Run `devtools::document()`, `test()` and `check()`, and fix any issues.

Tag v0.2.0.

## v0.3 (released as 0.3.0)

Second external review: importance-sampling support as intervals, one
finiteness check everywhere, `insufficient_info` precision status,
three-way documentation, converged joints for simulation,
single-analysis template, explicit snapshot hazard-ratio estimands,
historical conversions confined to reproduction.

Threshold searches
([`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)):
rank, sign, target and relative-change conclusions for one input at a
time.

Sampled backend for the joint distribution (Monte Carlo moment matching,
sampling, raking) beyond about 20 diseases.

Uncertainty given in the input tables (`dist`, `p1`-`p4`), with no
separate uncertainty files; three-way ratios can be uncertain.

Five-disease example that uses every feature
(`inst/extdata/five_diseases`, with `run_all_features.R`).

Run `devtools::document()`, `test()` and `check()` (0 errors, 0
warnings, 0 notes).

Tag v0.3.0.

## v0.4 (planned)

In this order:

**Test 0.3.0 in use.**

Run `run_all_features.R` end to end and review every printed result as a
new user would.

Have someone else work through the CSV workflow with their own data,
without help.

Optional: a short third review of the 0.3 additions (in-table
uncertainty,
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md),
the sampled backend).

`run_all_features.R`: optional `out_dir` that saves the printed output
and a PDF of all plots.

Update `run_all_features.R` after items 2-3 (no valuation section; the
published method only where it remains; the new features).

**Slim the interface to de-conflation itself.**

Demote the published method (eq. 16 of Rasmussen et al. 2022, and HR - 1
for hazard ratios):

Keep it in
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md),
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
(to reproduce the 2024 Monte Carlo) and the `reproduce_*()` functions
only.

Remove it from the method choices of
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
and
[`deconflate_hr()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate_hr.md),
from
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
and from the screens; simplify the pole handling that exists mainly for
it.

Documentation: use the simultaneous method when every pair has an
association estimate or a defensible independence assumption and impacts
are additive; otherwise the global method; the snapshot model for
hazards.

Remove valuation; users compute gaps and values from the adjusted
impacts themselves:

Remove
[`productivity_gap()`](https://rasmussenphilip.github.io/deconflate/reference/productivity_gap.md),
[`value_losses()`](https://rasmussenphilip.github.io/deconflate/reference/value_losses.md),
[`cm_mc_gap()`](https://rasmussenphilip.github.io/deconflate/reference/cm_mc_gap.md)
and the exported
[`uk_dairy_2022_economics()`](https://rasmussenphilip.github.io/deconflate/reference/uk_dairy_2022_economics.md).

Remove every `valuation` argument
([`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md),
[`sensitivity_oat()`](https://rasmussenphilip.github.io/deconflate/reference/sensitivity_oat.md),
the screens,
[`compare_scenarios()`](https://rasmussenphilip.github.io/deconflate/reference/compare_scenarios.md),
[`plot_burden()`](https://rasmussenphilip.github.io/deconflate/reference/plots.md));
sensitivity tools use the adjusted aggregate in the impacts’ own units.

Remove `unit_value` from
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md)
(keep the attributable risk itself: it needs the joint distribution).

Keep the adjusted aggregate and each disease’s contribution, including
the Shapley split of interaction effects.

Move the gap and valuation code inside
[`reproduce_rasmussen_2022()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
and
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md)
as internal helpers, so the published economic tables are still
reproduced.

A short vignette section showing the gap and value calculation in base R
(for percentage impacts, the loss is relative to the disease-free level:
observed / (1 - aggregate)).

**Several estimates per disease and per pair**, each with its own
`estimand`/`adjusted_for` and its own uncertainty.

Impact tables: several rows per disease (e.g. a crude and an adjusted
estimate from different studies).

Association tables: several rows per pair, converted to a common scale
and pooled (random effects, heterogeneity reported); adjusted odds
ratios enter as conditional associations (exactly in the global method,
approximately in the simultaneous method).

Mortality and culling tables: a `measure` column (`HR`, `rate_ratio`,
`RR`, `OR`, `RD`), crude or stratified, all mapped onto the snapshot
hazard model; the overall period risk is required when risk-based
measures are used (with hazard and rate ratios only, it is needed only
for the attributable risk), and risk-based measures must refer to the
same period.

Uncertainty for the mortality model: distributions for mortality
estimates (and optionally the overall risk);
[`cm_monte_carlo()`](https://rasmussenphilip.github.io/deconflate/reference/cm_monte_carlo.md)
runs the snapshot de-conflation and the attributable risk in every draw,
with intervals, rejections and stability checks as for additive
analyses.

**A shared estimator with an informativeness toggle**, used by the
simultaneous and global methods and the hazard model.

Fits all estimates jointly, each weighted by its uncertainty; the toggle
fixes the associations (`"associations"`, the default and current
behaviour), weights both kinds of evidence (`"balanced"`), or fixes the
impacts (`"impacts"`).

Diagnostics: a standardised residual for each estimate, and how far each
association moved in SD units, so that estimates that do not fit the
others stand out.

Optional constraints on adjusted impacts (e.g. non-negative) applied
probabilistically, reporting posterior means and intervals and how often
each constraint was violated without it; the best-fit version (binding
impacts exactly at the bound) only as a fast diagnostic of which inputs
conflict.

**Probabilistic machinery** (MCMC or directed importance sampling,
shared by the constraints in item 4).

Probabilistic thresholds: combine
[`cm_threshold()`](https://rasmussenphilip.github.io/deconflate/reference/cm_threshold.md)
with Monte Carlo runs, e.g. the input value at which a disease ranks
first (or the aggregate exceeds a target) in a given share of draws,
with Monte Carlo error.

**Revised 2024 inputs**, separate from
[`reproduce_rasmussen_2024()`](https://rasmussenphilip.github.io/deconflate/reference/reproduce.md),
which stays a faithful reproduction.

`example_global_dairy(inputs = "revised")`, using the features of items
3-4.

A vignette section documenting each correction: the SCK conversion from
kg to percent, the milk fever and ketosis estimates of Bareille et
al. (2003) as direct effects, the adjusted SCK estimates and odds ratios
of Raboisson et al. (2014), the clinical mastitis sources by estimand,
and the displaced abomasum sources. Unit and measurement-window
corrections are made here only, not in the main pipeline.

**Simpler, more intuitive input tables.**

Accept missing optional columns everywhere with sensible defaults, and
fewer required columns.

Disease labels (full names) next to the short ids, used in printed
tables and plots.

Clearer column names and messages; review the template and the example
folders accordingly (more important now that tables hold several
estimates).

Make uncertainty visible:
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
notes when the inputs have distributions, and an argument
(e.g. `n_draws`) runs the Monte Carlo analysis in the same call,
returning point estimates with intervals.

**Clearer output and summary tables.** The package assumes no species,
outcome or unit: tables show the user’s own labels and units, and the
reference population and time horizon from the diseases table when given
(otherwise “per average individual in the population”).

Consistent column names, order, rounding and layout across
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md),
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`contribution_table()`](https://rasmussenphilip.github.io/deconflate/reference/contribution_table.md)
and the Monte Carlo summaries; wide tables do not wrap into several
blocks.

Raw vs adjusted attributable risk for the mortality model (the same
calculation with the raw estimates), so that the double counting removed
is visible, as with the raw and adjusted aggregates of additive
analyses.

**Better plots.** The same rule as for tables (the user’s labels and
units only); a consistent look across plot methods, and plots for
uncertainty, thresholds and estimator diagnostics.

Documentation examples from more than dairy cattle (e.g. another
livestock system and a human-health example), so that the generality is
visible.

**Dependencies.** Keep the package on base R (see the Dependencies
section of the README); list any new dependency there with its reason
(the MCMC in item 5 is the most likely candidate).

Afterwards: a Zenodo DOI for each release, and possibly a CRAN
submission and a short software paper.

## To discuss

- Pairwise-only bounds as a complementary sensitivity output for large
  systems; pseudo-likelihood fitting when individual-level disease
  records are available.
- Higher-order and non-linear interaction models in the adjustment
  itself.

Dropped (October 2026): advanced component accounting, and lost
productive life as a culling outcome.
