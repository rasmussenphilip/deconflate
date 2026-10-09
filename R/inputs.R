#' Describe the diseases in a system
#'
#' @param id Character vector of unique disease identifiers (no `|`, `;` or
#'   `:`; `"all"` is reserved).
#' @param value Numeric vector of disease occurrence values. Interpreted
#'   according to `type`.
#' @param type Character, one per disease (or recycled):
#'   * `"prevalence"`: proportion of animals affected at a representative
#'     point in time (used as the probability directly).
#'   * `"probability"`: probability of at least one case within
#'     `time_horizon` (used directly).
#'   * `"incidence_rate"`: expected number of cases per animal within
#'     `time_horizon` (e.g. `0.48` for 48 cases per 100 lactations). Converted
#'     to a probability assuming Poisson-distributed events,
#'     `1 - exp(-value)`, as in Rasmussen et al. (2024), eq. 1.
#' @param time_horizon Character label of the period the probabilities refer
#'   to (e.g. `"lactation"`, `"year"`). Probabilities, associations and impacts
#'   should all refer to the same period.
#' @param reference_population Optional description of the population the
#'   values come from.
#' @param source Optional citation.
#'
#' @return A `cm_diseases` data frame with the converted probability in
#'   column `prob` and the inputs retained as metadata.
#' @export
#' @examples
#' cm_diseases(c("SCK", "PTB"), c(0.4789, 0.1001),
#'             type = c("incidence_rate", "prevalence"),
#'             time_horizon = "lactation")
cm_diseases <- function(id, value, type = "prevalence",
                        time_horizon = NA_character_,
                        reference_population = NA_character_,
                        source = NA_character_) {
  id <- as.character(id)
  n <- length(id)
  if (n < 1L) cm_abort("At least one disease is required.")
  check_ids(id)
  if (length(value) != n) cm_abort("`value` must have one entry per disease.")
  check_numeric(value, "value")
  if (any(value < 0)) cm_abort("`value` must be non-negative.")
  type <- as.character(recycle_arg(type, n, "type"))
  check_choices(type, c("prevalence", "probability", "incidence_rate"), "type")

  prob <- ifelse(type == "incidence_rate", 1 - exp(-value), value)
  if (any(prob <= 0 | prob >= 1)) {
    cm_abort(sprintf("Disease probabilities must lie strictly between 0 and 1 (check: %s).",
                     paste(id[prob <= 0 | prob >= 1], collapse = ", ")),
             class = "deconflate_infeasible")
  }
  horizon <- as.character(recycle_arg(time_horizon, n, "time_horizon"))
  if (length(unique(stats::na.omit(horizon))) > 1L) {
    cm_warn("Diseases refer to different time horizons; probabilities, associations and impacts should share one period.")
  }

  out <- data.frame(
    id = id, prob = prob, value = value, type = type,
    time_horizon = horizon,
    reference_population = as.character(recycle_arg(reference_population, n,
                                                    "reference_population")),
    source = as.character(recycle_arg(source, n, "source")),
    stringsAsFactors = FALSE
  )
  class(out) <- c("cm_diseases", "data.frame")
  out
}

#' Describe statistical associations between disease pairs
#'
#' Each row gives one association measure for one pair of diseases. A pair
#' without a row is unknown: the global model ([fit_joint()]) fills in its
#' association from the others (see [deconflate()]). To state that two
#' diseases are unrelated, give the pair an odds ratio of 1.
#'
#' Odds ratios (and the other measures) are applied to the modelled
#' population's own marginal probabilities ("transported"). This assumes the
#' measure is the same in the source study and in the modelled population.
#'
#' @param disease1,disease2 Character vectors of disease ids.
#' @param value Numeric association values (required).
#' @param measure Character, one per row (or recycled). Directional measures
#'   treat `disease1` as the outcome and `disease2` as the conditioning
#'   disease:
#'   * `"OR"` (default): odds ratio (symmetric). An odds ratio of 1 means the
#'     two diseases are independent.
#'   * `"RR"`: risk ratio `P(d1 | d2) / P(d1 | not d2)`.
#'   * `"RD"`: risk difference `P(d1 | d2) - P(d1 | not d2)`.
#'   * `"cond_prob"`: conditional probability `P(d1 | d2)`.
#'   * `"phi"`: binary (phi) correlation coefficient (symmetric).
#'
#'   For a study's 2x2 table, compute its odds ratio,
#'   `n11 * n00 / (n10 * n01)`, and enter it as an odds ratio.
#' @param adjusted Logical: was the measure adjusted for covariates (e.g. an
#'   odds ratio from multivariable logistic regression)? The 2x2 algebra needs
#'   marginal (crude) measures, so [cm_population()] rejects adjusted measures
#'   unless told to use them as marginal.
#' @param adjusted_for Optional description of the adjustment set.
#' @param source Optional citation.
#'
#' @return A `cm_associations` data frame.
#' @export
#' @examples
#' cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3), measure = "OR")
#' # Two diseases stated to be unrelated
#' cm_associations("d1", "d3", 1)
cm_associations <- function(disease1, disease2, value, measure = "OR",
                            adjusted = FALSE, adjusted_for = NA_character_,
                            source = NA_character_) {
  disease1 <- as.character(disease1)
  disease2 <- as.character(disease2)
  n <- length(disease1)
  if (length(disease2) != n) cm_abort("`disease1` and `disease2` must have equal length.")
  if (any(disease1 == disease2)) cm_abort("An association must involve two different diseases.")
  measure <- as.character(recycle_arg(measure, n, "measure"))
  old <- intersect(measure, names(retired_measures))
  if (length(old)) cm_abort(retired_measures[[old[1]]], class = "deconflate_unsupported")
  check_choices(measure, association_measures, "measure")
  if (missing(value)) cm_abort("Give each association a `value`.")
  value <- as.numeric(recycle_arg(value, n, "value"))
  if (anyNA(value)) {
    cm_abort("Association `value` is missing for some rows. Leave a pair out (no row) if its association is unknown.")
  }
  if (any(!is.finite(value))) cm_abort("Association values must be finite.")
  if (any(value[measure %in% c("OR", "RR")] <= 0)) {
    cm_abort("Odds ratios and risk ratios must be positive.", class = "deconflate_infeasible")
  }
  key <- pair_key(disease1, disease2)
  if (anyDuplicated(key)) {
    cm_abort(sprintf("Duplicate associations for pair(s): %s.",
                     paste(unique(key[duplicated(key)]), collapse = ", ")))
  }
  out <- data.frame(
    disease1 = disease1, disease2 = disease2, measure = measure, value = value,
    adjusted = as.logical(recycle_arg(adjusted, n, "adjusted")),
    adjusted_for = as.character(recycle_arg(adjusted_for, n, "adjusted_for")),
    source = as.character(recycle_arg(source, n, "source")),
    stringsAsFactors = FALSE
  )
  class(out) <- c("cm_associations", "data.frame")
  out
}

# Association measures, and the measures of earlier versions with what to do
# instead.
association_measures <- c("OR", "RR", "RD", "cond_prob", "phi")
retired_measures <- c(
  independent = "measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).",
  unknown = "measure = 'unknown' is no longer used: leave the pair out. Pairs without a row are unknown, and the global model fills in their association.",
  table = "measure = 'table' is no longer used: compute the odds ratio from the counts, n11 * n00 / (n10 * n01), and enter it with measure OR."
)

#' Describe three-way disease associations (scenarios)
#'
#' Pairwise associations do not determine how often three diseases occur
#' together. The global model ([fit_joint()]) fills this in by maximum
#' entropy, which assumes no three-way association: the odds ratio of a pair
#' is the same whether or not a third disease is present (given the other
#' diseases). A three-way term relaxes this for a chosen triple, while all
#' pairwise associations are still matched.
#'
#' `ratio` is the ratio of conditional odds ratios,
#' `OR(d1, d2 | d3 present) / OR(d1, d2 | d3 absent)` (given the other
#' diseases), which is symmetric in the three diseases. `ratio = 1` is the
#' maximum-entropy assumption. Pairwise evidence cannot identify `ratio`, so
#' three-way terms are sensitivity scenarios unless there is direct evidence.
#'
#' When every pair in the triple has an association, a three-way term leaves
#' the pairwise tables unchanged, so additive results without interactions
#' do not change; it affects only results that depend on the joint
#' distribution: additive impacts with interactions, and event impacts
#' (`event_model = TRUE` in [deconflate()]). When a pair is unknown (no
#' association given), the global fit determines its table, and a three-way
#' term can change it, and with it the additive results. A model with
#' three-way terms is always adjusted with the global method.
#'
#' @param disease1,disease2,disease3 Character vectors of disease ids.
#' @param ratio Positive ratio of conditional odds ratios (see Details).
#' @param source Optional citation or scenario label.
#' @return A `cm_three_way` data frame.
#' @export
#' @examples
#' cm_three_way("d1", "d2", "d3", ratio = 2)
cm_three_way <- function(disease1, disease2, disease3, ratio, source = NA_character_) {
  d <- cbind(as.character(disease1), as.character(disease2), as.character(disease3))
  n <- nrow(d)
  ratio <- as.numeric(recycle_arg(ratio, n, "ratio"))
  check_numeric(ratio, "ratio")
  if (any(ratio <= 0)) cm_abort("Three-way ratios must be positive.", class = "deconflate_infeasible")
  if (any(apply(d, 1, function(x) anyDuplicated(x) > 0))) {
    cm_abort("A three-way term needs three different diseases.")
  }
  key <- apply(d, 1, function(x) paste(sort(x), collapse = "|"))
  if (anyDuplicated(key)) cm_abort("Duplicate three-way terms.")
  out <- data.frame(disease1 = d[, 1], disease2 = d[, 2], disease3 = d[, 3],
                    ratio = ratio, source = as.character(recycle_arg(source, n, "source")),
                    stringsAsFactors = FALSE)
  class(out) <- c("cm_three_way", "data.frame")
  out
}

#' Describe the population: diseases and their associations
#'
#' A population holds the disease probabilities, the pairwise associations
#' and, optionally, three-way association scenarios. Combine it with one
#' impact table in [cm_model()].
#'
#' Pairs without an association are unknown: the global model fills in
#' their association from the others (see [fit_joint()]).
#'
#' @param diseases A [cm_diseases()] object.
#' @param associations Optional [cm_associations()] object. [deconflate()]
#'   needs at least one association; the sensitivity tools (e.g.
#'   [screen_associations()]) also run without any.
#' @param three_way Optional [cm_three_way()] object.
#' @param adjusted_associations Covariate-adjusted association measures are
#'   not marginal 2x2 associations. `"error"` (default) rejects them;
#'   `"use_as_marginal"` uses them as if they were marginal, which is an
#'   approximation, and records this.
#' @return A `cm_population` object.
#' @export
#' @examples
#' pop <- cm_population(
#'   cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
#'   cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
#' )
#' pop
cm_population <- function(diseases, associations = NULL, three_way = NULL,
                          adjusted_associations = c("error", "use_as_marginal")) {
  adjusted_associations <- match.arg(adjusted_associations)
  if (!inherits(diseases, "cm_diseases")) cm_abort("`diseases` must be created with cm_diseases().")
  ids <- diseases$id
  if (!is.null(associations)) {
    if (!inherits(associations, "cm_associations")) {
      cm_abort("`associations` must be created with cm_associations().")
    }
    unk <- setdiff(c(associations$disease1, associations$disease2), ids)
    if (length(unk)) cm_abort(sprintf("Associations refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
    adj <- associations$adjusted %in% TRUE
    if (any(adj) && adjusted_associations == "error") {
      cm_abort(sprintf(
        "Covariate-adjusted associations (%s) are not marginal 2x2 associations. Use crude measures, or set adjusted_associations = 'use_as_marginal' to use them as an approximation.",
        paste(associations$disease1[adj], associations$disease2[adj], sep = ":", collapse = ", ")),
        class = "deconflate_unsupported")
    }
  }
  if (!is.null(three_way)) {
    if (!inherits(three_way, "cm_three_way")) cm_abort("`three_way` must be created with cm_three_way().")
    unk <- setdiff(unlist(three_way[, c("disease1", "disease2", "disease3")]), ids)
    if (length(unk)) cm_abort(sprintf("Three-way terms refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
  }
  structure(list(diseases = diseases, associations = associations, three_way = three_way,
                 adjusted_associations = adjusted_associations),
            class = "cm_population")
}

#' Describe raw impact estimates
#'
#' One impact table holds one raw estimate per disease for one outcome. There
#' are two kinds:
#'
#' * **Additive impacts** (no `measure`): amounts in the outcome's own units,
#'   the same for every row (e.g. kg of milk, percent of yield, days, a
#'   welfare score). The package does not convert units; results come back
#'   in the units supplied. Adjust them with [deconflate()].
#' * **Event impacts** (`measure` given): comparisons of the risk of an event
#'   (e.g. death or culling) between animals with and without the disease,
#'   as hazard ratios, rate ratios, risk ratios, odds ratios or risk
#'   differences. Adjust them with `deconflate(..., event_model = TRUE)`,
#'   which uses the snapshot hazard model.
#'
#' For several outcomes, make one impact table each and adjust each in turn.
#'
#' @section Estimands of additive impacts:
#' * `"crude"` (default): the difference in the outcome between animals with
#'   and without the disease (unadjusted for other diseases).
#' * `"adjusted_linear"`: the coefficient of the disease in an additive
#'   (linear) regression of the outcome on the disease and the diseases in
#'   `adjusted_for`, in the same source population. The adjustment then uses
#'   the population projection of the omitted diseases (see [deconflate()]).
#'   `adjusted_for = "all"` means every other disease in the model.
#'
#' Other adjusted estimands (e.g. matched or propensity-score estimates) are
#' not supported. The probabilities and associations must describe the
#' population the estimates come from.
#'
#' @section Event impacts:
#' `measure` is one of `"HR"` (hazard ratio), `"rate_ratio"`, `"RR"` (risk
#' ratio over the period), `"OR"` (odds ratio of the event over the period)
#' or `"RD"` (risk difference over the period); measures can be mixed. The
#' estimand must be stated (there is no default, because it is an
#' assumption):
#' * `"snapshot_crude"`: the comparison of animals with and without the
#'   disease at the start of the period, in the population described by the
#'   probabilities and associations (unadjusted for other diseases);
#' * `"snapshot_stratified"`: the same comparison within strata of the
#'   diseases in `adjusted_for`, combined across strata with
#'   Mantel-Haenszel-type weights; with `adjusted_for = "all"` a hazard
#'   ratio is the disease's own hazard multiplier.
#'
#' A Cox hazard ratio estimated over follow-up is not exactly a snapshot
#' estimand (see [deconflate()], "Event impacts"); entering one as such is an
#' approximation and the user's assumption. Risk ratios, odds ratios and
#' risk differences refer to the period of the overall risk given to
#' [deconflate()].
#'
#' @param disease Character vector of disease ids (one row per disease).
#' @param value Numeric raw impacts. Every disease in the model needs a value
#'   (for no effect: 0 for additive impacts and risk differences, 1 for
#'   ratios).
#' @param estimand Additive impacts: `"crude"` (default) or
#'   `"adjusted_linear"`. Event impacts: `"snapshot_crude"` or
#'   `"snapshot_stratified"` (required). One per row or recycled.
#' @param adjusted_for For `"adjusted_linear"` and `"snapshot_stratified"`:
#'   the diseases the estimate was adjusted for, separated by `";"` (e.g.
#'   `"LAM; CM"`), or `"all"`.
#' @param source Optional citation.
#' @param label,units Optional table-level metadata (e.g.
#'   `label = "milk yield loss"`, `units = "% of yield"`), carried into the
#'   results.
#' @param measure `NULL` (default) for additive impacts; for event impacts,
#'   one measure per row or recycled (see Event impacts).
#'
#' @return A `cm_impacts` data frame (attributes `label`, `units` and
#'   `kind`, `"additive"` or `"event"`).
#' @export
#' @examples
#' cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%")
#' # Event impacts: culling, with mixed measures
#' cm_impacts(c("d1", "d2", "d3"), c(1.5, 1.3, 0.04), measure = c("HR", "RR", "RD"),
#'            estimand = "snapshot_crude", label = "culling")
cm_impacts <- function(disease, value, estimand = NULL, adjusted_for = NA_character_,
                       source = NA_character_, label = NULL, units = NULL, measure = NULL) {
  disease <- as.character(disease)
  n <- length(disease)
  if (!n) cm_abort("At least one impact is required.")
  if (length(value) != n) cm_abort("`value` must have one entry per disease.")
  check_numeric(value, "value")
  if (anyDuplicated(disease)) {
    cm_abort(sprintf("Each disease may have only one impact (duplicated: %s).",
                     paste(unique(disease[duplicated(disease)]), collapse = ", ")))
  }
  event <- !is.null(measure)
  if (event) {
    measure <- as.character(recycle_arg(measure, n, "measure"))
    check_choices(measure, event_measures, "measure")
    if (is.null(estimand)) {
      cm_abort("State the estimand of event impacts: estimand = \"snapshot_crude\" or \"snapshot_stratified\" (see ?cm_impacts).")
    }
    estimand <- as.character(recycle_arg(estimand, n, "estimand"))
    check_choices(estimand, c("snapshot_crude", "snapshot_stratified"), "estimand")
    bad <- measure != "RD" & value <= 0
    if (any(bad)) {
      cm_abort(sprintf("Ratios must be positive (check: %s).", paste(disease[bad], collapse = ", ")))
    }
    bad <- measure == "RD" & (value <= -1 | value >= 1)
    if (any(bad)) {
      cm_abort(sprintf("Risk differences must lie between -1 and 1 (check: %s).", paste(disease[bad], collapse = ", ")))
    }
    adj_estimand <- "snapshot_stratified"
  } else {
    estimand <- as.character(recycle_arg(estimand %||% "crude", n, "estimand"))
    if (any(estimand %in% c("snapshot_crude", "snapshot_stratified"))) {
      cm_abort("Snapshot estimands belong to event impacts: give the `measure` of each row (e.g. \"HR\") and use event_model = TRUE in deconflate().",
               class = "deconflate_unsupported")
    }
    check_choices(estimand, c("crude", "adjusted_linear"), "estimand")
    measure <- rep(NA_character_, n)
    adj_estimand <- "adjusted_linear"
  }
  adjusted_for <- as.character(recycle_arg(adjusted_for, n, "adjusted_for"))
  has_adj <- vapply(adjusted_for, function(x) length(split_ids(x)) > 0, logical(1))
  bad <- estimand != adj_estimand & has_adj
  if (any(bad)) {
    cm_abort(sprintf(
      "`adjusted_for` is given for %s estimates (%s). Set estimand = '%s' if these estimates are adjusted for those diseases; other adjusted estimands are not supported.",
      estimand[bad][1], paste(disease[bad], collapse = ", "), adj_estimand), class = "deconflate_unsupported")
  }
  bad <- estimand == adj_estimand & !has_adj
  if (any(bad)) {
    cm_abort(sprintf("estimand = '%s' needs `adjusted_for` (check: %s).", adj_estimand,
                     paste(disease[bad], collapse = ", ")))
  }
  out <- data.frame(disease = disease, value = value, estimand = estimand,
                    adjusted_for = adjusted_for, measure = measure,
                    source = as.character(recycle_arg(source, n, "source")),
                    stringsAsFactors = FALSE)
  attr(out, "label") <- label
  attr(out, "units") <- units
  attr(out, "kind") <- if (event) "event" else "additive"
  class(out) <- c("cm_impacts", "data.frame")
  out
}

# Measures of event impacts.
event_measures <- c("HR", "rate_ratio", "RR", "OR", "RD")

# Kind of an impact table: "event" when it has measures, else "additive".
impact_kind <- function(impacts) {
  if (is.null(impacts)) return("additive")
  k <- attr(impacts, "kind")
  if (!is.null(k)) return(k)
  if (!is.null(impacts$measure) && any(!is.na(impacts$measure))) "event" else "additive"
}

#' Describe pairwise impact interactions
#'
#' An interaction `delta` is the additional impact when both diseases are
#' present, in the same units as the impact vector: positive values are
#' synergistic (more than the sum), negative values antagonistic.
#' Interactions cannot be inferred from associations and must come from
#' evidence or explicit scenarios. A model with interactions is always
#' adjusted with the global method (see [deconflate()]). Interactions apply
#' to additive impacts only.
#'
#' @param disease1,disease2 Character vectors of disease ids.
#' @param value Numeric interaction values (same units as the impacts).
#' @param source Optional citation or scenario label.
#' @return A `cm_interactions` data frame.
#' @export
#' @examples
#' cm_interactions("d1", "d2", 0.5)
cm_interactions <- function(disease1, disease2, value, source = NA_character_) {
  disease1 <- as.character(disease1)
  disease2 <- as.character(disease2)
  n <- length(disease1)
  if (length(disease2) != n || length(value) != n) {
    cm_abort("`disease1`, `disease2` and `value` must have equal length.")
  }
  check_numeric(value, "value")
  if (any(disease1 == disease2)) cm_abort("An interaction must involve two different diseases.")
  key <- pair_key(disease1, disease2)
  if (anyDuplicated(key)) cm_abort("Duplicate interactions for the same pair.")
  out <- data.frame(disease1 = disease1, disease2 = disease2, value = value,
                    source = as.character(recycle_arg(source, n, "source")),
                    stringsAsFactors = FALSE)
  class(out) <- c("cm_interactions", "data.frame")
  out
}

#' Combine a population with an impact table
#'
#' @param population A [cm_population()]. For convenience, a [cm_diseases()]
#'   object can be given instead, together with `associations` and
#'   `three_way`.
#' @param impacts A [cm_impacts()] object with one value per disease.
#' @param interactions Optional [cm_interactions()] object (additive impacts
#'   only).
#' @param associations,three_way,adjusted_associations Used only when
#'   `population` is a [cm_diseases()] object; see [cm_population()].
#' @param distributions Optional named list of distributions ([distributions])
#'   for uncertain inputs, used when [deconflate()] runs draws. Names are
#'   keys: `"prob:<disease>"` (the disease's `value`, on the scale it was
#'   entered), `"assoc:<d1>:<d2>"` (the association, on its own measure),
#'   `"three:<d1>:<d2>:<d3>"`, `"impact:<disease>"` and `"inter:<d1>:<d2>"`.
#'   [cm_read_inputs()] fills this in from the `dist` columns of the tables;
#'   [cm_dist_table()] builds it from a table of keys. When `population` is a
#'   model and `distributions` is not given, its distributions are kept,
#'   except those of impacts or interactions that are replaced.
#' @return A `cm_model` object (which is also a `cm_population`).
#' @export
#' @examples
#' m <- cm_model(
#'   cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
#'   cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
#'   associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3)),
#'   distributions = list("impact:d1" = dist_normal(2.5, 0.5))
#' )
#' m
cm_model <- function(population, impacts = NULL, interactions = NULL,
                     associations = NULL, three_way = NULL,
                     adjusted_associations = c("error", "use_as_marginal"),
                     distributions = NULL) {
  if (inherits(population, "cm_diseases")) {
    population <- cm_population(population, associations, three_way,
                                adjusted_associations = match.arg(adjusted_associations))
  }
  if (!inherits(population, "cm_population")) {
    cm_abort("`population` must be created with cm_population() (or be a cm_diseases object).")
  }
  ids <- population$diseases$id
  if (!is.null(impacts)) {
    if (!inherits(impacts, "cm_impacts")) cm_abort("`impacts` must be created with cm_impacts().")
    unk <- setdiff(impacts$disease, ids)
    if (length(unk)) cm_abort(sprintf("Impacts refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
    kind <- impact_kind(impacts)
    miss <- setdiff(ids, impacts$disease)
    if (length(miss)) {
      cm_abort(sprintf("No impact for: %s. Every disease needs a value (for no effect: %s).",
                       paste(miss, collapse = ", "),
                       if (kind == "event") "1 for ratios, 0 for risk differences" else "0"))
    }
    for (r in seq_len(nrow(impacts))) {
      s <- split_ids(impacts$adjusted_for[r])
      if (length(s) == 1L && tolower(s) == "all") next
      unk <- setdiff(s, ids)
      if (length(unk)) {
        cm_abort(sprintf("`adjusted_for` of %s refers to unknown diseases: %s.",
                         impacts$disease[r], paste(unk, collapse = ", ")))
      }
    }
    lab <- attr(impacts, "label")
    un <- attr(impacts, "units")
    impacts <- impacts[match(ids, impacts$disease), , drop = FALSE]
    rownames(impacts) <- NULL
    attr(impacts, "label") <- lab
    attr(impacts, "units") <- un
    attr(impacts, "kind") <- kind
  }
  if (!is.null(interactions)) {
    if (!inherits(interactions, "cm_interactions")) {
      cm_abort("`interactions` must be created with cm_interactions().")
    }
    if (is.null(impacts)) cm_abort("Interactions need an impact table.")
    if (impact_kind(impacts) == "event") {
      cm_abort("Interactions apply to additive impacts only, not to event impacts.",
               class = "deconflate_unsupported")
    }
    unk <- setdiff(c(interactions$disease1, interactions$disease2), ids)
    if (length(unk)) cm_abort(sprintf("Interactions refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
  }
  out <- population
  out$impacts <- impacts
  out$interactions <- interactions
  class(out) <- c("cm_model", "cm_population")
  if (is.null(distributions)) {
    # Distributions carried over from a model belong to its tables: those of
    # replaced impacts or interactions are dropped.
    distributions <- population$distributions
    if (length(distributions) && inherits(population, "cm_model")) {
      type <- vapply(strsplit(names(distributions), ":", fixed = TRUE), `[`, character(1), 1)
      drop <- (type == "impact" & !identical(impacts, population$impacts)) |
        (type == "inter" & !identical(interactions, population$interactions))
      distributions <- distributions[!drop]
    }
  }
  out$distributions <- check_distributions(distributions, out)
  out
}

# Check the distributions of a model: a named list of cm_dist objects whose
# keys match the model's inputs. Returns the list in canonical key form
# (associations and interactions as listed in their tables, three-way terms
# in their table's order).
check_distributions <- function(d, model) {
  if (is.null(d) || !length(d)) return(NULL)
  if (!is.list(d) || is.null(names(d)) || any(!nzchar(names(d))) ||
      !all(vapply(d, inherits, logical(1), "cm_dist"))) {
    cm_abort("`distributions` must be a named list of distributions (see ?distributions).")
  }
  if (anyDuplicated(names(d))) {
    cm_abort(sprintf("`distributions` names an input more than once: %s.",
                     paste(unique(names(d)[duplicated(names(d))]), collapse = ", ")))
  }
  out <- list()
  for (k in names(d)) {
    ck <- canonical_key(k, model)
    if (!is.null(out[[ck]])) cm_abort(sprintf("The input '%s' is given more than once in `distributions`.", k))
    out[[ck]] <- d[[k]]
  }
  out
}

# Canonical form of an input key, checked against the model.
canonical_key <- function(key, model) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  type <- parts[1]
  ids <- model$diseases$id
  bad <- function(why) cm_abort(sprintf("Distribution '%s': %s.", key, why))
  pair_row <- function(tab, what) {
    if (length(parts) != 3L) bad(sprintf("use '%s:<d1>:<d2>'", type))
    if (is.null(tab) || !nrow(tab)) bad(sprintf("the model has no %s", what))
    r <- match(pair_key(parts[2], parts[3]), pair_key(tab$disease1, tab$disease2))
    if (is.na(r)) bad(sprintf("no %s row for this pair", what))
    r
  }
  switch(type,
    prob = {
      if (length(parts) != 2L || !(parts[2] %in% ids)) bad("unknown disease")
      key
    },
    impact = {
      if (length(parts) != 2L || !(parts[2] %in% ids)) bad("unknown disease")
      if (is.null(model$impacts)) bad("the model has no impacts")
      key
    },
    assoc = {
      r <- pair_row(model$associations, "association")
      paste0("assoc:", model$associations$disease1[r], ":", model$associations$disease2[r])
    },
    inter = {
      r <- pair_row(model$interactions, "interaction")
      paste0("inter:", model$interactions$disease1[r], ":", model$interactions$disease2[r])
    },
    three = {
      tw <- model$three_way
      if (length(parts) != 4L) bad("use 'three:<d1>:<d2>:<d3>'")
      if (is.null(tw) || !nrow(tw)) bad("the model has no three-way terms")
      tkey <- apply(tw[, c("disease1", "disease2", "disease3")], 1, function(x) paste(sort(x), collapse = "|"))
      r <- match(paste(sort(parts[2:4]), collapse = "|"), tkey)
      if (is.na(r)) bad("no three-way row for this triple")
      paste0("three:", tw$disease1[r], ":", tw$disease2[r], ":", tw$disease3[r])
    },
    bad("unknown input type (use prob:, assoc:, three:, impact: or inter:)")
  )
}

check_population <- function(x) {
  if (!inherits(x, "cm_population")) {
    cm_abort("Expected a cm_population() or cm_model() object.")
  }
  invisible(x)
}

check_model <- function(model) {
  if (!inherits(model, "cm_model")) cm_abort("`model` must be created with cm_model().")
  invisible(model)
}
