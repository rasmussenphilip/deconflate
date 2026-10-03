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
#' Each row gives one association measure for one pair. Pairs not listed are
#' handled by `missing_associations` in [cm_population()].
#'
#' Odds ratios (and the other measures) are applied to the modelled
#' population's own marginal probabilities ("transported"). This assumes the
#' measure is the same in the source study and in the modelled population.
#'
#' @param disease1,disease2 Character vectors of disease ids.
#' @param value Numeric association values (ignored for `"independent"`,
#'   `"unknown"` and `"table"`).
#' @param measure Character, one per row (or recycled). Directional measures
#'   treat `disease1` as the outcome and `disease2` as the conditioning
#'   disease:
#'   * `"OR"`: odds ratio (symmetric).
#'   * `"RR"`: risk ratio `P(d1 | d2) / P(d1 | not d2)`.
#'   * `"RD"`: risk difference `P(d1 | d2) - P(d1 | not d2)`.
#'   * `"cond_prob"`: conditional probability `P(d1 | d2)`.
#'   * `"phi"`: binary (phi) correlation coefficient (symmetric).
#'   * `"table"`: a study contingency table given in `n11`, `n10`, `n01` and
#'     `n00`. The table's odds ratio is used (transported to the modelled
#'     marginals), because the study's marginal frequencies generally differ
#'     from the modelled population's.
#'   * `"independent"`: independence is imposed (odds ratio 1).
#'   * `"unknown"`: the association is unknown. The pairwise methods reject
#'     unknown pairs; the global method leaves them unconstrained, so their
#'     association is implied by the maximum-entropy fit.
#' @param n11,n10,n01,n00 Counts for `measure = "table"`: both diseases, `d1`
#'   only, `d2` only and neither.
#' @param zero_cell For tables with a zero cell: `"haldane"` (default) adds
#'   0.5 to every cell (Haldane-Anscombe correction) and records it in column
#'   `corrected`; `"error"` rejects such tables.
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
cm_associations <- function(disease1, disease2, value = NA_real_,
                            measure = "OR",
                            n11 = NA_real_, n10 = NA_real_,
                            n01 = NA_real_, n00 = NA_real_,
                            zero_cell = c("haldane", "error"),
                            adjusted = FALSE, adjusted_for = NA_character_,
                            source = NA_character_) {
  zero_cell <- match.arg(zero_cell)
  disease1 <- as.character(disease1)
  disease2 <- as.character(disease2)
  n <- length(disease1)
  if (length(disease2) != n) cm_abort("`disease1` and `disease2` must have equal length.")
  if (any(disease1 == disease2)) cm_abort("An association must involve two different diseases.")
  measure <- as.character(recycle_arg(measure, n, "measure"))
  check_choices(measure, c("OR", "RR", "RD", "cond_prob", "phi", "table",
                           "independent", "unknown"), "measure")
  value <- as.numeric(recycle_arg(value, n, "value"))
  tab <- cbind(n11 = recycle_arg(n11, n, "n11"), n10 = recycle_arg(n10, n, "n10"),
               n01 = recycle_arg(n01, n, "n01"), n00 = recycle_arg(n00, n, "n00"))
  corrected <- rep(FALSE, n)

  is_tab <- measure == "table"
  if (any(is_tab)) {
    tt <- tab[is_tab, , drop = FALSE]
    if (anyNA(tt) || any(tt < 0)) cm_abort("Contingency tables need non-negative n11, n10, n01 and n00.")
    if (any(rowSums(tt) == 0)) cm_abort("A contingency table is empty (all counts are zero).")
    zero <- apply(tt == 0, 1, any)
    if (any(zero) && zero_cell == "error") {
      cm_abort("A contingency table has a zero cell; use zero_cell = 'haldane' to add 0.5 to every cell.")
    }
    tt[zero, ] <- tt[zero, ] + 0.5
    value[is_tab] <- (tt[, "n11"] * tt[, "n00"]) / (tt[, "n10"] * tt[, "n01"])
    corrected[is_tab] <- zero
  }
  value[measure == "independent"] <- 1
  needs_value <- !(measure %in% c("unknown", "independent", "table"))
  if (anyNA(value[needs_value])) cm_abort("Association `value` is missing for some rows.")
  if (any(!is.finite(value[needs_value]))) cm_abort("Association values must be finite.")
  if (any(value[measure %in% c("OR", "RR", "table")] <= 0)) {
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
    corrected = corrected,
    stringsAsFactors = FALSE
  )
  out <- cbind(out, as.data.frame(tab))
  class(out) <- c("cm_associations", "data.frame")
  out
}

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
#' When every pair in the triple is constrained (an association is given, or
#' the pair is independent by default), a three-way term leaves the pairwise
#' tables unchanged, so additive results without interactions do not change;
#' it affects only results that depend on the joint distribution: the global
#' method with interactions, the hazard-ratio snapshot model and
#' [attributable_risk()]. When a pair is unknown (`missing_associations =
#' "unknown"`), the global fit determines its table, and a three-way term can
#' change it, and with it the additive results.
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
#' A population holds everything that is shared by all impact analyses: the
#' disease probabilities, the pairwise associations and, optionally,
#' three-way association scenarios. Combine it with an impact vector in
#' [cm_model()], or with several in [cm_analyses()].
#'
#' @param diseases A [cm_diseases()] object.
#' @param associations Optional [cm_associations()] object.
#' @param three_way Optional [cm_three_way()] object (global model only).
#' @param missing_associations How to treat pairs without an association:
#'   `"independent"` (default; odds ratio 1, as in Rasmussen et al. 2022) or
#'   `"unknown"` (unconstrained). These are different assumptions.
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
                          missing_associations = c("independent", "unknown"),
                          adjusted_associations = c("error", "use_as_marginal")) {
  missing_associations <- match.arg(missing_associations)
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
                 missing_associations = missing_associations,
                 adjusted_associations = adjusted_associations),
            class = "cm_population")
}

#' Describe raw impact estimates (one impact vector)
#'
#' One analysis adjusts one set of compatible, additive impact estimates:
#' one value per disease, all in the same units (e.g. kg of milk per cow,
#' percent of yield, days, euros, welfare scores). The engine does not
#' convert units; results come back in the units supplied. For several
#' types of impact, use one impact vector each ([cm_analyses()]).
#'
#' @section Estimands:
#' Each value must be one of the supported estimands:
#' * `"crude"`: the difference in the outcome between animals with and
#'   without the disease (unadjusted for other diseases).
#' * `"adjusted_linear"`: the coefficient of the disease in an additive
#'   (linear) regression of the outcome on the disease and the diseases in
#'   `adjusted_for`, in the same source population. The adjustment then uses
#'   the population projection of the omitted diseases (see [deconflate()]).
#'   `adjusted_for = "all"` means every other disease in the model.
#'
#' `adjusted_for` is only used with `estimand = "adjusted_linear"`; it is never
#' used to infer the estimand. Other adjusted estimands (e.g. matched or
#' propensity-score estimates) are not supported. The probabilities and
#' associations must describe the population the estimates come from.
#'
#' @param disease Character vector of disease ids (one row per disease).
#' @param value Numeric raw impacts. Every disease in the model needs a value
#'   (use 0 for no impact).
#' @param estimand `"crude"` (default) or `"adjusted_linear"`, one per row or
#'   recycled.
#' @param adjusted_for For `"adjusted_linear"`: the diseases the estimate was
#'   adjusted for, separated by `";"` (e.g. `"LAM; CM"`), or `"all"`.
#' @param source Optional citation.
#' @param label,units Optional analysis-level metadata (e.g.
#'   `label = "milk yield loss"`, `units = "% of yield"`), carried into the
#'   results.
#'
#' @return A `cm_impacts` data frame (attributes `label` and `units`).
#' @export
#' @examples
#' cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%")
cm_impacts <- function(disease, value, estimand = "crude", adjusted_for = NA_character_,
                       source = NA_character_, label = NULL, units = NULL) {
  disease <- as.character(disease)
  n <- length(disease)
  if (!n) cm_abort("At least one impact is required.")
  if (length(value) != n) cm_abort("`value` must have one entry per disease.")
  check_numeric(value, "value")
  if (anyDuplicated(disease)) {
    cm_abort(sprintf("Each disease may have only one impact (duplicated: %s).",
                     paste(unique(disease[duplicated(disease)]), collapse = ", ")))
  }
  estimand <- as.character(recycle_arg(estimand, n, "estimand"))
  check_choices(estimand, c("crude", "adjusted_linear"), "estimand")
  adjusted_for <- as.character(recycle_arg(adjusted_for, n, "adjusted_for"))
  has_adj <- vapply(adjusted_for, function(x) length(split_ids(x)) > 0, logical(1))
  bad <- estimand == "crude" & has_adj
  if (any(bad)) {
    cm_abort(sprintf(
      "`adjusted_for` is given for crude estimates (%s). Set estimand = 'adjusted_linear' if these are coefficients from an additive regression adjusted for those diseases; other adjusted estimands are not supported.",
      paste(disease[bad], collapse = ", ")), class = "deconflate_unsupported")
  }
  bad <- estimand == "adjusted_linear" & !has_adj
  if (any(bad)) {
    cm_abort(sprintf("estimand = 'adjusted_linear' needs `adjusted_for` (check: %s).",
                     paste(disease[bad], collapse = ", ")))
  }
  out <- data.frame(disease = disease, value = value, estimand = estimand,
                    adjusted_for = adjusted_for,
                    source = as.character(recycle_arg(source, n, "source")),
                    stringsAsFactors = FALSE)
  attr(out, "label") <- label
  attr(out, "units") <- units
  class(out) <- c("cm_impacts", "data.frame")
  out
}

#' Describe pairwise impact interactions
#'
#' An interaction `delta` is the additional impact when both diseases are
#' present, in the same units as the impact vector: positive values are
#' synergistic (more than the sum), negative values antagonistic.
#' Interactions cannot be inferred from associations and must come from
#' evidence or explicit scenarios. They require `method = "global"` in
#' [deconflate()].
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

#' Combine a population with an impact vector
#'
#' @param population A [cm_population()]. For convenience, a [cm_diseases()]
#'   object can be given instead, together with `associations`,
#'   `missing_associations` and `three_way`.
#' @param impacts A [cm_impacts()] object with one value per disease.
#' @param interactions Optional [cm_interactions()] object.
#' @param associations,three_way,missing_associations,adjusted_associations
#'   Used only when `population` is a [cm_diseases()] object; see
#'   [cm_population()].
#' @return A `cm_model` object (which is also a `cm_population`).
#' @export
#' @examples
#' m <- cm_model(
#'   cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
#'   cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
#'   associations = cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
#' )
#' m
cm_model <- function(population, impacts = NULL, interactions = NULL,
                     associations = NULL, three_way = NULL,
                     missing_associations = c("independent", "unknown"),
                     adjusted_associations = c("error", "use_as_marginal")) {
  if (inherits(population, "cm_diseases")) {
    population <- cm_population(population, associations, three_way,
                                missing_associations = match.arg(missing_associations),
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
    miss <- setdiff(ids, impacts$disease)
    if (length(miss)) {
      cm_abort(sprintf("No impact for: %s. Every disease needs a value (use 0 for no impact).",
                       paste(miss, collapse = ", ")))
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
  }
  if (!is.null(interactions)) {
    if (!inherits(interactions, "cm_interactions")) {
      cm_abort("`interactions` must be created with cm_interactions().")
    }
    if (is.null(impacts)) cm_abort("Interactions need an impact vector.")
    unk <- setdiff(c(interactions$disease1, interactions$disease2), ids)
    if (length(unk)) cm_abort(sprintf("Interactions refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
  }
  out <- population
  out$impacts <- impacts
  out$interactions <- interactions
  class(out) <- c("cm_model", "cm_population")
  out
}

#' Several impact analyses on one population
#'
#' Bundles several impact vectors (e.g. milk yield, calving interval and
#' welfare) that share the same diseases and associations. Each is adjusted
#' separately; Monte Carlo runs on a batch sampler ([cm_batch_sampler()])
#' reuse the same disease and association draws for every analysis.
#'
#' @param population A [cm_population()] (or a [cm_model()], whose impacts
#'   are dropped).
#' @param ... Named [cm_impacts()] objects, one per analysis.
#' @param interactions Optional named list of [cm_interactions()] objects,
#'   with names matching the analyses.
#' @return A `cm_analyses` object with `population` and `models` (a named
#'   list of [cm_model()] objects).
#' @export
#' @examples
#' pop <- example_supplement()
#' a <- cm_analyses(pop,
#'   yield = cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), units = "%"),
#'   fertility = cm_impacts(c("d1", "d2", "d3"), c(1, 2, 0), units = "%"))
#' deconflate(a)
cm_analyses <- function(population, ..., interactions = list()) {
  if (!inherits(population, "cm_population")) cm_abort("`population` must come from cm_population().")
  imps <- list(...)
  if (!length(imps) || is.null(names(imps)) || any(!nzchar(names(imps)))) {
    cm_abort("Give each impact vector a name, e.g. cm_analyses(pop, yield = imp1, fertility = imp2).")
  }
  if (anyDuplicated(names(imps))) cm_abort("Analysis names must be unique.")
  if (length(interactions) && (is.null(names(interactions)) ||
                               length(setdiff(names(interactions), names(imps))))) {
    cm_abort("`interactions` must be a list named after the analyses.")
  }
  pop <- population
  pop$impacts <- NULL
  pop$interactions <- NULL
  class(pop) <- "cm_population"
  models <- lapply(names(imps), function(nm) {
    imp <- imps[[nm]]
    if (is.null(attr(imp, "label"))) attr(imp, "label") <- nm
    cm_model(pop, imp, interactions[[nm]])
  })
  names(models) <- names(imps)
  structure(list(population = pop, models = models), class = "cm_analyses")
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
