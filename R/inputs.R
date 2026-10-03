#' Describe the diseases in a system
#'
#' @param id Character vector of unique disease identifiers.
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
  if (anyNA(id) || any(!nzchar(id))) cm_abort("Disease ids must be non-empty.")
  if (anyDuplicated(id)) cm_abort("Disease ids must be unique.")
  if (any(grepl("[|;]", id))) cm_abort("Disease ids must not contain '|' or ';'.")
  if (length(value) != n) cm_abort("`value` must have one entry per disease.")
  check_numeric(value, "value")
  if (any(value < 0)) cm_abort("`value` must be non-negative.")
  type <- as.character(recycle_arg(type, n, "type"))
  check_choices(type, c("prevalence", "probability", "incidence_rate"), "type")

  prob <- ifelse(type == "incidence_rate", 1 - exp(-value), value)
  if (any(prob <= 0 | prob >= 1)) {
    cm_abort(sprintf("Disease probabilities must lie strictly between 0 and 1 (check: %s).",
                     paste(id[prob <= 0 | prob >= 1], collapse = ", ")))
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
#' handled by `missing_associations` in [cm_model()].
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
#'     `n00`. The table's odds ratio is used, because the study's marginal
#'     frequencies generally differ from the modelled population's.
#'   * `"independent"`: independence is imposed (odds ratio 1).
#'   * `"unknown"`: the association is unknown. The pairwise methods reject
#'     unknown pairs; the global method leaves them unconstrained, so their
#'     association is implied by the maximum-entropy fit.
#' @param n11,n10,n01,n00 Counts for `measure = "table"`: both diseases, `d1`
#'   only, `d2` only and neither.
#' @param adjusted Logical: was the measure adjusted for covariates (e.g. an
#'   odds ratio from multivariable logistic regression)? The 2x2 algebra
#'   assumes marginal measures; adjusted measures are used as if marginal and
#'   flagged in the model summary.
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
                            adjusted = FALSE, adjusted_for = NA_character_,
                            source = NA_character_) {
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

  is_tab <- measure == "table"
  if (any(is_tab)) {
    tt <- tab[is_tab, , drop = FALSE]
    if (anyNA(tt) || any(tt < 0)) cm_abort("Contingency tables need non-negative n11, n10, n01 and n00.")
    # Haldane-Anscombe correction when any cell is zero.
    tt[apply(tt == 0, 1, any), ] <- tt[apply(tt == 0, 1, any), ] + 0.5
    value[is_tab] <- (tt[, "n11"] * tt[, "n00"]) / (tt[, "n10"] * tt[, "n01"])
  }
  value[measure == "independent"] <- 1
  needs_value <- !(measure %in% c("unknown", "independent", "table"))
  if (anyNA(value[needs_value])) cm_abort("Association `value` is missing for some rows.")
  if (any(value[measure %in% c("OR", "RR", "table")] <= 0)) {
    cm_abort("Odds ratios and risk ratios must be positive.")
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
  out <- cbind(out, as.data.frame(tab))
  class(out) <- c("cm_associations", "data.frame")
  out
}

#' Describe raw (unadjusted) disease impact estimates
#'
#' @param disease Character vector of disease ids.
#' @param value Numeric raw impact estimates.
#' @param outcome Character outcome label(s), e.g. `"yield"`, `"fertility"`.
#'   Each disease needs exactly one impact per outcome (use `0` for no
#'   impact).
#' @param scale The scale of `value`, constant within an outcome:
#'   * `"proportion"`: proportional change relative to the disease-free value
#'     (e.g. `0.025`);
#'   * `"percent"`: converted to a proportion;
#'   * `"absolute"`: in `units` (e.g. excess culling risk);
#'   * `"hazard_ratio"`: a hazard ratio (e.g. for culling or mortality),
#'     adjusted on the log scale; see [deconflate()] and
#'     [attributable_risk()].
#' @param units Optional units label (required for `"absolute"`).
#' @param direction `"decrease"` if disease lowers the outcome (e.g. yield) or
#'   `"increase"` if it raises it (e.g. calving interval). Used by
#'   [productivity_gap()]. Set to `"increase"` for hazard ratios.
#' @param adjusted_for Diseases the raw estimate was already adjusted for,
#'   separated by `";"` (e.g. `"LAM; CM"`). Their conflation terms are
#'   removed from the adjustment for this impact.
#' @param source Optional citation.
#'
#' @return A `cm_impacts` data frame.
#' @export
#' @examples
#' cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), outcome = "yield",
#'            scale = "percent")
#' cm_impacts(c("d1", "d2", "d3"), c(1.5, 2.3, 1.1), outcome = "culling",
#'            scale = "hazard_ratio")
cm_impacts <- function(disease, value, outcome = "impact",
                       scale = "proportion", units = NA_character_,
                       direction = "decrease", adjusted_for = NA_character_,
                       source = NA_character_) {
  disease <- as.character(disease)
  n <- length(disease)
  if (length(value) != n) cm_abort("`value` must have one entry per row.")
  check_numeric(value, "value")
  outcome <- as.character(recycle_arg(outcome, n, "outcome"))
  scale <- as.character(recycle_arg(scale, n, "scale"))
  check_choices(scale, c("proportion", "percent", "absolute", "hazard_ratio"), "scale")
  direction <- as.character(recycle_arg(direction, n, "direction"))
  check_choices(direction, c("decrease", "increase"), "direction")
  units <- as.character(recycle_arg(units, n, "units"))
  if (any(scale == "absolute" & is.na(units))) {
    cm_abort("`units` is required for impacts on the absolute scale.")
  }
  is_hr <- scale == "hazard_ratio"
  if (any(is_hr & !(value > 0))) cm_abort("Hazard ratios must be positive.")
  direction[is_hr] <- "increase"
  units[is_hr & is.na(units)] <- "hazard ratio"
  input_scale <- scale
  value <- ifelse(scale == "percent", value / 100, value)
  scale[scale == "percent"] <- "proportion"

  out <- data.frame(
    disease = disease, outcome = outcome, value = value, scale = scale,
    units = units, direction = direction,
    adjusted_for = as.character(recycle_arg(adjusted_for, n, "adjusted_for")),
    source = as.character(recycle_arg(source, n, "source")),
    input_scale = input_scale,
    stringsAsFactors = FALSE
  )
  validate_impacts(out)
}

validate_impacts <- function(out) {
  if (anyDuplicated(paste(out$outcome, out$disease, sep = "|"))) {
    cm_abort("Each disease may have only one impact per outcome.")
  }
  for (o in unique(out$outcome)) {
    rows <- out[out$outcome == o, , drop = FALSE]
    if (length(unique(rows$scale)) > 1L || length(unique(rows$direction)) > 1L ||
        length(unique(rows$units)) > 1L) {
      cm_abort(sprintf("Outcome '%s' mixes scales, units or directions.", o))
    }
  }
  rownames(out) <- NULL
  class(out) <- c("cm_impacts", "data.frame")
  out
}

#' Combine impact tables
#'
#' Stacks several [cm_impacts()] objects (e.g. yield and fertility impacts
#' built separately, or culling impacts from [as_impacts()]) into one, and
#' re-validates the result.
#'
#' @param ... `cm_impacts` objects.
#' @return A `cm_impacts` data frame.
#' @export
combine_impacts <- function(...) {
  parts <- list(...)
  if (!length(parts) || !all(vapply(parts, inherits, logical(1), "cm_impacts"))) {
    cm_abort("All arguments must be cm_impacts objects.")
  }
  cols <- c("disease", "outcome", "value", "scale", "units", "direction",
            "adjusted_for", "source", "input_scale")
  out <- do.call(rbind, lapply(parts, function(p) {
    p <- as.data.frame(unclass(p), stringsAsFactors = FALSE)
    if (is.null(p$input_scale)) p$input_scale <- p$scale
    p[, cols]
  }))
  validate_impacts(out)
}

#' Describe pairwise impact interactions
#'
#' An interaction `delta` is the additional impact when both diseases are
#' present, on the same scale as the outcome's impacts: positive values are
#' synergistic (more loss than the sum), negative values antagonistic.
#' Interactions cannot be inferred from associations and must come from
#' evidence or explicit scenarios. They require `method = "global"` in
#' [deconflate()].
#'
#' @param disease1,disease2 Character vectors of disease ids.
#' @param value Numeric interaction values.
#' @param outcome Outcome label(s) matching [cm_impacts()].
#' @param scale `"proportion"` or `"percent"`; additive on that scale.
#' @param source Optional citation or scenario label.
#'
#' @return A `cm_interactions` data frame.
#' @export
cm_interactions <- function(disease1, disease2, value, outcome = "impact",
                            scale = "proportion", source = NA_character_) {
  disease1 <- as.character(disease1)
  disease2 <- as.character(disease2)
  n <- length(disease1)
  if (length(disease2) != n || length(value) != n) {
    cm_abort("`disease1`, `disease2` and `value` must have equal length.")
  }
  if (any(disease1 == disease2)) cm_abort("An interaction must involve two different diseases.")
  check_numeric(value, "value")
  scale <- as.character(recycle_arg(scale, n, "scale"))
  check_choices(scale, c("proportion", "percent"), "scale")
  input_scale <- scale
  value <- ifelse(scale == "percent", value / 100, value)
  outcome <- as.character(recycle_arg(outcome, n, "outcome"))
  key <- paste(outcome, pair_key(disease1, disease2))
  if (anyDuplicated(key)) cm_abort("Duplicate interactions for the same pair and outcome.")
  out <- data.frame(disease1 = disease1, disease2 = disease2, outcome = outcome,
                    value = value, scale = "proportion",
                    source = as.character(recycle_arg(source, n, "source")),
                    input_scale = input_scale,
                    stringsAsFactors = FALSE)
  class(out) <- c("cm_interactions", "data.frame")
  out
}

#' Combine inputs into a comorbidity model
#'
#' @param diseases A [cm_diseases()] object.
#' @param associations A [cm_associations()] object, or `NULL`.
#' @param impacts A [cm_impacts()] object, or `NULL` (e.g. when only the
#'   joint distribution or simulated impacts are needed).
#' @param interactions A [cm_interactions()] object, or `NULL`.
#' @param missing_associations How to treat disease pairs without a row in
#'   `associations`: `"independent"` imposes an odds ratio of 1 (as in
#'   Rasmussen et al. 2022, 2024); `"unknown"` leaves them unconstrained
#'   (only usable with the global method). These are different assumptions.
#'
#' @return A `cm_model` object.
#' @export
cm_model <- function(diseases, associations = NULL, impacts = NULL,
                     interactions = NULL,
                     missing_associations = c("independent", "unknown")) {
  missing_associations <- match.arg(missing_associations)
  if (!inherits(diseases, "cm_diseases")) cm_abort("`diseases` must be created with cm_diseases().")
  ids <- diseases$id
  if (!is.null(associations)) {
    if (!inherits(associations, "cm_associations")) cm_abort("`associations` must be created with cm_associations().")
    unk <- setdiff(c(associations$disease1, associations$disease2), ids)
    if (length(unk)) cm_abort(sprintf("Associations refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
  }
  if (!is.null(impacts)) {
    if (!inherits(impacts, "cm_impacts")) cm_abort("`impacts` must be created with cm_impacts().")
    unk <- setdiff(impacts$disease, ids)
    if (length(unk)) cm_abort(sprintf("Impacts refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
    for (o in unique(impacts$outcome)) {
      miss <- setdiff(ids, impacts$disease[impacts$outcome == o])
      if (length(miss)) {
        cm_abort(sprintf("Outcome '%s' has no impact for: %s. Use 0 for no impact.",
                         o, paste(miss, collapse = ", ")))
      }
    }
    adj <- unlist(lapply(impacts$adjusted_for, split_ids))
    unk <- setdiff(adj, ids)
    if (length(unk)) cm_abort(sprintf("`adjusted_for` refers to unknown diseases: %s.", paste(unk, collapse = ", ")))
  }
  if (!is.null(interactions)) {
    if (!inherits(interactions, "cm_interactions")) cm_abort("`interactions` must be created with cm_interactions().")
    unk <- setdiff(c(interactions$disease1, interactions$disease2), ids)
    if (length(unk)) cm_abort(sprintf("Interactions refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
    if (is.null(impacts) || length(setdiff(interactions$outcome, impacts$outcome))) {
      cm_abort("Every interaction outcome must also appear in `impacts`.")
    }
    if (any(impacts$scale[impacts$outcome %in% interactions$outcome] != "proportion")) {
      cm_abort("Interactions are only supported for outcomes on the proportion scale.")
    }
  }
  structure(list(diseases = diseases, associations = associations,
                 impacts = impacts, interactions = interactions,
                 missing_associations = missing_associations),
            class = "cm_model")
}

check_model <- function(model) {
  if (!inherits(model, "cm_model")) cm_abort("`model` must be created with cm_model().")
  invisible(model)
}
