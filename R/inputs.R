#' Describe the diseases in a system
#'
#' One row per disease: how common it is in the whole population. Rows about
#' two or more diseases (associations, combinations, prevalences within
#' subgroups) belong in [cm_associations()].
#'
#' @param disease Character vector of unique disease ids (no `|`, `;` or
#'   `:`; `"all"` and `"others"` are reserved).
#' @param value Numeric vector of disease occurrence values, interpreted
#'   according to `measure`.
#' @param measure Character, one per disease (or recycled):
#'   * `"prevalence"`: proportion of animals affected at a representative
#'     point in time (used as the probability directly).
#'   * `"probability"`: probability of at least one case within
#'     `time_horizon` (used directly).
#'   * `"incidence_rate"`: expected number of cases per animal within
#'     `time_horizon` (e.g. `0.48` for 48 cases per 100 lactations). Converted
#'     to a probability assuming Poisson-distributed events,
#'     `1 - exp(-value)`, as in Rasmussen et al. (2024), eq. 1.
#' @param label Optional full names of the diseases (e.g. `"lameness"`).
#' @param time_horizon Character label of the period the probabilities refer
#'   to (e.g. `"lactation"`, `"year"`). Probabilities, associations and impacts
#'   should all refer to the same period.
#' @param reference_population Optional description of the population the
#'   values come from.
#' @param source,note Optional citation and free-text note.
#' @param ... Not used (arguments of earlier versions give a message).
#'
#' @return A `cm_diseases` data frame with the converted probability in
#'   column `prob` and the inputs retained.
#' @export
#' @examples
#' cm_diseases(c("SCK", "PTB"), c(0.4789, 0.1001),
#'             measure = c("incidence_rate", "prevalence"),
#'             label = c("subclinical ketosis", "paratuberculosis"),
#'             time_horizon = "lactation")
cm_diseases <- function(disease, value, measure = "prevalence", label = NA_character_,
                        time_horizon = NA_character_, reference_population = NA_character_,
                        source = NA_character_, note = NA_character_, ...) {
  old_args(list(...), c(type = "`type` is now called `measure` (prevalence, probability or incidence_rate).",
                        id = "`id` is now called `disease`."))
  disease <- as.character(disease)
  n <- length(disease)
  if (n < 1L) cm_abort("At least one disease is required.", class = "deconflate_invalid_input")
  check_ids(disease)
  if (length(value) != n) cm_abort("`value` must have one entry per disease.", class = "deconflate_invalid_input")
  check_numeric(value, "value")
  if (any(value < 0)) cm_abort("`value` must be non-negative.", class = "deconflate_invalid_input")
  measure <- as.character(recycle_arg(measure, n, "measure"))
  check_choices(measure, disease_measures, "measure")

  prob <- ifelse(measure == "incidence_rate", 1 - exp(-value), value)
  if (any(prob <= 0 | prob >= 1)) {
    cm_abort(sprintf("Disease probabilities must lie strictly between 0 and 1 (check: %s).",
                     paste(disease[prob <= 0 | prob >= 1], collapse = ", ")),
             class = "deconflate_infeasible")
  }
  horizon <- as.character(recycle_arg(time_horizon, n, "time_horizon"))
  if (length(unique(stats::na.omit(horizon))) > 1L) {
    cm_warn("Diseases refer to different time horizons; probabilities, associations and impacts should share one period.")
  }
  out <- data.frame(
    disease = disease, prob = prob, value = value, measure = measure,
    label = as.character(recycle_arg(label, n, "label")),
    time_horizon = horizon,
    reference_population = as.character(recycle_arg(reference_population, n, "reference_population")),
    source = as.character(recycle_arg(source, n, "source")),
    note = as.character(recycle_arg(note, n, "note")),
    stringsAsFactors = FALSE
  )
  class(out) <- c("cm_diseases", "data.frame")
  out
}

disease_measures <- c("prevalence", "probability", "incidence_rate")

#' Describe associations between diseases
#'
#' Each row is one estimate from a study: an association between two
#' diseases, or between combinations of diseases, possibly within a subgroup
#' or adjusted for other diseases. Several rows can describe the same
#' association (they are pooled). Pairs that no row mentions are unknown: the
#' model fills in their association from the others (see [deconflate()]). To
#' state that two diseases are unrelated, give the pair an odds ratio of 1.
#'
#' @section What a row means:
#' `disease` (G) is the outcome and `associated_with` (H) the comparison:
#' animals with all of H against animals with none of H (for one disease,
#' with against without). Each can be one disease or a combination, with ids
#' separated by `";"` (e.g. `"MAS;SCK"`). `among` narrows the population to a
#' subgroup: ids that are present (`"SCK"`), absent (`"no SCK"`), both
#' (`"SCK;no MET"`), or `"no others"` (none of the diseases not named in the
#' row).
#'
#' * `"OR"` (default), `"RR"`, `"RD"`, `"phi"`: odds ratio, risk ratio
#'   `P(G | H) / P(G | no H)`, risk difference and phi correlation of the
#'   2x2 table of G by H, within `among`.
#' * `"prevalence"`: `P(G | H)` within `among`; without `associated_with`,
#'   the prevalence of G (a combination, or a disease within a subgroup).
#'   The prevalence of one disease in the whole population belongs in
#'   [cm_diseases()].
#' * `"OR_ratio"`: the odds ratio of G and H among `among` divided by their
#'   odds ratio among the complement of `among` (one disease, present or
#'   absent): how much a third disease changes an association, as a study
#'   that reports stratum-specific odds ratios would give it.
#'
#' `adjusted_for`: an odds ratio adjusted for other diseases of the model
#' (ids separated by `";"`, or `"all"` for every other disease), from a
#' logistic regression of G on H and those diseases, is an exact conditional
#' association and is used as such (G and H must be single diseases). Names
#' that are not diseases of the model (e.g. `"parity"`) are factors outside
#' the model: they are ignored, with a note, and a row adjusted only for them
#' is used as a crude association.
#'
#' Measures are applied to the modelled population's own probabilities
#' ("transported"): this assumes the measure is the same in the source study
#' and in the modelled population.
#'
#' @param disease,associated_with Disease ids or combinations (`"a;b"`).
#' @param value Numeric association values (required).
#' @param measure Character, one per row or recycled (see What a row means).
#' @param among Optional subgroups (see What a row means).
#' @param adjusted_for Optional adjustment sets (see What a row means).
#' @param study,source,note Optional study label (used in diagnostics),
#'   citation and free-text note.
#' @param ... Not used (arguments of earlier versions give a message). The
#'   arguments after it must be named.
#'
#' @return A `cm_associations` data frame (with the row number of each row
#'   in column `row`).
#' @export
#' @examples
#' cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
#' # An odds ratio within a subgroup, and how much d3 changes it
#' cm_associations(c("d1", "d1"), c("d2", "d2"), c(2.5, 1.8),
#'                 measure = c("OR", "OR_ratio"), among = c("d3", "d3"))
cm_associations <- function(disease, associated_with = NA_character_, value, measure = "OR", ...,
                            among = NA_character_, adjusted_for = NA_character_,
                            study = NA_character_, source = NA_character_,
                            note = NA_character_) {
  old_args(list(...), c(
    adjusted = "`adjusted` is no longer used: an association is adjusted when `adjusted_for` names diseases of the model.",
    disease1 = "`disease1` and `disease2` are now called `disease` and `associated_with`.",
    disease2 = "`disease1` and `disease2` are now called `disease` and `associated_with`."))
  disease <- as.character(disease)
  n <- length(disease)
  if (!n) cm_abort("At least one association row is required (or leave the associations out).", class = "deconflate_invalid_input")
  if (missing(value)) cm_abort("Give each association a `value`.", class = "deconflate_invalid_input")
  value <- as.numeric(recycle_arg(value, n, "value"))
  if (anyNA(value)) {
    cm_abort("Association `value` is missing for some rows. Leave a pair out (no row) if its association is unknown.", class = "deconflate_invalid_input")
  }
  if (any(!is.finite(value))) cm_abort("Association values must be finite.", class = "deconflate_invalid_input")
  measure <- as.character(recycle_arg(measure, n, "measure"))
  old <- intersect(measure, names(retired_measures))
  if (length(old)) cm_abort(retired_measures[[old[1]]], class = "deconflate_unsupported")
  check_choices(measure, association_measures, "measure")
  bad <- measure %in% c("OR", "RR", "OR_ratio") & value <= 0
  if (any(bad)) {
    cm_abort("Odds ratios, risk ratios and odds ratio ratios must be positive.", class = "deconflate_infeasible")
  }
  bad <- measure %in% c("RD", "phi") & (value < -1 | value > 1)
  if (any(bad)) cm_abort("Risk differences and phi must lie between -1 and 1.", class = "deconflate_infeasible")
  bad <- measure == "prevalence" & (value <= 0 | value >= 1)
  if (any(bad)) cm_abort("A prevalence must lie strictly between 0 and 1.", class = "deconflate_infeasible")
  out <- data.frame(
    row = seq_len(n), disease = disease,
    associated_with = clean_na(recycle_arg(associated_with, n, "associated_with")),
    among = clean_na(recycle_arg(among, n, "among")),
    adjusted_for = clean_na(recycle_arg(adjusted_for, n, "adjusted_for")),
    measure = measure, value = value,
    study = as.character(recycle_arg(study, n, "study")),
    source = as.character(recycle_arg(source, n, "source")),
    note = as.character(recycle_arg(note, n, "note")),
    stringsAsFactors = FALSE
  )
  class(out) <- c("cm_associations", "data.frame")
  out
}

# Association measures, and the measures of earlier versions with what to do
# instead.
association_measures <- c("OR", "RR", "RD", "phi", "prevalence", "OR_ratio")
retired_measures <- c(
  independent = "measure = 'independent' is no longer used: give the pair an odds ratio of 1 (value 1, measure OR).",
  unknown = "measure = 'unknown' is no longer used: leave the pair out. Pairs without a row are unknown, and the model fills in their association.",
  table = "measure = 'table' is no longer used: compute the odds ratio from the counts, n11 * n00 / (n10 * n01), and enter it with measure OR.",
  cond_prob = "measure = 'cond_prob' is now 'prevalence': P(disease | associated_with) is the prevalence of `disease` among animals with `associated_with`."
)

#' Describe the population: diseases and their associations
#'
#' A population holds the disease probabilities and the association rows.
#' Combine it with one impact table in [cm_model()]. Every row is checked
#' against the diseases (see [cm_associations()] for what a row means).
#'
#' @param diseases A [cm_diseases()] object.
#' @param associations Optional [cm_associations()] object. [deconflate()]
#'   needs at least one association; the sensitivity tools (e.g.
#'   [screen_associations()]) also run without any.
#' @param ... Not used (arguments of earlier versions give a message).
#' @return A `cm_population` object.
#' @export
#' @examples
#' pop <- cm_population(
#'   cm_diseases(c("d1", "d2", "d3"), c(0.10, 0.15, 0.20)),
#'   cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3))
#' )
#' pop
cm_population <- function(diseases, associations = NULL, ...) {
  old_args(list(...), c(three_way = three_way_message,
                        adjusted_associations = "`adjusted_associations` is no longer used: an association adjusted for diseases of the model is used exactly, and one adjusted only for other factors is used as crude (with a note)."))
  if (!inherits(diseases, "cm_diseases")) cm_abort("`diseases` must be created with cm_diseases().", class = "deconflate_invalid_input")
  ids <- diseases$disease
  if (!is.null(associations)) {
    if (!inherits(associations, "cm_associations")) {
      cm_abort("`associations` must be created with cm_associations().", class = "deconflate_invalid_input")
    }
    if (!nrow(associations)) associations <- NULL
  }
  if (!is.null(associations)) {
    if (is.null(associations$row)) associations$row <- seq_len(nrow(associations))
    for (r in seq_len(nrow(associations))) association_spec(associations[r, ], ids)
  }
  structure(list(diseases = diseases, associations = associations), class = "cm_population")
}

three_way_message <- "Three-way terms are now association rows: disease d1, associated_with d2, among d3, measure OR_ratio (the odds ratio of d1 and d2 among d3 divided by that among animals without d3)."
interactions_message <- "Interactions are now rows of the impact table: disease 'd1;d2' (the combination) with adjusted_for 'all' gives the interaction itself; a crude row compares animals with both diseases to animals with neither."

#' Describe raw impact estimates
#'
#' One impact table holds the estimates of one outcome, one row per
#' estimate. A disease can have several rows (e.g. from several studies, or
#' crude and adjusted), and a row can be about a combination of diseases or
#' a subgroup. There are two kinds of impact table:
#'
#' * **Additive impacts** (no `measure`): amounts in the outcome's own units,
#'   the same for every row (e.g. kg of milk, percent of yield, days, a
#'   welfare score). The package does not convert units; results come back
#'   in the units supplied. Adjust them with [deconflate()].
#' * **Event impacts** (`measure` given): comparisons of the risk of an event
#'   (e.g. death or culling), as hazard ratios (`"HR"`), rate ratios
#'   (`"rate_ratio"`), risk ratios (`"RR"`), odds ratios (`"OR"`) or risk
#'   differences (`"RD"`), which can be mixed. Adjust them with
#'   `deconflate(..., event_model = TRUE)` (snapshot hazard model).
#'
#' @section What a row means:
#' The model has an impact for each disease and an interaction for each
#' combination that appears in `disease`. A row compares animals with all
#' the diseases in `disease` against animals with none of them (for one
#' disease: with against without), within the subgroup `among` (as in
#' [cm_associations()], including `"no others"`, e.g. "only this disease
#' against healthy animals").
#'
#' * Crude (empty `adjusted_for`): the difference in the outcome between the
#'   two groups (event impacts: the ratio or difference of their hazards or
#'   risks).
#' * Adjusted (`adjusted_for` names diseases, separated by `";"`): additive
#'   impacts are coefficients of a linear regression of the outcome on the
#'   disease (for a combination: its members, their smaller combinations and
#'   the combination itself) and the diseases in `adjusted_for`; event
#'   impacts are comparisons within strata of those diseases (for a
#'   combination, the interaction contrast), combined with
#'   Mantel-Haenszel-type weights. `"all"` means every other disease and
#'   combination term of the model: the row is then the term's own value (a
#'   combination adjusted for `all` is its interaction; for event impacts
#'   this holds for hazard and rate ratios). Within a subgroup `among`
#'   naming present diseases, it is the sum of the term and the terms that
#'   combine it with those diseases (e.g. `a` among `c`, adjusted for all,
#'   is the impact of `a` plus the `a;c` interaction, when the model has
#'   one). Names that are not diseases of the model are factors outside the
#'   model, ignored with a note.
#'
#' Every disease needs at least one row of its own (for no effect: 0 for
#' additive impacts and risk differences, 1 for ratios). The rows must
#' determine every impact and interaction; [deconflate()] says which are not
#' determined otherwise. Several rows about the same impacts are fitted
#' together by weighted least squares (see [deconflate()]).
#'
#' A Cox hazard ratio estimated over follow-up is not exactly a snapshot
#' comparison (see [deconflate()], "Event impacts"); entering one is an
#' approximation and the user's assumption. Risk ratios, odds ratios and
#' risk differences refer to the period of the overall risk given to
#' [deconflate()].
#'
#' @param disease Character vector of disease ids or combinations (`"a;b"`).
#' @param value Numeric raw impacts.
#' @param among,adjusted_for Optional subgroups and adjustment sets (see What
#'   a row means).
#' @param measure `NULL` (default) for additive impacts; for event impacts,
#'   one measure per row or recycled.
#' @param study,source,note Optional study label, citation and note.
#' @param label,units Optional table-level metadata (e.g.
#'   `label = "milk yield loss"`, `units = "% of yield"`), carried into the
#'   results.
#' @param ... Not used (arguments of earlier versions give a message). The
#'   arguments after it must be named.
#'
#' @return A `cm_impacts` data frame (attributes `label`, `units` and
#'   `kind`, `"additive"` or `"event"`).
#' @export
#' @examples
#' cm_impacts(c("d1", "d2", "d3"), c(2.5, 5, 7.5), label = "yield", units = "%")
#' # Two estimates for d1, and the interaction of d1 and d2
#' cm_impacts(c("d1", "d1", "d2", "d3", "d1;d2"), c(2.5, 2.1, 5, 7.5, 1),
#'            adjusted_for = c(NA, "d2", NA, NA, "all"), units = "%")
#' # Event impacts: culling, with mixed measures
#' cm_impacts(c("d1", "d2", "d3"), c(1.5, 1.3, 0.04), measure = c("HR", "RR", "RD"),
#'            label = "culling")
cm_impacts <- function(disease, value, ..., among = NA_character_, adjusted_for = NA_character_,
                       measure = NULL, study = NA_character_, source = NA_character_,
                       note = NA_character_, label = NULL, units = NULL) {
  old_args(list(...), c(estimand = estimand_message))
  disease <- as.character(disease)
  n <- length(disease)
  if (!n) cm_abort("At least one impact is required.", class = "deconflate_invalid_input")
  if (length(value) != n) cm_abort("`value` must have one entry per row.", class = "deconflate_invalid_input")
  check_numeric(value, "value")
  event <- !is.null(measure)
  if (event) {
    measure <- as.character(recycle_arg(measure, n, "measure"))
    check_choices(measure, event_measures, "measure")
    bad <- measure != "RD" & value <= 0
    if (any(bad)) {
      cm_abort(sprintf("Ratios must be positive (check: %s).", paste(disease[bad], collapse = ", ")), class = "deconflate_invalid_input")
    }
    bad <- measure == "RD" & (value <= -1 | value >= 1)
    if (any(bad)) {
      cm_abort(sprintf("Risk differences must lie between -1 and 1 (check: %s).", paste(disease[bad], collapse = ", ")), class = "deconflate_invalid_input")
    }
  } else {
    measure <- rep(NA_character_, n)
  }
  out <- data.frame(row = seq_len(n), disease = disease,
                    among = clean_na(recycle_arg(among, n, "among")),
                    adjusted_for = clean_na(recycle_arg(adjusted_for, n, "adjusted_for")),
                    measure = measure, value = value,
                    study = as.character(recycle_arg(study, n, "study")),
                    source = as.character(recycle_arg(source, n, "source")),
                    note = as.character(recycle_arg(note, n, "note")),
                    stringsAsFactors = FALSE)
  attr(out, "label") <- label
  attr(out, "units") <- units
  attr(out, "kind") <- if (event) "event" else "additive"
  class(out) <- c("cm_impacts", "data.frame")
  out
}

estimand_message <- "`estimand` is no longer used: an estimate is crude when `adjusted_for` is empty and adjusted for the diseases it names otherwise."

# Measures of event impacts.
event_measures <- c("HR", "rate_ratio", "RR", "OR", "RD")

# Kind of an impact table: "event" when it has measures, else "additive".
impact_kind <- function(impacts) {
  if (is.null(impacts)) return("additive")
  k <- attr(impacts, "kind")
  if (!is.null(k)) return(k)
  if (!is.null(impacts$measure) && any(!is.na(impacts$measure))) "event" else "additive"
}

#' Combine a population with an impact table
#'
#' @param population A [cm_population()]. For convenience, a [cm_diseases()]
#'   object can be given instead, together with `associations`.
#' @param impacts A [cm_impacts()] object.
#' @param associations Used only when `population` is a [cm_diseases()]
#'   object; see [cm_population()].
#' @param distributions Optional named list of distributions ([distributions])
#'   for uncertain inputs, used when [deconflate()] runs draws. Names are
#'   keys: `"prob:<disease>"` (the disease's `value`, on the scale it was
#'   entered), `"associations:<row>"` and `"impacts:<row>"` (a row of those
#'   tables, by its number in `row`). For convenience, `"assoc:<d1>:<d2>"`
#'   names the crude row of a pair and `"impact:<id>"` the row of a disease
#'   or combination, when there is exactly one. [cm_read_inputs()] fills
#'   this in from the `dist` columns of the tables. When `population` is a
#'   model and `distributions` is not given, its distributions are kept,
#'   except those of impacts that are replaced.
#' @param ... Not used (arguments of earlier versions give a message).
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
cm_model <- function(population, impacts = NULL, associations = NULL, distributions = NULL, ...) {
  old_args(list(...), c(interactions = interactions_message, three_way = three_way_message,
                        adjusted_associations = "`adjusted_associations` is no longer used."))
  if (inherits(population, "cm_diseases")) population <- cm_population(population, associations)
  if (!inherits(population, "cm_population")) {
    cm_abort("`population` must be created with cm_population() (or be a cm_diseases object).", class = "deconflate_invalid_input")
  }
  ids <- population$diseases$disease
  if (!is.null(impacts)) {
    if (!inherits(impacts, "cm_impacts")) cm_abort("`impacts` must be created with cm_impacts().", class = "deconflate_invalid_input")
    if (is.null(impacts$row)) impacts$row <- seq_len(nrow(impacts))
    impacts <- check_impact_rows(impacts, ids)
  }
  out <- population
  out$impacts <- impacts
  class(out) <- c("cm_model", "cm_population")
  if (is.null(distributions)) {
    # Distributions carried over from a model belong to its tables: those of
    # replaced impacts are dropped.
    distributions <- population$distributions
    if (length(distributions) && inherits(population, "cm_model") &&
        !identical(impacts, population$impacts)) {
      distributions <- distributions[!startsWith(names(distributions), "impacts:")]
    }
  }
  out$distributions <- check_distributions(distributions, out)
  out
}

# Check the rows of an impact table against the diseases: every row parses,
# every disease has a row of its own. Returns the table.
check_impact_rows <- function(impacts, ids) {
  kind <- impact_kind(impacts)
  for (r in seq_len(nrow(impacts))) impact_spec(impacts[r, ], ids)
  own <- unique(impacts$disease[!grepl(";", impacts$disease, fixed = TRUE)])
  miss <- setdiff(ids, trimws(own))
  if (length(miss)) {
    cm_abort(sprintf("No impact for: %s. Every disease needs a row of its own (for no effect: %s).",
                     paste(miss, collapse = ", "),
                     if (kind == "event") "1 for ratios, 0 for risk differences" else "0"), class = "deconflate_invalid_input")
  }
  rownames(impacts) <- NULL
  attr(impacts, "kind") <- kind
  impacts
}

# Arguments of earlier versions, passed through `...`: stop with what to do
# instead.
old_args <- function(dots, messages) {
  if (!length(dots)) return(invisible(TRUE))
  nm <- names(dots) %||% rep("", length(dots))
  hit <- intersect(nm, names(messages))
  if (length(hit)) cm_abort(messages[[hit[1]]], class = "deconflate_unsupported")
  bad <- nm[nzchar(nm)]
  if (length(bad)) cm_abort(sprintf("Unknown argument%s: %s.", if (length(bad) > 1L) "s" else "", paste(bad, collapse = ", ")), class = "deconflate_invalid_input")
  cm_abort("Unexpected unnamed argument.", class = "deconflate_invalid_input")
}

clean_na <- function(x) {
  x <- trimws(as.character(x))
  x[!is.na(x) & !nzchar(x)] <- NA_character_
  x
}

# Check the distributions of a model: a named list of cm_dist objects whose
# keys match the model's inputs. Returns the list with canonical keys.
check_distributions <- function(d, model) {
  if (is.null(d) || !length(d)) return(NULL)
  if (!is.list(d) || is.null(names(d)) || any(!nzchar(names(d))) ||
      !all(vapply(d, inherits, logical(1), "cm_dist"))) {
    cm_abort("`distributions` must be a named list of distributions (see ?distributions).", class = "deconflate_invalid_input")
  }
  if (anyDuplicated(names(d))) {
    cm_abort(sprintf("`distributions` names an input more than once: %s.",
                     paste(unique(names(d)[duplicated(names(d))]), collapse = ", ")), class = "deconflate_invalid_input")
  }
  out <- list()
  for (k in names(d)) {
    ck <- canonical_key(k, model)
    if (!is.null(out[[ck]])) cm_abort(sprintf("The input '%s' is given more than once in `distributions`.", k), class = "deconflate_invalid_input")
    out[[ck]] <- d[[k]]
  }
  out
}

# Canonical form of an input key, checked against the model:
# "prob:<disease>", "associations:<row>", "impacts:<row>". Accepted aliases:
# "assoc:<d1>:<d2>" (the crude pair row of d1 and d2), "impact:<id>" (the row
# whose disease is <id>), "inter:<d1>:<d2>" (the row of the combination
# adjusted for all) and "three:<d1>:<d2>:<d3>" (the OR_ratio row of the
# triple), each when exactly one row matches.
canonical_key <- function(key, model) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  type <- parts[1]
  ids <- model$diseases$disease
  bad <- function(why) cm_abort(sprintf("Input '%s': %s.", key, why), class = "deconflate_invalid_input")
  by_row <- function(tab, what) {
    if (length(parts) != 2L) bad(sprintf("use '%s:<row>'", type))
    if (is.null(tab) || !nrow(tab)) bad(sprintf("the model has no %s", what))
    r <- suppressWarnings(as.integer(parts[2]))
    if (is.na(r) || !(r %in% tab$row)) bad(sprintf("no %s row %s", what, parts[2]))
    paste0(type, ":", r)
  }
  one <- function(hit, tab, what) {
    if (!sum(hit)) bad(sprintf("no %s row matches", what))
    if (sum(hit) > 1L) {
      bad(sprintf("%d %s rows match (rows %s); name one by its row, e.g. '%s:%d'", sum(hit), what,
                  paste(tab$row[hit], collapse = ", "), if (what == "impact") "impacts" else "associations",
                  tab$row[hit][1]))
    }
    tab$row[hit]
  }
  switch(type,
    prob = {
      if (length(parts) != 2L || !(parts[2] %in% ids)) bad("unknown disease")
      key
    },
    associations = by_row(model$associations, "association"),
    impacts = by_row(model$impacts, "impact"),
    assoc = {
      if (length(parts) != 3L) bad("use 'assoc:<d1>:<d2>'")
      paste0("associations:", one(pair_row_hits(model$associations, parts[2], parts[3], model$diseases$disease),
                                  model$associations, "association"))
    },
    impact = {
      if (length(parts) != 2L) bad("use 'impact:<disease>'")
      im <- model$impacts
      if (is.null(im)) bad("the model has no impacts")
      paste0("impacts:", one(same_set(im$disease, parts[2]), im, "impact"))
    },
    inter = {
      if (length(parts) < 3L) bad("use 'inter:<d1>:<d2>'")
      im <- model$impacts
      if (is.null(im)) bad("the model has no impacts")
      hit <- same_set(im$disease, paste(parts[-1], collapse = ";")) & is_all(im$adjusted_for) & is.na(im$among)
      paste0("impacts:", one(hit, im, "impact"))
    },
    three = {
      if (length(parts) != 4L) bad("use 'three:<d1>:<d2>:<d3>'")
      a <- model$associations
      if (is.null(a)) bad("the model has no associations")
      hit <- a$measure == "OR_ratio" & vapply(seq_len(nrow(a)), function(r) {
        setequal(c(a$disease[r], a$associated_with[r], sub("^no\\s+", "", a$among[r])), parts[2:4])
      }, logical(1))
      paste0("associations:", one(hit, a, "association"))
    },
    bad("unknown input type (use prob:, associations:, impacts:, assoc:, impact:, inter: or three:)")
  )
}

# Rows of an association table that are the crude pair row of d1 and d2.
# With the model's disease ids, a row adjusted only for factors outside the
# model (or for its own diseases) counts as crude, as in the model.
pair_row_hits <- function(a, d1, d2, ids = NULL) {
  if (is.null(a) || !nrow(a)) return(logical(0))
  adjusted <- vapply(a$adjusted_for, function(x) {
    x <- split_ids(x)
    if (!length(x)) return(FALSE)
    if (any(is_all(x)) || is.null(ids)) return(TRUE)
    any(x %in% setdiff(ids, c(d1, d2)))
  }, logical(1), USE.NAMES = FALSE)
  simple <- is.na(a$among) & !adjusted &
    a$measure %in% c("OR", "RR", "RD", "phi", "prevalence") & !is.na(a$associated_with) &
    !grepl(";", a$disease, fixed = TRUE) & !grepl(";", a$associated_with, fixed = TRUE)
  simple & ((trimws(a$disease) == d1 & trimws(a$associated_with) == d2) |
              (trimws(a$disease) == d2 & trimws(a$associated_with) == d1))
}

# Do the entries of x name the same set of ids as `target` ("a;b" = "b;a")?
same_set <- function(x, target) {
  t <- sort(split_ids(target))
  vapply(x, function(v) identical(sort(split_ids(v)), t), logical(1), USE.NAMES = FALSE)
}

is_all <- function(x) !is.na(x) & tolower(trimws(x)) == "all"

check_population <- function(x) {
  if (!inherits(x, "cm_population")) {
    cm_abort("Expected a cm_population() or cm_model() object.", class = "deconflate_invalid_input")
  }
  invisible(x)
}

check_model <- function(model) {
  if (!inherits(model, "cm_model")) cm_abort("`model` must be created with cm_model().", class = "deconflate_invalid_input")
  invisible(model)
}
