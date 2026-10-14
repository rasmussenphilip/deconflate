# Sensitivity and scenario tools ------------------------------------------------

# Metric used by the screening tools: the total and each disease's
# contribution. Additive impacts: the adjusted aggregate and the
# contributions (in the impacts' own units). Event impacts: the risk
# attributable to disease and its Shapley allocation.
burden_metric <- function(result) {
  if (inherits(result, "cm_event_result")) {
    ids <- result$adjusted$disease[result$adjusted$kind %||% rep("disease", nrow(result$adjusted)) == "disease"]
    by <- result$attributable$by_disease
    return(list(total = result$attributable$summary$attributable,
                by_disease = stats::setNames(if (is.null(by)) rep(NA_real_, length(ids)) else by$attributable, ids)))
  }
  list(total = result$totals$adjusted_total,
       by_disease = stats::setNames(result$contributions$total, result$contributions$disease))
}

metric_label <- function(event_model) if (isTRUE(event_model)) "attributable risk" else "aggregate"

# Common checks of the sensitivity tools: the model, the kind of impacts and
# the overall or baseline risk (its central value and kind).
sens_setup <- function(model, event_model, overall_risk, dots = list(), baseline_risk = NULL) {
  check_model(model)
  retired_args(dots)
  if (length(dots)) cm_abort(sprintf("Unknown argument%s: %s.", if (length(dots) > 1L) "s" else "",
                                     paste(names(dots), collapse = ", ")), class = "deconflate_invalid_input")
  check_event_model(model, event_model)
  risk <- check_overall_risk(overall_risk, event_model, baseline_risk)
  list(event_model = event_model, risk = risk$value, risk_kind = risk$kind %||% "overall")
}

compare_to_base <- function(m, base) {
  if (is.null(m)) {
    return(data.frame(total = NA_real_, change = NA_real_, rel_change = NA_real_,
                      max_rank_shift = NA_integer_, rank_corr = NA_real_))
  }
  rb <- rank(-base$by_disease, ties.method = "min")
  rm <- rank(-m$by_disease, ties.method = "min")
  data.frame(total = m$total, change = m$total - base$total,
             rel_change = if (isTRUE(abs(base$total) > 0)) (m$total - base$total) / abs(base$total) else NA_real_,
             max_rank_shift = as.integer(max(abs(rb - rm))),
             rank_corr = suppressWarnings(stats::cor(base$by_disease, m$by_disease,
                                                     method = "spearman")))
}

# Adjust a scenario as deconflate() would (point estimates), returning an
# empty list with the reason as an attribute when the scenario fails.
run_quiet <- function(model, st, joint = NULL) {
  tryCatch(
    withCallingHandlers({
      plan <- plan_model(model, st$event_model, has_joint = !is.null(joint), risk_kind = st$risk_kind)
      if (!is.null(joint)) plan$shortcut <- FALSE
      run_point(model, plan, risk = st$risk, joint = joint, warn = FALSE, raw = FALSE)
    }, deconflate_nonconvergence = function(w) {
      if (inherits(w, "warning")) invokeRestart("muffleWarning")
    }),
    deconflate_error = function(e) structure(list(), reason = condition_type(e))
  )
}
failed_run <- function(res) {
  is.null(res) || !(inherits(res, "cm_result") || inherits(res, "cm_event_result")) || !result_is_finite(res)
}

# Metric of a scenario run, or NULL with the reason the scenario failed.
scenario_metric <- function(res) {
  if (failed_run(res)) {
    reason <- if (inherits(res, "cm_result") || inherits(res, "cm_event_result")) {
      sprintf("undefined: non-finite adjusted values (%s)", nonfinite_diseases(res))
    } else {
      attr(res, "reason") %||% "failed"
    }
    return(list(met = NULL, reason = reason))
  }
  list(met = burden_metric(res), reason = NA_character_)
}

# One row of a screen: the scenario's labels, its comparison with the
# baseline and the reason it failed (NA if it ran).
screen_row <- function(labels, res, base) {
  sm <- scenario_metric(res)
  cbind(labels, compare_to_base(sm$met, base),
        data.frame(failed = sm$reason, stringsAsFactors = FALSE))
}

# A screen compares scenarios with the baseline, so the baseline itself must
# be a finite estimate.
check_baseline <- function(res) {
  if (!inherits(res, "cm_result") && !inherits(res, "cm_event_result")) {
    reason <- attr(res, "reason") %||% "failed"
    cm_abort(sprintf("The baseline could not be adjusted (%s).", reason), class = type_class(reason))
  }
  if (!result_is_finite(res)) {
    cm_abort(sprintf("The baseline result is undefined (non-finite adjusted impacts for %s).",
                     nonfinite_diseases(res)), class = "deconflate_nonfinite")
  }
  res
}

check_disease_ids <- function(ids, model) {
  unk <- setdiff(ids, model$diseases$disease)
  if (length(unk)) cm_abort(sprintf("Unknown diseases: %s.", paste(unk, collapse = ", ")), class = "deconflate_invalid_input")
  invisible(ids)
}

#' Set or replace one association
#'
#' Replaces the crude pair rows of a pair (one disease on each side, no
#' subgroup, not adjusted) by one row with the given value, or adds that row.
#' Other rows about the pair (e.g. within subgroups) are kept.
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param disease1,disease2 Disease ids.
#' @param value Association value.
#' @param measure Measure (see [cm_associations()]).
#' @return The modified object. Distributions of replaced rows are dropped
#'   (the scenario value is fixed).
#' @export
set_association <- function(model, disease1, disease2, value, measure = "OR") {
  check_population(model)
  check_disease_ids(c(disease1, disease2), model)
  new <- cm_associations(disease1, disease2, value, measure = measure, source = "scenario")
  replace_rows(model, "associations", pair_row_hits(model$associations, disease1, disease2, model$diseases$disease), new)
}

# Replace the rows `hit` of a table by `new` (one row; it takes the first
# replaced row's number, or the next free one), dropping their
# distributions.
replace_rows <- function(model, table, hit, new) {
  tab <- model[[table]]
  if (is.null(tab) || !nrow(tab)) {
    new$row <- 1L
    model[[table]] <- new
    return(model)
  }
  hit <- which(hit)
  new$row <- if (length(hit)) tab$row[hit[1]] else max(tab$row) + 1L
  if (length(model$distributions)) {
    drop <- names(model$distributions) %in% paste0(table, ":", tab$row[hit])
    model$distributions <- if (any(!drop)) model$distributions[!drop]
  }
  keep <- setdiff(seq_len(nrow(tab)), hit)
  at <- attributes(tab)
  at <- at[intersect(names(at), c("label", "units", "kind"))]
  cols <- names(tab)
  for (cc in setdiff(cols, names(new))) new[[cc]] <- NA
  out <- rbind(tab[keep, cols, drop = FALSE], new[, cols, drop = FALSE])
  out <- out[order(out$row), , drop = FALSE]
  rownames(out) <- NULL
  for (k in names(at)) attr(out, k) <- at[[k]]
  class(out) <- class(tab)
  model[[table]] <- out
  model
}

#' Set or replace one pairwise interaction
#'
#' Sets the interaction of two diseases to `value`: a row for the
#' combination adjusted for `all`, which is the interaction itself (see
#' [cm_impacts()]). Replaces such a row if present. For additive impacts the
#' value is in the units of the impacts (0: no interaction); for event
#' impacts it is a hazard ratio multiplier of the combination (1: no
#' interaction).
#'
#' @param model A [cm_model()].
#' @param disease1,disease2 Disease ids.
#' @param value Interaction value: in the units of the impacts, or for event
#'   impacts a positive hazard ratio multiplier.
#' @return The modified model. Distributions of replaced rows are dropped.
#' @export
set_interaction <- function(model, disease1, disease2, value) {
  check_model(model)
  event <- impact_kind(model$impacts) == "event"
  check_disease_ids(c(disease1, disease2), model)
  if (disease1 == disease2) cm_abort("An interaction needs two different diseases.", class = "deconflate_invalid_input")
  if (event && !(is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0)) {
    cm_abort("For event impacts, the interaction is a hazard ratio multiplier: a single positive number (1: no interaction).",
             class = "deconflate_invalid_input")
  }
  im <- model$impacts
  ids <- model$diseases$disease
  combo <- term_label(ids[ids %in% c(disease1, disease2)])
  hit <- same_set(im$disease, combo) & is_all(im$adjusted_for) & is.na(im$among)
  new <- if (event) {
    cm_impacts(combo, value, adjusted_for = "all", measure = "HR", source = "scenario")
  } else {
    cm_impacts(combo, value, adjusted_for = "all", source = "scenario")
  }
  out <- replace_rows(model, "impacts", hit, new)
  attr(out$impacts, "kind") <- if (event) "event" else "additive"
  out
}

#' Set or replace one three-way association
#'
#' Sets how much `disease3` changes the association of `disease1` and
#' `disease2`: an `OR_ratio` row (the odds ratio of disease1 and disease2
#' among animals with disease3, divided by that among animals without it;
#' see [cm_associations()]). Replaces such a row if present.
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param disease1,disease2,disease3 Disease ids.
#' @param ratio Positive ratio of the two odds ratios.
#' @return The modified object. Distributions of replaced rows are dropped.
#' @export
set_three_way <- function(model, disease1, disease2, disease3, ratio) {
  check_population(model)
  check_disease_ids(c(disease1, disease2, disease3), model)
  if (anyDuplicated(c(disease1, disease2, disease3))) cm_abort("A three-way association needs three different diseases.", class = "deconflate_invalid_input")
  a <- model$associations
  hit <- if (is.null(a)) logical(0) else a$measure == "OR_ratio" & vapply(seq_len(nrow(a)), function(r) {
    setequal(c(a$disease[r], a$associated_with[r], sub("^no\\s+", "", a$among[r])), c(disease1, disease2, disease3))
  }, logical(1))
  new <- cm_associations(disease1, disease2, ratio, measure = "OR_ratio", among = disease3, source = "scenario")
  replace_rows(model, "associations", hit, new)
}

#' Screen disease pairs for influential associations
#'
#' Re-runs the adjustment with each disease pair's association changed, one
#' pair at a time, and reports how much the total and the ranking of
#' diseases change. Pairs without an association (unknown) are set to each
#' of `or_values`; given associations of any measure are multiplied by each
#' of `multipliers`. Scenarios whose inputs are invalid or infeasible are
#' kept, with `NA` results and the reason, so that no scenario is skipped
#' silently. This identifies associations worth estimating, as in Rasmussen
#' et al. (2022), Fig. 3 and Rasmussen et al. (2024), Fig. 7C-D.
#'
#' The model needs no associations: with none, the baseline has every
#' disease independent (adjusted impacts equal the raw ones), and the screen
#' shows which pairs would matter if they were associated. Each scenario is
#' adjusted as [deconflate()] would (point estimates).
#'
#' The total is the adjusted aggregate for additive impacts, and the risk
#' attributable to disease for event impacts (`event_model = TRUE`); disease
#' rankings use the contributions to it.
#'
#' @param model A [cm_model()].
#' @param event_model,overall_risk,baseline_risk As in [deconflate()] (a distribution of
#'   the risk is used at its mean).
#' @param or_values Odds ratios tried for pairs without a crude pair row.
#' @param multipliers Factors applied to the crude pair rows of a pair (all
#'   of them, when a pair has several).
#' @param pairs Optional character vector of pairs (`"d1:d2"`); default all.
#' @param ... Not used.
#' @return A `cm_screen` data frame sorted by the absolute relative change in
#'   the total, with the scenario, total, change, relative change, the
#'   largest shift in any disease's rank, the Spearman correlation of disease
#'   contributions with the baseline, and `failed` (the reason a scenario
#'   could not be run).
#' @export
#' @examples
#' screen_associations(example_supplement(), or_values = c(0.5, 3))
#'
#' # Without any association estimates: which pairs would matter?
#' m <- example_supplement()
#' m$associations <- NULL
#' screen_associations(m, or_values = c(0.5, 2, 4))
screen_associations <- function(model, event_model = FALSE, overall_risk = NULL,
                                or_values = c(0.5, 2), multipliers = c(0.5, 2),
                                pairs = NULL, baseline_risk = NULL, ...) {
  st <- sens_setup(model, event_model, overall_risk, list(...), baseline_risk)
  base <- burden_metric(check_baseline(run_quiet(model, st)))
  pt <- pair_tables(model)
  if (!is.null(pairs)) {
    check_pair_ids(pairs, model$diseases$disease)
    ab <- strsplit(pairs, ":", fixed = TRUE)
    want <- vapply(ab, function(x) pair_key(x[1], x[2]), character(1))
    pt <- pt[pair_key(pt$disease1, pt$disease2) %in% want, , drop = FALSE]
  }
  key <- paste(pt$disease1, pt$disease2, sep = ":")
  a <- model$associations
  rows <- list()
  for (r in seq_len(nrow(pt))) {
    if (pt$status[r] == "specified") {
      hit <- which(pair_row_hits(a, pt$disease1[r], pt$disease2[r], model$diseases$disease))
      meas <- unique(a$measure[hit])
      labels <- sprintf("%s x %g", paste(meas, collapse = "/"), multipliers)
      setter <- function(x) {   # x: the multiplier
        m2 <- model
        for (h in hit) m2 <- set_inputs(m2, stats::setNames(a$value[h] * x, paste0("associations:", a$row[h])))
        m2
      }
      vals <- multipliers
    } else {
      labels <- sprintf("OR = %g", or_values)
      d1 <- pt$disease1[r]
      d2 <- pt$disease2[r]
      setter <- function(x) set_association(model, d1, d2, x, "OR")   # x: the odds ratio
      vals <- or_values
    }
    for (v in seq_along(vals)) {
      m2 <- tryCatch(setter(vals[v]),
                     deconflate_error = function(e) structure(list(), reason = "invalid value"))
      res <- if (inherits(m2, "cm_model")) run_quiet(m2, st) else m2
      rows[[length(rows) + 1L]] <- screen_row(
        data.frame(pair = key[r], status = pt$status[r], scenario = labels[v],
                   stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base, event_model)
}

# Pairs given as "d1:d2" must name two different diseases of the model.
check_pair_ids <- function(pairs, ids) {
  ab <- strsplit(as.character(pairs), ":", fixed = TRUE)
  bad <- !vapply(ab, function(x) length(x) == 2L && all(x %in% ids) && x[1] != x[2], logical(1))
  if (any(bad)) {
    cm_abort(sprintf("Unknown or malformed pairs (use 'd1:d2' with disease ids of the model): %s.",
                     paste(pairs[bad], collapse = ", ")), class = "deconflate_invalid_input")
  }
  invisible(pairs)
}

finish_screen <- function(rows, base, event_model = FALSE) {
  if (!length(rows)) {
    out <- data.frame(pair = character(0), status = character(0), scenario = character(0),
                      total = numeric(0), change = numeric(0), rel_change = numeric(0),
                      max_rank_shift = integer(0), rank_corr = numeric(0),
                      failed = character(0), stringsAsFactors = FALSE)
  } else {
    out <- do.call(rbind, rows)
    out <- out[order(-abs(out$rel_change), na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL
  }
  attr(out, "baseline_total") <- base$total
  attr(out, "metric") <- metric_label(event_model)
  class(out) <- c("cm_screen", "data.frame")
  out
}

#' Screen disease pairs for influential impact interactions
#'
#' Sets the interaction of each disease pair to each size in `values`, one
#' at a time (see [set_interaction()]), re-runs the adjustment and reports
#' the change in the aggregate (for event impacts, the attributable risk)
#' and rankings. The joint distribution is fitted once and reused.
#'
#' @param model A [cm_model()].
#' @param values Interaction sizes: for additive impacts in the units of the
#'   impacts (positive synergistic, negative antagonistic); for event impacts
#'   hazard ratio multipliers of the combination (above 1 synergistic, below
#'   1 antagonistic).
#' @param pairs Optional character vector of pairs (`"d1:d2"`).
#' @param event_model,overall_risk,baseline_risk As in [deconflate()].
#' @param ... Passed to [fit_joint()] (e.g. `backend = "sampled"` for many
#'   diseases).
#' @return A `cm_screen` data frame as in [screen_associations()].
#' @export
screen_interactions <- function(model, values, pairs = NULL, event_model = FALSE, overall_risk = NULL,
                                baseline_risk = NULL, ...) {
  check_model(model)
  st <- sens_setup(model, event_model, overall_risk, baseline_risk = baseline_risk)
  if (event_model && (!is.numeric(values) || any(!is.finite(values)) || any(values <= 0))) {
    cm_abort("For event impacts, `values` are hazard ratio multipliers: positive numbers (1: no interaction).",
             class = "deconflate_invalid_input")
  }
  # One joint distribution serves every scenario (interactions do not change
  # it); `...` goes to fit_joint(), e.g. backend = "sampled" (the default
  # beyond 20 diseases, as in deconflate()).
  fa <- list(...)
  if (nrow(model$diseases) > 20L && is.null(fa$backend)) fa$backend <- "sampled"
  joint <- do.call(fit_joint, c(list(as_population(model)), fa))
  base <- burden_metric(check_baseline(run_quiet(model, st, joint = joint)))
  ids <- model$diseases$disease
  pr <- if (!is.null(pairs)) {
    check_pair_ids(pairs, ids)
    matrix(unlist(strsplit(pairs, ":", fixed = TRUE)), ncol = 2, byrow = TRUE)
  } else if (length(ids) >= 2L) {
    t(utils::combn(ids, 2))
  } else {
    matrix(character(0), 0, 2)
  }
  rows <- list()
  for (r in seq_len(nrow(pr))) {
    for (v in values) {
      res <- run_quiet(set_interaction(model, pr[r, 1], pr[r, 2], v), st, joint = joint)
      rows[[length(rows) + 1L]] <- screen_row(
        data.frame(pair = paste(pr[r, 1], pr[r, 2], sep = ":"), status = "interaction",
                   scenario = sprintf(if (event_model) "interaction HR = %g" else "delta = %g", v),
                   stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base, event_model)
}

#' Screen three-way association scenarios
#'
#' Pairwise associations do not identify how often three diseases occur
#' together; the model assumes as little as possible (maximum entropy). This
#' screen sets, for each triple in turn, how much the third disease changes
#' the association of the first two (an `OR_ratio` row, see
#' [set_three_way()]), refits the joint distribution (the other rows are
#' still matched) and reports the change in the aggregate and rankings.
#'
#' When the pairs of a triple all have crude pair rows, such a row keeps
#' their tables fixed, so it changes additive results only through
#' interactions or rows about combinations or subgroups. When a pair is
#' unknown, its fitted table changes, and additive results can change even
#' without interactions. Three-way associations also matter for event
#' impacts (`event_model = TRUE`).
#'
#' @param model A [cm_model()].
#' @param event_model,overall_risk,baseline_risk As in [deconflate()].
#' @param ratios Odds ratio ratios to try.
#' @param triples Optional list of character vectors of three disease ids
#'   (the association of the first two, among the third); default all
#'   triples (which can be many).
#' @return A `cm_screen` data frame as in [screen_associations()].
#' @export
screen_three_way <- function(model, ratios = c(0.5, 2), triples = NULL, event_model = FALSE,
                             overall_risk = NULL, baseline_risk = NULL) {
  st <- sens_setup(model, event_model, overall_risk, baseline_risk = baseline_risk)
  base <- burden_metric(check_baseline(run_quiet(model, st)))
  ids <- model$diseases$disease
  if (is.null(triples)) {
    triples <- if (length(ids) >= 3L) utils::combn(ids, 3, simplify = FALSE) else list()
  } else {
    if (is.character(triples)) triples <- list(triples)
    for (tr in triples) {
      if (length(tr) != 3L || anyDuplicated(tr) || !all(tr %in% ids)) {
        cm_abort(sprintf("Each triple must name three different diseases of the model (check: %s).",
                         paste(tr, collapse = ", ")), class = "deconflate_invalid_input")
      }
    }
  }
  rows <- list()
  for (tr in triples) {
    for (rt in ratios) {
      m2 <- set_three_way(model, tr[1], tr[2], tr[3], rt)
      res <- run_quiet(m2, st)
      rows[[length(rows) + 1L]] <- screen_row(
        data.frame(pair = paste(tr, collapse = ":"), status = "three-way",
                   scenario = sprintf("ratio = %g", rt), stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base, event_model)
}

#' One-at-a-time sensitivity analysis
#'
#' Varies each input by `variation` (e.g. +/- 20%) with everything else
#' fixed, as in Rasmussen et al. (2024), Fig. 7, and reports the resulting
#' range of the total: the adjusted aggregate for additive impacts, the risk
#' attributable to disease for event impacts.
#'
#' @param model A [cm_model()].
#' @param event_model,overall_risk,baseline_risk As in [deconflate()].
#' @param variation Relative change applied downwards and upwards.
#' @param inputs Which inputs to vary: `"prob"` (disease probabilities),
#'   `"assoc"` (association rows), `"impact"` (impact rows), and for event
#'   impacts `"risk"` (the overall or baseline risk, whichever is given).
#' @param ... Not used.
#' @return A `cm_oat` data frame with one row per input: the input key
#'   (`"prob:<disease>"`, `"associations:<row>"`, `"impacts:<row>"` or
#'   `"risk"`), totals at the low and high values and the swing, sorted by
#'   swing.
#' @export
#' @examples
#' sensitivity_oat(example_supplement())
sensitivity_oat <- function(model, event_model = FALSE, overall_risk = NULL,
                            variation = 0.2, inputs = c("prob", "assoc", "impact", "risk"),
                            baseline_risk = NULL, ...) {
  st <- sens_setup(model, event_model, overall_risk, list(...), baseline_risk)
  inputs <- cm_match_arg(inputs, several_ok = TRUE)
  base <- burden_metric(check_baseline(run_quiet(model, st)))
  vals <- model_inputs(model)
  if (event_model) vals <- c(vals, risk = st$risk)
  keys <- names(vals)
  type <- c(prob = "prob", associations = "assoc", impacts = "impact", risk = "risk")[sub(":.*$", "", keys)]
  keep <- type %in% inputs & vals != 0
  keys <- keys[keep]
  rows <- lapply(keys, function(k) {
    x <- vals[[k]]
    tot <- vapply(c(1 - variation, 1 + variation), function(f) {
      st2 <- st
      m2 <- model
      if (k == "risk") {
        if (!(x * f > 0 && x * f < 1)) return(NA_real_)
        st2$risk <- x * f
      } else {
        m2 <- tryCatch(set_inputs(model, stats::setNames(x * f, k)), deconflate_error = function(e) NULL)
        if (is.null(m2)) return(NA_real_)
      }
      sm <- scenario_metric(run_quiet(m2, st2))
      if (is.null(sm$met)) NA_real_ else sm$met$total
    }, numeric(1))
    data.frame(input = k, value = x, total_low = tot[1], total_high = tot[2],
               swing = abs(tot[2] - tot[1]),
               rel_swing = if (isTRUE(abs(base$total) > 0)) abs(tot[2] - tot[1]) / abs(base$total) else NA_real_,
               stringsAsFactors = FALSE)
  })
  if (!length(rows)) {
    rows <- list(data.frame(input = character(0), value = numeric(0), total_low = numeric(0),
                            total_high = numeric(0), swing = numeric(0), rel_swing = numeric(0),
                            stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out <- out[order(-out$swing, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "baseline_total") <- base$total
  attr(out, "metric") <- metric_label(event_model)
  class(out) <- c("cm_oat", "data.frame")
  out
}

# The numeric inputs of a model, keyed as in cm_model(distributions = ).
model_inputs <- function(model) {
  v <- stats::setNames(model$diseases$value, paste0("prob:", model$diseases$disease))
  a <- model$associations
  if (!is.null(a) && nrow(a)) v <- c(v, stats::setNames(a$value, paste0("associations:", a$row)))
  i <- model$impacts
  if (!is.null(i) && nrow(i)) v <- c(v, stats::setNames(i$value, paste0("impacts:", i$row)))
  v
}
