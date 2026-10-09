# Sensitivity and scenario tools ------------------------------------------------

# Metric used by the screening tools: the total and each disease's
# contribution. Additive impacts: the adjusted aggregate and the
# contributions (in the impacts' own units). Event impacts: the risk
# attributable to disease and its Shapley allocation.
burden_metric <- function(result) {
  ids <- result$adjusted$disease
  if (inherits(result, "cm_event_result")) {
    by <- result$attributable$by_disease
    return(list(total = result$attributable$summary$attributable,
                by_disease = stats::setNames(if (is.null(by)) rep(NA_real_, length(ids)) else by$attributable, ids)))
  }
  list(total = result$totals$adjusted_total,
       by_disease = stats::setNames(result$contributions$total, ids))
}

metric_label <- function(event_model) if (isTRUE(event_model)) "attributable risk" else "aggregate"

# Common checks of the sensitivity tools: the model, the kind of impacts, the
# method and the overall risk (its central value).
sens_setup <- function(model, method, event_model, overall_risk, fun) {
  check_model(model)
  check_event_model(model, event_model)
  risk <- check_overall_risk(overall_risk, event_model)$value
  method <- public_method(method, c("auto", "simultaneous", "global"), fun)
  list(method = method, event_model = event_model, risk = risk)
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

# Adjust a scenario as deconflate() would (point estimates; the method is
# chosen for the scenario's own inputs), returning an empty list with the
# reason as an attribute when the scenario fails.
run_quiet <- function(model, st, joint = NULL) {
  tryCatch(
    withCallingHandlers({
      plan <- plan_method(model, st$method, st$event_model, given = FALSE)
      run_point(model, plan, risk = st$risk, joint = joint, warn = FALSE)
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
    cm_abort(sprintf("The baseline could not be adjusted (%s).", attr(res, "reason") %||% "failed"))
  }
  if (!result_is_finite(res)) {
    cm_abort(sprintf("The baseline result is undefined (non-finite adjusted impacts for %s); use another method.",
                     nonfinite_diseases(res)), class = "deconflate_nonfinite")
  }
  res
}

check_disease_ids <- function(ids, model) {
  unk <- setdiff(ids, model$diseases$id)
  if (length(unk)) cm_abort(sprintf("Unknown diseases: %s.", paste(unk, collapse = ", ")))
  invisible(ids)
}

#' Set or replace one association
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param disease1,disease2 Disease ids.
#' @param value Association value.
#' @param measure Measure (see [cm_associations()]).
#' @return The modified object. A distribution of the pair's association in
#'   the model's `distributions` is dropped (the scenario value is fixed).
#' @export
set_association <- function(model, disease1, disease2, value, measure = "OR") {
  check_population(model)
  check_disease_ids(c(disease1, disease2), model)
  a <- model$associations
  new <- cm_associations(disease1, disease2, value, measure = measure, source = "scenario")
  if (!is.null(model$distributions)) {
    k <- vapply(strsplit(names(model$distributions), ":", fixed = TRUE), function(p) {
      p[1] == "assoc" && pair_key(p[2], p[3]) == pair_key(disease1, disease2)
    }, logical(1))
    model$distributions <- if (any(!k)) model$distributions[!k]
  }
  if (is.null(a)) {
    model$associations <- new
    return(model)
  }
  hit <- match(pair_key(disease1, disease2), pair_key(a$disease1, a$disease2))
  if (is.na(hit)) {
    model$associations <- rbind(a, new)
  } else {
    a[hit, names(new)] <- new[1, ]
    model$associations <- a
  }
  model
}

#' Set or replace one pairwise interaction
#'
#' @param model A [cm_model()].
#' @param disease1,disease2 Disease ids.
#' @param value Interaction value (same units as the impacts).
#' @return The modified model. A distribution of the pair's interaction in
#'   the model's `distributions` is dropped (the scenario value is fixed).
#' @export
set_interaction <- function(model, disease1, disease2, value) {
  check_model(model)
  if (impact_kind(model$impacts) == "event") {
    cm_abort("Interactions apply to additive impacts only.", class = "deconflate_unsupported")
  }
  check_disease_ids(c(disease1, disease2), model)
  it <- model$interactions
  new <- cm_interactions(disease1, disease2, value, source = "scenario")
  if (!is.null(it)) {
    hit <- which(pair_key(it$disease1, it$disease2) == pair_key(disease1, disease2))
    if (length(hit)) it <- it[-hit, , drop = FALSE]
    new <- rbind(it, new)
  }
  model$interactions <- new
  if (!is.null(model$distributions)) {
    k <- vapply(strsplit(names(model$distributions), ":", fixed = TRUE), function(p) {
      p[1] == "inter" && pair_key(p[2], p[3]) == pair_key(disease1, disease2)
    }, logical(1))
    model$distributions <- if (any(!k)) model$distributions[!k]
  }
  model
}

#' Set or replace one three-way association term
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param disease1,disease2,disease3 Disease ids.
#' @param ratio Ratio of conditional odds ratios (see [cm_three_way()]).
#' @return The modified object. A distribution of the triple's ratio in the
#'   model's `distributions` is dropped (the scenario value is fixed).
#' @export
set_three_way <- function(model, disease1, disease2, disease3, ratio) {
  check_population(model)
  new <- cm_three_way(disease1, disease2, disease3, ratio, source = "scenario")
  check_disease_ids(c(disease1, disease2, disease3), model)
  tw <- model$three_way
  if (!is.null(tw) && nrow(tw)) {
    key <- apply(tw[, c("disease1", "disease2", "disease3")], 1, function(x) paste(sort(x), collapse = "|"))
    tw <- tw[key != paste(sort(c(disease1, disease2, disease3)), collapse = "|"), , drop = FALSE]
    new <- rbind(tw, new)
  }
  model$three_way <- new
  if (!is.null(model$distributions)) {
    target <- paste(sort(c(disease1, disease2, disease3)), collapse = "|")
    k <- vapply(strsplit(names(model$distributions), ":", fixed = TRUE), function(p) {
      p[1] == "three" && length(p) == 4L && paste(sort(p[2:4]), collapse = "|") == target
    }, logical(1))
    model$distributions <- if (any(!k)) model$distributions[!k]
  }
  model
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
#' adjusted as [deconflate()] would (point estimates; the method is chosen
#' for the scenario's own inputs, see `method`).
#'
#' The total is the adjusted aggregate for additive impacts, and the risk
#' attributable to disease for event impacts (`event_model = TRUE`); disease
#' rankings use the contributions to it.
#'
#' @param model A [cm_model()].
#' @param method As in [deconflate()].
#' @param event_model,overall_risk As in [deconflate()] (a distribution of
#'   the overall risk is used at its mean).
#' @param or_values Odds ratios tried for pairs without an association.
#' @param multipliers Factors applied to given associations.
#' @param pairs Optional character vector of pairs (`"d1:d2"`); default all.
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
screen_associations <- function(model, method = "auto", event_model = FALSE, overall_risk = NULL,
                                or_values = c(0.5, 2), multipliers = c(0.5, 2),
                                pairs = NULL) {
  st <- sens_setup(model, method, event_model, overall_risk, "screen_associations()")
  base <- burden_metric(check_baseline(run_quiet(model, st)))
  pt <- pair_tables(model)
  if (!is.null(pairs)) {
    check_pair_ids(pairs, model$diseases$id)
    ab <- strsplit(pairs, ":", fixed = TRUE)
    want <- vapply(ab, function(x) pair_key(x[1], x[2]), character(1))
    pt <- pt[pair_key(pt$disease1, pt$disease2) %in% want, , drop = FALSE]
  }
  key <- paste(pt$disease1, pt$disease2, sep = ":")
  a <- model$associations
  rows <- list()
  for (r in seq_len(nrow(pt))) {
    if (pt$status[r] == "specified") {
      hit <- match(pair_key(pt$disease1[r], pt$disease2[r]), pair_key(a$disease1, a$disease2))
      meas <- a$measure[hit]
      vals <- a$value[hit] * multipliers
      labels <- sprintf("%s x %g", meas, multipliers)
      d1 <- a$disease1[hit]
      d2 <- a$disease2[hit]
    } else {
      meas <- "OR"
      vals <- or_values
      labels <- sprintf("OR = %g", or_values)
      d1 <- pt$disease1[r]
      d2 <- pt$disease2[r]
    }
    for (v in seq_along(vals)) {
      m2 <- tryCatch(set_association(model, d1, d2, vals[v], meas),
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
                     paste(pairs[bad], collapse = ", ")))
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
#' Adds a pairwise interaction of each size in `values` to each disease pair,
#' one at a time, re-runs the global adjustment and reports the change in
#' the aggregate and rankings. The joint distribution is fitted once and
#' reused. Interactions apply to additive impacts only.
#'
#' @param model A [cm_model()].
#' @param values Interaction sizes (same units as the impacts; positive
#'   synergistic, negative antagonistic).
#' @param pairs Optional character vector of pairs (`"d1:d2"`).
#' @param ... Passed to [fit_joint()] (e.g. `backend = "sampled"` for many
#'   diseases).
#' @return A `cm_screen` data frame as in [screen_associations()].
#' @export
screen_interactions <- function(model, values, pairs = NULL, ...) {
  check_model(model)
  if (!is.null(model$impacts) && impact_kind(model$impacts) == "event") {
    cm_abort("Interactions apply to additive impacts only; screen_interactions() does not apply to event impacts.",
             class = "deconflate_unsupported")
  }
  st <- sens_setup(model, "global", FALSE, NULL, "screen_interactions()")
  # One joint distribution serves every scenario (interactions do not change
  # it); `...` goes to fit_joint(), e.g. backend = "sampled" (the default
  # beyond 20 diseases, as in deconflate()).
  fa <- list(...)
  if (nrow(model$diseases) > 20L && is.null(fa$backend)) fa$backend <- "sampled"
  joint <- do.call(fit_joint, c(list(as_population(model)), fa))
  base <- burden_metric(check_baseline(run_quiet(model, st, joint = joint)))
  ids <- model$diseases$id
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
                   scenario = sprintf("delta = %g", v), stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base, FALSE)
}

#' Screen three-way association scenarios
#'
#' Pairwise associations do not identify how often three diseases occur
#' together; the global model assumes no three-way association. This screen
#' sets a three-way term (see [cm_three_way()]) for each triple in turn,
#' refits the joint distribution (all pairwise associations are still
#' matched) and reports the change in the aggregate and rankings.
#'
#' When the pairs of a triple are all constrained, a three-way term keeps
#' their tables fixed, so it changes additive results only through
#' interactions: without interactions such a scenario reproduces the
#' baseline. When a pair is unknown (unconstrained), the fitted pairwise
#' table changes with the three-way term, and additive results can change
#' even without interactions. Three-way terms also matter for event impacts
#' (`event_model = TRUE`).
#'
#' @param model A [cm_model()].
#' @param event_model,overall_risk As in [deconflate()].
#' @param ratios Ratios of conditional odds ratios to try.
#' @param triples Optional list of character vectors of three disease ids;
#'   default all triples (which can be many).
#' @return A `cm_screen` data frame as in [screen_associations()].
#' @export
screen_three_way <- function(model, ratios = c(0.5, 2), triples = NULL, event_model = FALSE,
                             overall_risk = NULL) {
  st <- sens_setup(model, "global", event_model, overall_risk, "screen_three_way()")
  base <- burden_metric(check_baseline(run_quiet(model, st)))
  ids <- model$diseases$id
  if (is.null(triples)) {
    triples <- if (length(ids) >= 3L) utils::combn(ids, 3, simplify = FALSE) else list()
  } else {
    if (is.character(triples)) triples <- list(triples)
    for (tr in triples) {
      if (length(tr) != 3L || anyDuplicated(tr) || !all(tr %in% ids)) {
        cm_abort(sprintf("Each triple must name three different diseases of the model (check: %s).",
                         paste(tr, collapse = ", ")))
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
#' @param method,event_model,overall_risk As in [deconflate()].
#' @param variation Relative change applied downwards and upwards.
#' @param inputs Which inputs to vary: `"prob"`, `"assoc"`, `"impact"`, and
#'   for event impacts `"risk"` (the overall risk).
#' @return A `cm_oat` data frame with one row per input: the input key,
#'   totals at the low and high values and the swing, sorted by swing.
#' @export
#' @examples
#' sensitivity_oat(example_supplement())
sensitivity_oat <- function(model, method = "auto", event_model = FALSE, overall_risk = NULL,
                            variation = 0.2, inputs = c("prob", "assoc", "impact", "risk")) {
  st <- sens_setup(model, method, event_model, overall_risk, "sensitivity_oat()")
  inputs <- match.arg(inputs, several.ok = TRUE)
  base <- burden_metric(check_baseline(run_quiet(model, st)))
  vals <- model_inputs(model)
  if (event_model) vals <- c(vals, risk = st$risk)
  keys <- names(vals)
  type <- sub(":.*$", "", keys)
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
  v <- stats::setNames(model$diseases$value, paste0("prob:", model$diseases$id))
  a <- model$associations
  if (!is.null(a) && nrow(a)) {
    v <- c(v, stats::setNames(a$value, paste0("assoc:", a$disease1, ":", a$disease2)))
  }
  i <- model$impacts
  if (!is.null(i) && nrow(i)) v <- c(v, stats::setNames(i$value, paste0("impact:", i$disease)))
  x <- model$interactions
  if (!is.null(x) && nrow(x)) {
    v <- c(v, stats::setNames(x$value, paste0("inter:", x$disease1, ":", x$disease2)))
  }
  tw <- model$three_way
  if (!is.null(tw) && nrow(tw)) {
    v <- c(v, stats::setNames(tw$ratio, paste0("three:", tw$disease1, ":", tw$disease2, ":", tw$disease3)))
  }
  v
}
