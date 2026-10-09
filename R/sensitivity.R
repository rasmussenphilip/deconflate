# Sensitivity and scenario tools ------------------------------------------------

# Metric used by the screening tools: the adjusted aggregate and the
# contributions (in the impacts' own units).
burden_metric <- function(result) {
  ids <- result$adjusted$disease
  list(total = result$totals$adjusted_total,
       by_disease = stats::setNames(result$contributions$total, ids))
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

# Adjust, returning NULL (with the reason as an attribute) when the scenario
# fails.
run_quiet <- function(model, method, joint = NULL) {
  tryCatch(
    withCallingHandlers(deconflate(model, method = method, joint = joint, warn = FALSE),
                        deconflate_nonconvergence = function(w) {
                          if (inherits(w, "warning")) invokeRestart("muffleWarning")
                        }),
    deconflate_error = function(e) structure(list(), reason = condition_type(e))
  )
}
failed_run <- function(res) is.null(res) || !inherits(res, "cm_result") || !result_is_finite(res)

# Metric of a scenario run, or NULL with the reason the scenario failed.
scenario_metric <- function(res) {
  if (failed_run(res)) {
    reason <- if (inherits(res, "cm_result")) {
      sprintf("undefined: non-finite adjusted impacts (%s)", nonfinite_diseases(res))
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
#' @return The modified object.
#' @export
set_association <- function(model, disease1, disease2, value, measure = "OR") {
  check_population(model)
  check_disease_ids(c(disease1, disease2), model)
  a <- model$associations
  new <- cm_associations(disease1, disease2, value, measure = measure, source = "scenario")
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
#' @return The modified model.
#' @export
set_interaction <- function(model, disease1, disease2, value) {
  check_model(model)
  check_disease_ids(c(disease1, disease2), model)
  it <- model$interactions
  new <- cm_interactions(disease1, disease2, value, source = "scenario")
  if (!is.null(it)) {
    hit <- which(pair_key(it$disease1, it$disease2) == pair_key(disease1, disease2))
    if (length(hit)) it <- it[-hit, , drop = FALSE]
    new <- rbind(it, new)
  }
  model$interactions <- new
  model
}

#' Set or replace one three-way association term
#'
#' @param model A [cm_population()] or [cm_model()].
#' @param disease1,disease2,disease3 Disease ids.
#' @param ratio Ratio of conditional odds ratios (see [cm_three_way()]).
#' @return The modified object.
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
  model
}

#' Screen disease pairs for influential associations
#'
#' Re-runs the adjustment with each disease pair's association changed, one
#' pair at a time, and reports how much the aggregate and the ranking of
#' diseases change. Pairs without a specified association (independent by
#' default, or unknown) are set to each of `or_values`; specified
#' associations of any measure are multiplied by each of `multipliers`
#' (contingency tables via their odds ratio). Scenarios whose inputs are
#' invalid or infeasible are kept, with `NA` results and the reason, so that
#' no measure is skipped silently. This identifies associations worth
#' estimating, as in Rasmussen et al. (2022), Fig. 3 and Rasmussen et al.
#' (2024), Fig. 7C-D.
#'
#' @param model A [cm_model()].
#' @param method Adjustment method.
#' @param or_values Odds ratios tried for unspecified pairs.
#' @param multipliers Factors applied to specified associations.
#' @param pairs Optional character vector of pairs (`"d1:d2"`); default all.
#' @return A `cm_screen` data frame sorted by the absolute relative change in
#'   the total, with the scenario, total, change, relative change, the
#'   largest shift in any disease's rank, the Spearman correlation of disease
#'   contributions with the baseline, and `failed` (the reason a scenario
#'   could not be run).
#' @export
#' @examples
#' screen_associations(example_supplement(), or_values = c(0.5, 3))
screen_associations <- function(model, method = "simultaneous",
                                or_values = c(0.5, 2), multipliers = c(0.5, 2),
                                pairs = NULL) {
  check_model(model)
  base_res <- deconflate(model, method = method, warn = FALSE)
  base <- burden_metric(check_baseline(base_res))
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
    specified <- pt$status[r] == "specified" && !(pt$measure[r] %in% c("independent", "unknown"))
    if (specified) {
      hit <- match(pair_key(pt$disease1[r], pt$disease2[r]), pair_key(a$disease1, a$disease2))
      meas <- if (a$measure[hit] == "table") "OR" else a$measure[hit]
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
      res <- if (inherits(m2, "cm_model")) run_quiet(m2, method) else m2
      rows[[length(rows) + 1L]] <- screen_row(
        data.frame(pair = key[r], status = pt$status[r], scenario = labels[v],
                   stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base)
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

finish_screen <- function(rows, base) {
  if (!length(rows)) {
    out <- data.frame(pair = character(0), status = character(0), scenario = character(0),
                      total = numeric(0), change = numeric(0), rel_change = numeric(0),
                      max_rank_shift = integer(0), rank_corr = numeric(0),
                      failed = character(0), stringsAsFactors = FALSE)
    attr(out, "baseline_total") <- base$total
    class(out) <- c("cm_screen", "data.frame")
    return(out)
  }
  out <- do.call(rbind, rows)
  out <- out[order(-abs(out$rel_change), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "baseline_total") <- base$total
  class(out) <- c("cm_screen", "data.frame")
  out
}

#' Screen disease pairs for influential impact interactions
#'
#' Adds a pairwise interaction of each size in `values` to each disease pair,
#' one at a time, re-runs the global adjustment and reports the change in
#' the aggregate and rankings. The joint distribution is fitted once and
#' reused.
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
  # One joint distribution serves every scenario (interactions do not change
  # it); `...` goes to fit_joint(), e.g. backend = "sampled".
  joint <- fit_joint(model, ...)
  base_res <- deconflate(model, method = "global", joint = joint, warn = FALSE)
  base <- burden_metric(check_baseline(base_res))
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
      res <- run_quiet(set_interaction(model, pr[r, 1], pr[r, 2], v), "global", joint = joint)
      rows[[length(rows) + 1L]] <- screen_row(
        data.frame(pair = paste(pr[r, 1], pr[r, 2], sep = ":"), status = "interaction",
                   scenario = sprintf("delta = %g", v), stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base)
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
#' even without interactions. Three-way terms also matter for hazard ratios
#' ([deconflate_hr()]) and [attributable_risk()].
#'
#' @param model A [cm_model()].
#' @param ratios Ratios of conditional odds ratios to try.
#' @param triples Optional list of character vectors of three disease ids;
#'   default all triples (which can be many).
#' @return A `cm_screen` data frame as in [screen_associations()].
#' @export
screen_three_way <- function(model, ratios = c(0.5, 2), triples = NULL) {
  check_model(model)
  base_res <- deconflate(model, method = "global", warn = FALSE)
  base <- burden_metric(check_baseline(base_res))
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
      res <- run_quiet(m2, "global")
      rows[[length(rows) + 1L]] <- screen_row(
        data.frame(pair = paste(tr, collapse = ":"), status = "three-way",
                   scenario = sprintf("ratio = %g", rt), stringsAsFactors = FALSE),
        res, base)
    }
  }
  finish_screen(rows, base)
}

#' One-at-a-time sensitivity analysis
#'
#' Varies each input by `variation` (e.g. +/- 20%) with everything else
#' fixed, as in Rasmussen et al. (2024), Fig. 7, and reports the resulting
#' range of the adjusted aggregate.
#'
#' @param model A [cm_model()].
#' @param method Adjustment method.
#' @param variation Relative change applied downwards and upwards.
#' @param inputs Which inputs to vary: `"prob"`, `"assoc"`, `"impact"`.
#' @return A `cm_oat` data frame with one row per input: the input key,
#'   totals at the low and high values and the swing, sorted by swing.
#' @export
#' @examples
#' sensitivity_oat(example_supplement())
sensitivity_oat <- function(model, method = "simultaneous", variation = 0.2,
                            inputs = c("prob", "assoc", "impact")) {
  check_model(model)
  inputs <- match.arg(inputs, several.ok = TRUE)
  base <- burden_metric(check_baseline(deconflate(model, method = method, warn = FALSE)))
  vals <- flatten_model(model)
  keys <- names(vals)
  type <- sub(":.*$", "", keys)
  keep <- type %in% inputs
  if ("assoc" %in% inputs && !is.null(model$associations)) {
    a <- model$associations
    num <- paste0("assoc:", a$disease1, ":", a$disease2)[!(a$measure %in% c("independent", "unknown", "table"))]
    keep <- keep & (type != "assoc" | keys %in% num)
  }
  keep <- keep & unlist(vals[1, ]) != 0
  keys <- keys[keep]
  rows <- lapply(keys, function(k) {
    x <- vals[[k]]
    tot <- vapply(c(1 - variation, 1 + variation), function(f) {
      m2 <- tryCatch(set_input(model, k, x * f), deconflate_error = function(e) NULL)
      if (is.null(m2)) return(NA_real_)
      sm <- scenario_metric(run_quiet(m2, method))
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
  class(out) <- c("cm_oat", "data.frame")
  out
}

# Set one input (key as in cm_sampler params) to a value, via a fixed-value sampler.
set_input <- function(model, key, value) {
  parts <- strsplit(key, ":", fixed = TRUE)[[1]]
  d <- list(dist_fixed(value))
  s <- switch(parts[1],
    prob = cm_sampler(model, diseases = stats::setNames(d, parts[2])),
    assoc = cm_sampler(model, associations = stats::setNames(d, paste(parts[2], parts[3], sep = ":"))),
    impact = cm_sampler(model, impacts = stats::setNames(d, parts[2])),
    inter = cm_sampler(model, interactions = stats::setNames(d, paste(parts[2], parts[3], sep = ":"))),
    three = cm_sampler(model, three_way = stats::setNames(d, paste(parts[2], parts[3], parts[4], sep = ":"))),
    cm_abort(sprintf("Unknown input key '%s'.", key))
  )
  # Pass the value directly, so that no random number is drawn.
  s(1, values = stats::setNames(value, names(attr(s, "specs"))))
}

#' Compare scenarios
#'
#' Adjusts several versions of a model (e.g. with different association,
#' three-way or interaction assumptions) and tabulates the aggregate and
#' disease rankings side by side.
#'
#' @param ... Named [cm_model()] objects.
#' @param method Adjustment method (`"global"` is used automatically for
#'   models with interactions or three-way terms).
#' @return A list with `totals` (one row per scenario, with the reason in
#'   `failed` when a scenario could not be adjusted or gave an undefined
#'   result) and `by_disease` (contribution and rank per disease and
#'   scenario).
#' @export
#' @examples
#' base <- example_supplement()
#' strong <- set_association(base, "d1", "d3", 3)
#' compare_scenarios(base = base, strong_d1_d3 = strong)
compare_scenarios <- function(..., method = "simultaneous") {
  method <- public_method(method, c("simultaneous", "global"), "compare_scenarios()")
  models <- list(...)
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    cm_abort("Scenarios must be named, e.g. compare_scenarios(base = m1, alt = m2).")
  }
  not_model <- names(models)[!vapply(models, inherits, logical(1), "cm_model")]
  if (length(not_model)) {
    cm_abort(sprintf("Scenarios must be cm_model() objects (not: %s).", paste(not_model, collapse = ", ")))
  }
  sms <- lapply(models, function(m) {
    needs_global <- (!is.null(m$interactions) && nrow(m$interactions)) ||
      (!is.null(m$three_way) && nrow(m$three_way))
    scenario_metric(run_quiet(m, if (needs_global) "global" else method))
  })
  mets <- lapply(sms, `[[`, "met")
  totals <- data.frame(scenario = names(models),
                       total = vapply(mets, function(x) if (is.null(x)) NA_real_ else x$total, numeric(1)),
                       failed = vapply(sms, function(x) x$reason, character(1)),
                       stringsAsFactors = FALSE)
  totals$rel_to_first <- if (isTRUE(abs(totals$total[1]) > 0)) totals$total / totals$total[1] - 1 else NA_real_
  by <- do.call(rbind, lapply(names(mets), function(nm) {
    x <- mets[[nm]]
    if (is.null(x)) return(NULL)
    data.frame(scenario = nm, disease = names(x$by_disease), contribution = unname(x$by_disease),
               rank = rank(-x$by_disease, ties.method = "min"), stringsAsFactors = FALSE)
  }))
  rownames(totals) <- NULL
  if (!is.null(by)) rownames(by) <- NULL
  list(totals = totals, by_disease = by)
}
