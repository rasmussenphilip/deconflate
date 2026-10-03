# Sensitivity and scenario tools ------------------------------------------------

# Metric used by the screening tools: a total and a per-disease vector.
# With economics: monetary losses summed over the outcomes in
# economics$observed. Without: the expected loss L of one outcome.
burden_metric <- function(result, outcome = NULL, economics = NULL) {
  if (!is.null(economics)) {
    ev <- evaluate_economics(result, economics, allocate = TRUE)
    ids <- result$model$diseases$id
    by <- tapply(ev$by_disease$value, factor(ev$by_disease$disease, levels = ids), sum)
    by[is.na(by)] <- 0
    return(list(total = sum(ev$by_outcome$value), by_disease = stats::setNames(as.vector(by), ids)))
  }
  b <- attribute_burden(result)
  outcome <- outcome %||% b$outcome[1]
  b <- b[b$outcome == outcome, , drop = FALSE]
  if (!nrow(b)) cm_abort(sprintf("Outcome '%s' is not in the result.", outcome))
  list(total = sum(b$total), by_disease = stats::setNames(b$total, b$disease))
}

compare_to_base <- function(m, base) {
  if (is.null(m)) {
    return(data.frame(total = NA_real_, change = NA_real_, rel_change = NA_real_,
                      max_rank_shift = NA_integer_, rank_corr = NA_real_))
  }
  rb <- rank(-base$by_disease, ties.method = "min")
  rm <- rank(-m$by_disease, ties.method = "min")
  data.frame(total = m$total, change = m$total - base$total,
             rel_change = (m$total - base$total) / base$total,
             max_rank_shift = as.integer(max(abs(rb - rm))),
             rank_corr = suppressWarnings(stats::cor(base$by_disease, m$by_disease,
                                                     method = "spearman")))
}

run_quiet <- function(model, method, joint = NULL) {
  tryCatch(
    withCallingHandlers(deconflate(model, method = method, joint = joint, warn = FALSE),
                        deconflate_nonconvergence = function(w) invokeRestart("muffleWarning")),
    deconflate_infeasible = function(e) NULL
  )
}

#' Set or replace one association in a model
#'
#' @param model A [cm_model()].
#' @param disease1,disease2 Disease ids.
#' @param value Association value.
#' @param measure Measure (see [cm_associations()]).
#' @return The modified model.
#' @export
set_association <- function(model, disease1, disease2, value, measure = "OR") {
  check_model(model)
  a <- model$associations
  new <- cm_associations(disease1, disease2, value, measure = measure,
                         source = "scenario")
  if (is.null(a)) {
    model$associations <- new
    return(model)
  }
  hit <- match(pair_key(disease1, disease2), pair_key(a$disease1, a$disease2))
  if (is.na(hit)) {
    model$associations <- rbind(a, new)
  } else {
    a$disease1[hit] <- disease1
    a$disease2[hit] <- disease2
    a$value[hit] <- new$value
    a$measure[hit] <- measure
    a$source[hit] <- "scenario"
    model$associations <- a
  }
  model
}

#' Set or replace one pairwise interaction in a model
#'
#' @param model A [cm_model()].
#' @param disease1,disease2 Disease ids.
#' @param value Interaction value on the outcome's proportion scale.
#' @param outcome Outcome label.
#' @return The modified model.
#' @export
set_interaction <- function(model, disease1, disease2, value, outcome) {
  check_model(model)
  it <- model$interactions
  new <- cm_interactions(disease1, disease2, value, outcome = outcome, source = "scenario")
  if (!is.null(it)) {
    hit <- which(it$outcome == outcome &
                   pair_key(it$disease1, it$disease2) == pair_key(disease1, disease2))
    if (length(hit)) it <- it[-hit, , drop = FALSE]
    new <- rbind(it, new)
  }
  cm_model(model$diseases, model$associations, model$impacts, new,
           missing_associations = model$missing_associations)
}

#' Screen disease pairs for influential associations
#'
#' Re-runs the adjustment with each disease pair's association changed, one
#' pair at a time, and reports how much the total burden and the ranking of
#' diseases change. Pairs without a specified association (independent by
#' default, or unknown) are set to each of `or_values`; specified odds ratios
#' and risk ratios are multiplied by each of `multipliers`. This identifies
#' associations worth estimating, as in Rasmussen et al. (2022), Fig. 3 and
#' Rasmussen et al. (2024), Fig. 7C-D.
#'
#' @param model A [cm_model()].
#' @param method Adjustment method.
#' @param or_values Odds ratios tried for unspecified pairs.
#' @param multipliers Factors applied to specified odds ratios or risk ratios.
#' @param pairs Optional character vector of pairs (`"d1:d2"`) to screen;
#'   default all pairs.
#' @param outcome Outcome whose burden is the metric (ignored if `economics`
#'   is given); default the first outcome.
#' @param economics Optional list (`observed`, `unit_value`, `additional`) to
#'   use total monetary losses as the metric (see [value_losses()]).
#' @return A data frame sorted by the absolute relative change in the total,
#'   with the scenario, total, change, relative change, the largest shift in
#'   any disease's rank and the Spearman correlation of disease burdens with
#'   the baseline. Infeasible scenarios have `NA` totals.
#' @export
#' @examples
#' screen_associations(example_supplement(), or_values = c(0.5, 3))
screen_associations <- function(model, method = "simultaneous",
                                or_values = c(0.5, 2), multipliers = c(0.5, 2),
                                pairs = NULL, outcome = NULL, economics = NULL) {
  check_model(model)
  base_res <- deconflate(model, method = method, warn = FALSE)
  base <- burden_metric(base_res, outcome, economics)
  pt <- pair_tables(model)
  key <- paste(pt$disease1, pt$disease2, sep = ":")
  if (!is.null(pairs)) {
    ab <- strsplit(pairs, ":", fixed = TRUE)
    want <- vapply(ab, function(x) pair_key(x[1], x[2]), character(1))
    pt <- pt[pair_key(pt$disease1, pt$disease2) %in% want, , drop = FALSE]
    key <- paste(pt$disease1, pt$disease2, sep = ":")
  }
  rows <- list()
  for (r in seq_len(nrow(pt))) {
    specified <- pt$status[r] == "specified" &&
      !(pt$measure[r] %in% c("independent", "unknown"))
    if (specified) {
      if (!(pt$measure[r] %in% c("OR", "RR", "table"))) next
      a <- model$associations
      hit <- match(pair_key(pt$disease1[r], pt$disease2[r]), pair_key(a$disease1, a$disease2))
      meas <- if (a$measure[hit] == "table") "OR" else a$measure[hit]
      vals <- a$value[hit] * multipliers
      labels <- sprintf("%s x %g", meas, multipliers)
    } else {
      meas <- "OR"
      vals <- or_values
      labels <- sprintf("OR = %g", or_values)
    }
    for (v in seq_along(vals)) {
      m2 <- if (specified) {
        a <- model$associations
        set_association(model, a$disease1[hit], a$disease2[hit], vals[v], meas)
      } else {
        set_association(model, pt$disease1[r], pt$disease2[r], vals[v], meas)
      }
      res <- run_quiet(m2, method)
      met <- if (is.null(res)) NULL else burden_metric(res, outcome, economics)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(pair = key[r], status = pt$status[r], scenario = labels[v],
                   stringsAsFactors = FALSE),
        compare_to_base(met, base))
    }
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
#' total burden and rankings. The joint distribution is fitted once and
#' reused. Use it to see which assumed interactions would matter, before
#' looking for evidence on them.
#'
#' @param model A [cm_model()].
#' @param outcome Outcome to add interactions to.
#' @param values Interaction sizes (proportion scale; positive synergistic,
#'   negative antagonistic).
#' @param pairs Optional character vector of pairs (`"d1:d2"`).
#' @param economics Optional economics list; otherwise the burden of
#'   `outcome` is the metric.
#' @return A `cm_screen` data frame as in [screen_associations()].
#' @export
screen_interactions <- function(model, outcome, values = c(-0.01, 0.01),
                                pairs = NULL, economics = NULL) {
  check_model(model)
  joint <- fit_joint(model)
  base_res <- deconflate(model, method = "global", joint = joint, warn = FALSE)
  base <- burden_metric(base_res, outcome, economics)
  ids <- model$diseases$id
  pr <- if (is.null(pairs)) t(utils::combn(ids, 2)) else
    do.call(rbind, strsplit(pairs, ":", fixed = TRUE))
  rows <- list()
  for (r in seq_len(nrow(pr))) {
    for (v in values) {
      m2 <- set_interaction(model, pr[r, 1], pr[r, 2], v, outcome)
      res <- run_quiet(m2, "global", joint = joint)
      met <- if (is.null(res)) NULL else burden_metric(res, outcome, economics)
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(pair = paste(pr[r, 1], pr[r, 2], sep = ":"), status = "interaction",
                   scenario = sprintf("delta = %g", v), stringsAsFactors = FALSE),
        compare_to_base(met, base))
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(-abs(out$rel_change), na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "baseline_total") <- base$total
  class(out) <- c("cm_screen", "data.frame")
  out
}

#' One-at-a-time sensitivity analysis
#'
#' Varies each input by `variation` (e.g. +/- 20%) with everything else
#' fixed, as in Rasmussen et al. (2024), Fig. 7, and reports the resulting
#' range of the total burden.
#'
#' @param model A [cm_model()].
#' @param method Adjustment method.
#' @param variation Relative change applied downwards and upwards.
#' @param inputs Which inputs to vary: `"prob"`, `"assoc"`, `"impact"`.
#' @param outcome,economics Metric, as in [screen_associations()].
#' @return A `cm_screen`-like data frame with one row per input: the input
#'   key, totals at the low and high values and the swing, sorted by swing.
#' @export
#' @examples
#' sensitivity_oat(example_supplement())
sensitivity_oat <- function(model, method = "simultaneous", variation = 0.2,
                            inputs = c("prob", "assoc", "impact"),
                            outcome = NULL, economics = NULL) {
  check_model(model)
  inputs <- match.arg(inputs, several.ok = TRUE)
  base <- burden_metric(deconflate(model, method = method, warn = FALSE), outcome, economics)
  vals <- flatten_model(model)
  keys <- names(vals)
  type <- sub(":.*$", "", keys)
  keep <- type %in% inputs
  if ("assoc" %in% inputs && !is.null(model$associations)) {
    a <- model$associations
    num <- paste0("assoc:", a$disease1, ":", a$disease2)[a$measure %in% c("OR", "RR")]
    keep <- keep & (type != "assoc" | keys %in% num)
  }
  keep <- keep & unlist(vals[1, ]) != 0
  keys <- keys[keep]
  rows <- lapply(keys, function(k) {
    x <- vals[[k]]
    tot <- vapply(c(1 - variation, 1 + variation), function(f) {
      m2 <- tryCatch(set_input(model, k, x * f), deconflate_infeasible = function(e) NULL)
      if (is.null(m2)) return(NA_real_)
      res <- run_quiet(m2, method)
      if (is.null(res)) NA_real_ else burden_metric(res, outcome, economics)$total
    }, numeric(1))
    data.frame(input = k, value = x, total_low = tot[1], total_high = tot[2],
               swing = abs(tot[2] - tot[1]), rel_swing = abs(tot[2] - tot[1]) / base$total,
               stringsAsFactors = FALSE)
  })
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
    impact = cm_sampler(model, impacts = stats::setNames(d, paste(parts[2], parts[3], sep = ":"))),
    inter = cm_sampler(model, interactions = stats::setNames(d, paste(parts[2:4], collapse = ":"))),
    cm_abort(sprintf("Unknown input key '%s'.", key))
  )
  s(1)
}

#' Compare scenarios
#'
#' Adjusts several versions of a model (e.g. with different association or
#' interaction assumptions) and tabulates total burden and disease rankings
#' side by side.
#'
#' @param ... Named [cm_model()] objects.
#' @param method Adjustment method (`"global"` is used automatically for
#'   models with interactions).
#' @param outcome,economics Metric, as in [screen_associations()].
#' @return A list with `totals` (one row per scenario) and `by_disease`
#'   (burden and rank per disease and scenario).
#' @export
#' @examples
#' base <- example_supplement()
#' strong <- set_association(base, "d1", "d3", 3)
#' compare_scenarios(base = base, strong_d1_d3 = strong)
compare_scenarios <- function(..., method = "simultaneous", outcome = NULL,
                              economics = NULL) {
  models <- list(...)
  if (is.null(names(models)) || any(!nzchar(names(models)))) {
    cm_abort("Scenarios must be named, e.g. compare_scenarios(base = m1, alt = m2).")
  }
  mets <- lapply(models, function(m) {
    meth <- if (!is.null(m$interactions) && nrow(m$interactions)) "global" else method
    res <- run_quiet(m, meth)
    if (is.null(res)) NULL else burden_metric(res, outcome, economics)
  })
  totals <- data.frame(scenario = names(models),
                       total = vapply(mets, function(x) if (is.null(x)) NA_real_ else x$total, numeric(1)),
                       stringsAsFactors = FALSE)
  totals$rel_to_first <- totals$total / totals$total[1] - 1
  by <- do.call(rbind, lapply(names(mets), function(nm) {
    x <- mets[[nm]]
    if (is.null(x)) return(NULL)
    data.frame(scenario = nm, disease = names(x$by_disease), burden = unname(x$by_disease),
               rank = rank(-x$by_disease, ties.method = "min"), stringsAsFactors = FALSE)
  }))
  rownames(totals) <- NULL
  if (!is.null(by)) rownames(by) <- NULL
  list(totals = totals, by_disease = by)
}
