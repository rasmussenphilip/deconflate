# Sensitivity of the results to the inputs ------------------------------------
#
# cm_sensitivity(): design note section 17; reference values:
# inst/validation/reference_v06_sensitivity.py.

# Metric used by the screening tools: the total and each disease's
# contribution. Additive impacts: the adjusted aggregate and the
# contributions (in the impacts' own units). Event impacts: the risk
# attributable to disease and its Shapley allocation.
burden_metric <- function(result) {
  if (inherits(result, "cm_event_result")) {
    ids <- result$adjusted$disease[result$adjusted$kind %||% rep("disease", nrow(result$adjusted)) == "disease"]
    by <- result$contributions
    return(list(total = result$attributable$summary$attributable,
                by_disease = stats::setNames(if (is.null(by)) rep(NA_real_, length(ids)) else by$contribution, ids)))
  }
  list(total = result$totals$adjusted_total,
       by_disease = stats::setNames(result$contributions$contribution, result$contributions$disease))
}

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
#' @noRd
set_association <- function(model, disease1, disease2, value, measure = "OR") {
  check_population(model)
  check_disease_ids(c(disease1, disease2), model)
  new <- cm_associations(disease1, disease2, value, measure = measure, note = "scenario")
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
#' @noRd
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
    cm_impacts(combo, value, adjusted_for = "all", measure = "HR", note = "scenario")
  } else {
    cm_impacts(combo, value, adjusted_for = "all", note = "scenario")
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
#' @noRd
set_three_way <- function(model, disease1, disease2, disease3, ratio) {
  check_population(model)
  check_disease_ids(c(disease1, disease2, disease3), model)
  if (anyDuplicated(c(disease1, disease2, disease3))) cm_abort("A three-way association needs three different diseases.", class = "deconflate_invalid_input")
  a <- model$associations
  hit <- if (is.null(a)) logical(0) else a$measure == "OR_ratio" & vapply(seq_len(nrow(a)), function(r) {
    setequal(c(a$disease[r], a$associated_with[r], sub("^no\\s+", "", a$among[r])), c(disease1, disease2, disease3))
  }, logical(1))
  new <- cm_associations(disease1, disease2, ratio, measure = "OR_ratio", among = disease3, note = "scenario")
  replace_rows(model, "associations", hit, new)
}


# ---- Shared engine: rerunning a central result -------------------------------
#
# cm_sensitivity(), cm_threshold() and the diagnosis rerun the model of a
# central result with one input changed. `ctx` holds what that needs: the
# model, the plan, the risk (value, distribution, kind), the arguments of
# the joint fit and the central result.

# The context of a deconflate() result.
result_context <- function(res) {
  s <- res$settings
  if (is.null(s)) {
    cm_abort("This result does not carry the settings of deconflate(); run deconflate() again with this version of the package.",
             class = "deconflate_invalid_input")
  }
  plan <- s$plan
  dots <- s$dots %||% list()
  # With the sampled backend, every rerun uses the same seed, so that the
  # differences between runs are not Monte Carlo noise.
  if (identical(plan$fit_args$backend %||% dots$backend, "sampled") && is.null(plan$fit_args$seed) &&
      is.null(dots$seed)) {
    plan$fit_args$seed <- 1L
  }
  list(model = res$model, plan = plan, risk = s$risk, dots = dots, res = res,
       event = inherits(res, "cm_event_result"))
}

# Adjust the context's model again, as the central result was: with some
# inputs changed (`values`, by canonical key, "risk" for the overall or
# baseline risk), or a changed model (`model`, e.g. with an added row; it is
# planned again). The pooling of association rows and the impact-row SDs
# stay at the central estimate (`frozen`; NULL recomputes them, the same
# for rows that did not change). `same_joint`: the change leaves the joint
# distribution as it is (impact rows only), so the central one is reused.
# Returns the result, or the condition if the adjustment fails.
rerun_result <- function(ctx, values = list(), model = NULL, allocate = TRUE, frozen = "central",
                         same_joint = NULL) {
  res <- ctx$res
  plan <- ctx$plan
  if (identical(frozen, "central")) {
    frozen <- if (is.null(model)) list(associations = res$pooling$frozen, impacts = res$row_se)
  }
  tryCatch(withCallingHandlers({
    m <- if (is.null(model)) {
      set_inputs(ctx$model, values[names(values) != "risk"])
    } else {
      plan <- plan_model(model, ctx$event, ctx$dots, pool = ctx$plan$pool %||% "random",
                         risk_kind = ctx$risk$kind %||% "overall")
      if (!is.null(ctx$plan$fit_args$seed)) plan$fit_args$seed <- ctx$plan$fit_args$seed
      model
    }
    rk <- if ("risk" %in% names(values)) values[["risk"]] else ctx$risk$value
    if (!is.null(rk) && !(is.finite(rk) && rk > 0 && rk < 1)) {
      cm_abort(sprintf("The %s must lie between 0 and 1.", risk_label(ctx$risk$kind)), class = "deconflate_infeasible")
    }
    same_joint <- same_joint %||% (is.null(model) && length(values) > 0L &&
                                     all(startsWith(names(values), "impacts:") | names(values) == "risk"))
    jt <- if (same_joint && !isTRUE(plan$shortcut)) res$joint
    joint_start <- if (!is.null(res$joint) && identical(res$joint$backend, "exact")) res$joint
    beta_start <- if (ctx$event && is.null(model)) log(res$adjusted$adjusted)
    run_point(m, plan, risk = rk, joint = jt, warn = FALSE, dots = ctx$dots, joint_start = joint_start,
              start = beta_start, frozen = frozen, allocate = allocate, raw = FALSE)
  }, deconflate_nonconvergence = function(w) {
    if (inherits(w, "warning")) invokeRestart("muffleWarning")
  }), deconflate_error = function(e) e)
}

# The value of an input with distribution `d` at `shift` SDs from `x0`, on the
# normal score u = qnorm(F(x)) (clamped to [-5, 5]).
shifted_value <- function(d, x0, shift) {
  tryCatch({
    u0 <- stats::qnorm(d$p(x0))
    u0 <- if (is.finite(u0)) min(max(u0, -5), 5) else if (is.na(u0)) 0 else sign(u0) * 5
    d$q(stats::pnorm(u0 + shift))
  }, error = function(e) NA_real_)
}

# The inputs of a model that have a distribution, by canonical key (with
# "risk" for the overall or baseline risk); fixed distributions are left out.
input_dists <- function(model, risk_dist = NULL) {
  dists <- model$distributions %||% list()
  if (length(dists)) {
    dists <- dists[!startsWith(names(dists), "between:")]
    keys <- vapply(names(dists), function(k) {
      if (grepl("^(prob|associations|impacts):", k)) k else
        tryCatch(canonical_key(k, model), deconflate_error = function(e) NA_character_)
    }, character(1), USE.NAMES = FALSE)
    names(dists) <- keys
    dists <- dists[!is.na(keys)]
  }
  if (!is.null(risk_dist)) dists[["risk"]] <- risk_dist
  dists[!vapply(dists, function(d) identical(d$type, "fixed"), logical(1))]
}

# ---- cm_sensitivity() -------------------------------------------------------

sens_unknown_pair <- c(0.5, 3)
sens_unknown_or_ratio <- c(0.5, 2)
sens_unknown_hr <- c(0.67, 1.5)
sens_unknown_share <- 0.5

#' Sensitivity of the results to entered and unknown inputs
#'
#' Reruns the adjustment with one input at a time at a low and a high value,
#' everything else at its central value, and reports how much the total
#' (the adjusted aggregate; for event impacts, the risk attributable to
#' disease), each disease's adjusted impact and the ranking of diseases
#' change. Inputs you entered and inputs you did not enter (unknown
#' associations, interactions, three-way associations) are varied in the
#' same table, so you can see which missing estimates would matter most
#' next to the uncertainty of the ones you have.
#'
#' @section Inputs and their ranges:
#' Entered inputs (`inputs = "entered"`), named as in the input tables:
#' `prob:<disease>`, `associations:<row>`, `impacts:<row>`, and for event
#' impacts `risk` (the overall or baseline risk, whichever was given).
#' * With a distribution: one SD down and up, on the normal score of the
#'   distribution (the 16% and 84% quantiles when the value is the median).
#' * Without one: the value times `1 - change` and `1 + change`
#'   (probabilities stay inside (0, 1), risk differences and phi inside
#'   (-1, 1)). A value of zero has no relative change: give it a range with
#'   `ranges`.
#'
#' Unknown inputs (`inputs = "unknown"`), which are not in the tables:
#' * `pair:<d1>:<d2>`: a pair that no association row mentions, as a crude
#'   odds ratio from 0.5 to 3 (a pair mentioned only by other rows, e.g.
#'   within a subgroup or in an `OR_ratio` row, is not varied);
#' * `interaction:<d1>;<d2>`: an interaction of two diseases that have no
#'   combination row, as a row for the combination adjusted for all (see
#'   [cm_impacts()]): for additive impacts from minus to plus half the
#'   smaller of the two adjusted impacts (in absolute value), for event
#'   impacts a hazard ratio from 0.67 to 1.5;
#' * `or_ratio:<d1>:<d2>|<d3>`: only on request (`three_way`), how much
#'   disease 3 changes the association of diseases 1 and 2 (an `OR_ratio`
#'   row), from 0.5 to 2. Triples that already have an `OR_ratio` row are
#'   entered inputs (`associations:<row>`).
#'
#' Any range can be changed with `ranges`. With crude impact rows, a
#' positive interaction lowers the aggregate: the rows already contain it,
#' once in each disease's row, so the model takes it out twice and adds it
#' back once.
#'
#' @section Variance shares:
#' When the result has draws (`n_draws` in [deconflate()]), each uncertain
#' input's share of the variance of the total is estimated from them: the
#' total is regressed on the normal scores of all drawn inputs, and the share
#' of an input is its squared standardised coefficient. `r_squared` says how
#' much of the variance the linear fit explains; below 0.8 the shares are
#' rough (the response is far from linear, or inputs interact). Only inputs
#' with a distribution have a share.
#'
#' @param x A [deconflate()] result, or a [cm_model()]: a model is adjusted
#'   first (point estimates; it may have no association rows, to see which
#'   associations would matter).
#' @param inputs `"entered"`, `"unknown"` or both (default).
#' @param change Relative change of entered inputs without a distribution
#'   (default 0.2: 20% down and up).
#' @param ranges Optional data frame with columns `input` (a key as above),
#'   `low` and `high`, replacing the default range of those inputs. A key
#'   `or_ratio:<d1>:<d2>|<d3>` adds that three-way input.
#' @param three_way Three-way inputs to add: a list of character vectors of
#'   three disease ids (the association of the first two, among the third),
#'   or `"all"` (every triple, which can be many). Default none.
#' @param keys Optional input keys: vary only these (to save time).
#' @param event_model,overall_risk,baseline_risk With a model: as in
#'   [deconflate()]. With a result they are taken from it.
#' @param ... With a model: `pool`, `joint`, or arguments of [fit_joint()].
#' @return A `cm_sensitivity` object:
#' * `inputs`: one row per input, sorted by swing: `input` (key), `label`,
#'   `kind` (`"entered"` or `"unknown"`), `range` (how the range was set),
#'   `value` (entered inputs), `low`, `high`, `total_low`, `total_high`,
#'   `swing` (`|total_high - total_low|`), `rank_shift` (the most places any
#'   disease moves in the ranking by contribution, at low or high),
#'   `first` (the disease ranked first at low or high, when it is not the
#'   central one), `share` (of the variance over the draws, a fraction; `NA`
#'   without draws) and `failed` (why a run could not be used);
#' * `by_disease`: each input's effect on each disease (`adjusted_low`,
#'   `adjusted_high`; hazard ratios for event impacts; `contribution_low`,
#'   `contribution_high`);
#' * `variance`: the variance shares of all drawn inputs (with
#'   `r_squared`), or `NULL` without draws;
#' * `central`: the central total, contributions and ranking.
#' @seealso [cm_threshold()] for the value of one input at which a
#'   conclusion changes; `$ranks` of a [deconflate()] result for how often
#'   each disease ranks first over the draws.
#' @export
#' @examples
#' res <- deconflate(example_supplement(), n_draws = 0)
#' s <- cm_sensitivity(res)
#' s
#' print(s, n = Inf)
#' # An unknown association from 1 to 5 instead of 0.5 to 3
#' m <- example_supplement()
#' m$associations <- m$associations[m$associations$associated_with != "d3" |
#'                                    m$associations$disease != "d1", ]
#' cm_sensitivity(deconflate(m, n_draws = 0),
#'                ranges = data.frame(input = "pair:d1:d3", low = 1, high = 5))
cm_sensitivity <- function(x, inputs = c("entered", "unknown"), change = 0.2, ranges = NULL,
                           three_way = NULL, keys = NULL, event_model = FALSE, overall_risk = NULL,
                           baseline_risk = NULL, ...) {
  dots <- list(...)
  retired_args(dots)
  inputs <- cm_match_arg(inputs, several_ok = TRUE)
  if (!is.numeric(change) || length(change) != 1L || !is.finite(change) || change <= 0 || change >= 1) {
    cm_abort("`change` must be a single number between 0 and 1 (e.g. 0.2 for 20%).", class = "deconflate_invalid_input")
  }
  if (inherits(x, "cm_result") || inherits(x, "cm_event_result")) {
    if (!isFALSE(event_model) || !is.null(overall_risk) || !is.null(baseline_risk) || length(dots)) {
      cm_abort("With a deconflate() result, the settings come from the result: leave out event_model, overall_risk, baseline_risk and other arguments.",
               class = "deconflate_invalid_input")
    }
    res <- x
  } else {
    res <- sens_central(x, event_model, overall_risk, baseline_risk, dots)
  }
  ctx <- result_context(res)
  model <- ctx$model
  ids <- model$diseases$disease
  base <- burden_metric(res)
  if (!result_is_finite(res) || !all(is.finite(base$by_disease))) {
    cm_abort("The central result is not finite, so there is nothing to compare with.", class = "deconflate_nonfinite")
  }
  rank0 <- rank(-base$by_disease, ties.method = "first")
  first0 <- names(base$by_disease)[which.min(rank0)]

  # The inputs and their ranges.
  todo <- list()
  if ("entered" %in% inputs) todo <- sens_entered(ctx, change)
  if ("unknown" %in% inputs) todo <- c(todo, sens_unknowns(ctx, res, three_way))
  todo <- sens_apply_ranges(todo, ranges, ctx)
  if (!is.null(keys)) {
    if (!is.character(keys)) cm_abort("`keys` must be a character vector of input keys.", class = "deconflate_invalid_input")
    have <- vapply(todo, `[[`, character(1), "key")
    unk <- setdiff(keys, have)
    if (length(unk)) {
      cm_abort(sprintf("Unknown input keys: %s (see ?cm_sensitivity for the keys).", paste(unk, collapse = ", ")),
               class = "deconflate_invalid_input")
    }
    todo <- todo[have %in% keys]
  }

  # One input at a time.
  term_ids <- res$adjusted$disease[res$adjusted$kind == "disease"]
  rows <- list()
  long <- list()
  for (inp in todo) {
    outs <- lapply(c(inp$low, inp$high), function(v) sens_run(ctx, inp, v))
    tot <- vapply(outs, function(o) if (inherits(o, "condition")) NA_real_ else burden_metric(o)$total, numeric(1))
    bad <- vapply(outs, function(o) inherits(o, "condition") || !result_is_finite(o), logical(1))
    tot[bad] <- NA_real_
    shifts <- vapply(seq_along(outs), function(j) {
      if (bad[j]) return(NA_real_)
      bd <- burden_metric(outs[[j]])$by_disease[names(base$by_disease)]
      max(abs(rank(-bd, ties.method = "first") - rank0))
    }, numeric(1))
    firsts <- vapply(seq_along(outs), function(j) {
      if (bad[j]) return(NA_character_)
      bd <- burden_metric(outs[[j]])$by_disease[names(base$by_disease)]
      names(bd)[which.max(bd)]
    }, character(1))
    changed <- !is.na(firsts) & firsts != first0
    first <- if (any(changed)) paste(sprintf("%s (%s)", firsts[changed], c("low", "high")[changed]), collapse = ", ") else ""
    reason <- vapply(seq_along(outs), function(j) {
      o <- outs[[j]]
      if (!bad[j]) return(NA_character_)
      sprintf("%s: %s", c("low", "high")[j], if (inherits(o, "condition")) condition_type(o) else "non-finite result")
    }, character(1))
    rows[[length(rows) + 1L]] <- data.frame(
      input = inp$key, label = inp$label, kind = inp$kind, range = inp$how,
      value = inp$value %||% NA_real_, low = inp$low, high = inp$high,
      total_low = tot[1], total_high = tot[2], swing = abs(tot[2] - tot[1]),
      rank_shift = if (all(is.na(shifts))) NA_integer_ else as.integer(max(shifts, na.rm = TRUE)),
      first = first, share = NA_real_,
      failed = if (all(is.na(reason))) NA_character_ else paste(reason[!is.na(reason)], collapse = "; "),
      stringsAsFactors = FALSE)
    adj_of <- function(j) {
      if (bad[j]) return(stats::setNames(rep(NA_real_, length(term_ids)), term_ids))
      a <- outs[[j]]$adjusted
      stats::setNames(a$adjusted[match(term_ids, a$disease)], term_ids)
    }
    con_of <- function(j) {
      if (bad[j]) return(stats::setNames(rep(NA_real_, length(ids)), ids))
      burden_metric(outs[[j]])$by_disease[ids]
    }
    long[[length(long) + 1L]] <- data.frame(
      input = inp$key, disease = term_ids, adjusted_low = unname(adj_of(1)), adjusted_high = unname(adj_of(2)),
      contribution_low = unname(con_of(1)[term_ids]), contribution_high = unname(con_of(2)[term_ids]),
      stringsAsFactors = FALSE)
  }
  tab <- if (length(rows)) do.call(rbind, rows) else
    data.frame(input = character(0), label = character(0), kind = character(0), range = character(0),
               value = numeric(0), low = numeric(0), high = numeric(0), total_low = numeric(0),
               total_high = numeric(0), swing = numeric(0), rank_shift = integer(0), first = character(0),
               share = numeric(0), failed = character(0), stringsAsFactors = FALSE)
  variance <- sens_variance(res, model)
  if (!is.null(variance)) tab$share <- variance$share[match(tab$input, variance$input)]
  tab <- tab[order(-tab$swing, na.last = TRUE), , drop = FALSE]
  rownames(tab) <- NULL
  by_disease <- if (length(long)) do.call(rbind, long) else NULL
  if (!is.null(by_disease)) {
    by_disease <- by_disease[order(match(by_disease$input, tab$input)), , drop = FALSE]
    rownames(by_disease) <- NULL
  }
  structure(list(inputs = tab, by_disease = by_disease, variance = variance,
                 r_squared = attr(variance, "r_squared"),
                 central = list(total = base$total, contributions = base$by_disease, first = first0,
                                adjusted = stats::setNames(res$adjusted$adjusted, res$adjusted$disease)),
                 metric = if (ctx$event) "attributable risk" else "aggregate",
                 units = res$units, label = res$label, event = ctx$event,
                 settings = list(change = change, inputs = inputs)),
            class = "cm_sensitivity")
}

# The central result for cm_sensitivity(model): as deconflate() (point
# estimates), but a model without association rows is allowed (the diseases
# are then unrelated).
sens_central <- function(model, event_model, overall_risk, baseline_risk, dots) {
  check_model(model)
  check_event_model(model, event_model)
  risk <- check_overall_risk(overall_risk, event_model, baseline_risk)
  pool <- dots$pool %||% "random"
  if (!(identical(pool, "random") || identical(pool, "fixed"))) {
    cm_abort("`pool` must be \"random\" or \"fixed\".", class = "deconflate_invalid_input")
  }
  joint <- dots$joint
  dots$pool <- NULL
  dots$joint <- NULL
  unknown <- setdiff(names(dots), joint_arg_names)
  if (length(unknown)) {
    cm_abort(sprintf("Unknown argument%s: %s.", if (length(unknown) > 1L) "s" else "", paste(unknown, collapse = ", ")),
             class = "deconflate_invalid_input")
  }
  plan <- plan_model(model, event_model, dots, has_joint = !is.null(joint), pool = pool,
                     risk_kind = risk$kind %||% "overall")
  res <- run_point(model, plan, risk = risk$value, joint = joint, warn = FALSE, dots = dots, central = TRUE)
  res <- finish_result(res, model, plan)
  res$settings <- list(plan = plan, risk = risk, dots = dots)
  res
}

# Entered inputs with their ranges.
sens_entered <- function(ctx, change) {
  model <- ctx$model
  dists <- input_dists(model, ctx$risk$dist)
  keys <- c(paste0("prob:", model$diseases$disease),
            if (!is.null(model$associations) && nrow(model$associations)) paste0("associations:", model$associations$row),
            paste0("impacts:", model$impacts$row),
            if (ctx$event) "risk")
  labels <- input_labels(model, keys, ctx$risk$kind %||% "overall")
  out <- vector("list", length(keys))
  for (i in seq_along(keys)) {
    k <- keys[i]
    v <- input_value(model, k, ctx$risk$value)
    d <- dists[[k]]
    if (!is.null(d)) {
      lo <- shifted_value(d, v, -1)
      hi <- shifted_value(d, v, 1)
      how <- "1 SD"
    } else {
      lo <- v * (1 - change)
      hi <- v * (1 + change)
      meas <- input_measure(model, k)
      if (k == "risk" || meas %in% c("prevalence", "probability")) {
        lo <- max(lo, 1e-6)
        hi <- min(hi, 1 - 1e-6)
      }
      if (meas %in% c("RD", "phi")) {
        lo <- max(lo, -1 + 1e-6)
        hi <- min(hi, 1 - 1e-6)
      }
      how <- sprintf("%g%%", 100 * change)
      # A value of zero has no relative change: it needs a range.
      if (v == 0) {
        lo <- hi <- NA_real_
        how <- "value 0: give a range"
      }
    }
    out[[i]] <- list(key = k, label = labels[i], kind = "entered", how = how, value = v,
                     low = if (anyNA(c(lo, hi))) NA_real_ else min(lo, hi),
                     high = if (anyNA(c(lo, hi))) NA_real_ else max(lo, hi))
  }
  out
}

# Unknown inputs: pairs that no row mentions, interactions of pairs without
# a combination term, and the requested three-way associations.
sens_unknowns <- function(ctx, res, three_way = NULL) {
  model <- ctx$model
  ids <- model$diseases$disease
  out <- list()
  pt <- pair_tables(as_population(model))
  for (r in which(pt$status == "unknown")) {
    a <- pt$disease1[r]
    b <- pt$disease2[r]
    out[[length(out) + 1L]] <- list(key = sprintf("pair:%s:%s", a, b), label = sprintf("%s-%s OR (unknown)", a, b),
                                    kind = "unknown", how = sprintf("OR %g-%g", sens_unknown_pair[1], sens_unknown_pair[2]),
                                    value = NA_real_, low = sens_unknown_pair[1], high = sens_unknown_pair[2])
  }
  terms <- res$adjusted$disease
  adj <- stats::setNames(res$adjusted$adjusted, res$adjusted$disease)
  if (length(ids) >= 2L) {
    pr <- utils::combn(ids, 2)
    for (j in seq_len(ncol(pr))) {
      a <- pr[1, j]
      b <- pr[2, j]
      if (term_label(c(a, b)) %in% terms) next
      if (ctx$event) {
        lo <- sens_unknown_hr[1]
        hi <- sens_unknown_hr[2]
        how <- sprintf("HR %g-%g", lo, hi)
      } else {
        half <- sens_unknown_share * min(abs(adj[[a]]), abs(adj[[b]]))
        lo <- -half
        hi <- half
        how <- sprintf("%g%% of smaller impact", 100 * sens_unknown_share)
      }
      out[[length(out) + 1L]] <- list(key = sprintf("interaction:%s;%s", a, b),
                                      label = sprintf("%s;%s interaction (unknown)", a, b),
                                      kind = "unknown", how = how, value = NA_real_, low = lo, high = hi)
    }
  }
  if (!is.null(three_way)) {
    if (identical(three_way, "all")) {
      three_way <- if (length(ids) >= 3L) utils::combn(ids, 3, simplify = FALSE) else list()
    }
    if (is.character(three_way)) three_way <- list(three_way)
    for (tr in three_way) {
      if (!is.character(tr) || length(tr) != 3L || anyDuplicated(tr) || !all(tr %in% ids)) {
        cm_abort(sprintf("Each element of `three_way` must name three different diseases of the model (check: %s).",
                         paste(tr, collapse = ", ")), class = "deconflate_invalid_input")
      }
      if (!or_ratio_entered(model, tr)) out[[length(out) + 1L]] <- sens_or_ratio_input(tr)
    }
  }
  out
}

# Is there an OR_ratio row about these three diseases (in any order)? It is
# then an entered input (associations:<row>), not an unknown one.
or_ratio_entered <- function(model, tr) {
  a <- model$associations
  if (is.null(a) || !nrow(a)) return(FALSE)
  any(a$measure == "OR_ratio" & vapply(seq_len(nrow(a)), function(r) {
    setequal(c(a$disease[r], a$associated_with[r], sub("^no\\s+", "", a$among[r])), tr)
  }, logical(1)))
}

sens_or_ratio_input <- function(tr) {
  list(key = sprintf("or_ratio:%s:%s|%s", tr[1], tr[2], tr[3]),
       label = sprintf("%s-%s OR ratio among %s (unknown)", tr[1], tr[2], tr[3]),
       kind = "unknown", how = sprintf("ratio %g-%g", sens_unknown_or_ratio[1], sens_unknown_or_ratio[2]),
       value = NA_real_, low = sens_unknown_or_ratio[1], high = sens_unknown_or_ratio[2])
}

# Ranges given by the user replace the defaults (an or_ratio key adds that
# input).
sens_apply_ranges <- function(todo, ranges, ctx) {
  if (is.null(ranges)) return(todo)
  if (!is.data.frame(ranges) || !all(c("input", "low", "high") %in% names(ranges))) {
    cm_abort("`ranges` must be a data frame with columns input, low and high.", class = "deconflate_invalid_input")
  }
  have <- vapply(todo, `[[`, character(1), "key")
  ids <- ctx$model$diseases$disease
  for (i in seq_len(nrow(ranges))) {
    k <- as.character(ranges$input[i])
    lo <- ranges$low[i]
    hi <- ranges$high[i]
    if (!is.numeric(lo) || !is.numeric(hi) || !is.finite(lo) || !is.finite(hi) || lo >= hi) {
      cm_abort(sprintf("The range of '%s' must be two finite numbers, low below high.", k), class = "deconflate_invalid_input")
    }
    j <- match(k, have)
    if (is.na(j) && grepl("^(assoc|impact):", k)) {
      # Keys of the distributions' older form name an entered row.
      ck <- tryCatch(canonical_key(k, ctx$model), deconflate_error = function(e) NA_character_)
      j <- match(ck, have)
    }
    if (is.na(j) && startsWith(k, "or_ratio:")) {
      tr <- strsplit(sub("^or_ratio:", "", k), "[:|]")[[1]]
      if (length(tr) != 3L || anyDuplicated(tr) || !all(tr %in% ids)) {
        cm_abort(sprintf("'%s' must be 'or_ratio:<d1>:<d2>|<d3>' with three diseases of the model.", k),
                 class = "deconflate_invalid_input")
      }
      if (or_ratio_entered(ctx$model, tr)) {
        cm_abort(sprintf("The model has an OR_ratio row about %s: vary it as 'associations:<row>'.",
                         paste(tr, collapse = ", ")), class = "deconflate_invalid_input")
      }
      same <- vapply(todo, function(t) startsWith(t$key, "or_ratio:") &&
                       setequal(strsplit(sub("^or_ratio:", "", t$key), "[:|]")[[1]], tr), logical(1))
      if (any(same)) {
        j <- which(same)[1]
      } else {
        todo[[length(todo) + 1L]] <- sens_or_ratio_input(tr)
        have <- c(have, k)
        j <- length(todo)
      }
    }
    if (is.na(j)) {
      cm_abort(sprintf("No input '%s' to set a range for (see ?cm_sensitivity for the keys).", k),
               class = "deconflate_invalid_input")
    }
    todo[[j]]$low <- lo
    todo[[j]]$high <- hi
    todo[[j]]$how <- "given"
  }
  todo
}

# One run of cm_sensitivity(): the result, or the condition.
sens_run <- function(ctx, inp, v) {
  k <- inp$key
  if (!is.finite(v)) {
    return(errorCondition("The range of this input could not be set.",
                          class = c("deconflate_invalid_input", "deconflate_error")))
  }
  if (inp$kind == "entered") return(rerun_result(ctx, stats::setNames(list(v), k)))
  model <- ctx$model
  m2 <- tryCatch({
    if (startsWith(k, "pair:")) {
      p <- strsplit(sub("^pair:", "", k), ":", fixed = TRUE)[[1]]
      set_association(model, p[1], p[2], v, "OR")
    } else if (startsWith(k, "interaction:")) {
      p <- strsplit(sub("^interaction:", "", k), ";", fixed = TRUE)[[1]]
      set_interaction(model, p[1], p[2], v)
    } else {
      p <- strsplit(sub("^or_ratio:", "", k), "[:|]")[[1]]
      set_three_way(model, p[1], p[2], p[3], v)
    }
  }, deconflate_error = function(e) e)
  if (inherits(m2, "condition")) return(m2)
  rerun_result(ctx, model = m2, frozen = NULL, same_joint = startsWith(k, "interaction:"))
}

# Variance shares of the drawn inputs: the total regressed on their normal
# scores over the accepted draws (squared standardised coefficients).
sens_variance <- function(res, model) {
  dr <- res$draws
  if (is.null(dr) || is.null(dr$values) || is.null(dr$uniforms) || !("total" %in% colnames(dr$values))) return(NULL)
  y <- dr$values[, "total"]
  U <- dr$uniforms
  fixed <- vapply(colnames(U), function(k) identical(dr$specs[[k]]$type, "fixed"), logical(1))
  U <- U[, !fixed, drop = FALSE]
  if (!ncol(U)) return(NULL)
  ok <- is.finite(y) & apply(U, 1, function(u) all(is.finite(u)))
  y <- y[ok]
  U <- U[ok, , drop = FALSE]
  if (length(y) < ncol(U) + 3L || stats::sd(y) == 0) return(NULL)
  Z <- stats::qnorm(pmin(pmax(U, 1e-12), 1 - 1e-12))
  vary <- apply(Z, 2, stats::sd) > 0
  Z <- Z[, vary, drop = FALSE]
  fit <- stats::lm.fit(cbind(1, Z), y)
  beta <- fit$coefficients[-1]
  beta[is.na(beta)] <- 0
  share <- (beta * apply(Z, 2, stats::sd) / stats::sd(y))^2
  r2 <- 1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
  keys <- vapply(colnames(Z), function(k) {
    if (k == "risk" || startsWith(k, "between:") || grepl("^(prob|associations|impacts):", k)) k else
      tryCatch(canonical_key(k, model), deconflate_error = function(e) k)
  }, character(1), USE.NAMES = FALSE)
  lab <- vapply(keys, function(k) {
    if (startsWith(k, "between:")) sprintf("between-study variation (%s)", sub("^between:", "", k)) else
      tryCatch(input_labels(model, k, res$settings$risk$kind %||% "overall"), error = function(e) k)
  }, character(1), USE.NAMES = FALSE)
  out <- data.frame(input = keys, label = lab, share = unname(share), stringsAsFactors = FALSE)
  out <- out[order(-out$share), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "r_squared") <- r2
  out
}

#' @export
print.cm_sensitivity <- function(x, n = 10, digits = 3, ...) {
  tab <- x$inputs
  units <- if (!is.null(x$units) && !is.na(x$units) && !x$event) sprintf(" [%s]", x$units) else ""
  cat_wrapped(sprintf("Sensitivity of the %s%s: central value %s (ranked first: %s)", x$metric, units,
                      fmt_num(x$central$total, digits), x$central$first))
  if (!nrow(tab)) {
    cat("No inputs were varied.\n")
    return(invisible(x))
  }
  show <- utils::head(tab, if (is.finite(n)) n else nrow(tab))
  cat(sprintf("%d inputs varied one at a time (%d entered, %d unknown)%s:\n", nrow(tab), sum(tab$kind == "entered"),
              sum(tab$kind == "unknown"),
              if (nrow(show) < nrow(tab)) sprintf("; the %d largest swings", nrow(show)) else ""))
  f <- function(v) fmt_num(v, digits)
  out <- data.frame(input = show$label,
                    low = f(show$low), high = f(show$high), `at low` = f(show$total_low),
                    `at high` = f(show$total_high), swing = f(show$swing),
                    rank = ifelse(is.na(show$rank_shift), "", show$rank_shift),
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (any(is.finite(show$share))) {
    out$share <- f(show$share)
  }
  failed <- !is.na(show$failed)
  out$input[failed] <- paste0(out$input[failed], " *")
  # The labels take the width the numbers leave (at least 16 characters).
  rest <- sum(vapply(names(out)[-1], function(nm) max(nchar(c(nm, out[[nm]]))) + 1L, numeric(1)))
  room <- max(16L, getOption("width", 80L) - rest - 2L)
  long <- nchar(out$input) > room
  out$input[long] <- paste0(substr(out$input[long], 1L, room - 3L), "...")
  w <- max(nchar(c("input", out$input)))
  out$input <- formatC(out$input, width = -w)
  names(out)[1] <- formatC("input", width = -w)
  print(out, row.names = FALSE)
  if (nrow(tab) > nrow(show)) cat(sprintf("... %d more: print(x, n = Inf) shows all; keys and ranges in x$inputs.\n",
                                          nrow(tab) - nrow(show)))
  notes <- c(sprintf("low, high: entered inputs move 1 SD (with a distribution) or by %g%%; unknown inputs over their default ranges (see ?cm_sensitivity) or the ranges given. at low, at high: the %s. rank: the most places any disease moves in the ranking by contribution.",
                     100 * x$settings$change, x$metric))
  if (any(failed)) {
    notes <- c(notes, sprintf("* could not be run at low or high: %s.",
                              paste(sprintf("%s (%s)", show$label[failed], show$failed[failed]), collapse = "; ")))
  }
  ch <- nzchar(show$first)
  if (any(ch)) {
    notes <- c(notes, sprintf("First place changes: %s.",
                              paste(sprintf("%s: %s", trimws(show$label[ch]), show$first[ch]), collapse = "; ")))
  }
  if (!is.null(x$variance)) {
    r2 <- x$r_squared
    notes <- c(notes, sprintf("share: share of the variance of the %s over the draws (linear fit R^2 = %.2f%s).",
                              x$metric, r2, if (is.finite(r2) && r2 < 0.8) "; the shares are rough" else ""))
  }
  for (nt in notes) cat(strwrap(nt, width = getOption("width", 80L)), sep = "\n")
  invisible(x)
}

#' @rdname plots
#' @param n For `cm_sensitivity` plots: the number of inputs to show.
#' @param main For `cm_sensitivity` plots: the title.
#' @export
plot.cm_sensitivity <- function(x, n = 10, main = "Sensitivity (one input at a time)", ...) {
  d <- x$inputs[is.finite(x$inputs$swing), , drop = FALSE]
  if (!nrow(d)) cm_abort("No input has a swing to plot.", class = "deconflate_invalid_input")
  d <- utils::head(d, n)
  d <- d[rev(seq_len(nrow(d))), , drop = FALSE]
  base <- x$central$total
  lo <- pmin(d$total_low, d$total_high)
  hi <- pmax(d$total_low, d$total_high)
  k <- nrow(d)
  op <- graphics::par(mar = c(4, 14, 3, 1))
  on.exit(graphics::par(op))
  graphics::plot(base, 1, type = "n", xlim = range(c(lo, hi, base)), ylim = c(0.5, k + 0.5),
                 yaxt = "n", xlab = sprintf("%s%s", tools_title(x$metric),
                                            if (!is.null(x$units) && !is.na(x$units) && !x$event) sprintf(" (%s)", x$units) else ""),
                 ylab = "", main = main, ...)
  cols <- ifelse(d$kind == "unknown", "#DDCC77", "#88CCEE")
  graphics::rect(lo, seq_len(k) - 0.4, hi, seq_len(k) + 0.4, col = cols, border = NA)
  graphics::axis(2, at = seq_len(k), labels = d$label, las = 1, cex.axis = 0.7)
  graphics::abline(v = base)
  kinds <- intersect(c("entered", "unknown"), d$kind)
  graphics::legend("bottomright", legend = kinds, fill = c(entered = "#88CCEE", unknown = "#DDCC77")[kinds],
                   border = NA, bty = "n", cex = 0.8)
  invisible(x)
}

tools_title <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))
