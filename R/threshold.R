#' Threshold search: where does a conclusion change?
#'
#' Varies one input over a range, with everything else fixed at the model's
#' values, and finds the input values at which a conclusion changes:
#' * `"rank"`: two diseases swap places, ranked by their contribution to the
#'   total (including interaction shares; `by = "adjusted"` ranks by
#'   adjusted impact instead);
#' * `"sign"`: a disease's adjusted impact crosses zero (for event impacts,
#'   its adjusted hazard ratio crosses 1). This is a sign change implied by
#'   the model and its inputs, not by itself evidence of a protective effect;
#' * `"total"`: the loss of a set of diseases (the sum of their
#'   contributions; all diseases when `diseases` is left out, which is the
#'   adjusted aggregate) crosses `target`. The target is a value in the units
#'   of the impacts, or, with `relative_to`, a fraction of the set's adjusted
#'   value (`"adjusted"`, at the model's own inputs) or of its raw value
#'   (`"raw"`, the diseases unrelated): `target = 0.8, relative_to =
#'   "adjusted"` finds where the loss is 20% below its current value.
#'
#' For event impacts (`event_model = TRUE`) the total is the risk
#' attributable to disease, and the contributions are its Shapley
#' allocation. A set's loss includes its shares of the interactions with
#' diseases outside the set; it is the set's part of the current total, not
#' the loss there would be if only those diseases existed.
#'
#' Without `input`, the inputs are chosen automatically: the `top` inputs
#' (3 by default) whose low and high values in [cm_sensitivity()] move the
#' compared quantity most (the set's loss for `"total"`, the difference
#' between two diseases for `"rank"`, a disease's adjusted value for
#' `"sign"`), each searched over its default range (see Search ranges).
#' Several inputs can also be named. The result is then a
#' `cm_threshold_set`: one table of the inputs and one of the crossings.
#'
#' @section Search ranges:
#' When `range` (or a row of `ranges`) is not given, each input is searched
#' over a default range:
#' * entered inputs with a distribution: its 0.5% and 99.5% quantiles;
#' * entered inputs without one: the value minus and plus 50%
#'   (probabilities stay inside (0, 1), risk differences and phi inside
#'   (-1, 1)); a value of zero needs a range;
#' * a pair that no row mentions (`pair:`): odds ratios 0.2 to 5;
#' * an interaction without a combination row (`interaction:`): minus to
#'   plus the smaller of the two adjusted impacts (in absolute value); for
#'   event impacts a hazard ratio from 0.5 to 2;
#' * a three-way association without an `OR_ratio` row (`or_ratio:`): 0.25
#'   to 4.
#'
#' @section Method:
#' The response need not be monotonic, so the range is first scanned on a
#' grid. Every grid point is adjusted and classified: `"ok"`, or the reason
#' it could not be used (`"infeasible"` inputs, a `"singular"` system, an
#' `"undefined"` (non-finite) result, an `"unresolved"` joint fit, an
#' `"unsupported"` combination or another `"error"`); an input value that is
#' itself invalid (e.g. a non-positive odds ratio or a probability of 1) is
#' `"infeasible"`. Each pair of neighbouring usable grid points whose
#' conclusion differs brackets a crossing, which is refined by bisection. A
#' grid point where the compared quantity is exactly zero (or a stretch of
#' such points) is a crossing only if the quantity has opposite signs at the
#' usable points on either side; the threshold is then the first zero point.
#' A refined crossing is reported as a `"threshold"` only if the quantity
#' being compared is close to zero on both sides of the final bracket
#' (continuity); a jump across a pole (e.g. where the system becomes
#' singular) is reported as a `"discontinuity"`, never as a threshold; so is
#' a bracket whose bisection reaches a singular point. If a bisection step
#' lands on any other unusable point, the
#' crossing is `"unresolved"`. All crossings in the range are reported;
#' stretches of unusable grid points are listed in `regions`, so that "no
#' crossing found" can be told apart from "part of the range could not be
#' evaluated". A sign change between usable points separated by unusable
#' ones is not reported as a crossing; such stretches appear in `regions`.
#'
#' Each point is adjusted as [deconflate()] would (point estimates). When
#' the input is an impact row, an interaction or the overall risk (which do
#' not change the joint distribution) and the adjustment uses the joint
#' distribution, it is fitted once and reused (or the `joint` passed through
#' `...` is used).
#'
#' Thresholds are deterministic: they are for the model's input values. For
#' the probability that a conclusion holds when the other inputs are
#' uncertain, use [cm_threshold_prob()].
#'
#' @param model A [deconflate()] result (its settings are used), or a
#'   [cm_model()].
#' @param input The input to vary, named as in [cm_sensitivity()]:
#'   `"prob:<disease>"` (the disease's `value`, on the scale it was
#'   entered), `"associations:<row>"` or `"impacts:<row>"` (any row, by its
#'   number), `"risk"` (event impacts: the overall or baseline risk,
#'   whichever is given), or an input that is not entered:
#'   `"pair:<d1>:<d2>"` (a pair's crude pair row on its measure; a pair
#'   without one is varied as an odds ratio, so unknown pairs can be
#'   explored), `"interaction:<d1>;<d2>"` (an interaction: a row for the
#'   combination adjusted for all; for event impacts a hazard ratio) or
#'   `"or_ratio:<d1>:<d2>|<d3>"` (how much disease 3 changes the
#'   association of diseases 1 and 2). Several keys search each in turn;
#'   left out, the inputs are chosen automatically (see `top`).
#' @param range Two numbers: the range of the input to search (one input
#'   only). Left out, the default range (see Search ranges).
#' @param conclusion `"rank"`, `"sign"` or `"total"`.
#' @param target For `"total"`: the value of the set's loss (for event
#'   impacts, its attributable risk), or with `relative_to` a fraction of
#'   its adjusted or raw value (e.g. `0.8`).
#' @param diseases Disease ids: for `"total"`, the set whose loss is
#'   compared; for `"rank"` the pairs among them; for `"sign"` the diseases
#'   checked. Left out: all diseases.
#' @param relative_to For `"total"`: `NULL` (default; `target` is a value),
#'   `"adjusted"` or `"raw"` (`target` is a fraction of the set's adjusted
#'   or raw value at the model's own inputs).
#' @param by For `"rank"`: `"contribution"` (default) or `"adjusted"`.
#' @param top Without `input`: how many inputs to search (default 3).
#' @param inputs Without `input`: the candidates, `"entered"`, `"unknown"`
#'   or both (default), as in [cm_sensitivity()].
#' @param ranges Optional data frame with columns `input`, `low` and `high`:
#'   the search ranges of those inputs (instead of the defaults).
#' @param sensitivity Without `input`: a [cm_sensitivity()] result of the
#'   same model to choose from (it is computed otherwise).
#' @param three_way Without `input`: three-way inputs to consider, as in
#'   [cm_sensitivity()].
#' @param event_model,overall_risk,baseline_risk With a model: as in
#'   [deconflate()] (a distribution of the risk is used at its mean). With a
#'   result they are taken from it.
#' @param n_grid Number of grid points.
#' @param log_scale Use a log-spaced grid? Default `TRUE` for associations
#'   measured as ratios and three-way ratios with a positive range.
#' @param tol Relative tolerance on the threshold. (The bracket is also
#'   narrowed to at most 1e-6 of the grid spacing, so that the continuity
#'   test is reliable.) Because of this argument, the `tol` of [fit_joint()]
#'   cannot be passed through `...`.
#' @param ... `joint`, or arguments of [fit_joint()] (e.g.
#'   `backend = "sampled"`).
#' @return For one input, a `cm_threshold` object: `thresholds` (one row per
#'   crossing: conclusion, item, status, threshold, bracket `lower`/`upper`,
#'   the compared quantity just below and above, and a description),
#'   `details` (for each crossing, the adjusted results at the two ends of
#'   its bracket), `scan` (each grid point with its status, the aggregate and
#'   the condition number), `regions` (stretches of unusable grid points),
#'   `summary` (per item: number of thresholds, and whether part of the
#'   range could not be evaluated), and the settings (`input` and its
#'   `label`, `range`, `baseline` value of the input, `method` (the
#'   computations used), `conclusion`, `target`, `level` (the value the loss
#'   is compared with), `diseases`, `relative_to`, `fixed`, and the impact
#'   table's `outcome` and `units`).
#'
#'   For several inputs (or none named), a `cm_threshold_set`: `inputs` (one
#'   row per input: key, label, kind, `effect` (how much it moved the
#'   compared quantity in [cm_sensitivity()], when chosen automatically),
#'   the range searched and how it was set, the baseline value, the number
#'   of thresholds and the result), `thresholds` (the crossings of all
#'   inputs, with an `input` column) and `searches` (each input's
#'   `cm_threshold`).
#' @export
#' @examples
#' m <- example_supplement()
#' # At what odds ratio between d1 and d2 does the ranking change?
#' th <- cm_threshold(m, "pair:d1:d2", c(0.2, 20), conclusion = "rank")
#' th
#' # What odds ratio between d1 and d3 (an odds ratio of 1 in the example)
#' # would reduce the aggregate by 10%?
#' cm_threshold(m, "pair:d1:d3", c(1, 50), conclusion = "total", target = 0.9,
#'              relative_to = "adjusted")$thresholds
#' # The loss of d2 and d3 together, from a result: its settings are used
#' res <- deconflate(m, n_draws = 0)
#' cm_threshold(res, "pair:d2:d3", conclusion = "total", diseases = c("d2", "d3"), target = 1.8)
#' # The three inputs that matter most for the ranking of d1 and d2
#' cm_threshold(res, conclusion = "rank", diseases = c("d1", "d2"), top = 3)
cm_threshold <- function(model, input = NULL, range = NULL, conclusion = c("rank", "sign", "total"),
                         target = NULL, diseases = NULL, relative_to = NULL, by = c("contribution", "adjusted"),
                         top = 3, inputs = c("entered", "unknown"), ranges = NULL, sensitivity = NULL,
                         three_way = NULL, event_model = FALSE, overall_risk = NULL,
                         n_grid = 101L, log_scale = NULL, tol = 1e-8, baseline_risk = NULL, ...) {
  dots <- list(...)
  retired_args(dots)
  if (identical(conclusion, "change")) {
    cm_abort("conclusion = 'change' has been replaced: use conclusion = 'total', relative_to = 'adjusted' and target = 1 + the change (e.g. target = 0.9 for 10% below the current total).",
             class = "deconflate_unsupported")
  }
  th <- threshold_setup(model, event_model, overall_risk, baseline_risk, dots)
  spec <- threshold_spec(conclusion, target, diseases, relative_to, by, th$model)
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 3) {
    cm_abort("`n_grid` must be a single number, at least 3.", class = "deconflate_invalid_input")
  }
  n_grid <- as.integer(n_grid)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    cm_abort("`tol` must be a single positive number.", class = "deconflate_invalid_input")
  }
  if (!is.null(log_scale) && (!is.logical(log_scale) || length(log_scale) != 1L || is.na(log_scale))) {
    cm_abort("`log_scale` must be NULL, TRUE or FALSE.", class = "deconflate_invalid_input")
  }
  if (!is.null(input) && (!is.character(input) || !length(input) || anyNA(input))) {
    cm_abort("`input` must be an input key such as \"pair:d1:d2\", several keys, or left out.",
             class = "deconflate_invalid_input")
  }
  sel <- threshold_inputs(th, spec, input, range, top, inputs, ranges, sensitivity, three_way)
  if (!is.null(input) && length(input) == 1L) {
    return(threshold_search(th, sel$inputs$input[1], c(sel$inputs$low[1], sel$inputs$high[1]), spec,
                            n_grid, log_scale, tol))
  }
  searches <- lapply(seq_len(nrow(sel$inputs)), function(i) {
    if (!is.finite(sel$inputs$low[i])) return(NULL)
    tryCatch(threshold_search(th, sel$inputs$input[i], c(sel$inputs$low[i], sel$inputs$high[i]), spec,
                              n_grid, log_scale, tol),
             deconflate_error = function(e) e)
  })
  names(searches) <- sel$inputs$input
  threshold_set(sel, searches, spec)
}

# Settings of a threshold search: the model, the risk and how it is given,
# the fit arguments, and the central result when the model came as one (or
# a function that computes it when needed).
threshold_setup <- function(model, event_model, overall_risk, baseline_risk, dots) {
  if (inherits(model, "cm_result") || inherits(model, "cm_event_result")) {
    if (!isFALSE(event_model) || !is.null(overall_risk) || !is.null(baseline_risk)) {
      cm_abort("With a deconflate() result, the settings come from the result: leave out event_model, overall_risk and baseline_risk.",
               class = "deconflate_invalid_input")
    }
    ctx <- result_context(model)
    st <- list(event_model = ctx$event, risk = ctx$risk$value, risk_kind = ctx$risk$kind %||% "overall",
               pool = ctx$plan$pool %||% "random")
    return(list(model = ctx$model, st = st, dots = utils::modifyList(ctx$dots, dots), risk_dist = ctx$risk$dist,
                res = model))
  }
  st <- sens_setup(model, event_model, overall_risk, baseline_risk = baseline_risk)
  st$pool <- "random"
  risk <- check_overall_risk(overall_risk, event_model, baseline_risk)
  list(model = model, st = st, dots = dots, risk_dist = risk$dist, res = NULL,
       args = list(event_model = event_model, overall_risk = overall_risk, baseline_risk = baseline_risk))
}

# The central result of a threshold setup (computed when the model came
# without one).
threshold_central <- function(th) {
  if (!is.null(th$res)) return(th$res)
  dots <- th$dots[intersect(names(th$dots), c("pool", "joint", joint_arg_names))]
  sens_central(th$model, th$args$event_model, th$args$overall_risk, th$args$baseline_risk, dots)
}

# What is compared: the conclusion with its target, set of diseases and
# reference, checked against the model.
threshold_spec <- function(conclusion, target, diseases, relative_to, by, model) {
  conclusion <- cm_match_arg(conclusion, c("rank", "sign", "total"))
  by <- cm_match_arg(by, c("contribution", "adjusted"))
  ids <- model$diseases$disease
  if (!is.null(diseases)) {
    if (!is.character(diseases) || !length(diseases) || anyNA(diseases)) {
      cm_abort("`diseases` must be a character vector of disease ids.", class = "deconflate_invalid_input")
    }
    check_disease_ids(diseases, model)
    diseases <- unique(diseases)
  }
  if (conclusion == "rank" && length(diseases %||% ids) < 2L) {
    cm_abort("conclusion = 'rank' needs at least two diseases.", class = "deconflate_invalid_input")
  }
  if (conclusion == "total") {
    if (is.null(target) || !is.numeric(target) || length(target) != 1L || !is.finite(target)) {
      cm_abort("conclusion = 'total' needs a single finite `target`.", class = "deconflate_invalid_input")
    }
    if (!is.null(relative_to)) relative_to <- cm_match_arg(relative_to, c("adjusted", "raw"))
  } else {
    if (!is.null(target)) cm_abort(sprintf("`target` is used only with conclusion = 'total' (not '%s').", conclusion),
                                    class = "deconflate_invalid_input")
    if (!is.null(relative_to)) cm_abort("`relative_to` is used only with conclusion = 'total'.",
                                        class = "deconflate_invalid_input")
  }
  # A set of all diseases is the aggregate itself.
  set <- if (conclusion == "total" && !is.null(diseases) && length(diseases) < length(ids)) ids[ids %in% diseases]
  list(conclusion = conclusion, target = target, diseases = diseases, set = set, relative_to = relative_to,
       by = by, ids = ids)
}

# The loss of a set of diseases in a result: the sum of their contributions
# (NULL set: the total).
set_total <- function(res, set) {
  bm <- burden_metric(res)
  if (is.null(set)) bm$total else sum(bm$by_disease[set])
}

# The level the loss is compared with: the target, or the target times the
# set's adjusted or raw value at the model's own inputs (`base`: the central
# result of the model as given).
threshold_level <- function(spec, model, st, dots, base) {
  if (spec$conclusion != "total") return(list(level = NA_real_, reference = NA_real_))
  if (is.null(spec$relative_to)) return(list(level = spec$target, reference = NA_real_))
  if (is.null(base)) {
    cm_abort("The model as given cannot be adjusted, so a relative target has no reference value.",
             class = "deconflate_nonfinite")
  }
  ref <- if (spec$relative_to == "adjusted") set_total(base, spec$set) else {
    plan <- plan_model(model, st$event_model, dots = dots[intersect(names(dots), joint_arg_names)],
                       pool = st$pool, risk_kind = st$risk_kind)
    rt <- raw_total(model, plan, st$risk, base$row_se, by_disease = !is.null(spec$set))
    if (is.null(spec$set)) rt else sum(rt[spec$set])
  }
  if (!is.finite(ref) || ref <= 0) {
    cm_abort(sprintf("The %s value of the %s is %s, so a target relative to it is undefined; give the target as a value.",
                     spec$relative_to, if (is.null(spec$set)) "total" else "set",
                     if (!is.finite(ref)) "not available" else if (ref == 0) "zero" else sprintf("negative (%s)", fmt_num(ref))),
             class = "deconflate_nonfinite")
  }
  list(level = spec$target * ref, reference = ref)
}

# The items compared for a conclusion: one per disease ("sign"), per pair
# ("rank") or one ("total"), with what they compare.
threshold_items <- function(spec, level) {
  use <- spec$diseases %||% spec$ids
  switch(spec$conclusion,
    total = data.frame(item = if (is.null(spec$set)) "total" else sprintf("total of %s", paste(spec$set, collapse = ", ")),
                       a = NA_character_, b = NA_character_, level = level, stringsAsFactors = FALSE),
    sign = data.frame(item = use, a = use, b = NA_character_, level = 0, stringsAsFactors = FALSE),
    rank = {
      if (length(use) < 2L) cm_abort("conclusion = 'rank' needs at least two diseases.", class = "deconflate_invalid_input")
      pr <- utils::combn(use, 2)
      data.frame(item = paste(pr[1, ], pr[2, ], sep = " vs "), a = pr[1, ], b = pr[2, ],
                 level = 0, stringsAsFactors = FALSE)
    })
}

# The compared quantities of a result: each crosses zero where the
# conclusion changes.
threshold_metric <- function(res, spec, items) {
  event <- inherits(res, "cm_event_result")
  switch(spec$conclusion,
    total = stats::setNames(set_total(res, spec$set) - items$level, items$item),
    sign = {
      a <- res$adjusted
      v <- a$adjusted[match(items$a, a$disease)]
      stats::setNames(if (event) log(v) else v, items$item)
    },
    rank = {
      sc <- if (spec$by == "contribution") burden_metric(res)$by_disease else
        stats::setNames(res$adjusted$adjusted, res$adjusted$disease)
      stats::setNames(unname(sc[items$a] - sc[items$b]), items$item)
    })
}

# The words for what is compared with the target: "the aggregate", "the
# attributable risk of LAM and SCK", ...
total_words <- function(spec, event) {
  what <- if (event) "attributable risk" else "aggregate"
  if (is.null(spec$set)) sprintf("the %s", what) else
    sprintf("the %s of %s", if (event) "attributable risk" else "loss", join_and(spec$set))
}

# The target in words: "2.00", "0.8 times its adjusted value (1.64)".
target_words <- function(spec, level, reference) {
  if (is.null(spec$relative_to)) return(fmt_num(level))
  sprintf("%s times its %s value (%s)", format(signif(spec$target, 3)), spec$relative_to, fmt_num(level))
}

# One threshold search over `range` for one input.
threshold_search <- function(th, input, range, spec, n_grid, log_scale, tol) {
  model <- th$model
  st <- th$st
  dots <- th$dots
  check_numeric(range, "range")
  if (length(range) != 2L || range[1] >= range[2]) cm_abort("`range` must be two increasing numbers.", class = "deconflate_invalid_input")
  setter <- threshold_setter(model, input, st)
  if (setter$type == "inter" && isTRUE(st$event_model) && range[1] <= 0) {
    cm_abort("An interaction of event impacts is a hazard ratio: `range` must be positive.",
             class = "deconflate_invalid_input")
  }
  log_scale <- log_scale %||% (setter$ratio_scale && range[1] > 0)
  if (log_scale && range[1] <= 0) cm_abort("A log-spaced grid needs a positive range.", class = "deconflate_invalid_input")
  xs <- if (log_scale) exp(seq(log(range[1]), log(range[2]), length.out = n_grid)) else
    seq(range[1], range[2], length.out = n_grid)
  xs[c(1L, n_grid)] <- range   # exact ends (exp(log(x)) can differ from x)

  # Impacts, interactions and the overall risk do not change the joint
  # distribution, so for those inputs it is fitted once (when it is used).
  joint <- dots[["joint"]]
  dots[["joint"]] <- NULL
  fj <- dots[intersect(names(dots), joint_arg_names)]
  fixed_joint <- setter$type %in% c("impact", "inter", "risk")
  if (!is.null(joint) && !fixed_joint) {
    cm_abort("A fixed `joint` can be used only when varying an impact row, an interaction or the overall risk; other inputs change the joint distribution.", class = "deconflate_unsupported")
  }
  # (Planned on the model at the start of the range: varying an interaction
  # can add a combination row.)
  m_plan <- tryCatch(setter$set(range[1])$model, deconflate_error = function(e) model)
  base_plan <- plan_model(m_plan, st$event_model, dots = fj, pool = st$pool, risk_kind = st$risk_kind)
  uses_joint <- isTRUE(base_plan$event) || !base_plan$shortcut
  if (is.null(joint) && fixed_joint && uses_joint) {
    fa <- utils::modifyList(fj, base_plan$fit_args)
    joint <- tryCatch(
      withCallingHandlers(do.call(fit_joint, c(list(as_population(model)), fa)),
                          deconflate_nonconvergence = function(w) {
                            if (inherits(w, "warning")) invokeRestart("muffleWarning")
                          }),
      deconflate_error = function(e) NULL)   # each point then reports its own failure
  }

  evaluate <- function(x) {
    m2 <- tryCatch(setter$set(x), deconflate_error = function(e) e)
    if (inherits(m2, "condition")) {
      # The setter fails only for an invalid value of the input itself.
      s2 <- threshold_status(m2)
      return(list(status = if (s2 == "error") "infeasible" else s2, res = NULL))
    }
    evaluate_model(m2$model, m2$st)
  }
  used <- character(0)
  evaluate_model <- function(m2, st2) {
    res <- tryCatch(
      withCallingHandlers({
        plan <- plan_model(m2, st2$event_model, dots = fj, has_joint = fixed_joint && !is.null(joint),
                           pool = st2$pool, risk_kind = st2$risk_kind)
        if (fixed_joint && !is.null(joint)) plan$shortcut <- FALSE
        r <- run_point(m2, plan, risk = st2$risk, joint = if (fixed_joint) joint, warn = FALSE, dots = fj,
                       raw = FALSE)
        used <<- union(used, r$diagnostics$method)
        r
      }, deconflate_nonconvergence = function(w) {
        if (inherits(w, "warning")) invokeRestart("muffleWarning")
      }),
      deconflate_error = function(e) e)
    if (inherits(res, "condition")) return(list(status = threshold_status(res), res = NULL))
    if (!result_is_finite(res)) return(list(status = "undefined", res = res))
    list(status = "ok", res = res)
  }
  total_of <- function(res) burden_metric(res)$total

  # The compared quantities: each must cross zero where the conclusion
  # changes. A relative target refers to the model as given (the input at
  # its own value).
  base_eval <- if (!is.null(spec$relative_to)) evaluate_model(model, st)
  lv <- threshold_level(spec, model, st, dots, if (identical(base_eval$status, "ok")) base_eval$res)
  items <- threshold_items(spec, lv$level)
  metric <- function(res) threshold_metric(res, spec, items)

  evals <- lapply(xs, evaluate)
  status <- vapply(evals, `[[`, character(1), "status")
  ok <- status == "ok"
  vals <- matrix(NA_real_, n_grid, nrow(items), dimnames = list(NULL, items$item))
  for (g in which(ok)) vals[g, ] <- metric(evals[[g]]$res)
  scan <- data.frame(
    value = xs, status = status,
    total = vapply(evals, function(e) if (identical(e$status, "ok")) total_of(e$res) else NA_real_,
                   numeric(1)),
    condition_number = vapply(evals, function(e) {
      if (is.null(e$res)) NA_real_ else e$res$diagnostics$condition_number
    }, numeric(1)),
    stringsAsFactors = FALSE)

  # Scale of each compared quantity, for the continuity test.
  scale_f <- apply(vals, 2, function(v) {
    v <- v[is.finite(v)]
    if (length(v) && max(abs(v)) > 0) max(abs(v)) else 1
  })
  width_tol <- function(lo, hi) tol * max(abs(lo), abs(hi), 1)
  event <- isTRUE(base_plan$event)

  rows <- list()
  details <- list()
  for (k in seq_len(nrow(items))) {
    f <- vals[, k]
    br <- threshold_brackets(f, ok)
    for (b in seq_len(nrow(br))) {
      g <- br[b, 1L]
      # A crossing through exact zeros on the grid (points g + 1, ..., j - 1,
      # one or more): the bracket ends at the last zero and the threshold is
      # the first; no bisection is needed.
      zero_run <- br[b, 2L] > g + 1L
      h <- if (zero_run) br[b, 2L] - 1L else g + 1L
      lo <- xs[g]; hi <- xs[h]
      flo <- f[g]; fhi <- f[h]
      rlo <- evals[[g]]$res; rhi <- evals[[h]]$res
      stt <- "threshold"
      if (flo != 0 && fhi != 0) {
        w_stop <- 1e-6 * (hi - lo)
        for (it in seq_len(200L)) {
          if (hi - lo <= min(width_tol(lo, hi), w_stop)) break
          mid <- if (log_scale) sqrt(lo * hi) else (lo + hi) / 2
          e <- evaluate(mid)
          if (!identical(e$status, "ok")) {
            # A singular point between values of opposite sign is a pole of
            # the system (not a threshold); any other failure leaves the
            # crossing unresolved.
            stt <- if (identical(e$status, "singular")) "discontinuity" else "unresolved"
            break
          }
          fm <- metric(e$res)[[k]]
          if (fm == 0) {
            lo <- hi <- mid; flo <- fhi <- 0; rlo <- rhi <- e$res
            break
          }
          if (sign(fm) == sign(flo)) {
            lo <- mid; flo <- fm; rlo <- e$res
          } else {
            hi <- mid; fhi <- fm; rhi <- e$res
          }
        }
      }
      if (stt == "threshold" && flo != 0 && fhi != 0) {
        # Continuity: a genuine crossing has a small compared quantity on
        # both sides; a pole has large values of opposite sign.
        eps <- 1e-4 * scale_f[k]
        kap <- max(rlo$diagnostics$condition_number, rhi$diagnostics$condition_number)
        if (max(abs(flo), abs(fhi)) > eps || !is.finite(kap) || kap > 1e12) stt <- "discontinuity"
      }
      x0 <- if (zero_run) xs[g + 1L] else if (flo == 0) lo else if (fhi == 0) hi else
        if (log_scale) sqrt(lo * hi) else (lo + hi) / 2
      rows[[length(rows) + 1L]] <- data.frame(
        conclusion = spec$conclusion, item = items$item[k], status = stt,
        threshold = if (stt == "threshold") x0 else NA_real_,
        lower = lo, upper = hi, below = flo, above = fhi,
        description = threshold_description(spec, items[k, ], flo, fhi, stt, lv, event),
        stringsAsFactors = FALSE)
      details[[length(details) + 1L]] <- list(below = rlo$adjusted, above = rhi$adjusted,
                                              totals_below = threshold_totals(rlo),
                                              totals_above = threshold_totals(rhi))
    }
  }
  thresholds <- if (length(rows)) do.call(rbind, rows) else
    data.frame(conclusion = character(0), item = character(0), status = character(0),
               threshold = numeric(0), lower = numeric(0), upper = numeric(0),
               below = numeric(0), above = numeric(0), description = character(0),
               stringsAsFactors = FALSE)
  if (nrow(thresholds)) {
    o <- order(thresholds$lower)
    thresholds <- thresholds[o, , drop = FALSE]
    details <- details[o]
  }
  rownames(thresholds) <- NULL
  regions <- threshold_regions(xs, status)
  summary <- data.frame(
    item = items$item,
    n_thresholds = vapply(items$item, function(it) sum(thresholds$item == it & thresholds$status == "threshold"),
                          integer(1), USE.NAMES = FALSE),
    n_other = vapply(items$item, function(it) sum(thresholds$item == it & thresholds$status != "threshold"),
                     integer(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE)
  summary$result <- ifelse(summary$n_thresholds > 0, "threshold(s) found",
                    ifelse(summary$n_other > 0, "no threshold; see discontinuities or unresolved crossings",
                    ifelse(nrow(regions) > 0, "none found where the range could be evaluated",
                           "none found in the range")))
  rownames(summary) <- NULL
  structure(list(thresholds = thresholds, details = details, scan = scan, regions = regions,
                 summary = summary, input = input, label = threshold_label(input, model, st), range = range,
                 baseline = setter$baseline, outcome = impact_meta(model$impacts)$outcome,
                 units = impact_meta(model$impacts)$units,
                 method = paste(if (length(used)) used else "none", collapse = ", "),
                 event_model = st$event_model, conclusion = spec$conclusion, target = spec$target,
                 level = lv$level, reference = lv$reference, diseases = spec$diseases,
                 relative_to = spec$relative_to, by = spec$by, log_scale = log_scale,
                 fixed = sprintf("all inputs other than %s at the model's values", input),
                 values = vals),
            class = "cm_threshold")
}

# The totals of a result, for the details of a crossing.
threshold_totals <- function(res) {
  if (inherits(res, "cm_event_result")) res$attributable$summary else res$totals
}

# Status label of a failed evaluation.
threshold_status <- function(e) {
  switch(condition_type(e),
         infeasible = "infeasible", singular = "singular", nonconvergence = "unresolved",
         unsupported = "unsupported", nonfinite = "undefined", "error")
}

threshold_description <- function(spec, item, flo, fhi, status, lv, event = FALSE) {
  if (status != "threshold") {
    return(switch(status,
      discontinuity = "sign change across a pole or a nearly singular system (not a threshold)",
      unresolved = "the crossing could not be refined: an unusable point lies inside the bracket"))
  }
  up <- fhi > flo
  switch(spec$conclusion,
    total = sprintf("%s %s %s", total_words(spec, event), if (up) "rises above" else "falls below",
                    target_words(spec, lv$level, lv$reference)),
    sign = sprintf("the adjusted value of %s moves %s no effect (a model-implied sign change)",
                   item$a, if (up) "above" else "below"),
    rank = sprintf("%s %s %s (by %s)", item$a, if (up) "moves above" else "moves below", item$b, spec$by))
}

# Default search ranges (see ?cm_threshold, Search ranges).
thr_unknown_pair <- c(0.2, 5)
thr_unknown_hr <- c(0.5, 2)
thr_unknown_or_ratio <- c(0.25, 4)

# The inputs to search, with their ranges: the named ones, or the `top`
# chosen from cm_sensitivity(). Returns a data frame (input, label, kind,
# effect, low, high, how, value).
threshold_inputs <- function(th, spec, input, range, top, inputs, ranges, sensitivity, three_way) {
  model <- th$model
  if (!is.null(range) && (is.null(input) || length(input) != 1L)) {
    cm_abort("`range` is for one input; give several ranges with `ranges` (a data frame input, low, high).",
             class = "deconflate_invalid_input")
  }
  if (!is.null(ranges) && (!is.data.frame(ranges) || !all(c("input", "low", "high") %in% names(ranges)))) {
    cm_abort("`ranges` must be a data frame with columns input, low and high.", class = "deconflate_invalid_input")
  }
  central <- NULL
  get_central <- function() {
    if (is.null(central)) central <<- threshold_central(th)
    central
  }
  if (is.null(input)) {
    if (!is.numeric(top) || length(top) != 1L || !is.finite(top) || top < 1 || top != round(top)) {
      cm_abort("`top` must be a whole number, 1 or more.", class = "deconflate_invalid_input")
    }
    s <- sensitivity
    if (is.null(s)) {
      s <- cm_sensitivity(get_central(), inputs = inputs, three_way = three_way)
    } else if (!inherits(s, "cm_sensitivity")) {
      cm_abort("`sensitivity` must be a cm_sensitivity() result.", class = "deconflate_invalid_input")
    } else if (!all(unique(s$by_disease$disease) %in% model$diseases$disease)) {
      cm_abort("`sensitivity` is for another model (its diseases differ).", class = "deconflate_invalid_input")
    }
    pick <- threshold_select(s, spec, top, isTRUE(th$st$event_model))
    if (!nrow(pick)) {
      cm_abort("No input moved the compared quantity in cm_sensitivity(), so there is nothing to search.",
               class = "deconflate_invalid_input")
    }
    keys <- pick$input
  } else {
    keys <- unique(input)
    # Each key is checked as the search would check it (retired keys say
    # what to write instead).
    for (k in keys) threshold_setter(model, k, th$st)
    pick <- data.frame(input = keys, label = NA_character_, kind = NA_character_, effect = NA_real_,
                       stringsAsFactors = FALSE)
  }
  rows <- lapply(seq_along(keys), function(i) {
    k <- keys[i]
    given <- if (!is.null(range)) range else if (!is.null(ranges) && k %in% ranges$input) {
      r <- ranges[match(k, ranges$input), ]
      c(r$low, r$high)
    }
    d <- if (!is.null(given)) {
      list(low = given[1], high = given[2], how = "given", value = NA_real_)
    } else {
      threshold_default_range(k, th, get_central)
    }
    data.frame(low = d$low, high = d$high, how = d$how, value = d$value, stringsAsFactors = FALSE)
  })
  out <- cbind(pick, do.call(rbind, rows))
  # Labels and kinds of named inputs.
  lab <- vapply(out$input, function(k) threshold_label(k, model, th$st), character(1))
  out$label <- ifelse(is.na(out$label), lab, out$label)
  out$kind <- ifelse(is.na(out$kind),
                     vapply(out$input, function(k) if (is.null(entered_key(k, model))) "unknown" else "entered",
                            character(1)), out$kind)
  rownames(out) <- NULL
  if (!is.null(input) && length(input) == 1L && !is.finite(out$low[1])) {
    cm_abort(sprintf("No default range for '%s' (%s): give `range`.", input, out$how[1]), class = "deconflate_invalid_input")
  }
  list(inputs = out, auto = is.null(input), sensitivity = if (is.null(input)) s)
}

# A short label of an input key.
threshold_label <- function(k, model, st) {
  if (grepl("^(prob|associations|impacts):", k) || k == "risk") {
    return(tryCatch(input_labels(model, k, st$risk_kind %||% "overall"), error = function(e) k))
  }
  if (startsWith(k, "pair:")) {
    p <- strsplit(sub("^pair:", "", k), ":", fixed = TRUE)[[1]]
    return(sprintf("%s-%s OR", p[1], p[2]))
  }
  if (startsWith(k, "interaction:")) return(sprintf("%s interaction", sub("^interaction:", "", k)))
  if (startsWith(k, "or_ratio:")) {
    p <- strsplit(sub("^or_ratio:", "", k), "[:|]")[[1]]
    return(sprintf("%s-%s OR ratio among %s", p[1], p[2], p[3]))
  }
  k
}

# The entered row a key names ("pair:a:b" with one crude pair row,
# "interaction:a;b" with a combination row adjusted for all, "or_ratio:..."
# with an OR_ratio row), the key itself for entered keys, or NULL.
entered_key <- function(k, model) {
  if (grepl("^(prob|associations|impacts):", k) || k == "risk") return(k)
  parts <- strsplit(k, ":", fixed = TRUE)[[1]]
  if (parts[1] == "pair" && length(parts) == 3L && !is.null(model$associations)) {
    hit <- pair_row_hits(model$associations, parts[2], parts[3], model$diseases$disease)
    if (sum(hit) == 1L) return(paste0("associations:", model$associations$row[hit]))
  } else if (parts[1] == "interaction" && length(parts) == 2L) {
    ab <- strsplit(parts[2], ";", fixed = TRUE)[[1]]
    im <- model$impacts
    hit <- same_set(im$disease, paste(ab, collapse = ";")) & is_all(im$adjusted_for) & is.na(im$among)
    if (any(hit)) return(paste0("impacts:", im$row[which(hit)[1]]))
  } else if (parts[1] == "or_ratio") {
    tr <- strsplit(sub("^or_ratio:", "", k), "[:|]")[[1]]
    a <- model$associations
    if (length(tr) == 3L && !is.null(a) && or_ratio_entered(model, tr)) {
      hit <- a$measure == "OR_ratio" & vapply(seq_len(nrow(a)), function(r) {
        setequal(c(a$disease[r], a$associated_with[r], sub("^no\\s+", "", a$among[r])), tr)
      }, logical(1))
      return(paste0("associations:", a$row[which(hit)[1]]))
    }
  }
  NULL
}

# The default search range of one input (see ?cm_threshold, Search ranges):
# list(low, high, how, value). `get_central()` gives the central result
# (for the adjusted impacts that scale an unknown interaction).
threshold_default_range <- function(k, th, get_central) {
  model <- th$model
  st <- th$st
  none <- function(why) list(low = NA_real_, high = NA_real_, how = why, value = NA_real_)
  parts <- strsplit(k, ":", fixed = TRUE)[[1]]
  ek <- entered_key(k, model)
  if (!is.null(ek)) {
    ck <- if (ek == "risk") "risk" else tryCatch(canonical_key(ek, model), deconflate_error = function(e) NA_character_)
    if (is.na(ck)) return(none("not an input of the model"))
    if (ck == "risk" && !isTRUE(st$event_model)) return(none("not an event model"))
    v <- input_value(model, ck, st$risk)
    d <- input_dists(model, th$risk_dist)[[ck]]
    # An OR_ratio row among "no <d3>" holds the inverse of the ratio that
    # the key "or_ratio:..." varies.
    inv <- parts[1] == "or_ratio" && grepl("^no\\s+",
                                         model$associations$among[match(as.integer(sub("^associations:", "", ck)),
                                                                        model$associations$row)])
    if (!is.null(d)) {
      lo <- d$q(0.005)
      hi <- d$q(0.995)
      if (inv) return(list(low = 1 / hi, high = 1 / lo, how = "99% of its distribution", value = 1 / v))
      return(list(low = lo, high = hi, how = "99% of its distribution", value = v))
    }
    if (inv) v <- 1 / v
    if (!is.finite(v) || v == 0) return(none("value 0: give a range"))
    lo <- min(v * 0.5, v * 1.5)
    hi <- max(v * 0.5, v * 1.5)
    meas <- input_measure(model, ck)
    if (ck == "risk" || isTRUE(meas %in% c("prevalence", "probability"))) {
      lo <- max(lo, 1e-6)
      hi <- min(hi, 1 - 1e-6)
    }
    if (isTRUE(meas %in% c("RD", "phi"))) {
      lo <- max(lo, -1 + 1e-6)
      hi <- min(hi, 1 - 1e-6)
    }
    return(list(low = lo, high = hi, how = "value -50% to +50%", value = v))
  }
  switch(parts[1],
    pair = list(low = thr_unknown_pair[1], high = thr_unknown_pair[2],
                how = sprintf("OR %g-%g", thr_unknown_pair[1], thr_unknown_pair[2]), value = NA_real_),
    or_ratio = list(low = thr_unknown_or_ratio[1], high = thr_unknown_or_ratio[2],
                    how = sprintf("ratio %g-%g", thr_unknown_or_ratio[1], thr_unknown_or_ratio[2]), value = NA_real_),
    interaction = {
      if (isTRUE(st$event_model)) {
        list(low = thr_unknown_hr[1], high = thr_unknown_hr[2],
             how = sprintf("HR %g-%g", thr_unknown_hr[1], thr_unknown_hr[2]), value = NA_real_)
      } else {
        ab <- strsplit(parts[2], ";", fixed = TRUE)[[1]]
        res <- get_central()
        adj <- stats::setNames(res$adjusted$adjusted, res$adjusted$disease)
        if (length(ab) != 2L || !all(ab %in% names(adj))) return(none("unknown diseases"))
        w <- min(abs(adj[[ab[1]]]), abs(adj[[ab[2]]]))
        if (!is.finite(w) || w == 0) return(none("adjusted impact 0: give a range"))
        list(low = -w, high = w, how = "+/- the smaller adjusted impact", value = 0)
      }
    },
    none("unknown input"))
}

# The `top` inputs of a cm_sensitivity() result by their effect on the
# compared quantity (data frame input, label, kind, effect).
threshold_select <- function(s, spec, top, event) {
  tab <- s$inputs
  bd <- s$by_disease
  eff <- vapply(tab$input, function(k) {
    b <- bd[bd$input == k, , drop = FALSE]
    val <- function(col, ds) stats::setNames(b[[col]][match(ds, b$disease)], ds)
    switch(spec$conclusion,
      total = if (is.null(spec$set)) abs(tab$total_high[tab$input == k] - tab$total_low[tab$input == k]) else
        abs(sum(val("contribution_high", spec$set)) - sum(val("contribution_low", spec$set))),
      sign = {
        ds <- spec$diseases %||% spec$ids
        lo <- val("adjusted_low", ds)
        hi <- val("adjusted_high", ds)
        if (event) {
          lo <- log(lo)
          hi <- log(hi)
        }
        d <- abs(hi - lo)
        if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE)
      },
      rank = {
        ds <- spec$diseases %||% spec$ids
        col <- if (spec$by == "contribution") "contribution" else "adjusted"
        lo <- val(paste0(col, "_low"), ds)
        hi <- val(paste0(col, "_high"), ds)
        pr <- utils::combn(ds, 2)
        d <- abs((hi[pr[1, ]] - hi[pr[2, ]]) - (lo[pr[1, ]] - lo[pr[2, ]]))
        if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE)
      })
  }, numeric(1), USE.NAMES = FALSE)
  keep <- is.finite(eff) & eff > 0
  out <- data.frame(input = tab$input[keep], label = tab$label[keep], kind = tab$kind[keep], effect = eff[keep],
                    stringsAsFactors = FALSE)
  out <- out[order(-out$effect), , drop = FALSE]
  out <- utils::head(out, top)
  rownames(out) <- NULL
  out
}

# Combine the searches of several inputs into a cm_threshold_set.
threshold_set <- function(sel, searches, spec) {
  inp <- sel$inputs
  nth <- vapply(seq_len(nrow(inp)), function(i) {
    s <- searches[[i]]
    if (inherits(s, "cm_threshold")) sum(s$thresholds$status == "threshold") else NA_integer_
  }, integer(1))
  result <- vapply(seq_len(nrow(inp)), function(i) {
    s <- searches[[i]]
    if (is.null(s)) return(sprintf("not searched (%s)", inp$how[i]))
    if (inherits(s, "condition")) return(sprintf("could not be searched: %s", conditionMessage(s)))
    r <- s$summary$result
    if (any(r == "threshold(s) found")) "threshold(s) found" else r[1]
  }, character(1))
  inp$baseline <- vapply(seq_len(nrow(inp)), function(i) {
    s <- searches[[i]]
    if (inherits(s, "cm_threshold") && !is.null(s$baseline)) s$baseline else inp$value[i]
  }, numeric(1))
  inp$value <- NULL
  inp$n_thresholds <- nth
  inp$result <- result
  ths <- lapply(seq_len(nrow(inp)), function(i) {
    s <- searches[[i]]
    if (!inherits(s, "cm_threshold") || !nrow(s$thresholds)) return(NULL)
    cbind(data.frame(input = inp$input[i], label = inp$label[i], stringsAsFactors = FALSE), s$thresholds)
  })
  ths <- do.call(rbind, ths[!vapply(ths, is.null, logical(1))])
  if (is.null(ths)) {
    ths <- data.frame(input = character(0), label = character(0), conclusion = character(0), item = character(0),
                      status = character(0), threshold = numeric(0), lower = numeric(0), upper = numeric(0),
                      below = numeric(0), above = numeric(0), description = character(0), stringsAsFactors = FALSE)
  }
  rownames(ths) <- NULL
  first <- Find(function(s) inherits(s, "cm_threshold"), searches)
  structure(list(inputs = inp, thresholds = ths, searches = searches, auto = sel$auto,
                 conclusion = spec$conclusion, target = spec$target, diseases = spec$diseases,
                 relative_to = spec$relative_to, by = spec$by,
                 level = if (!is.null(first)) first$level else NA_real_,
                 reference = if (!is.null(first)) first$reference else NA_real_,
                 event_model = if (!is.null(first)) first$event_model else NA,
                 outcome = first$outcome, units = first$units),
            class = "cm_threshold_set")
}

# What a threshold search compares, in words (the first line of print()).
threshold_question <- function(x) {
  spec <- list(conclusion = x$conclusion, set = if (x$conclusion == "total" && length(x$diseases)) x$diseases,
               relative_to = x$relative_to, target = x$target)
  event <- isTRUE(x$event_model)
  switch(x$conclusion,
    total = sprintf("where %s crosses %s", total_words(spec, event),
                    if (is.finite(x$level)) target_words(spec, x$level, x$reference) else format(x$target)),
    sign = sprintf("where the adjusted value of %s crosses no effect",
                   if (length(x$diseases)) join_and(x$diseases) else "a disease"),
    rank = sprintf("where %s swap places (by %s)",
                   if (length(x$diseases)) join_and(x$diseases) else "two diseases", x$by))
}

#' @export
print.cm_threshold_set <- function(x, digits = 3, ...) {
  f <- function(v) fmt_num(v, digits, na = "NA")
  inp <- x$inputs
  cat_wrapped(sprintf("<cm_threshold_set> %s; %d input%s%s", threshold_question(x), nrow(inp),
                      if (nrow(inp) == 1L) "" else "s",
                      if (isTRUE(x$auto)) ", chosen by their effect in cm_sensitivity()" else ""))
  cat("\n")
  tab <- data.frame(input = inp$label, now = f(inp$baseline), from = f(inp$low), to = f(inp$high),
                    stringsAsFactors = FALSE)
  tab$now[!is.finite(inp$baseline)] <- "-"
  if (isTRUE(x$auto)) tab$effect <- f(inp$effect)
  tab$thresholds <- ifelse(is.na(inp$n_thresholds), "-", inp$n_thresholds)
  rest <- sum(vapply(names(tab)[-1], function(nm) max(nchar(c(nm, tab[[nm]]))) + 1L, numeric(1)))
  room <- max(16L, getOption("width", 80L) - rest - 2L)
  long <- nchar(tab$input) > room
  tab$input[long] <- paste0(substr(tab$input[long], 1L, room - 3L), "...")
  print_table(tab, left = "input")
  th <- x$thresholds
  cat("\n")
  if (nrow(th)) {
    cat("Crossings:\n")
    for (i in seq_len(nrow(th))) {
      lo <- f(th$lower[i])
      hi <- f(th$upper[i])
      at <- if (identical(th$status[i], "threshold")) sprintf("at %s", f(th$threshold[i])) else
        sprintf("%s between %s and %s", th$status[i], lo, hi)
      item <- if (x$conclusion == "total") "" else paste0(th$item[i], ", ")
      cat(strwrap(sprintf("%s: %s%s; %s", th$label[i], item, at, th$description[i]),
                  width = getOption("width", 80L), indent = 2L, exdent = 4L), sep = "\n")
    }
  }
  none <- inp$result != "threshold(s) found"
  if (any(none)) {
    cat(if (nrow(th)) "\n" else "", "No threshold:\n", sep = "")
    for (i in which(none)) {
      cat(strwrap(sprintf("%s: %s", inp$label[i], inp$result[i]), width = getOption("width", 80L),
                  indent = 2L, exdent = 4L), sep = "\n")
    }
  }
  cat_wrapped("now: the input's value in the model; from, to: the range searched. Each search is in $searches.",
              exdent = 0L)
  invisible(x)
}

# Contiguous stretches of unusable grid points.
threshold_regions <- function(xs, status) {
  bad <- status != "ok"
  if (!any(bad)) {
    return(data.frame(from = numeric(0), to = numeric(0), status = character(0),
                      n_points = integer(0), stringsAsFactors = FALSE))
  }
  r <- rle(paste(bad, status))
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  keep <- bad[starts]
  data.frame(from = xs[starts[keep]], to = xs[ends[keep]], status = status[starts[keep]],
             n_points = r$lengths[keep], stringsAsFactors = FALSE)
}

# A function that sets the input to a value (returning the model and the
# settings), its baseline value, its type and whether it is a ratio (log
# grid).
threshold_setter <- function(model, input, st) {
  parts <- strsplit(input, ":", fixed = TRUE)[[1]]
  type <- if (length(parts)) parts[1] else ""
  old <- c(assoc = "'pair:<d1>:<d2>' (or 'associations:<row>')", inter = "'interaction:<d1>;<d2>'",
           three = "'or_ratio:<d1>:<d2>|<d3>'", impact = "'impacts:<row>'")
  if (type %in% names(old)) {
    cm_abort(sprintf("'%s' is a key of earlier versions: write %s.", input, old[[type]]), class = "deconflate_invalid_input")
  }
  # Inputs that are not entered: parsed into the internal types.
  if (type == "pair") {
    type <- parts[1] <- "assoc"
  } else if (type == "interaction") {
    if (length(parts) != 2L) cm_abort(sprintf("`input` '%s' must be 'interaction:<d1>;<d2>'.", input), class = "deconflate_invalid_input")
    parts <- c("inter", strsplit(parts[2], ";", fixed = TRUE)[[1]])
    type <- "inter"
  } else if (type == "or_ratio") {
    parts <- c("three", strsplit(sub("^or_ratio:", "", input), "[:|]")[[1]])
    type <- "three"
  } else if (!(type %in% c("prob", "risk", "associations", "impacts"))) {
    cm_abort(sprintf("Unknown input '%s' (use prob:, associations:, impacts:, risk, pair:, interaction: or or_ratio:).", input),
             class = "deconflate_invalid_input")
  }
  ids <- model$diseases$disease
  need <- function(k) {
    if (length(parts) != k + 1L) cm_abort(sprintf("`input` '%s' has the wrong number of parts.", input), class = "deconflate_invalid_input")
    check_disease_ids(parts[-1], model)
    if (anyDuplicated(parts[-1])) cm_abort(sprintf("`input` '%s' names the same disease twice.", input), class = "deconflate_invalid_input")
  }
  with_st <- function(m) list(model = m, st = st)
  by_row <- function(key) {
    ck <- canonical_key(key, model)
    p2 <- strsplit(ck, ":", fixed = TRUE)[[1]]
    tab <- model[[p2[1]]]
    r <- match(as.integer(p2[2]), tab$row)
    meas <- tab$measure[r]
    list(type = if (p2[1] == "impacts") "impact" else "assoc", baseline = tab$value[r], baseline_set = tab$value[r],
         ratio_scale = !is.na(meas) && meas %in% c("OR", "RR", "OR_ratio", "HR", "rate_ratio"),
         set = function(x) with_st(set_inputs(model, stats::setNames(x, ck))))
  }
  switch(type,
    associations = , impacts = by_row(input),
    assoc = {
      need(2)
      hit <- pair_row_hits(model$associations, parts[2], parts[3], model$diseases$disease)
      if (sum(hit) > 1L) {
        cm_abort(sprintf("The pair %s:%s has %d crude pair rows (rows %s); vary one of them, e.g. 'associations:%d'.",
                         parts[2], parts[3], sum(hit), paste(model$associations$row[hit], collapse = ", "),
                         model$associations$row[hit][1]), class = "deconflate_invalid_input")
      }
      if (sum(hit) == 1L) return(by_row(paste0("associations:", model$associations$row[hit])))
      # A pair without a crude pair row is varied as an odds ratio.
      list(type = type, baseline = NA_real_, baseline_set = NULL, ratio_scale = TRUE,
           set = function(x) {
             if (is.null(x)) return(with_st(model))
             with_st(set_association(model, parts[2], parts[3], x, measure = "OR"))
           })
    },
    inter = {
      need(2)
      # Additive impacts: an interaction in the impacts' units (0: none);
      # event impacts: a hazard ratio multiplier (1: none).
      ev <- isTRUE(st$event_model)
      im <- model$impacts
      hit <- same_set(im$disease, paste(parts[2:3], collapse = ";")) & is_all(im$adjusted_for) & is.na(im$among)
      base <- if (any(hit)) im$value[hit][1] else if (ev) 1 else 0
      list(type = type, baseline = base, baseline_set = if (any(hit)) base else NULL, ratio_scale = ev,
           set = function(x) {
             if (is.null(x)) return(with_st(model))
             with_st(set_interaction(model, parts[2], parts[3], x))
           })
    },
    three = {
      need(3)
      a <- model$associations
      hit <- if (is.null(a)) logical(0) else a$measure == "OR_ratio" & vapply(seq_len(nrow(a)), function(r) {
        setequal(c(a$disease[r], a$associated_with[r], sub("^no\\s+", "", a$among[r])), parts[2:4])
      }, logical(1))
      # (A row among "no S" holds the inverse ratio.)
      base <- if (any(hit)) {
        h1 <- which(hit)[1]
        if (grepl("^no\\s+", a$among[h1])) 1 / a$value[h1] else a$value[h1]
      } else NA_real_
      list(type = type, baseline = base, baseline_set = if (any(hit)) base else NULL, ratio_scale = TRUE,
           set = function(x) {
             if (is.null(x)) return(with_st(model))
             with_st(set_three_way(model, parts[2], parts[3], parts[4], x))
           })
    },
    prob = {
      need(1)
      b <- model$diseases$value[match(parts[2], ids)]
      list(type = type, baseline = b, baseline_set = b, ratio_scale = FALSE,
           set = function(x) with_st(set_inputs(model, stats::setNames(x, input))))
    },
    risk = {
      if (length(parts) != 1L) cm_abort("`input` 'risk' has no further parts.", class = "deconflate_invalid_input")
      if (!isTRUE(st$event_model)) cm_abort("The overall or baseline risk can be varied only with event_model = TRUE.", class = "deconflate_unsupported")
      list(type = type, baseline = st$risk, baseline_set = st$risk, ratio_scale = FALSE,
           set = function(x) {
             if (!(x > 0 && x < 1)) {
               cm_abort(sprintf("A%s %s of %g is not a proportion.", if (identical(st$risk_kind, "baseline")) "" else "n",
                                risk_label(st$risk_kind), x), class = "deconflate_infeasible")
             }
             st2 <- st
             st2$risk <- x
             list(model = model, st = st2)
           })
    }
  )
}

# Grid brackets of the sign changes of `f` (rows: first and last index),
# using only usable points (`ok`). An exact zero, or a stretch of zeros, is a
# crossing only if the usable points on either side have opposite signs (it
# only touches zero otherwise); its bracket runs from the point before the
# zeros to the first non-zero point after them.
threshold_brackets <- function(f, ok) {
  n <- length(f)
  out <- list()
  for (g in seq_len(n - 1L)) {
    if (!(ok[g] && ok[g + 1L]) || f[g] == 0) next
    j <- g + 1L
    while (j <= n && ok[j] && f[j] == 0) j <- j + 1L
    if (j > n || !ok[j]) next
    if (sign(f[j]) != sign(f[g])) out[[length(out) + 1L]] <- c(g, j)
  }
  if (!length(out)) return(matrix(integer(0), 0L, 2L))
  do.call(rbind, out)
}

#' @export
print.cm_threshold <- function(x, digits = 3, ...) {
  f <- function(v) fmt_num(v, digits)
  cat_wrapped(c(sprintf("<cm_threshold> %s, as %s varies from %s to %s (now %s; method %s)",
                        threshold_question(x), x$input, format(signif(x$range[1], 3)), format(signif(x$range[2], 3)),
                        if (is.null(x$baseline) || is.na(x$baseline)) "not entered" else format(signif(x$baseline, 3)),
                        x$method),
                sprintf("Fixed: %s.", x$fixed)))
  n_bad <- sum(x$scan$status != "ok")
  if (n_bad) {
    cat(sprintf("%d of %d grid points could not be used:\n", n_bad, nrow(x$scan)))
    print_table(x$regions, digits = digits)
  }
  th <- x$thresholds
  if (!nrow(th)) {
    cat(if (n_bad) "\nNo crossing found where the range could be evaluated.\n" else "\nNo crossing found in the range.\n")
  } else {
    cat("\nCrossings:\n")
    for (i in seq_len(nrow(th))) {
      lo <- f(th$lower[i])
      hi <- f(th$upper[i])
      at <- if (identical(th$status[i], "threshold")) {
        sprintf("at %s%s", f(th$threshold[i]), if (lo != hi) sprintf(" (between %s and %s)", lo, hi) else "")
      } else {
        sprintf("%s between %s and %s", th$status[i], lo, hi)
      }
      cat(strwrap(sprintf("%s: %s; %s", th$item[i], at, th$description[i]), width = getOption("width", 80L),
                  indent = 2L, exdent = 4L), sep = "\n")
    }
  }
  invisible(x)
}

