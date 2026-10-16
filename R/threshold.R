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
#' * `"total"`: the total crosses `target`;
#' * `"change"`: the total departs from its value at the model's own input
#'   by the relative amount `target` (e.g. `0.1` for +10%), which answers
#'   questions such as "what association strength would increase the total
#'   by 10%?".
#'
#' The total is the adjusted aggregate (in the units of the impacts) for
#' additive impacts, and the risk attributable to disease for event impacts
#' (`event_model = TRUE`), whose contributions are its Shapley allocation.
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
#' Thresholds are deterministic: they are for the model's input values. To
#' see first which inputs matter, use [cm_sensitivity()].
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
#'   association of diseases 1 and 2).
#' @param range Two numbers: the range of the input to search.
#' @param conclusion `"rank"`, `"sign"`, `"total"` or `"change"`.
#' @param target For `"total"`, the value of the aggregate (for event
#'   impacts, of the attributable risk); for `"change"`, the relative change
#'   from the baseline (e.g. `0.1`, or `-0.1`).
#' @param diseases Optional disease ids to restrict `"rank"` (pairs among
#'   them) and `"sign"`.
#' @param by For `"rank"`: `"contribution"` (default) or `"adjusted"`.
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
#' @return A `cm_threshold` object: `thresholds` (one row per crossing:
#'   conclusion, item, status, threshold, bracket `lower`/`upper`, the
#'   compared quantity just below and above, and a description), `details`
#'   (for each crossing, the adjusted results at the two ends of its
#'   bracket), `scan` (each grid point with its status, the aggregate and the
#'   condition number), `regions` (stretches of unusable grid points),
#'   `summary` (per item: number of thresholds, and whether part of the range
#'   could not be evaluated), and the settings (`input`, `range`, `baseline`
#'   value of the input, `method` (the computations used), `conclusion`,
#'   `target`, `fixed`).
#' @export
#' @examples
#' m <- example_supplement()
#' # At what odds ratio between d1 and d2 does the ranking change?
#' th <- cm_threshold(m, "pair:d1:d2", c(0.2, 20), conclusion = "rank")
#' th
#' # What odds ratio between d1 and d3 (an odds ratio of 1 in the example)
#' # would reduce the aggregate by 10%?
#' cm_threshold(m, "pair:d1:d3", c(1, 50), conclusion = "change", target = -0.1)$thresholds
#' # From a result: its settings are used
#' res <- deconflate(m, n_draws = 0)
#' cm_threshold(res, "interaction:d2;d3", c(-2, 2), conclusion = "rank")$summary
cm_threshold <- function(model, input, range, conclusion = c("rank", "sign", "total", "change"),
                         target = NULL, diseases = NULL, by = c("contribution", "adjusted"),
                         event_model = FALSE, overall_risk = NULL,
                         n_grid = 101L, log_scale = NULL, tol = 1e-8, baseline_risk = NULL, ...) {
  dots <- list(...)
  retired_args(dots)
  if (inherits(model, "cm_result") || inherits(model, "cm_event_result")) {
    if (!isFALSE(event_model) || !is.null(overall_risk) || !is.null(baseline_risk)) {
      cm_abort("With a deconflate() result, the settings come from the result: leave out event_model, overall_risk and baseline_risk.",
               class = "deconflate_invalid_input")
    }
    ctx <- result_context(model)
    st <- list(event_model = ctx$event, risk = ctx$risk$value, risk_kind = ctx$risk$kind %||% "overall",
               pool = ctx$plan$pool %||% "random")
    dots <- utils::modifyList(ctx$dots, dots)
    model <- ctx$model
    event_model <- ctx$event
  } else {
    st <- sens_setup(model, event_model, overall_risk, baseline_risk = baseline_risk)
    st$pool <- "random"
  }
  conclusion <- cm_match_arg(conclusion)
  by <- cm_match_arg(by)
  if (!is.character(input) || length(input) != 1L || is.na(input)) {
    cm_abort("`input` must be a single string such as \"pair:d1:d2\".", class = "deconflate_invalid_input")
  }
  check_numeric(range, "range")
  if (length(range) != 2L || range[1] >= range[2]) cm_abort("`range` must be two increasing numbers.", class = "deconflate_invalid_input")
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
  if (conclusion %in% c("total", "change")) {
    if (is.null(target) || !is.numeric(target) || length(target) != 1L || !is.finite(target)) {
      cm_abort(sprintf("conclusion = '%s' needs a single finite `target`.", conclusion), class = "deconflate_invalid_input")
    }
  }
  ids <- model$diseases$disease
  if (!is.null(diseases)) {
    if (!is.character(diseases) || !length(diseases)) cm_abort("`diseases` must be a character vector of disease ids.", class = "deconflate_invalid_input")
    check_disease_ids(diseases, model)
    diseases <- unique(diseases)
  }
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
  # changes.
  base_eval <- evaluate_model(model, st)
  base_total <- if (identical(base_eval$status, "ok")) total_of(base_eval$res) else NA_real_
  if (conclusion == "change" && !is.finite(base_total)) {
    # (The baseline is the model as given, with the input at its own value.)
    cm_abort(sprintf("The baseline result is not usable (%s), so a relative change cannot be computed.",
                     base_eval$status),
             class = type_class(switch(base_eval$status, unresolved = "nonconvergence", undefined = "nonfinite",
                                       base_eval$status)))
  }
  if (conclusion == "change" && base_total == 0) {
    cm_abort("The baseline aggregate is zero, so a relative change is undefined; use conclusion = 'total'.", class = "deconflate_nonfinite")
  }
  items <- threshold_items(conclusion, ids, diseases, target, base_total)
  metric <- function(res) {
    event <- inherits(res, "cm_event_result")
    switch(conclusion,
      total = , change = stats::setNames(total_of(res) - items$level, items$item),
      sign = {
        a <- res$adjusted
        v <- a$adjusted[match(items$a, a$disease)]
        stats::setNames(if (event) log(v) else v, items$item)
      },
      rank = {
        sc <- if (by == "contribution") burden_metric(res)$by_disease else
          stats::setNames(res$adjusted$adjusted, res$adjusted$disease)
        stats::setNames(unname(sc[items$a] - sc[items$b]), items$item)
      })
  }

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
      st <- "threshold"
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
            st <- if (identical(e$status, "singular")) "discontinuity" else "unresolved"
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
      if (st == "threshold" && flo != 0 && fhi != 0) {
        # Continuity: a genuine crossing has a small compared quantity on
        # both sides; a pole has large values of opposite sign.
        eps <- 1e-4 * scale_f[k]
        kap <- max(rlo$diagnostics$condition_number, rhi$diagnostics$condition_number)
        if (max(abs(flo), abs(fhi)) > eps || !is.finite(kap) || kap > 1e12) st <- "discontinuity"
      }
      x0 <- if (zero_run) xs[g + 1L] else if (flo == 0) lo else if (fhi == 0) hi else
        if (log_scale) sqrt(lo * hi) else (lo + hi) / 2
      rows[[length(rows) + 1L]] <- data.frame(
        conclusion = conclusion, item = items$item[k], status = st,
        threshold = if (st == "threshold") x0 else NA_real_,
        lower = lo, upper = hi, below = flo, above = fhi,
        description = threshold_description(conclusion, items[k, ], flo, fhi, by, st, target,
                                            isTRUE(base_plan$event)),
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
                 summary = summary, input = input, range = range, baseline = setter$baseline,
                 method = paste(if (length(used)) used else "none", collapse = ", "),
                 event_model = event_model, conclusion = conclusion, target = target, by = by,
                 log_scale = log_scale,
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

# Items compared for a conclusion, with what they compare.
threshold_items <- function(conclusion, ids, diseases, target, base_total) {
  use <- diseases %||% ids
  switch(conclusion,
    total = data.frame(item = "total", a = NA_character_, b = NA_character_, level = target,
                       stringsAsFactors = FALSE),
    change = data.frame(item = "total", a = NA_character_, b = NA_character_,
                        level = base_total + target * abs(base_total), stringsAsFactors = FALSE),
    sign = data.frame(item = use, a = use, b = NA_character_, level = 0, stringsAsFactors = FALSE),
    rank = {
      if (length(use) < 2L) cm_abort("conclusion = 'rank' needs at least two diseases.", class = "deconflate_invalid_input")
      pr <- utils::combn(use, 2)
      data.frame(item = paste(pr[1, ], pr[2, ], sep = " vs "), a = pr[1, ], b = pr[2, ],
                 level = 0, stringsAsFactors = FALSE)
    })
}

threshold_description <- function(conclusion, item, flo, fhi, by, status, target = NULL,
                                  event = FALSE) {
  if (status != "threshold") {
    return(switch(status,
      discontinuity = "sign change across a pole or a nearly singular system (not a threshold)",
      unresolved = "the crossing could not be refined: an unusable point lies inside the bracket"))
  }
  up <- fhi > flo
  what <- if (event) "the attributable risk" else "the aggregate"
  switch(conclusion,
    total = sprintf("%s %s the target", what, if (up) "rises above" else "falls below"),
    change = sprintf("%s %s the baseline %s %g%%", what, if (up) "rises above" else "falls below",
                     if (isTRUE(target >= 0)) "+" else "-", abs(100 * target)),
    sign = sprintf("the adjusted value of %s moves %s no effect (a model-implied sign change)",
                   item$a, if (up) "above" else "below"),
    rank = sprintf("%s %s %s (by %s)", item$a, if (up) "moves above" else "moves below", item$b, by))
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
  cat_wrapped(c(sprintf("<cm_threshold> %s over %s in [%s, %s] (baseline %s; method %s)", x$conclusion, x$input,
                        sprintf("%g", x$range[1]), sprintf("%g", x$range[2]),
                        if (is.null(x$baseline) || is.na(x$baseline)) "not entered" else sprintf("%g", x$baseline), x$method),
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

#' @rdname plots
#' @param items For `cm_threshold` plots: the items to show (default all).
#' @export
plot.cm_threshold <- function(x, items = NULL, ...) {
  v <- x$values
  items <- items %||% colnames(v)
  unk <- setdiff(items, colnames(v))
  if (length(unk)) {
    cm_abort(sprintf("Unknown items: %s (available: %s).", paste(unk, collapse = ", "),
                     paste(colnames(v), collapse = ", ")), class = "deconflate_invalid_input")
  }
  v <- v[, items, drop = FALSE]
  xs <- x$scan$value
  cols <- grDevices::hcl.colors(max(2L, ncol(v)), "Dark 3")[seq_len(ncol(v))]
  ylim <- range(v, 0, finite = TRUE)
  graphics::matplot(xs, v, type = "l", lty = 1, col = cols, log = if (isTRUE(x$log_scale)) "x" else "",
                    xlab = x$input, ylab = switch(x$conclusion, rank = "difference in score",
                                                  sign = if (isTRUE(x$event_model)) "log adjusted hazard ratio" else "adjusted impact",
                                                  "total minus target"),
                    main = sprintf("Threshold search: %s", x$conclusion), ylim = ylim, ...)
  graphics::abline(h = 0, col = "grey50")
  if (isTRUE(is.finite(x$baseline)) && (!isTRUE(x$log_scale) || x$baseline > 0)) {
    graphics::abline(v = x$baseline, lty = 3, col = "grey40")
  }
  bad <- x$scan$status != "ok"
  if (any(bad)) graphics::rug(xs[bad], col = "#CC6677", lwd = 2)
  th <- x$thresholds[x$thresholds$status == "threshold" & x$thresholds$item %in% items, , drop = FALSE]
  if (nrow(th)) graphics::abline(v = th$threshold, lty = 2, col = cols[match(th$item, items)])
  if (ncol(v) <= 8) {
    graphics::legend("topright", legend = items, col = cols, lty = 1, bty = "n", cex = 0.8)
  }
  invisible(x)
}
