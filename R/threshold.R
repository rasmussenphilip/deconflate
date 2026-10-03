#' Threshold search: where does a conclusion change?
#'
#' Varies one input over a range, with everything else fixed at the model's
#' values, and finds the input values at which a conclusion changes:
#' * `"rank"`: two diseases swap places, ranked by their contribution to the
#'   aggregate (including interaction shares; `by = "adjusted"` ranks by
#'   adjusted impact instead);
#' * `"sign"`: a disease's adjusted impact crosses zero. This is a sign
#'   change implied by the model and its inputs, not by itself evidence of a
#'   protective effect;
#' * `"total"`: the aggregate adjusted burden crosses `target` (in the units
#'   of the impacts);
#' * `"change"`: the aggregate departs from its value at the model's own
#'   input by the relative amount `target` (e.g. `0.1` for +10%), which
#'   answers questions such as "what association strength would increase the
#'   total by 10%?".
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
#' (continuity); a jump across a pole of the published approximation, or
#' across a nearly singular system, is reported as a `"discontinuity"`, never
#' as a threshold. If a bisection step lands on an unusable point, the
#' crossing is `"unresolved"`. All crossings in the range are reported;
#' stretches of unusable grid points are listed in `regions`, so that "no
#' crossing found" can be told apart from "part of the range could not be
#' evaluated". A sign change between usable points separated by unusable
#' ones is not reported as a crossing; such stretches appear in `regions`.
#'
#' With the global method, when the input is an impact or an interaction
#' (which do not change the joint distribution), the joint distribution is
#' fitted once and reused (or the `joint` passed through `...` is used).
#'
#' Thresholds are deterministic: they are for the model's input values.
#' Probabilistic statements (e.g. the probability that a disease ranks
#' first) need the Monte Carlo tools.
#'
#' @param model A [cm_model()].
#' @param input The input to vary, keyed as in [cm_sampler()]:
#'   `"assoc:<d1>:<d2>"` (the pair's association; a pair without a numeric
#'   measure is varied as an odds ratio, so unknown or unlisted pairs can be
#'   explored), `"inter:<d1>:<d2>"` (an interaction), `"prob:<disease>"` (the
#'   disease's `value`, on the scale it was entered), `"impact:<disease>"`
#'   (a raw impact), or `"three:<d1>:<d2>:<d3>"` (a three-way ratio).
#' @param range Two numbers: the range of the input to search.
#' @param conclusion `"rank"`, `"sign"`, `"total"` or `"change"`.
#' @param target For `"total"`, the value of the aggregate; for `"change"`,
#'   the relative change from the baseline aggregate (e.g. `0.1`, or `-0.1`).
#' @param diseases Optional disease ids to restrict `"rank"` (pairs among
#'   them) and `"sign"`.
#' @param by For `"rank"`: `"contribution"` (default) or `"adjusted"`.
#' @param method Adjustment method; by default `"global"` for models with
#'   interactions or three-way terms (and for `"inter:"` and `"three:"`
#'   inputs), otherwise `"simultaneous"`. The pairwise methods are rejected
#'   (class `deconflate_unsupported`) for `"inter:"` and `"three:"` inputs and
#'   for models with interactions.
#' @param n_grid Number of grid points.
#' @param log_scale Use a log-spaced grid? Default `TRUE` for associations
#'   measured as ratios and three-way ratios with a positive range.
#' @param tol Relative tolerance on the threshold. (The bracket is also
#'   narrowed to at most 1e-6 of the grid spacing, so that the continuity
#'   test is reliable.) Because of this argument, the `tol` of [fit_joint()]
#'   cannot be passed through `...`.
#' @param ... Passed to [deconflate()] (e.g. `feasibility`, `joint`, or
#'   arguments of [fit_joint()] for the global method).
#' @return A `cm_threshold` object: `thresholds` (one row per crossing:
#'   conclusion, item, status, threshold, bracket `lower`/`upper`, the
#'   compared quantity just below and above, and a description), `details`
#'   (for each crossing, the adjusted results at the two ends of its
#'   bracket), `scan` (each grid point with its status, the aggregate and the
#'   condition number), `regions` (stretches of unusable grid points),
#'   `summary` (per item: number of thresholds, and whether part of the range
#'   could not be evaluated), and the settings (`input`, `range`, `baseline`
#'   value of the input, `method`, `conclusion`, `target`, `fixed`).
#' @export
#' @examples
#' m <- example_supplement()
#' # At what odds ratio between d1 and d2 does the ranking change?
#' th <- cm_threshold(m, "assoc:d1:d2", c(0.2, 20), conclusion = "rank")
#' th
#' # What odds ratio between d1 and d3 (independent in the example) would
#' # reduce the aggregate by 10%?
#' cm_threshold(m, "assoc:d1:d3", c(1, 50), conclusion = "change", target = -0.1)$thresholds
cm_threshold <- function(model, input, range, conclusion = c("rank", "sign", "total", "change"),
                         target = NULL, diseases = NULL, by = c("contribution", "adjusted"),
                         method = NULL, n_grid = 101L, log_scale = NULL, tol = 1e-8, ...) {
  check_model(model)
  conclusion <- match.arg(conclusion)
  by <- match.arg(by)
  if (!is.character(input) || length(input) != 1L || is.na(input)) {
    cm_abort("`input` must be a single string such as \"assoc:d1:d2\".")
  }
  check_numeric(range, "range")
  if (length(range) != 2L || range[1] >= range[2]) cm_abort("`range` must be two increasing numbers.")
  if (!is.numeric(n_grid) || length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 3) {
    cm_abort("`n_grid` must be a single number, at least 3.")
  }
  n_grid <- as.integer(n_grid)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    cm_abort("`tol` must be a single positive number.")
  }
  if (!is.null(log_scale) && (!is.logical(log_scale) || length(log_scale) != 1L || is.na(log_scale))) {
    cm_abort("`log_scale` must be NULL, TRUE or FALSE.")
  }
  if (conclusion %in% c("total", "change")) {
    if (is.null(target) || !is.numeric(target) || length(target) != 1L || !is.finite(target)) {
      cm_abort(sprintf("conclusion = '%s' needs a single finite `target`.", conclusion))
    }
  }
  if (is.null(model$impacts)) cm_abort("The model has no impacts to adjust.")
  ids <- model$diseases$id
  if (!is.null(diseases)) {
    if (!is.character(diseases) || !length(diseases)) cm_abort("`diseases` must be a character vector of disease ids.")
    check_disease_ids(diseases, model)
    diseases <- unique(diseases)
  }
  setter <- threshold_setter(model, input)
  has_inter <- !is.null(model$interactions) && nrow(model$interactions) > 0L
  if (is.null(method)) {
    needs_global <- has_inter || (!is.null(model$three_way) && nrow(model$three_way) > 0L) ||
      setter$type %in% c("inter", "three")
    method <- if (needs_global) "global" else "simultaneous"
  }
  method <- match.arg(method, c("simultaneous", "published", "global"))
  if (method != "global" && (has_inter || setter$type %in% c("inter", "three"))) {
    # The pairwise methods reject interactions and ignore three-way terms,
    # so the search would fail at every point or find nothing.
    cm_abort(sprintf("%s needs method = 'global'.",
                     if (setter$type %in% c("inter", "three")) sprintf("Varying '%s'", input) else "A model with interactions"),
             class = "deconflate_unsupported")
  }
  log_scale <- log_scale %||% (setter$ratio_scale && range[1] > 0)
  if (log_scale && range[1] <= 0) cm_abort("A log-spaced grid needs a positive range.")
  xs <- if (log_scale) exp(seq(log(range[1]), log(range[2]), length.out = n_grid)) else
    seq(range[1], range[2], length.out = n_grid)
  xs[c(1L, n_grid)] <- range   # exact ends (exp(log(x)) can differ from x)

  # Arguments for deconflate(). Impacts and interactions do not change the
  # joint distribution, so for those inputs the global method fits it once.
  dots <- list(...)
  joint <- dots[["joint"]]
  dots[["joint"]] <- NULL
  if (!is.null(joint) && !(setter$type %in% c("impact", "inter"))) {
    cm_abort("A fixed `joint` can be used only when varying an impact or an interaction; other inputs change the joint distribution.")
  }
  if (is.null(joint) && method == "global" && setter$type %in% c("impact", "inter")) {
    fj <- dots[intersect(names(dots), joint_arg_names)]
    joint <- tryCatch(
      withCallingHandlers(do.call(fit_joint, c(list(model), fj)),
                          deconflate_nonconvergence = function(w) {
                            if (inherits(w, "warning")) invokeRestart("muffleWarning")
                          }),
      deconflate_error = function(e) NULL)   # each point then reports its own failure
  }

  evaluate <- function(x) {
    m2 <- tryCatch(setter$set(x), deconflate_error = function(e) e)
    if (inherits(m2, "condition")) {
      # The setter fails only for an invalid value of the input itself.
      st <- threshold_status(m2)
      return(list(status = if (st == "error") "infeasible" else st, res = NULL))
    }
    evaluate_model(m2)
  }
  evaluate_model <- function(m2) {
    res <- tryCatch(
      withCallingHandlers(do.call(deconflate, c(list(m2, method = method, joint = joint, warn = FALSE), dots)),
                          deconflate_nonconvergence = function(w) {
                            if (inherits(w, "warning")) invokeRestart("muffleWarning")
                          }),
      deconflate_error = function(e) e)
    if (inherits(res, "condition")) return(list(status = threshold_status(res), res = NULL))
    if (!result_is_finite(res)) return(list(status = "undefined", res = res))
    list(status = "ok", res = res)
  }

  # The compared quantities: each must cross zero where the conclusion
  # changes.
  base_eval <- evaluate_model(model)
  base_total <- if (identical(base_eval$status, "ok")) base_eval$res$totals$adjusted_total else NA_real_
  if (conclusion == "change" && !is.finite(base_total)) {
    # (The baseline is the model as given, with the input at its own value.)
    cm_abort(sprintf("The baseline result is not usable (%s), so a relative change cannot be computed. For a model with unknown pairs, use method = \"global\".",
                     base_eval$status))
  }
  if (conclusion == "change" && base_total == 0) {
    cm_abort("The baseline aggregate is zero, so a relative change is undefined; use conclusion = 'total'.")
  }
  items <- threshold_items(conclusion, ids, diseases, target, base_total)
  metric <- function(res) {
    switch(conclusion,
      total = , change = stats::setNames(res$totals$adjusted_total - items$level, items$item),
      sign = {
        a <- res$adjusted
        stats::setNames(a$adjusted[match(items$a, a$disease)], items$item)
      },
      rank = {
        sc <- if (by == "contribution") res$contributions$total else res$adjusted$adjusted
        d <- if (by == "contribution") res$contributions$disease else res$adjusted$disease
        stats::setNames(sc[match(items$a, d)] - sc[match(items$b, d)], items$item)
      })
  }

  evals <- lapply(xs, evaluate)
  status <- vapply(evals, `[[`, character(1), "status")
  ok <- status == "ok"
  vals <- matrix(NA_real_, n_grid, nrow(items), dimnames = list(NULL, items$item))
  for (g in which(ok)) vals[g, ] <- metric(evals[[g]]$res)
  scan <- data.frame(
    value = xs, status = status,
    total = vapply(evals, function(e) if (identical(e$status, "ok")) e$res$totals$adjusted_total else NA_real_,
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
            st <- "unresolved"
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
        description = threshold_description(conclusion, items[k, ], flo, fhi, by, st),
        stringsAsFactors = FALSE)
      details[[length(details) + 1L]] <- list(below = rlo$adjusted, above = rhi$adjusted,
                                              totals_below = rlo$totals, totals_above = rhi$totals)
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
                 method = method, conclusion = conclusion, target = target, by = by,
                 log_scale = log_scale,
                 fixed = sprintf("all inputs other than %s at the model's values", input),
                 values = vals),
            class = "cm_threshold")
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
      if (length(use) < 2L) cm_abort("conclusion = 'rank' needs at least two diseases.")
      pr <- utils::combn(use, 2)
      data.frame(item = paste(pr[1, ], pr[2, ], sep = " vs "), a = pr[1, ], b = pr[2, ],
                 level = 0, stringsAsFactors = FALSE)
    })
}

threshold_description <- function(conclusion, item, flo, fhi, by, status) {
  if (status != "threshold") {
    return(switch(status,
      discontinuity = "sign change across a pole or a nearly singular system (not a threshold)",
      unresolved = "the crossing could not be refined: an unusable point lies inside the bracket"))
  }
  up <- fhi > flo
  switch(conclusion,
    total = sprintf("the aggregate %s the target", if (up) "rises above" else "falls below"),
    change = sprintf("the aggregate %s the baseline by the target change", if (up) "rises above" else "falls below"),
    sign = sprintf("the adjusted impact of %s becomes %s (a model-implied sign change)",
                   item$a, if (up) "positive" else "negative"),
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

# A function that sets the input to a value, its baseline value, its type and
# whether it is a ratio (log grid).
threshold_setter <- function(model, input) {
  parts <- strsplit(input, ":", fixed = TRUE)[[1]]
  type <- if (length(parts)) parts[1] else ""
  if (!(type %in% c("assoc", "inter", "three", "prob", "impact"))) {
    cm_abort(sprintf("Unknown input '%s' (use assoc:, inter:, prob:, impact: or three:).", input))
  }
  ids <- model$diseases$id
  need <- function(k) {
    if (length(parts) != k + 1L) cm_abort(sprintf("`input` '%s' has the wrong number of parts.", input))
    check_disease_ids(parts[-1], model)
    if (anyDuplicated(parts[-1])) cm_abort(sprintf("`input` '%s' names the same disease twice.", input))
  }
  # Disease probabilities and impacts are set through a sampler with a fixed
  # value (which validates the value, e.g. a probability in (0, 1)); passing
  # the value directly draws no random numbers.
  fixed_setter <- function(key, sampler) {
    function(x) sampler(1L, values = stats::setNames(x, key))
  }
  switch(type,
    assoc = {
      need(2)
      a <- model$associations
      hit <- if (is.null(a)) NA_integer_ else match(pair_key(parts[2], parts[3]), pair_key(a$disease1, a$disease2))
      numeric_measure <- !is.na(hit) && a$measure[hit] %in% c("OR", "RR", "RD", "cond_prob", "phi")
      meas <- if (numeric_measure) a$measure[hit] else "OR"
      d1 <- if (is.na(hit)) parts[2] else a$disease1[hit]
      d2 <- if (is.na(hit)) parts[3] else a$disease2[hit]
      # The model's own value of the input: its measure, the odds ratio of a
      # table, 1 for an independent or unlisted pair (when missing pairs are
      # independent), and NA for an unknown pair.
      base <- if (numeric_measure || (!is.na(hit) && a$measure[hit] == "table")) {
        a$value[hit]
      } else if ((!is.na(hit) && a$measure[hit] == "unknown") ||
                 (is.na(hit) && identical(model$missing_associations, "unknown"))) {
        NA_real_
      } else 1
      list(type = type, baseline = base, ratio_scale = meas %in% c("OR", "RR"),
           set = function(x) set_association(model, d1, d2, x, measure = meas))
    },
    inter = {
      need(2)
      it <- model$interactions
      hit <- if (is.null(it)) NA_integer_ else match(pair_key(parts[2], parts[3]), pair_key(it$disease1, it$disease2))
      list(type = type, baseline = if (is.na(hit)) 0 else it$value[hit], ratio_scale = FALSE,
           set = function(x) set_interaction(model, parts[2], parts[3], x))
    },
    three = {
      need(3)
      tw <- model$three_way
      key <- paste(sort(parts[2:4]), collapse = "|")
      hit <- if (is.null(tw) || !nrow(tw)) NA_integer_ else
        match(key, apply(tw[, c("disease1", "disease2", "disease3")], 1, function(x) paste(sort(x), collapse = "|")))
      list(type = type, baseline = if (is.na(hit)) 1 else tw$ratio[hit], ratio_scale = TRUE,
           set = function(x) set_three_way(model, parts[2], parts[3], parts[4], x))
    },
    prob = {
      need(1)
      base <- model$diseases$value[match(parts[2], ids)]
      s <- cm_sampler(model, diseases = stats::setNames(list(dist_fixed(base)), parts[2]))
      list(type = type, baseline = base, ratio_scale = FALSE,
           set = fixed_setter(paste0("prob:", parts[2]), s))
    },
    impact = {
      need(1)
      base <- model$impacts$value[match(parts[2], model$impacts$disease)]
      s <- cm_sampler(model, impacts = stats::setNames(list(dist_fixed(base)), parts[2]))
      list(type = type, baseline = base, ratio_scale = FALSE,
           set = fixed_setter(paste0("impact:", parts[2]), s))
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
print.cm_threshold <- function(x, digits = 4, ...) {
  cat(sprintf("<cm_threshold> %s over %s in [%g, %g] (baseline %g; method %s)\n",
              x$conclusion, x$input, x$range[1], x$range[2], x$baseline, x$method))
  cat(sprintf("Fixed: %s.\n", x$fixed))
  n_bad <- sum(x$scan$status != "ok")
  if (n_bad) {
    cat(sprintf("%d of %d grid points could not be used:\n", n_bad, nrow(x$scan)))
    print(x$regions, row.names = FALSE, digits = digits)
  }
  th <- x$thresholds
  if (!nrow(th)) {
    cat(if (n_bad) "\nNo crossing found where the range could be evaluated.\n" else "\nNo crossing found in the range.\n")
  } else {
    cat("\nCrossings:\n")
    print(th[, c("item", "status", "threshold", "lower", "upper", "description")],
          row.names = FALSE, digits = digits)
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
                     paste(colnames(v), collapse = ", ")))
  }
  v <- v[, items, drop = FALSE]
  xs <- x$scan$value
  cols <- grDevices::hcl.colors(max(2L, ncol(v)), "Dark 3")[seq_len(ncol(v))]
  ylim <- range(v, 0, finite = TRUE)
  graphics::matplot(xs, v, type = "l", lty = 1, col = cols, log = if (isTRUE(x$log_scale)) "x" else "",
                    xlab = x$input, ylab = switch(x$conclusion, rank = "difference in score",
                                                  sign = "adjusted impact", "aggregate minus target"),
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
