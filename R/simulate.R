#' Simulate the raw impacts that studies would report
#'
#' For validation: given a population whose disease combinations follow the
#' maximum-entropy distribution of `model` (see [fit_joint()]), and known true
#' additive impacts (and optional interactions), computes what a study would
#' estimate for each disease: the crude difference between animals with and
#' without the disease, or the coefficient of an additive regression adjusted
#' for other diseases. With `n = NULL` the exact population values are
#' returned; with `n` set, `n` animals are sampled, which adds sampling error.
#'
#' Adjusting the returned raw impacts with [deconflate()] should recover
#' `true_impacts` exactly (for `n = NULL`): with the simultaneous or global
#' method without interactions, and with the global method with them.
#'
#' @param model A [cm_population()] or [cm_model()] (its impacts are ignored).
#' @param true_impacts Named numeric vector of true impacts, one per disease
#'   (any units).
#' @param interactions Optional [cm_interactions()] with true interaction
#'   values (same units).
#' @param estimand `"crude"` or `"adjusted_linear"`, one per disease or
#'   recycled.
#' @param adjusted_for For `"adjusted_linear"`: adjustment sets (ids
#'   separated by `";"`, or `"all"`), one per disease or recycled.
#' @param n Optional number of animals to sample.
#' @param joint Optional pre-computed [fit_joint()] result.
#' @param units Optional units label for the returned impacts.
#'
#' @return A [cm_impacts()] object with the simulated raw impacts.
#' @export
#' @examples
#' raw <- simulate_raw_impacts(example_supplement(), c(d1 = 2, d2 = 4, d3 = 6))
#' raw$value
#' deconflate(cm_model(example_supplement(), raw))$adjusted$adjusted
simulate_raw_impacts <- function(model, true_impacts, interactions = NULL,
                                 estimand = "crude", adjusted_for = NA_character_,
                                 n = NULL, joint = NULL, units = NULL) {
  check_population(model)
  ids <- model$diseases$id
  if (is.null(names(true_impacts)) || !setequal(names(true_impacts), ids)) {
    cm_abort("`true_impacts` must be named, with one value per disease.")
  }
  k <- length(ids)
  estimand <- recycle_arg(estimand, k, "estimand")
  adjusted_for <- recycle_arg(adjusted_for, k, "adjusted_for")
  if (is.null(joint)) joint <- fit_joint(model) else validate_joint(joint, model)
  cells <- joint$cells
  y <- as.vector(cells %*% true_impacts[ids])
  if (!is.null(interactions)) {
    if (!inherits(interactions, "cm_interactions")) {
      cm_abort("`interactions` must be created with cm_interactions().")
    }
    unk <- setdiff(c(interactions$disease1, interactions$disease2), ids)
    if (length(unk)) cm_abort(sprintf("Interactions refer to unknown diseases: %s.", paste(unk, collapse = ", ")))
    for (r in seq_len(nrow(interactions))) {
      y <- y + interactions$value[r] * cells[, interactions$disease1[r]] *
        cells[, interactions$disease2[r]]
    }
  }
  if (is.null(n)) {
    X0 <- cells
    w <- joint$prob
    yy <- y
  } else {
    draw <- sample.int(nrow(cells), n, replace = TRUE, prob = joint$prob)
    X0 <- cells[draw, , drop = FALSE]
    w <- rep(1, n)
    yy <- y[draw]
  }
  raw <- vapply(seq_len(k), function(i) {
    S <- if (estimand[i] == "adjusted_linear") {
      match(resolve_adjusted_for(adjusted_for[i], ids[i], ids), ids)
    } else integer(0)
    X <- cbind(1, X0[, c(i, S), drop = FALSE])
    XtW <- t(X * w)
    solve(XtW %*% X, XtW %*% yy)[2]
  }, numeric(1))
  cm_impacts(ids, raw, estimand = estimand,
             adjusted_for = ifelse(estimand == "adjusted_linear", adjusted_for, NA_character_),
             source = if (is.null(n)) "simulated (expected)" else sprintf("simulated (n = %d)", n),
             units = units)
}
