#' Simulate the raw impacts a single-disease study would report
#'
#' For validation: given a population whose disease combinations follow the
#' maximum-entropy distribution implied by `model` (see [fit_joint()]), and
#' known true additive impacts (and optional interactions), computes the crude
#' difference in outcome between animals with and without each disease. With
#' `n = NULL` the exact expectation is returned. With `n` set, `n` animals
#' are sampled and the empirical crude differences are returned, which adds
#' sampling error.
#'
#' Adjusting the returned raw impacts with [deconflate()] should recover
#' `true_impacts` (exactly for `"simultaneous"`/`"global"` without
#' interactions, and for `"global"` with interactions, when `n = NULL`).
#'
#' @param model A [cm_model()] (its impacts, if any, are ignored).
#' @param true_impacts Named numeric vector of true proportional impacts, one
#'   per disease.
#' @param interactions Optional [cm_interactions()] with true interaction
#'   values (only rows for `outcome` are used).
#' @param outcome Outcome label for the returned impacts.
#' @param n Optional number of animals to sample.
#' @param joint Optional pre-computed [fit_joint()] result.
#'
#' @return A [cm_impacts()] object with the simulated raw impacts.
#' @export
#' @examples
#' m <- example_supplement()
#' raw <- simulate_raw_impacts(m, c(d1 = 0.02, d2 = 0.04, d3 = 0.06))
#' raw$value
simulate_raw_impacts <- function(model, true_impacts, interactions = NULL,
                                 outcome = "impact", n = NULL, joint = NULL) {
  check_model(model)
  ids <- model$diseases$id
  if (is.null(names(true_impacts)) || !setequal(names(true_impacts), ids)) {
    cm_abort("`true_impacts` must be named, with one value per disease.")
  }
  joint <- joint %||% fit_joint(model)
  cells <- joint$cells
  loss <- as.vector(cells %*% true_impacts[ids])
  if (!is.null(interactions)) {
    int <- interactions[interactions$outcome == outcome, , drop = FALSE]
    for (r in seq_len(nrow(int))) {
      loss <- loss + int$value[r] * cells[, int$disease1[r]] * cells[, int$disease2[r]]
    }
  }
  if (is.null(n)) {
    raw <- crude_difference(cells, joint$prob, loss)
  } else {
    draw <- sample.int(nrow(cells), n, replace = TRUE, prob = joint$prob)
    raw <- crude_difference(cells[draw, , drop = FALSE], rep(1, n), loss[draw])
  }
  cm_impacts(ids, raw, outcome = outcome, scale = "proportion",
             source = if (is.null(n)) "simulated (expected)" else sprintf("simulated (n = %d)", n))
}
