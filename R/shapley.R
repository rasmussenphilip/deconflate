#' Shapley attribution of a general loss function over disease combinations
#'
#' Computes each disease's expected Shapley value for an arbitrary loss
#' function of the diseases an animal has. The calculation runs cell by cell
#' over the joint distribution ([fit_joint()]): within a combination, only
#' the diseases present are enumerated, so the cost depends on how many
#' diseases co-occur rather than on the total number of diseases. Shapley
#' values are linear in the game, so the expected Shapley value equals the
#' probability-weighted sum over cells.
#'
#' For additive loss functions the result has a closed form: each term is
#' split equally among the diseases it involves ([loss_additive()]). The
#' cell-wise calculation is needed for non-additive losses such as
#' [loss_multiplicative()].
#'
#' @param joint A [fit_joint()] result.
#' @param loss A function of a named 0/1 vector (one element per disease)
#'   returning the loss for an animal with those diseases; a loss of zero for
#'   no disease is assumed. See [loss_additive()] and [loss_multiplicative()].
#' @param max_present Cells with more diseases present than this are skipped;
#'   the skipped probability mass is reported.
#' @return A data frame with each disease's expected Shapley value and share,
#'   with attributes `total` (expected loss over the cells evaluated) and
#'   `skipped_mass`.
#' @export
#' @examples
#' j <- fit_joint(example_supplement())
#' m <- c(d1 = 0.02, d2 = 0.04, d3 = 0.06)
#' shapley_by_cell(j, loss_multiplicative(m))
shapley_by_cell <- function(joint, loss, max_present = 10L) {
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().")
  ids <- joint$diseases
  n <- length(ids)
  phi <- stats::setNames(numeric(n), ids)
  total <- 0
  skipped <- 0
  wcache <- list()
  for (r in seq_len(nrow(joint$cells))) {
    pr <- joint$prob[r]
    if (pr <= 0) next
    present <- which(joint$cells[r, ] == 1)
    k <- length(present)
    if (k == 0L) next
    if (k > max_present) {
      skipped <- skipped + pr
      next
    }
    masks <- 0:(2^k - 1)
    vals <- vapply(masks, function(mask) {
      x <- stats::setNames(numeric(n), ids)
      x[present[bitwAnd(mask, 2^(0:(k - 1))) > 0]] <- 1
      loss(x)
    }, numeric(1))
    size <- vapply(masks, function(mask) sum(bitwAnd(mask, 2^(0:(k - 1))) > 0), numeric(1))
    key <- as.character(k)
    if (is.null(wcache[[key]])) {
      wcache[[key]] <- factorial(0:(k - 1)) * factorial(k - 1 - 0:(k - 1)) / factorial(k)
    }
    w <- wcache[[key]]
    for (j in seq_len(k)) {
      bit <- 2^(j - 1)
      without <- masks[bitwAnd(masks, bit) == 0]
      contrib <- vals[without + bit + 1] - vals[without + 1]
      phi[present[j]] <- phi[present[j]] + pr * sum(w[size[without + 1] + 1] * contrib)
    }
    total <- total + pr * vals[2^k]
  }
  out <- data.frame(disease = ids, shapley = unname(phi),
                    share = unname(phi) / sum(phi), stringsAsFactors = FALSE)
  attr(out, "total") <- total
  attr(out, "skipped_mass") <- skipped
  out
}

#' Loss functions for Shapley attribution
#'
#' Builders for the `loss` argument of [shapley_by_cell()].
#'
#' * `loss_additive()`: `sum_i m[i] x[i]` plus higher-order terms, each
#'   applying only when all its diseases are present.
#' * `loss_multiplicative()`: `1 - prod_i (1 - m[i])^x[i]`, i.e.
#'   proportional impacts that compound.
#'
#' @param main Named numeric vector of per-disease impacts.
#' @param terms Optional list of terms, each a list with `diseases`
#'   (character vector of two or more ids) and `value`.
#' @return A function of a named 0/1 vector.
#' @name loss_functions
#' @examples
#' f <- loss_additive(c(a = 0.02, b = 0.03),
#'                    terms = list(list(diseases = c("a", "b"), value = 0.01)))
#' f(c(a = 1, b = 1))
NULL

#' @rdname loss_functions
#' @export
loss_additive <- function(main, terms = NULL) {
  if (is.null(names(main))) cm_abort("`main` must be named by disease id.")
  for (t in terms) {
    if (is.null(t$diseases) || length(t$diseases) < 2L || is.null(t$value)) {
      cm_abort("Each term needs `diseases` (two or more ids) and `value`.")
    }
  }
  function(x) {
    out <- sum(main * x[names(main)])
    for (t in terms) out <- out + t$value * prod(x[t$diseases])
    out
  }
}

#' @rdname loss_functions
#' @export
loss_multiplicative <- function(main) {
  if (is.null(names(main))) cm_abort("`main` must be named by disease id.")
  function(x) 1 - prod((1 - main)^x[names(main)])
}
