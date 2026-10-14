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
#' split equally among the diseases it involves. The cell-wise calculation
#' is needed for non-additive losses, such as the attributable risk of
#' event impacts.
#'
#' @param joint A [fit_joint()] result.
#' @param loss A function of a named 0/1 vector (one element per disease)
#'   returning the loss for an animal with those diseases; a loss of zero for
#'   no disease is assumed.
#' @param max_present Cells with more diseases present than this are skipped
#'   (the cost grows as 2^k for a cell with k diseases). The default skips
#'   nothing. Skipped cells are reported: their probability (`skipped_mass`)
#'   and their expected loss (`skipped_loss`), which is the part of the total
#'   that is not allocated. The warning (class
#'   `deconflate_incomplete_allocation`) has the same values as fields
#'   (`skipped_mass`, `skipped_loss`, `total`).
#' @return A data frame with each disease's expected Shapley value and share
#'   (of the allocated amount, which is the total only when no cell is
#'   skipped), with attributes `total` (expected loss over all cells),
#'   `allocated` (the sum of the Shapley values), `skipped_mass` and
#'   `skipped_loss`.
#' @noRd
shapley_by_cell <- function(joint, loss, max_present = Inf) {
  if (!inherits(joint, "cm_joint")) cm_abort("`joint` must come from fit_joint().", class = "deconflate_invalid_input")
  ids <- joint$diseases
  n <- length(ids)
  phi <- stats::setNames(numeric(n), ids)
  total <- 0
  skipped <- 0
  skipped_loss <- 0
  wcache <- list()
  # Every subset of the diseases present in a combination is itself a
  # combination, so the loss is evaluated at most once per combination
  # (cached by its binary code) rather than once per subset and cell.
  if (n > 52L) cm_abort("Cell-wise Shapley allocation supports at most 52 diseases.", class = "deconflate_unsupported")
  bits <- 2^(seq_len(n) - 1)
  # Small n: a vector indexed by code; large n (sampled joints): a hash.
  if (n <= 20L) {
    vcache <- rep(NA_real_, 2^n)
    loss_of <- function(codes) {
      need <- unique(codes[is.na(vcache[codes + 1])])
      for (cd in need) {
        x <- stats::setNames(as.numeric(floor(cd / bits) %% 2 == 1), ids)
        vcache[cd + 1] <<- loss(x)
      }
      vcache[codes + 1]
    }
  } else {
    cache <- new.env(hash = TRUE, parent = emptyenv())
    loss_of <- function(codes) {
      vapply(codes, function(cd) {
        key <- sprintf("%.0f", cd)
        v <- cache[[key]]
        if (is.null(v)) {
          x <- stats::setNames(as.numeric(floor(cd / bits) %% 2 == 1), ids)
          v <- loss(x)
          assign(key, v, envir = cache)
        }
        v
      }, numeric(1))
    }
  }
  for (r in seq_len(nrow(joint$cells))) {
    pr <- joint$prob[r]
    if (pr <= 0) next
    present <- which(joint$cells[r, ] == 1)
    k <- length(present)
    if (k == 0L) next
    if (k > max_present) {
      lv <- loss_of(sum(bits[present]))
      skipped <- skipped + pr
      skipped_loss <- skipped_loss + pr * lv
      total <- total + pr * lv
      next
    }
    masks <- 0:(2^k - 1)
    # Membership of each present disease in each subset (rows = masks).
    member <- vapply(seq_len(k), function(j) bitwAnd(masks, 2^(j - 1)) > 0, logical(2^k))
    member <- matrix(member, nrow = 2^k)
    vals <- loss_of(as.vector(member %*% bits[present]))
    size <- rowSums(member)
    key <- as.character(k)
    if (is.null(wcache[[key]])) {
      wcache[[key]] <- factorial(0:(k - 1)) * factorial(k - 1 - 0:(k - 1)) / factorial(k)
    }
    w <- wcache[[key]]
    for (j in seq_len(k)) {
      bit <- 2^(j - 1)
      without <- masks[!member[, j]]
      contrib <- vals[without + bit + 1] - vals[without + 1]
      phi[present[j]] <- phi[present[j]] + pr * sum(w[size[without + 1] + 1] * contrib)
    }
    total <- total + pr * vals[2^k]
  }
  tot_phi <- sum(phi)
  out <- data.frame(disease = ids, shapley = unname(phi),
                    share = if (abs(tot_phi) > 0) unname(phi) / tot_phi else NA_real_,
                    stringsAsFactors = FALSE)
  attr(out, "total") <- total
  attr(out, "allocated") <- tot_phi
  attr(out, "skipped_mass") <- skipped
  attr(out, "skipped_loss") <- skipped_loss
  if (skipped > 0) {
    cm_warn(sprintf("Shapley allocation skipped combinations with more than %g diseases (probability %.3g); %.3g of the total %.3g is not allocated.",
                    max_present, skipped, skipped_loss, total), class = "deconflate_incomplete_allocation",
            skipped_mass = skipped, skipped_loss = skipped_loss, total = total)
  }
  out
}
