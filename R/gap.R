#' Attribute the aggregate burden to diseases (Shapley allocation)
#'
#' For outcomes on the proportion scale, the expected proportional loss is
#'
#' `L = sum_i m[i] P(i) + sum_{j < k} delta[j, k] P(j and k)`.
#'
#' Removing any disease involved in a term removes that term. Under this
#' accounting convention, the Shapley value of each disease is its own term
#' plus an equal share of every interaction term it is involved in:
#' `s[i] = m[i] P(i) + 1/2 sum_k delta[i, k] P(i and k)`. The shares sum to `L`.
#'
#' Outcomes on the hazard-ratio scale are not additive burdens and are
#' skipped; use [attributable_risk()] for them.
#'
#' @param result A [deconflate()] result.
#' @return A data frame with, per outcome and disease, the main-effect
#'   burden, the disease's share of interaction burden, the total and the
#'   fraction of `L`.
#' @export
attribute_burden <- function(result) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  ids <- result$model$diseases$id
  P <- stats::setNames(result$model$diseases$prob, ids)
  J <- result$joint_pairs
  adj <- result$adjusted[result$adjusted$scale != "hazard_ratio", , drop = FALSE]
  if (!nrow(adj)) {
    return(data.frame(outcome = character(0), disease = character(0), main = numeric(0),
                      interaction = numeric(0), total = numeric(0), share = numeric(0),
                      stringsAsFactors = FALSE))
  }
  out <- lapply(split(adj, adj$outcome), function(a) {
    o <- a$outcome[1]
    m <- stats::setNames(a$adjusted, a$disease)[ids]
    D <- result$interactions[[o]][ids, ids]
    main <- m * P
    inter <- 0.5 * rowSums(D * J[ids, ids])
    total <- main + inter
    data.frame(outcome = o, disease = ids, main = main, interaction = inter,
               total = total, share = total / sum(total),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, unname(out))
  rownames(out) <- NULL
  out
}

#' Productivity gaps and their attribution
#'
#' Generalises Rasmussen et al. (2022), eqs. 17-22. For an outcome with
#' observed mean `x` and expected loss `L` (see [attribute_burden()]):
#' * proportion scale: the disease-free value is `x / (1 - L)` for outcomes
#'   that disease decreases (e.g. yield), or `x / (1 + L)` for outcomes that
#'   disease increases (e.g. calving interval);
#' * absolute scale: the disease-free value is `x + L` (decrease) or `x - L`
#'   (increase), with `x` in the impacts' units (e.g. a culling risk of
#'   `0.27` with excess-risk impacts from [as_impacts()]).
#'
#' The gap is attributed to diseases in proportion to their Shapley shares,
#' which reduces to eq. 22 when there are no interactions.
#'
#' @param result A [deconflate()] result.
#' @param observed Named numeric vector of observed mean values, one per
#'   outcome to evaluate (names = outcome labels).
#' @return A list with `summary` (observed, disease-free value, gap and `L`
#'   per outcome) and `attribution` (gap attributed to each disease, split
#'   into main and interaction parts).
#' @export
#' @examples
#' res <- deconflate(example_supplement(), method = "published")
#' productivity_gap(res, c(yield = 10000))
productivity_gap <- function(result, observed) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  if (is.null(names(observed)) || any(!nzchar(names(observed)))) {
    cm_abort("`observed` must be a named vector (names = outcomes).")
  }
  check_numeric(observed, "observed")
  burden <- attribute_burden(result)
  summ <- list()
  attr_rows <- list()
  for (o in names(observed)) {
    a <- result$adjusted[result$adjusted$outcome == o, , drop = FALSE]
    if (!nrow(a)) cm_abort(sprintf("Outcome '%s' is not in the result.", o))
    if (a$scale[1] == "hazard_ratio") {
      cm_abort(sprintf("Outcome '%s' is on the hazard-ratio scale; use attributable_risk() instead of productivity_gap().", o))
    }
    b <- burden[burden$outcome == o, , drop = FALSE]
    L <- sum(b$total)
    x <- observed[[o]]
    if (a$scale[1] == "absolute") {
      xh <- if (a$direction[1] == "decrease") x + L else x - L
      gap <- L
    } else if (a$direction[1] == "decrease") {
      if (L >= 1) cm_abort(sprintf("Outcome '%s': total proportional loss is >= 1.", o))
      xh <- x / (1 - L)
      gap <- xh - x
    } else {
      xh <- x / (1 + L)
      gap <- x - xh
    }
    summ[[o]] <- data.frame(outcome = o, observed = x, disease_free = xh,
                            gap = gap, total_loss_fraction = L,
                            stringsAsFactors = FALSE)
    attr_rows[[o]] <- data.frame(outcome = o, disease = b$disease,
                                 gap = b$total / L * gap,
                                 gap_main = b$main / L * gap,
                                 gap_interaction = b$interaction / L * gap,
                                 stringsAsFactors = FALSE)
  }
  list(summary = do.call(rbind, unname(summ)),
       attribution = do.call(rbind, unname(attr_rows)))
}

#' Value productivity gaps in monetary terms
#'
#' Multiplies each outcome's gap (from [productivity_gap()]) by a unit value
#' and returns disease-specific and total losses, as in Rasmussen et al.
#' (2022), Tables 9-10.
#'
#' @param gap A [productivity_gap()] result.
#' @param unit_value Named numeric vector: value per unit of gap for each
#'   outcome. For example, with yield in kg the milk price per kg; with
#'   calving interval in days the value of a day's lost milk; with culling
#'   in percentage points the replacement cost per percentage point.
#' @param additional Optional named numeric vector of costs added to the
#'   total as lump sums (e.g. veterinary expenditure per animal).
#' @return A list with `by_disease` (outcome, disease, gap, value),
#'   `by_outcome` (gap and value per outcome), `by_disease_total` (value per
#'   disease summed over outcomes) and `total` (including `additional`).
#' @export
#' @examples
#' res <- deconflate(example_supplement(), method = "published")
#' value_losses(productivity_gap(res, c(yield = 10000)), c(yield = 0.30))
value_losses <- function(gap, unit_value, additional = 0) {
  if (!is.list(gap) || is.null(gap$attribution)) cm_abort("`gap` must come from productivity_gap().")
  if (is.null(names(unit_value))) cm_abort("`unit_value` must be named by outcome.")
  att <- gap$attribution
  miss <- setdiff(unique(att$outcome), names(unit_value))
  if (length(miss)) cm_abort(sprintf("No unit value for outcome(s): %s.", paste(miss, collapse = ", ")))
  by_disease <- data.frame(outcome = att$outcome, disease = att$disease, gap = att$gap,
                           value = att$gap * unname(unit_value[att$outcome]),
                           stringsAsFactors = FALSE)
  rownames(by_disease) <- NULL
  summ <- gap$summary
  by_outcome <- data.frame(outcome = summ$outcome, gap = summ$gap,
                           value = summ$gap * unname(unit_value[summ$outcome]),
                           stringsAsFactors = FALSE)
  rownames(by_outcome) <- NULL
  per_disease <- stats::aggregate(value ~ disease, data = by_disease, FUN = sum)
  list(by_disease = by_disease, by_outcome = by_outcome,
       by_disease_total = per_disease[order(-per_disease$value), , drop = FALSE],
       total = sum(by_outcome$value) + sum(additional))
}
