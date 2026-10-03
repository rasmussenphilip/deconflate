#' Productivity gap of one adjusted impact vector (optional helper)
#'
#' Converts the adjusted aggregate impact of one analysis into the gap
#' between the observed mean of an outcome and its disease-free value, and
#' attributes the gap to diseases. This generalises Rasmussen et al. (2022),
#' eqs. 17-22. It is a separate step from the adjustment, because it needs to
#' know what the impacts mean.
#'
#' With aggregate `L` and observed mean `x`:
#' * `effect = "proportion"` (or `"percent"`, divided by 100 first): impacts
#'   are proportional changes relative to the disease-free value. For an
#'   outcome that disease decreases (e.g. yield) the disease-free value is
#'   `x / (1 - L)`; for one it increases (e.g. calving interval), `x / (1 + L)`.
#' * `effect = "absolute"`: impacts are in the outcome's units; the
#'   disease-free value is `x + L` (decrease) or `x - L` (increase).
#'
#' Each disease's part of the gap is its contribution (see
#' [attribute_burden()]) times the same factor, computed directly rather
#' than by dividing by `L`, so it stays defined when contributions cancel.
#'
#' @param result A [deconflate()] result.
#' @param observed Observed mean of the outcome.
#' @param direction `"decrease"` or `"increase"`: what disease does to the
#'   outcome.
#' @param effect `"proportion"`, `"percent"` or `"absolute"`: how the impacts
#'   relate to the outcome.
#' @return A list with `summary` (observed, disease-free value, gap and
#'   aggregate) and `attribution` (gap attributed to each disease, split into
#'   main and interaction parts).
#' @export
#' @examples
#' res <- deconflate(example_supplement(), method = "published")
#' productivity_gap(res, observed = 10000, direction = "decrease", effect = "percent")
productivity_gap <- function(result, observed, direction = c("decrease", "increase"),
                             effect = c("proportion", "percent", "absolute")) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  direction <- match.arg(direction)
  effect <- match.arg(effect)
  check_numeric(observed, "observed")
  if (length(observed) != 1L) cm_abort("`observed` must be a single number.")
  ct <- result$contributions
  scale <- if (effect == "percent") 100 else 1
  L <- sum(ct$total) / scale
  x <- observed
  if (effect == "absolute") {
    factor <- 1
    xh <- if (direction == "decrease") x + L else x - L
  } else if (direction == "decrease") {
    if (!(L < 1)) cm_abort(sprintf("The aggregate proportional loss is %g; it must be below 1.", L))
    factor <- x / (1 - L)
    xh <- x / (1 - L)
  } else {
    if (!(L > -1)) cm_abort(sprintf("The aggregate proportional increase is %g; it must be above -1.", L))
    factor <- x / (1 + L)
    xh <- x / (1 + L)
  }
  gap <- if (direction == "decrease") xh - x else x - xh
  list(
    summary = data.frame(observed = x, disease_free = xh, gap = gap, aggregate = L,
                         direction = direction, effect = effect, stringsAsFactors = FALSE),
    attribution = data.frame(disease = ct$disease,
                             gap = factor * ct$total / scale,
                             gap_main = factor * ct$main / scale,
                             gap_interaction = factor * ct$interaction / scale,
                             stringsAsFactors = FALSE)
  )
}

#' Value a productivity gap in monetary terms (optional helper)
#'
#' Multiplies a gap (from [productivity_gap()]) by a unit value, as in
#' Rasmussen et al. (2022), Tables 9-10. Lump-sum costs that are not part of
#' any adjustment (e.g. veterinary expenditure) can be added with
#' `additional`; they are reported separately and included in `total`.
#'
#' @param gap A [productivity_gap()] result.
#' @param unit_value Value per unit of gap (e.g. milk price per kg).
#' @param additional Optional named numeric vector of lump-sum costs.
#' @return A list with `by_disease` (gap and value), `value` (the gap's
#'   value), `additional` and `total`.
#' @export
#' @examples
#' res <- deconflate(example_supplement(), method = "published")
#' value_losses(productivity_gap(res, 10000, "decrease", "percent"), unit_value = 0.30)
value_losses <- function(gap, unit_value, additional = 0) {
  if (!is.list(gap) || is.null(gap$attribution)) cm_abort("`gap` must come from productivity_gap().")
  check_numeric(unit_value, "unit_value")
  if (length(unit_value) != 1L) cm_abort("`unit_value` must be a single number.")
  check_numeric(additional, "additional")
  att <- gap$attribution
  by_disease <- data.frame(disease = att$disease, gap = att$gap, value = att$gap * unit_value,
                           stringsAsFactors = FALSE)
  value <- gap$summary$gap * unit_value
  list(by_disease = by_disease, value = value, additional = additional,
       total = value + sum(additional))
}
