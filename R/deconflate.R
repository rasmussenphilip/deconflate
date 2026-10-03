#' Adjust raw impact estimates for comorbidity
#'
#' Raw impact estimates (comparisons of animals with and without a disease)
#' are treated as conflations of the disease's own impact and the impacts of
#' associated diseases. Under additive impacts, the raw impact of disease `i`
#' is
#'
#' `m_raw[i] = m[i] + sum_k E[k, i] * m[k]`
#'
#' and, with pairwise interactions `delta`, additionally
#'
#' `+ sum_k delta[i, k] * P(k | i) + sum_{j < k; j, k != i} delta[j, k] * (P(j, k | i) - P(j, k | not i))`.
#'
#' @param model A [cm_model()] with impacts.
#' @param method
#'   * `"simultaneous"` (default): solves the additive equations exactly
#'     (`(I + t(E)) m = m_raw`) using the pairwise 2x2 tables. No joint
#'     distribution is needed.
#'   * `"published"`: the proportional approximation of Rasmussen et al.
#'     (2022), eq. 16: `m[i] = m_raw[i]^2 / (m_raw[i] + sum_k E[k, i] * m_raw[k])`.
#'     Provided for reproduction and comparison; it can mask incompatible
#'     inputs because it cannot return negative impacts.
#'   * `"global"`: fits the maximum-entropy joint distribution
#'     ([fit_joint()]) and solves the full equations, including any
#'     interactions in the model. Unknown pairs are left unconstrained.
#'     Without interactions and unknown pairs, it equals `"simultaneous"`.
#' @param joint Optional pre-computed [fit_joint()] result for
#'   `method = "global"`.
#' @param warn Logical: warn when adjusted impacts change sign or the method
#'   produces non-finite values?
#' @param ... Passed to [fit_joint()].
#'
#' @details Impacts flagged with `adjusted_for` in [cm_impacts()] have the
#'   corresponding conflation terms set to zero.
#'
#'   The diagnostics report, per outcome: the maximum absolute difference
#'   between the supplied raw impacts and those reconstructed from the
#'   adjusted impacts under the additive (or interaction) model; the
#'   number of adjusted impacts whose sign differs from the raw impact; and
#'   the condition number of the conflation matrix. For the published method
#'   a non-zero reconstruction residual is expected: it measures the
#'   approximation error.
#'
#' @return A `cm_result` with elements `adjusted` (long data frame of raw and
#'   adjusted impacts), `diagnostics`, `conflation` (the matrix `A` per
#'   outcome, with `m_raw = A m + offset`), `joint_pairs` (matrix of
#'   `P(j and k)`), the `model`, the `joint` fit (global only) and `method`.
#' @export
#' @examples
#' res <- deconflate(example_supplement())
#' res$adjusted
deconflate <- function(model, method = c("simultaneous", "published", "global"),
                       joint = NULL, warn = TRUE, ...) {
  check_model(model)
  method <- match.arg(method)
  if (is.null(model$impacts)) cm_abort("The model has no impacts to adjust.")
  ids <- model$diseases$id
  n <- length(ids)
  has_int <- !is.null(model$interactions) && nrow(model$interactions) > 0L

  if (method == "global") {
    joint <- joint %||% fit_joint(model, ...)
    if (!joint$converged) {
      cm_abort("The joint distribution did not converge; the pairwise associations may be jointly infeasible (see check_feasibility()).",
               class = "deconflate_infeasible")
    }
    cells <- joint$cells
    base_A <- vapply(seq_len(n), function(k) crude_difference(cells, joint$prob, cells[, k]),
                     numeric(n))
    base_A <- matrix(base_A, n, n, dimnames = list(ids, ids))
    J <- crossprod(cells * joint$prob, cells)
    dimnames(J) <- list(ids, ids)
  } else {
    if (has_int) {
      cm_abort("Interactions require method = 'global' (they need probabilities of disease triples).")
    }
    E <- excess_matrix(model)
    base_A <- diag(n) + t(E)
    dimnames(base_A) <- list(ids, ids)
    J <- pairwise_joint_matrix(model)
    joint <- NULL
  }

  outcomes <- unique(model$impacts$outcome)
  adjusted <- list()
  diagnostics <- list()
  conflation <- list()
  deltas <- list()
  for (o in outcomes) {
    imp <- model$impacts[model$impacts$outcome == o, , drop = FALSE]
    imp <- imp[match(ids, imp$disease), , drop = FALSE]
    m_raw <- imp$value
    A <- base_A
    for (i in seq_len(n)) {
      adj <- setdiff(split_ids(imp$adjusted_for[i]), ids[i])
      if (length(adj)) A[i, adj] <- 0
    }

    D <- matrix(0, n, n, dimnames = list(ids, ids))
    offset <- rep(0, n)
    if (method == "global" && has_int) {
      int <- model$interactions[model$interactions$outcome == o, , drop = FALSE]
      for (r in seq_len(nrow(int))) {
        D[int$disease1[r], int$disease2[r]] <- int$value[r]
        D[int$disease2[r], int$disease1[r]] <- int$value[r]
      }
      # Interaction burden of each combination: sum_{j<k} delta_jk D_j D_k.
      burden <- rowSums((cells %*% D) * cells) / 2
      offset <- crude_difference(cells, joint$prob, burden)
    }

    m_adj <- switch(method,
      published = {
        conf <- as.vector((A - diag(diag(A))) %*% m_raw)
        ifelse(m_raw == 0, 0, m_raw^2 / (m_raw + conf))
      },
      solve(A, m_raw - offset)
    )
    m_adj <- as.vector(m_adj)
    recon <- as.vector(A %*% m_adj) + offset
    sign_change <- is.finite(m_adj) & abs(m_raw) > 1e-12 & abs(m_adj) > 1e-12 &
      sign(m_adj) != sign(m_raw)

    if (warn && any(sign_change)) {
      cm_warn(sprintf(
        "Outcome '%s': adjusted impacts change sign for %s. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check whether these estimates were already adjusted for co-diseases or come from populations with different comorbidity patterns.",
        o, paste(ids[sign_change], collapse = ", ")), class = "deconflate_sign_change")
    }
    if (warn && any(!is.finite(m_adj))) {
      cm_warn(sprintf("Outcome '%s': non-finite adjusted impacts for %s.", o,
                      paste(ids[!is.finite(m_adj)], collapse = ", ")))
    }

    adjusted[[o]] <- data.frame(
      outcome = o, disease = ids, raw = m_raw, adjusted = m_adj,
      change = ifelse(m_raw != 0, m_adj / m_raw - 1, NA_real_),
      scale = imp$scale, units = imp$units, direction = imp$direction,
      stringsAsFactors = FALSE
    )
    diagnostics[[o]] <- data.frame(
      outcome = o, method = method,
      max_reconstruction_residual = max(abs(recon - m_raw)),
      n_sign_changes = sum(sign_change),
      sign_changes = paste(ids[sign_change], collapse = ", "),
      condition_number = kappa(A, exact = TRUE),
      stringsAsFactors = FALSE
    )
    conflation[[o]] <- list(A = A, offset = stats::setNames(offset, ids))
    deltas[[o]] <- D
  }

  structure(list(
    adjusted = do.call(rbind, unname(adjusted)),
    diagnostics = do.call(rbind, unname(diagnostics)),
    conflation = conflation, interactions = deltas, joint_pairs = J,
    model = model, joint = joint, method = method
  ), class = "cm_result")
}

#' Compare adjustment methods
#'
#' Runs [deconflate()] with several methods and returns adjusted impacts side
#' by side.
#'
#' @param model A [cm_model()].
#' @param methods Methods to compare.
#' @param ... Passed to [deconflate()].
#' @return A data frame with one column of adjusted impacts per method.
#' @export
compare_methods <- function(model, methods = c("published", "simultaneous", "global"), ...) {
  res <- lapply(methods, function(m) deconflate(model, method = m, warn = FALSE, ...)$adjusted)
  out <- res[[1]][, c("outcome", "disease", "raw")]
  for (k in seq_along(methods)) out[[methods[k]]] <- res[[k]]$adjusted
  out
}
