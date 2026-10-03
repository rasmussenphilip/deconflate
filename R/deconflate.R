#' Adjust raw impact estimates for comorbidity
#'
#' Raw impact estimates (comparisons of animals with and without a disease)
#' are treated as conflations of the disease's own impact and the impacts of
#' associated diseases. For one vector of additive impacts `b` (in any
#' units), the raw estimates satisfy
#'
#' `raw = A %*% b + offset`,
#'
#' where `A` depends on the disease probabilities and associations (and on
#' each estimate's estimand), and `offset` holds the contribution of any
#' pairwise interactions. Results are in the units of the impacts supplied.
#'
#' @section Estimands and the conflation matrix:
#' For a crude estimate of disease `i`, `A[i, k] = P(k | i) - P(k | not i)`
#' (the excess probability; Rasmussen et al. 2022, eq. 14). For a coefficient
#' from an additive regression adjusted for the diseases `S`
#' (`estimand = "adjusted_linear"` in [cm_impacts()]), `A[i, k]` is the
#' coefficient of `D_i` in the population linear projection of `D_k` on
#' `D_i` and `D_S`, and `A[i, S] = 0`. These coefficients depend only on the
#' pairwise tables. They assume the probabilities and associations describe
#' the population the estimates come from. An adjustment set that is
#' collinear with the disease is not identifiable and is rejected.
#'
#' With interactions `delta` (global method), the offset of disease `i` is
#' the corresponding regression coefficient of the interaction burden
#' `sum_{j<k} delta[j, k] D_j D_k`, which needs the joint distribution of
#' disease combinations.
#'
#' @section Aggregate and contributions:
#' The expected aggregate impact per animal is
#' `sum_i p_i b_i + sum_{j<k} P(j and k) delta[j, k]`. Each disease's
#' contribution is its own term plus half of each interaction term it is
#' involved in (the closed-form Shapley value); contributions add up to the
#' aggregate. Shares are `NA` when the aggregate is zero.
#'
#' @section Feasibility:
#' Pairwise tables can each be valid while no population has all of them
#' (an invertible `A` does not mean the inputs are feasible). By default the
#' pairwise methods screen every triple of diseases (a necessary condition;
#' see [check_feasibility()]); `feasibility = "lp"` runs the exact check. The
#' global method fits the joint distribution, and stops if it does not
#' converge.
#'
#' @param model A [cm_model()], or a [cm_analyses()] object (each analysis is
#'   adjusted in turn).
#' @param method
#'   * `"simultaneous"` (default): solves `raw = A b` exactly using the
#'     pairwise 2x2 tables. No joint distribution is needed.
#'   * `"published"`: the proportional approximation of Rasmussen et al.
#'     (2022), eq. 16: `b[i] = raw[i]^2 / (raw[i] + sum_{k != i} A[i, k] raw[k])`,
#'     for crude estimates only. Provided for reproduction and comparison;
#'     it can mask incompatible inputs, and its value is undefined when the
#'     denominator is zero.
#'   * `"global"`: fits the maximum-entropy joint distribution
#'     ([fit_joint()]; the "iterative" model) and solves the equations
#'     including any interactions. Unknown pairs are left unconstrained.
#'     Without interactions and unknown pairs, it equals `"simultaneous"`.
#' @param joint Optional [fit_joint()] result for `method = "global"`. It is
#'   checked against the model (diseases, probabilities, associations and
#'   three-way terms); impact-only changes do not require a refit.
#' @param warn Logical: warn when adjusted impacts change sign or are not
#'   finite?
#' @param feasibility For the pairwise methods: `"screen"` (default; triple
#'   screen), `"lp"` (exact check, needs `lpSolve`) or `"none"`.
#' @param ... Passed to [fit_joint()] (method `"global"`), or to the method
#'   for each analysis.
#'
#' @return A `cm_result` with elements:
#'   * `adjusted`: raw and adjusted impacts, relative change and estimand;
#'   * `totals`: the sum of `p_i * raw_i` (the naive aggregate), the adjusted
#'     aggregate and its interaction part;
#'   * `contributions`: per disease, main and interaction contributions,
#'     their total and share;
#'   * `diagnostics`: reconstruction residual, rank and condition number of
#'     `A`, sign changes and the feasibility check;
#'   * `conflation` (`A` and `offset`), `interactions` (matrix of `delta`),
#'     `joint_pairs` (matrix of `P(j and k)`), `model`, `joint`, `method`,
#'     `label` and `units`.
#'   For a [cm_analyses()] object, a named list of results (class
#'   `cm_results`).
#' @export
#' @examples
#' res <- deconflate(example_supplement())
#' res$adjusted
#' res$totals
#'
#' # A coefficient adjusted for one other disease
#' m <- example_supplement()
#' m$impacts <- cm_impacts(c("d1", "d2", "d3"), c(2.2, 5, 7.5),
#'                         estimand = c("adjusted_linear", "crude", "crude"),
#'                         adjusted_for = c("d2", NA, NA), units = "%")
#' deconflate(m)$adjusted
deconflate <- function(model, ...) UseMethod("deconflate")

#' @export
deconflate.default <- function(model, ...) {
  cm_abort("`model` must be created with cm_model() or cm_analyses().")
}

#' @rdname deconflate
#' @export
deconflate.cm_model <- function(model, method = c("simultaneous", "published", "global"),
                                joint = NULL, warn = TRUE,
                                feasibility = c("screen", "lp", "none"), ...) {
  method <- match.arg(method)
  feasibility <- match.arg(feasibility)
  imp <- model$impacts
  if (is.null(imp)) cm_abort("The model has no impacts to adjust.")
  ids <- model$diseases$id
  if (!identical(imp$disease, ids)) {
    # Impacts replaced after cm_model() (e.g. `m$impacts <- ...`): validate
    # them and put them in disease order.
    model <- cm_model(model, imp, model$interactions)
    imp <- model$impacts
  }
  n <- length(ids)
  P <- stats::setNames(model$diseases$prob, ids)
  m_raw <- imp$value
  has_int <- !is.null(model$interactions) && nrow(model$interactions) > 0L
  if (method == "published" && any(imp$estimand != "crude")) {
    cm_abort("The published approximation is defined for crude estimates only; use method = 'simultaneous' or 'global' for adjusted estimands.",
             class = "deconflate_unsupported")
  }
  if (has_int && method != "global") {
    cm_abort("Interactions require method = 'global' (they need probabilities of disease triples).",
             class = "deconflate_unsupported")
  }

  cells <- NULL
  if (method == "global") {
    if (is.null(joint)) {
      joint <- withCallingHandlers(fit_joint(model, ...),
                                   deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
    } else {
      validate_joint(joint, model)
    }
    if (!joint$converged) {
      cm_abort(sprintf("The joint distribution did not converge (max residual %.2e). The pairwise associations may be jointly infeasible; see check_feasibility(method = 'lp').",
                       joint$max_residual), class = "deconflate_nonconvergence")
    }
    cells <- joint$cells
    J <- crossprod(cells * joint$prob, cells)
    dimnames(J) <- list(ids, ids)
    feas <- sprintf("joint distribution fitted (max residual %.1e)", joint$max_residual)
  } else {
    pt <- pair_tables(model)
    if (anyNA(pt$p11)) {
      unk <- paste(pt$disease1, pt$disease2, sep = "-")[is.na(pt$p11)]
      cm_abort(sprintf(
        "Unknown associations for %s. Specify them, set missing_associations = 'independent', or use method = 'global'.",
        paste(unk, collapse = ", ")), class = c("deconflate_unknown_pairs", "deconflate_unsupported"))
    }
    J <- diag(P, nrow = n)
    dimnames(J) <- list(ids, ids)
    J[cbind(pt$disease1, pt$disease2)] <- pt$p11
    J[cbind(pt$disease2, pt$disease1)] <- pt$p11
    feas <- switch(feasibility,
      screen = {
        tri <- triple_screen(pt, P)
        if (nrow(tri)) {
          cm_abort(sprintf("The pairwise associations are jointly infeasible: no population has all of them (conflicting triples: %s). See check_feasibility().",
                           paste(utils::head(paste(tri$disease1, tri$disease2, tri$disease3, sep = "-"), 5),
                                 collapse = ", ")), class = "deconflate_infeasible")
        }
        "triple screen passed (necessary condition only)"
      },
      lp = {
        cf <- check_feasibility(model, method = "lp")
        if (!isTRUE(cf$feasible)) {
          cm_abort("The pairwise associations are jointly infeasible (LP check); see check_feasibility().",
                   class = "deconflate_infeasible")
        }
        "LP check passed (jointly feasible)"
      },
      none = "unchecked"
    )
  }

  A <- projection_matrix(pair_covariance(J), imp, ids)
  offset <- rep(0, n)
  D <- matrix(0, n, n, dimnames = list(ids, ids))
  if (has_int) {
    int <- model$interactions
    D[cbind(int$disease1, int$disease2)] <- int$value
    D[cbind(int$disease2, int$disease1)] <- int$value
    # Interaction burden of each combination: sum_{j<k} delta_jk D_j D_k.
    burden <- rowSums((cells %*% D) * cells) / 2
    offset <- projection_offset(cells, joint$prob, burden, imp, ids)
  }

  if (method == "published") {
    conf <- as.vector((A - diag(n)) %*% m_raw)
    m_adj <- ifelse(m_raw == 0, 0, m_raw^2 / (m_raw + conf))
    rk <- qr(A)$rank
  } else {
    rk <- qr(A)$rank
    if (rk < n) {
      cm_abort(sprintf("The conflation matrix is singular (rank %d of %d): the impacts are not identifiable from these inputs.",
                       rk, n), class = "deconflate_singular")
    }
    m_adj <- tryCatch(as.vector(solve(A, m_raw - offset)),
                      error = function(e) cm_abort(sprintf("Could not solve the adjustment equations: %s", conditionMessage(e)),
                                                   class = "deconflate_singular"))
  }
  m_adj <- unname(m_adj)
  recon <- as.vector(A %*% m_adj) + offset
  resid <- if (all(is.finite(recon))) max(abs(recon - m_raw)) else NA_real_
  sign_change <- is.finite(m_adj) & abs(m_raw) > 1e-12 & abs(m_adj) > 1e-12 &
    sign(m_adj) != sign(m_raw)

  if (warn && any(sign_change)) {
    cm_warn(sprintf(
      "Adjusted impacts change sign for %s. The raw impacts are smaller than the associated diseases alone would produce under the additive model; check the estimands, and whether the estimates come from populations with different comorbidity patterns.",
      paste(ids[sign_change], collapse = ", ")), class = "deconflate_sign_change")
  }
  if (warn && any(!is.finite(m_adj))) {
    cm_warn(sprintf("Non-finite adjusted impacts for %s.", paste(ids[!is.finite(m_adj)], collapse = ", ")),
            class = "deconflate_nonfinite_warning")
  }

  main <- unname(P * m_adj)
  inter <- unname(0.5 * rowSums(D * J))
  contrib <- main + inter
  total <- sum(contrib)
  share <- if (is.finite(total) && abs(total) > 0) contrib / total else rep(NA_real_, n)

  structure(list(
    adjusted = data.frame(disease = ids, raw = m_raw, adjusted = m_adj,
                          change = ifelse(m_raw != 0, m_adj / m_raw - 1, NA_real_),
                          estimand = imp$estimand, adjusted_for = imp$adjusted_for,
                          stringsAsFactors = FALSE),
    totals = data.frame(raw_sum = sum(P * m_raw), adjusted_total = total,
                        interaction_total = sum(inter)),
    contributions = data.frame(disease = ids, main = main, interaction = inter,
                               total = contrib, share = share, stringsAsFactors = FALSE),
    diagnostics = data.frame(method = method, max_reconstruction_residual = resid,
                             rank = rk, condition_number = kappa(A, exact = TRUE),
                             n_sign_changes = sum(sign_change),
                             sign_changes = paste(ids[sign_change], collapse = ", "),
                             feasibility = feas, stringsAsFactors = FALSE),
    conflation = list(A = A, offset = stats::setNames(offset, ids)),
    interactions = D, joint_pairs = J, model = model, joint = joint, method = method,
    label = attr(imp, "label"), units = attr(imp, "units")
  ), class = "cm_result")
}

#' @rdname deconflate
#' @export
deconflate.cm_analyses <- function(model, method = c("simultaneous", "published", "global"),
                                   joint = NULL, ...) {
  method <- match.arg(method)
  if (method == "global" && is.null(joint)) {
    # One joint distribution serves every analysis (it depends on the
    # population only).
    dots <- list(...)
    fj <- dots[intersect(names(dots), c("tol", "max_iter", "max_diseases"))]
    joint <- withCallingHandlers(do.call(fit_joint, c(list(model$population), fj)),
                                 deconflate_nonconvergence = function(w) invokeRestart("muffleWarning"))
  }
  out <- lapply(model$models, deconflate, method = method, joint = joint, ...)
  structure(out, class = "cm_results")
}

# Regression coefficient of D_i for the interaction burden g, given the
# estimand of each impact (crude: difference in means; adjusted_linear: the
# coefficient in a population regression on D_i and D_S).
projection_offset <- function(cells, prob, g, impacts, ids) {
  n <- length(ids)
  p <- colSums(cells * prob)
  J <- crossprod(cells * prob, cells)
  Sigma <- J - outer(p, p)
  diag(Sigma) <- p * (1 - p)
  cov_g <- as.vector(crossprod(cells, prob * g)) - p * sum(prob * g)
  vapply(seq_len(n), function(i) {
    S <- if (impacts$estimand[i] == "adjusted_linear") {
      match(resolve_adjusted_for(impacts$adjusted_for[i], ids[i], ids), ids)
    } else integer(0)
    X <- c(i, S)
    tryCatch(solve(Sigma[X, X, drop = FALSE], cov_g[X])[1],
             error = function(e) cm_abort(sprintf("The interaction offset of %s is not identifiable: %s",
                                                  ids[i], conditionMessage(e)),
                                          class = "deconflate_singular"))
  }, numeric(1))
}

#' Attribute the aggregate to diseases (Shapley allocation)
#'
#' The expected aggregate impact is
#' `sum_i p_i b_i + sum_{j < k} delta[j, k] P(j and k)`. Removing any disease
#' involved in a term removes that term, so each disease's Shapley value is
#' its own term plus an equal share of every interaction term it is involved
#' in. The shares add up to the aggregate. This is the `contributions` table
#' of [deconflate()].
#'
#' @param result A [deconflate()] result.
#' @return A data frame with, per disease, the main contribution, the share
#'   of interaction terms, the total and the share of the aggregate (`NA`
#'   when the aggregate is zero).
#' @export
attribute_burden <- function(result) {
  if (!inherits(result, "cm_result")) cm_abort("`result` must come from deconflate().")
  result$contributions
}

#' @export
print.cm_results <- function(x, ...) {
  cat(sprintf("<cm_results> %d analyses: %s\n\n", length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("== %s ==\n", nm))
    print(x[[nm]])
    cat("\n")
  }
  invisible(x)
}
