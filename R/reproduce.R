#' Reproduce the published analyses
#'
#' Recompute the tables of Rasmussen et al. (2022) and Rasmussen et al.
#' (2024) with the published (eq. 16) method and the conversions used in the
#' papers, and set them beside the printed values.
#'
#' These functions exist to document and check the published numbers. Two of
#' the conversions they use are not offered for new analyses:
#' * 2022: culling hazard ratios were converted to excess annual culling
#'   risks by treating them as odds ratios, and adjusted hazard ratios were
#'   recovered by rescaling (eq. 23);
#' * 2024: hazard ratios minus 1 were adjusted as additive impacts.
#'
#' For new analyses of hazard ratios use [deconflate_hr()] and
#' [attributable_risk()].
#'
#' @section 2022 (Tables 8-10):
#' `reproduce_rasmussen_2022()` adjusts yield and calving interval for the UK
#' example ([example_uk_dairy_2022()]) and the historical culling analysis,
#' computes the productivity gaps and their values with the paper's economic
#' inputs (Table 1, and its culling valuation), adds veterinary expenditure,
#' and back-converts the culling impacts to hazard ratios. (The gap and
#' valuation steps exist only here: the package itself does not value
#' impacts.)
#' Fertility and the culling hazard ratios reproduce the paper; yield does not
#' reproduce exactly from the printed Table 4 (see
#' `vignette("reproducing-published")`).
#'
#' @section 2024 (Table 5):
#' `reproduce_rasmussen_2024()` runs [cm_monte_carlo()] with the samplers of
#' [sampler_global_dairy()] plus the historical culling analysis (HR - 1,
#' named `culling_hr_minus_1`) and the published method, and
#' reports the means beside Table 5. Culling is reported on the hazard-ratio
#' scale (1 + adjusted HR - 1). Table 5 used 50,000 draws; use at least
#' several thousand for stable means. Some fertility means are unstable for
#' the published method (see the `stability` column).
#'
#' @param yield_sck Raw yield impact of subclinical ketosis in percent (see
#'   [example_uk_dairy_2022()]).
#' @param n_draws Number of Monte Carlo draws.
#' @param seed Random seed.
#' @param inputs `"analysis"` or `"tables"` (see [example_global_dairy()]).
#' @return A `cm_reproduction` object: a list with `adjusted` (adjusted
#'   values per disease), `comparison` (package beside the published values),
#'   and for 2022 `gaps`, `values` and `total`, or for 2024 the Monte Carlo
#'   run `mc`.
#' @name reproduce
#' @examples
#' r22 <- reproduce_rasmussen_2022()
#' r22$comparison
#' \donttest{
#' r24 <- reproduce_rasmussen_2024(n_draws = 500)
#' r24$comparison[r24$comparison$analysis == "yield", ]
#' }
NULL

#' @rdname reproduce
#' @export
reproduce_rasmussen_2022 <- function(yield_sck = 3.05) {
  uk <- uk_dairy_2022_analyses(yield_sck, culling = TRUE)
  eco <- uk_dairy_2022_economics()
  eco$valuation$culling <- uk_dairy_2022_culling_valuation()
  res <- adjust_analyses(uk, method = "published", warn = FALSE)
  ids <- uk$population$diseases$id
  nms <- names(uk$models)

  adjusted <- data.frame(disease = ids, stringsAsFactors = FALSE)
  for (nm in nms) adjusted[[nm]] <- res[[nm]]$adjusted$adjusted
  hr <- attr(uk, "hazard_ratios")
  cu <- res$culling$adjusted
  adjusted$culling_hr <- legacy_eq23(hr$value[match(ids, hr$disease)], cu$raw, cu$adjusted)

  gaps <- lapply(nms, function(nm) {
    v <- eco$valuation[[nm]]
    legacy_gap(res[[nm]], v$observed, v$direction, v$effect)
  })
  names(gaps) <- nms
  values <- lapply(nms, function(nm) legacy_value(gaps[[nm]], eco$valuation[[nm]]$unit_value))
  names(values) <- nms

  summary_tab <- do.call(rbind, lapply(nms, function(nm) {
    s <- gaps[[nm]]$summary
    data.frame(analysis = nm, observed = s$observed, disease_free = s$disease_free,
               gap = s$gap, value = values[[nm]]$value, stringsAsFactors = FALSE)
  }))
  by_disease <- data.frame(disease = ids, stringsAsFactors = FALSE)
  for (nm in nms) by_disease[[nm]] <- values[[nm]]$by_disease$value
  by_disease$total <- rowSums(by_disease[nms])
  total <- sum(summary_tab$value) + sum(eco$additional)

  published <- c(yield_disease_free = 9306.32, yield_value = 172.05,
                 fertility_disease_free = 375.09, fertility_value = 101.79,
                 culling_disease_free = 22.56, culling_value = 59.27,
                 total_with_veterinary = 404.21)
  pkg <- c(summary_tab$disease_free[1], summary_tab$value[1],
           summary_tab$disease_free[2], summary_tab$value[2],
           summary_tab$disease_free[3], summary_tab$value[3], total)
  comparison <- data.frame(quantity = names(published), table9 = unname(published),
                           package = pkg, stringsAsFactors = FALSE)
  table8_hr <- c(CO = 1, DA = 2.68, DYS = 1.60, FAS = 1, GIN = 1, LAM = 3.00, MAS = 2.22,
                 MET = 1.55, MF = 1.89, NEO = 1.60, PTB = 1.64, RP = 1, SCK = 1.29)
  hr_comparison <- data.frame(disease = ids, hr = hr$value[match(ids, hr$disease)],
                              table8 = unname(table8_hr[ids]),
                              package = adjusted$culling_hr, stringsAsFactors = FALSE)

  structure(list(study = "Rasmussen et al. (2022)", adjusted = adjusted,
                 gaps = summary_tab, values = by_disease,
                 additional = eco$additional, total = total,
                 comparison = comparison, hr_comparison = hr_comparison,
                 results = res),
            class = "cm_reproduction")
}

#' @rdname reproduce
#' @export
reproduce_rasmussen_2024 <- function(n_draws = 5000, seed = 2024,
                                     inputs = c("analysis", "tables")) {
  inputs <- match.arg(inputs)
  mc <- cm_monte_carlo(global_dairy_sampler(inputs, culling = TRUE), n_draws,
                       method = "published", seed = seed)
  s <- summary(mc, diagnose = FALSE)
  ids <- c("CK", "CM", "DA", "DYS", "LAM", "MET", "MF", "OC", "PTB", "RP", "SCK", "SCM")
  table5 <- list(
    yield = c(0.03, 1.36, 1.18, 3.48, 2.62, 2.87, 0.07, 2.59, 3.37, 2.30, 7.11, 5.58),
    fertility = c(0.34, 6.09, 0.78, 1.11, 1.86, 11.22, 1.06, 9.03, 4.23, 3.74, 0.39, 0.04),
    culling_hr_minus_1 = c(1.18, 1.90, 2.75, 1.18, 1.40, 1.03, 2.64, 1.51, 2.07, 1.29,
                           1.67, 1.25))
  comparison <- do.call(rbind, lapply(names(table5), function(nm) {
    sm <- s[s$analysis == nm, , drop = FALSE]
    sm <- sm[match(ids, sm$disease), , drop = FALSE]
    add <- if (nm == "culling_hr_minus_1") 1 else 0
    data.frame(analysis = nm, disease = ids, table5 = table5[[nm]],
               mean = sm$mean + add, median = sm$q0.5 + add, mcse = sm$mcse,
               stability = sm$stability, stringsAsFactors = FALSE)
  }))
  comparison$analysis[comparison$analysis == "culling_hr_minus_1"] <- "culling_hr"
  rownames(comparison) <- NULL
  adjusted <- adjust_analyses(global_dairy_analyses(inputs, culling = TRUE),
                              method = "published", warn = FALSE)
  structure(list(study = "Rasmussen et al. (2024)", comparison = comparison,
                 central = adjusted, mc = mc, n_draws = n_draws, inputs = inputs),
            class = "cm_reproduction")
}

#' @export
print.cm_reproduction <- function(x, digits = 4, ...) {
  cat(sprintf("<cm_reproduction> %s\n", x$study))
  rnd <- function(d) {
    num <- vapply(d, is.numeric, logical(1))
    d[num] <- lapply(d[num], signif, digits = digits)
    d
  }
  if (!is.null(x$gaps)) {
    cat("\nGaps and values:\n")
    print(rnd(x$gaps), row.names = FALSE)
    cat(sprintf("\nTotal including %s: %.2f\n", paste(names(x$additional), collapse = ", "),
                x$total))
    cat("\nComparison with Table 9:\n")
    print(rnd(x$comparison), row.names = FALSE)
    cat("\nAdjusted culling hazard ratios (eq. 23) vs Table 8:\n")
    print(rnd(x$hr_comparison), row.names = FALSE)
  } else {
    cat(sprintf("Monte Carlo: %d draws, inputs = '%s'\n", x$n_draws, x$inputs))
    cat("\nMeans beside Table 5 (yield and fertility in %, culling as hazard ratios):\n")
    print(rnd(x$comparison), row.names = FALSE)
  }
  invisible(x)
}
