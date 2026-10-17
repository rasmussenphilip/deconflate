# Groups of diseases reported together (e.g. the subclinical and clinical
# forms of one disease). The diseases stay separate in the model; a group's
# contribution is the sum of its diseases' contributions (an interaction
# between two diseases of a group is then counted wholly in the group), and
# its interval comes from the sum within each draw.

# A checked list of groups: named, each a set of disease ids of the model,
# no disease in two groups, no group named after a disease.
check_groups <- function(groups, ids) {
  if (is.null(groups)) return(NULL)
  bad <- function(msg) cm_abort(msg, class = "deconflate_invalid_input")
  if (!is.list(groups) || !length(groups)) {
    bad('`groups` must be a named list of disease ids, e.g. list(ketosis = c("SCK", "CK")).')
  }
  nm <- names(groups) %||% rep("", length(groups))
  if (any(is.na(nm) | !nzchar(trimws(nm)))) bad("Every group in `groups` needs a name.")
  nm <- trimws(nm)
  if (anyDuplicated(nm)) bad(sprintf("Group names must be unique (check: %s).", paste(unique(nm[duplicated(nm)]), collapse = ", ")))
  clash <- intersect(nm, ids)
  if (length(clash)) {
    bad(sprintf("A group cannot have the name of a disease (check: %s).", paste(clash, collapse = ", ")))
  }
  groups <- lapply(groups, function(g) {
    if (!is.character(g) || !length(g) || anyNA(g)) bad("Each group in `groups` must be a character vector of disease ids.")
    unique(trimws(g))
  })
  names(groups) <- nm
  members <- unlist(groups, use.names = FALSE)
  unknown <- setdiff(members, ids)
  if (length(unknown)) {
    bad(sprintf("Groups name diseases that are not in the model: %s (diseases: %s).",
                paste(unknown, collapse = ", "), paste(ids, collapse = ", ")))
  }
  twice <- unique(members[duplicated(members)])
  if (length(twice)) bad(sprintf("A disease can be in one group only (check: %s).", paste(twice, collapse = ", ")))
  groups
}

# The units reported: each group, then each disease not in a group (so that
# the rows add up to the sum of all contributions).
group_units <- function(groups, ids) {
  rest <- setdiff(ids, unlist(groups, use.names = FALSE))
  c(groups, stats::setNames(as.list(rest), rest))
}

# Group results of a deconflate() result: `table` (one row per group and per
# disease not in a group) and, with draws, `ranks`.
group_results <- function(res, groups) {
  ct <- res$contributions
  if (is.null(ct)) return(NULL)
  units <- group_units(groups, ct$disease)
  add <- function(col) {
    v <- ct[[col]]
    if (is.null(v)) return(rep(NA_real_, length(units)))
    vapply(units, function(u) sum(v[match(u, ct$disease)]), numeric(1))
  }
  tab <- data.frame(group = names(units), diseases = vapply(units, paste, character(1), collapse = ", "),
                    raw_contribution = unname(add("raw_contribution")), contribution = unname(add("contribution")),
                    lower = NA_real_, upper = NA_real_, main = unname(add("main")),
                    interaction = unname(add("interaction")), share = unname(add("share")),
                    stringsAsFactors = FALSE)
  ranks <- NULL
  v <- res$draws$values
  if (!is.null(v)) {
    G <- group_draw_values(v, units)
    if (!is.null(G)) {
      s <- summarise_draws(list(values = G, block = res$draws$block, n_blocks = res$draws$n_blocks))
      if (!is.null(s)) {
        tab$lower <- s$lower[match(names(units), s$quantity)]
        tab$upper <- s$upper[match(names(units), s$quantity)]
      }
      ranks <- rank_table(G, names(units))
      if (!is.null(ranks)) names(ranks)[1] <- "group"
    }
  }
  rownames(tab) <- NULL
  list(table = tab, ranks = ranks)
}

# The contribution of each group in each draw (columns named by group).
group_draw_values <- function(v, units) {
  cols <- lapply(units, function(u) paste0("contribution:", u))
  if (!all(unlist(cols) %in% colnames(v))) return(NULL)
  G <- vapply(cols, function(cc) rowSums(v[, cc, drop = FALSE]), numeric(nrow(v)))
  G <- matrix(G, nrow = nrow(v), dimnames = list(NULL, names(units)))
  G
}

# The group of each disease (NA when it is in none).
disease_group <- function(ids, groups) {
  out <- rep(NA_character_, length(ids))
  for (g in names(groups)) out[ids %in% groups[[g]]] <- g
  out
}
