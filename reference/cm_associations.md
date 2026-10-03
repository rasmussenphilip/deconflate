# Describe statistical associations between disease pairs

Each row gives one association measure for one pair. Pairs not listed
are handled by `missing_associations` in
[`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md).

## Usage

``` r
cm_associations(
  disease1,
  disease2,
  value = NA_real_,
  measure = "OR",
  n11 = NA_real_,
  n10 = NA_real_,
  n01 = NA_real_,
  n00 = NA_real_,
  zero_cell = c("haldane", "error"),
  adjusted = FALSE,
  adjusted_for = NA_character_,
  source = NA_character_
)
```

## Arguments

- disease1, disease2:

  Character vectors of disease ids.

- value:

  Numeric association values (ignored for `"independent"`, `"unknown"`
  and `"table"`).

- measure:

  Character, one per row (or recycled). Directional measures treat
  `disease1` as the outcome and `disease2` as the conditioning disease:

  - `"OR"`: odds ratio (symmetric).

  - `"RR"`: risk ratio `P(d1 | d2) / P(d1 | not d2)`.

  - `"RD"`: risk difference `P(d1 | d2) - P(d1 | not d2)`.

  - `"cond_prob"`: conditional probability `P(d1 | d2)`.

  - `"phi"`: binary (phi) correlation coefficient (symmetric).

  - `"table"`: a study contingency table given in `n11`, `n10`, `n01`
    and `n00`. The table's odds ratio is used (transported to the
    modelled marginals), because the study's marginal frequencies
    generally differ from the modelled population's.

  - `"independent"`: independence is imposed (odds ratio 1).

  - `"unknown"`: the association is unknown. The pairwise methods reject
    unknown pairs; the global method leaves them unconstrained, so their
    association is implied by the maximum-entropy fit.

- n11, n10, n01, n00:

  Counts for `measure = "table"`: both diseases, `d1` only, `d2` only
  and neither.

- zero_cell:

  For tables with a zero cell: `"haldane"` (default) adds 0.5 to every
  cell (Haldane-Anscombe correction) and records it in column
  `corrected`; `"error"` rejects such tables.

- adjusted:

  Logical: was the measure adjusted for covariates (e.g. an odds ratio
  from multivariable logistic regression)? The 2x2 algebra needs
  marginal (crude) measures, so
  [`cm_population()`](https://rasmussenphilip.github.io/deconflate/reference/cm_population.md)
  rejects adjusted measures unless told to use them as marginal.

- adjusted_for:

  Optional description of the adjustment set.

- source:

  Optional citation.

## Value

A `cm_associations` data frame.

## Details

Odds ratios (and the other measures) are applied to the modelled
population's own marginal probabilities ("transported"). This assumes
the measure is the same in the source study and in the modelled
population.

## Examples

``` r
cm_associations(c("d1", "d2"), c("d2", "d3"), c(2, 3), measure = "OR")
#>   disease1 disease2 measure value adjusted adjusted_for source corrected n11
#> 1       d1       d2      OR     2    FALSE         <NA>   <NA>     FALSE  NA
#> 2       d2       d3      OR     3    FALSE         <NA>   <NA>     FALSE  NA
#>   n10 n01 n00
#> 1  NA  NA  NA
#> 2  NA  NA  NA
```
