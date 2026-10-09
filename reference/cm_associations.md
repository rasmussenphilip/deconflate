# Describe statistical associations between disease pairs

Each row gives one association measure for one pair of diseases. A pair
without a row is unknown: the global model
([`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md))
fills in its association from the others (see
[`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)).
To state that two diseases are unrelated, give the pair an odds ratio of
1.

## Usage

``` r
cm_associations(
  disease1,
  disease2,
  value,
  measure = "OR",
  adjusted = FALSE,
  adjusted_for = NA_character_,
  source = NA_character_
)
```

## Arguments

- disease1, disease2:

  Character vectors of disease ids.

- value:

  Numeric association values (required).

- measure:

  Character, one per row (or recycled). Directional measures treat
  `disease1` as the outcome and `disease2` as the conditioning disease:

  - `"OR"` (default): odds ratio (symmetric). An odds ratio of 1 means
    the two diseases are independent.

  - `"RR"`: risk ratio `P(d1 | d2) / P(d1 | not d2)`.

  - `"RD"`: risk difference `P(d1 | d2) - P(d1 | not d2)`.

  - `"cond_prob"`: conditional probability `P(d1 | d2)`.

  - `"phi"`: binary (phi) correlation coefficient (symmetric).

  For a study's 2x2 table, compute its odds ratio,
  `n11 * n00 / (n10 * n01)`, and enter it as an odds ratio.

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
#>   disease1 disease2 measure value adjusted adjusted_for source
#> 1       d1       d2      OR     2    FALSE         <NA>   <NA>
#> 2       d2       d3      OR     3    FALSE         <NA>   <NA>
# Two diseases stated to be unrelated
cm_associations("d1", "d3", 1)
#>   disease1 disease2 measure value adjusted adjusted_for source
#> 1       d1       d3      OR     1    FALSE         <NA>   <NA>
```
