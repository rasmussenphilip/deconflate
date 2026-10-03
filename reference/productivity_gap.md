# Productivity gaps and their attribution

Generalises Rasmussen et al. (2022), eqs. 17-22. For an outcome with
observed mean `x` and expected loss `L` (see
[`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md)):

- proportion scale: the disease-free value is `x / (1 - L)` for outcomes
  that disease decreases (e.g. yield), or `x / (1 + L)` for outcomes
  that disease increases (e.g. calving interval);

- absolute scale: the disease-free value is `x + L` (decrease) or
  `x - L` (increase), with `x` in the impacts' units (e.g. a culling
  risk of `0.27` with excess-risk impacts from
  [`as_impacts()`](https://rasmussenphilip.github.io/deconflate/reference/as_impacts.md)).

## Usage

``` r
productivity_gap(result, observed)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

- observed:

  Named numeric vector of observed mean values, one per outcome to
  evaluate (names = outcome labels).

## Value

A list with `summary` (observed, disease-free value, gap and `L` per
outcome) and `attribution` (gap attributed to each disease, split into
main and interaction parts).

## Details

The gap is attributed to diseases in proportion to their Shapley shares,
which reduces to eq. 22 when there are no interactions.

## Examples

``` r
res <- deconflate(example_supplement(), method = "published")
productivity_gap(res, c(yield = 10000))
#> $summary
#>   outcome observed disease_free      gap total_loss_fraction
#> 1   yield    10000     10215.66 215.6617          0.02111089
#> 
#> $attribution
#>   outcome disease       gap  gap_main gap_interaction
#> 1   yield      d1  21.09535  21.09535               0
#> 2   yield      d2  56.68610  56.68610               0
#> 3   yield      d3 137.88024 137.88024               0
#> 
```
