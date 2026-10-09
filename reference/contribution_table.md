# Contribution table

One row per disease with the raw estimate, the adjusted value and the
disease's contribution to the total, with 95% intervals when the result
has draws.

## Usage

``` r
contribution_table(result)
```

## Arguments

- result:

  A
  [`deconflate()`](https://rasmussenphilip.github.io/deconflate/reference/deconflate.md)
  result.

## Value

A data frame.

## Details

- Additive impacts: the contribution to the aggregate (its Shapley
  share, including any interaction share), in the units of the impacts.

- Event impacts: the adjusted hazard ratio and the disease's share of
  the risk attributable to disease (Shapley allocation), as a proportion
  of animals.

## Examples

``` r
contribution_table(deconflate(example_supplement()))
#>   disease raw adjusted      change contribution      main interaction     share
#> 1      d1 2.5 2.143250 -0.14269981    0.2143250 0.2143250           0 0.1016130
#> 2      d2 5.0 3.387080 -0.32258408    0.5080619 0.5080619           0 0.2408757
#> 3      d3 7.5 6.934209 -0.07543874    1.3868419 1.3868419           0 0.6575113
```
