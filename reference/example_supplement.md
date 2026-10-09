# Worked example from the Supplementary File of Rasmussen et al. (2022)

Three hypothetical diseases with prevalences 0.10, 0.15 and 0.20, yield
impacts of 2.5%, 5% and 7.5% (entered in percent, so results are in
percent), and odds ratios of 2 (d1:d2), 1 (d1:d3) and 3 (d2:d3). The
observed mean yield in the example is 10,000 units per animal per year.

## Usage

``` r
example_supplement()
```

## Value

A
[`cm_model()`](https://rasmussenphilip.github.io/deconflate/reference/cm_model.md).

## Details

Note: the published productivity gap (10,225; 20, 61 and 143 units)
rounds the adjusted impacts to 2%, 4% and 7% before computing the gap.
Without rounding, the published method gives 10,215.7 and 21.1, 56.7 and
137.9 units (see
[`vignette("reproducing-published")`](https://rasmussenphilip.github.io/deconflate/articles/reproducing-published.md)).

## Examples

``` r
deconflate(example_supplement())
#> <cm_result> method: simultaneous; yield [%]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    2.143 -0.14270
#>       d2 5.0    3.387 -0.32258
#>       d3 7.5    6.934 -0.07544
#> 
#> Raw sum: 2.5; adjusted total: 2.109
#> Diagnostics: residual 8.88e-16, condition number 1.53, sign changes 0
#> 
#> Notes:
#> * No input has a distribution, so no draws were run: the results are point
#>   estimates.
compare_methods(example_supplement())
#> <cm_comparison> methods: published, simultaneous, global
#> Impacts: yield
#> Units: %
#> 
#> Adjusted values:
#>  disease raw published simultaneous global
#>       d1 2.5      2.07         2.14   2.14
#>       d2 5.0      3.70         3.39   3.39
#>       d3 7.5      6.75         6.93   6.93
#> 
#> Totals:
#>        method raw_sum adjusted_total
#>     published     2.5          2.111
#>  simultaneous     2.5          2.109
#>        global     2.5          2.109
```
