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
137.9 units.

## Examples

``` r
deconflate(example_supplement(), method = "published")
#> <cm_result> method: published; yield [%]
#> 
#>  disease raw adjusted  change
#>       d1 2.5    2.065 -0.1740
#>       d2 5.0    3.699 -0.2601
#>       d3 7.5    6.748 -0.1002
#> 
#> Raw sum: 2.5; adjusted total: 2.111
#> Diagnostics: residual 2.67e-01, condition number 1.53, sign changes 0
```
