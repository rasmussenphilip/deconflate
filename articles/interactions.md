# Disease combinations, interactions and attribution

``` r

library(deconflate)
```

## The distribution of disease combinations

[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
fits the distribution of all 2^n combinations of n diseases, using
iterative proportional fitting from independence. The result is the
maximum-entropy distribution consistent with the disease probabilities
and the pairwise associations. That is a log-linear model with all
two-way terms and no explicit three-way or higher terms. Three-disease
co-occurrence still differs from independence, through the pairwise
terms.

``` r

j <- fit_joint(example_supplement())
j
#> <cm_joint> 3 diseases, 8 combinations
#>   Converged: TRUE after 9 sweeps (max residual 2.50e-11)
#>   Constrained pairs: 3
head(combination_probs(j))
#>   d1 d2 d3 n_diseases       prob
#> 1  0  0  0          0 0.64262256
#> 2  0  0  1          1 0.13185683
#> 3  0  1  0          1 0.07737744
#> 4  1  0  0          1 0.06410444
#> 5  0  1  1          2 0.04814317
#> 6  1  1  0          2 0.01589556
```

With three or more diseases the constraints do not determine a unique
distribution, which is why the maximum-entropy criterion is needed. If
the pairwise associations cannot hold together, no distribution exists.
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)
detects this and identifies the conflicting associations:

``` r

bad <- cm_model(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
                cm_associations(c("a", "a", "b"), c("b", "c", "c"), c(20, 20, 0.05)))
check_feasibility(bad, method = "triples")
#> <cm_feasibility> NOT jointly feasible [triples]
#> 
#> Conflicting triples:
#>  disease1 disease2 disease3   gap
#>         a        b        c 0.226
```

## Unknown associations are not independence

Setting an association to an odds ratio of 1 imposes independence.
Leaving it unknown lets the maximum-entropy fit imply an association
through the other diseases:

``` r

d <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
a <- cm_associations(c("a", "b"), c("b", "c"), c(3, 2))
imp <- cm_impacts(c("a", "b", "c"), c(0.03, 0.04, 0.05))
ind <- cm_model(d, a, imp, missing_associations = "independent")
unk <- cm_model(d, a, imp, missing_associations = "unknown")
deconflate(ind)$adjusted$adjusted
#> [1] 0.02263264 0.02946945 0.04548702
deconflate(unk, method = "global")$adjusted$adjusted
#> [1] 0.02100000 0.02987262 0.04481272
```

## Impact interactions

An interaction is the extra impact when two diseases occur together:
positive values are synergistic and negative values antagonistic.
Interactions cannot be inferred from associations; they need evidence or
explicit scenarios. They are adjusted with the global method:

``` r

m <- example_supplement()
mi <- set_interaction(m, "d1", "d2", 0.01, outcome = "yield")
res <- deconflate(mi, method = "global")
res$adjusted[, c("disease", "raw", "adjusted")]
#>   disease   raw   adjusted
#> 1      d1 0.025 0.01913880
#> 2      d2 0.050 0.03240641
#> 3      d3 0.075 0.06935621
attribute_burden(res)
#>   outcome disease        main interaction       total      share
#> 1   yield      d1 0.001913880 0.000122397 0.002036277 0.09747208
#> 2   yield      d2 0.004860961 0.000122397 0.004983358 0.23854231
#> 3   yield      d3 0.013871243 0.000000000 0.013871243 0.66398561
```

## Attribution

With additive terms, each disease’s Shapley value is its own term plus
an equal share of every interaction term it is involved in, which is
what
[`attribute_burden()`](https://rasmussenphilip.github.io/deconflate/reference/attribute_burden.md)
reports. For general loss functions,
[`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
computes Shapley values cell by cell over the joint distribution. For
example, with proportional impacts that compound:

``` r

shapley_by_cell(fit_joint(m), loss_multiplicative(c(d1 = 0.02, d2 = 0.04, d3 = 0.06)))
#>   disease     shapley      share
#> 1      d1 0.001978346 0.09981218
#> 2      d2 0.005922273 0.29879258
#> 3      d3 0.011920065 0.60139525
```

## Which interactions would matter?

``` r

screen_interactions(m, outcome = "yield", values = c(-0.01, 0.01))
#> <cm_screen> 6 scenarios; baseline total 0.02109
#>   pair      status      scenario  total    change rel_change max_rank_shift
#>  d2:d3 interaction  delta = 0.01 0.0207 -0.000383   -0.01815              0
#>  d2:d3 interaction delta = -0.01 0.0215  0.000383    0.01815              0
#>  d1:d3 interaction  delta = 0.01 0.0209 -0.000209   -0.00990              0
#>  d1:d3 interaction delta = -0.01 0.0213  0.000209    0.00990              0
#>  d1:d2 interaction  delta = 0.01 0.0209 -0.000201   -0.00955              0
#>  d1:d2 interaction delta = -0.01 0.0213  0.000201    0.00955              0
#>  rank_corr
#>          1
#>          1
#>          1
#>          1
#>          1
#>          1
```
