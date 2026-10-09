# Disease combinations, interactions and attribution

``` r

library(deconflate)
```

## The distribution of disease combinations

[`fit_joint()`](https://rasmussenphilip.github.io/deconflate/reference/fit_joint.md)
fits the distribution of all 2^n combinations of n diseases, using
iterative proportional fitting from independence. The result is the
maximum-entropy distribution consistent with the disease probabilities
and the pairwise associations.

``` r

m <- example_supplement()
j <- fit_joint(m)
j
#> <cm_joint> 3 diseases, 8 combinations (exact backend)
#>   Converged: TRUE after 8 sweeps (max residual 2.48e-11)
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

Enumerating 2^n combinations limits this exact fit to about 20 diseases.
For more, `fit_joint(backend = "sampled")` fits the same model by Monte
Carlo (Gibbs sampling with calibrated parameters, followed by raking to
the pairwise tables); the functions below that use the joint
distribution accept a fit from either backend. See
[`vignette("thresholds-and-scaling")`](https://rasmussenphilip.github.io/deconflate/articles/thresholds-and-scaling.md)
for how it works, its diagnostics and a 24-disease example.

## What maximum entropy assumes

With three or more diseases, the probabilities and pairwise associations
do not determine a unique distribution: there are 2^n - 1 free
probabilities and n(n + 1) / 2 constraints. Maximum entropy picks one
solution. It is an assumption, not something the pairwise evidence
identifies:

- the fitted distribution is a log-linear model with all two-way terms
  and no three-way or higher terms. The odds ratio of a pair is the same
  whether or not a third disease is present (given the other diseases);
- three diseases still occur together more (or less) often than under
  independence, through the pairwise terms;
- a pair whose association is unknown has no term of its own. Its
  association is implied by the associations both diseases have with
  other diseases.

Additive results without interactions do not depend on this assumption:
they need the pairwise tables only. It matters for interactions, for
three-way scenarios, and for the results that use the whole distribution
(the snapshot hazard-ratio model and
[`attributable_risk()`](https://rasmussenphilip.github.io/deconflate/reference/attributable_risk.md),
see
[`vignette("culling-hazard-ratios")`](https://rasmussenphilip.github.io/deconflate/articles/culling-hazard-ratios.md)).

If the pairwise associations cannot hold together, no distribution
exists.
[`check_feasibility()`](https://rasmussenphilip.github.io/deconflate/reference/check_feasibility.md)
detects this and identifies the conflicting associations:

``` r

bad <- cm_population(cm_diseases(c("a", "b", "c"), c(0.5, 0.5, 0.5)),
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
through the other diseases. Use the simultaneous method when every pair
has an association estimate or a defensible independence assumption;
with unknown pairs (or interactions, below), use the global method, the
only one that accepts them.
[`compare_methods()`](https://rasmussenphilip.github.io/deconflate/reference/compare_methods.md)
reports why the others fail (the `"published"` approximation is listed
for comparison only):

``` r

d <- cm_diseases(c("a", "b", "c"), c(0.2, 0.3, 0.25))
a <- cm_associations(c("a", "b"), c("b", "c"), c(3, 2))
imp <- cm_impacts(c("a", "b", "c"), c(3, 4, 5), units = "%")
ind <- cm_model(d, imp, associations = a, missing_associations = "independent")
unk <- cm_model(d, imp, associations = a, missing_associations = "unknown")
deconflate(ind)$adjusted$adjusted
#> [1] 2.263264 2.946945 4.548702
compare_methods(unk)
#> <cm_comparison> methods: global
#> Units: %
#> 
#> Adjusted values:
#>  disease raw global
#>        a   3   2.10
#>        b   4   2.99
#>        c   5   4.48
#> 
#> Totals:
#>  method raw_sum adjusted_total
#>  global    3.05          2.436
#> 
#> Not run:
#>   published: Unknown associations for a-c. Specify them, set missing_associations = 'independent', or use method = 'global'.
#>   simultaneous: Unknown associations for a-c. Specify them, set missing_associations = 'independent', or use method = 'global'.
```

The odds ratio of `a` and `c` implied by the fitted distribution:

``` r

ju <- fit_joint(unk)
p_ac <- sum(ju$prob[ju$cells[, "a"] == 1 & ju$cells[, "c"] == 1])
joint_to_or(p_ac, 0.2, 0.25)
#> [1] 1.194529
```

## Impact interactions

An interaction is the extra impact when two diseases occur together, in
the units of the impact vector: positive values are synergistic and
negative values antagonistic. Interactions cannot be inferred from
associations; they need evidence or explicit scenarios.

The raw estimate of disease `i` then includes, besides the conflated
main impacts, the part of the interaction burden that goes with `i`. Its
expected value depends on how often pairs occur together with each
disease, that is, on probabilities of disease triples, so interactions
need the joint distribution and the global method:

``` r

mi <- set_interaction(m, "d1", "d2", 1)   # 1 percentage point
tryCatch(deconflate(mi), deconflate_unsupported = function(e) conditionMessage(e))
#> [1] "Interactions require method = 'global' (they need probabilities of disease triples)."
res <- deconflate(mi, method = "global")
res
#> <cm_result> method: global; yield [%]
#> 
#>  disease raw adjusted   change
#>       d1 2.5    1.914 -0.23445
#>       d2 5.0    3.241 -0.35187
#>       d3 7.5    6.936 -0.07525
#> 
#> Raw sum: 2.5; adjusted total: 2.089 (interactions: 0.02448)
#> Diagnostics: residual 0.00e+00, condition number 1.53, sign changes 0
```

The raw impacts satisfy `raw = A %*% adjusted + offset`, where `offset`
holds the interaction part:

``` r

res$conflation$offset
#>         d1         d2         d3 
#> 0.24479393 0.16319595 0.02304974
as.vector(res$conflation$A %*% res$adjusted$adjusted + res$conflation$offset)
#> [1] 2.5 5.0 7.5
res$adjusted$raw
#> [1] 2.5 5.0 7.5
```

A simulated herd with known impacts (2, 4 and 6 percent) and
interactions (1 point for d1 with d2, 1.5 for d2 with d3) checks the
method:

``` r

int <- cm_interactions(c("d1", "d2"), c("d2", "d3"), c(1, 1.5))
raw <- simulate_raw_impacts(m, c(d1 = 2, d2 = 4, d3 = 6), interactions = int)
raw$value
#> [1] 2.714619 6.136904 7.116677
deconflate(cm_model(m, raw, int), method = "global")$adjusted
#>   disease      raw adjusted     change estimand adjusted_for
#> 1      d1 2.714619        2 -0.2632484    crude         <NA>
#> 2      d2 6.136904        4 -0.3482055    crude         <NA>
#> 3      d3 7.116677        6 -0.1569099    crude         <NA>
```

In files, interactions go in `interactions_<analysis>.csv` (see
[`vignette("own-data")`](https://rasmussenphilip.github.io/deconflate/articles/own-data.md)).

## Attribution

The expected aggregate is
`sum_i p_i b_i + sum_{j<k} P(j and k) delta_jk`. Removing any disease
involved in a term removes that term, so each disease’s Shapley value is
its own term plus an equal share of every interaction term it is
involved in. This is the `contributions` table:

``` r

attribute_burden(res)
#>   disease      main interaction     total      share
#> 1      d1 0.1913880   0.0122397 0.2036277 0.09747208
#> 2      d2 0.4860961   0.0122397 0.4983358 0.23854231
#> 3      d3 1.3871243   0.0000000 1.3871243 0.66398561
```

For general loss functions,
[`shapley_by_cell()`](https://rasmussenphilip.github.io/deconflate/reference/shapley_by_cell.md)
computes Shapley values cell by cell over the joint distribution. For
example, with proportional impacts that compound (here as proportions):

``` r

shapley_by_cell(j, loss_multiplicative(c(d1 = 0.02, d2 = 0.04, d3 = 0.06)))
#>   disease     shapley      share
#> 1      d1 0.001978346 0.09981218
#> 2      d2 0.005922273 0.29879258
#> 3      d3 0.011920065 0.60139525
```

## Three-way scenarios

Pairwise evidence cannot identify how often three diseases occur
together.
[`cm_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/cm_three_way.md)
(or
[`set_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/set_three_way.md))
sets a three-way term, the ratio of conditional odds ratios
`OR(d1, d2 | d3 present) / OR(d1, d2 | d3 absent)`. A ratio of 1 is the
maximum-entropy assumption. All pairwise associations are still matched:

``` r

m3 <- set_three_way(m, "d1", "d2", "d3", ratio = 2)
m3$three_way
#>   disease1 disease2 disease3 ratio   source
#> 1       d1       d2       d3     2 scenario
j3 <- fit_joint(m3)
c(maximum_entropy = sum(j$prob[rowSums(j$cells) == 3]),
  ratio_2 = sum(j3$prob[rowSums(j3$cells) == 3]))
#> maximum_entropy         ratio_2 
#>     0.008583836     0.010703834
```

Three-way terms are sensitivity scenarios. When all pairs of a triple
are constrained, they keep the pairwise tables fixed and change additive
results only through interactions: without interactions every scenario
reproduces the baseline, and with an interaction they change the offset.
When a pair is unknown, its fitted table moves with the three-way term,
so additive results can change even without interactions.
[`screen_three_way()`](https://rasmussenphilip.github.io/deconflate/reference/screen_three_way.md)
sets each triple in turn:

``` r

screen_three_way(m, ratios = c(0.5, 2))
#> <cm_screen> 2 scenarios; baseline total 2.109
#>      pair    status    scenario total    change rel_change max_rank_shift
#>  d1:d2:d3 three-way ratio = 0.5  2.11 -8.55e-11  -4.06e-11              0
#>  d1:d2:d3 three-way   ratio = 2  2.11  5.96e-11   2.82e-11              0
#>  rank_corr failed
#>          1   <NA>
#>          1   <NA>
screen_three_way(mi, ratios = c(0.5, 2))
#> <cm_screen> 2 scenarios; baseline total 2.089
#>      pair    status    scenario total   change rel_change max_rank_shift
#>  d1:d2:d3 three-way   ratio = 2  2.09 -0.00234   -0.00112              0
#>  d1:d2:d3 three-way ratio = 0.5  2.09  0.00231    0.00111              0
#>  rank_corr failed
#>          1   <NA>
#>          1   <NA>
```

They also change the joint-based results, such as snapshot hazard
ratios:

``` r

hr <- cm_hazard_ratios(c("d1", "d2", "d3"), c(1.5, 2.0, 1.3), estimand = "snapshot_crude")
rbind(maximum_entropy = deconflate_hr(cm_hr_model(m, hr))$adjusted$adjusted,
      ratio_2 = deconflate_hr(cm_hr_model(m3, hr))$adjusted$adjusted)
#>                     [,1]     [,2]     [,3]
#> maximum_entropy 1.383153 1.890932 1.143987
#> ratio_2         1.379893 1.891328 1.139917
```

In files, three-way terms go in `three_way.csv`.

## Which interactions would matter?

[`screen_interactions()`](https://rasmussenphilip.github.io/deconflate/reference/screen_interactions.md)
adds an interaction of each size to each pair in turn and reports the
change in the aggregate and in the ranking of diseases:

``` r

screen_interactions(m, values = c(-1, 1))
#> <cm_screen> 6 scenarios; baseline total 2.109
#>   pair      status   scenario total  change rel_change max_rank_shift rank_corr
#>  d2:d3 interaction delta = -1  2.15  0.0383    0.01815              0         1
#>  d2:d3 interaction  delta = 1  2.07 -0.0383   -0.01815              0         1
#>  d1:d3 interaction delta = -1  2.13  0.0209    0.00990              0         1
#>  d1:d3 interaction  delta = 1  2.09 -0.0209   -0.00990              0         1
#>  d1:d2 interaction delta = -1  2.13  0.0201    0.00955              0         1
#>  d1:d2 interaction  delta = 1  2.09 -0.0201   -0.00955              0         1
#>  failed
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
#>    <NA>
```
