# Describe the diseases in a system

Describe the diseases in a system

## Usage

``` r
cm_diseases(
  id,
  value,
  type = "prevalence",
  time_horizon = NA_character_,
  reference_population = NA_character_,
  source = NA_character_
)
```

## Arguments

- id:

  Character vector of unique disease identifiers (no `|`, `;` or `:`;
  `"all"` is reserved).

- value:

  Numeric vector of disease occurrence values. Interpreted according to
  `type`.

- type:

  Character, one per disease (or recycled):

  - `"prevalence"`: proportion of animals affected at a representative
    point in time (used as the probability directly).

  - `"probability"`: probability of at least one case within
    `time_horizon` (used directly).

  - `"incidence_rate"`: expected number of cases per animal within
    `time_horizon` (e.g. `0.48` for 48 cases per 100 lactations).
    Converted to a probability assuming Poisson-distributed events,
    `1 - exp(-value)`, as in Rasmussen et al. (2024), eq. 1.

- time_horizon:

  Character label of the period the probabilities refer to (e.g.
  `"lactation"`, `"year"`). Probabilities, associations and impacts
  should all refer to the same period.

- reference_population:

  Optional description of the population the values come from.

- source:

  Optional citation.

## Value

A `cm_diseases` data frame with the converted probability in column
`prob` and the inputs retained as metadata.

## Examples

``` r
cm_diseases(c("SCK", "PTB"), c(0.4789, 0.1001),
            type = c("incidence_rate", "prevalence"),
            time_horizon = "lactation")
#>    id      prob  value           type time_horizon reference_population source
#> 1 SCK 0.3805356 0.4789 incidence_rate    lactation                 <NA>   <NA>
#> 2 PTB 0.1001000 0.1001     prevalence    lactation                 <NA>   <NA>
```
