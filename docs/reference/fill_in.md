# Fill-in Data Below Limit Of Detection (LOD)

This function fills in data below LOD using the fill-in method by Helsel
(1990). It is designed to work with data that should have a normal
distribution.

## Usage

``` r
fill_in(var_to_fill, lod, loq = NULL)
```

## Arguments

- var_to_fill:

  A numeric vector of values to fill.

- lod:

  A numeric vector of the limit of detection (LOD), should be of the
  same length as `var_to_fill`.

- loq:

  A numeric vector of the limit of quantification (LOQ), should be of
  the same length as `var_to_fill`.

## Value

A numeric vector with the original values and where those are below LOD
replaced by fill-in values below the limit of detection.

## Details

The function flags values that are to be imputed and computes
distribution parameters using the
[`NADA::cenros`](https://rdrr.io/pkg/NADA/man/ros.html). It then
computes fill-in values by performing a random sample between 0 and
`lod` for a normal distribution with the previously computed parameters.
Finally, the function replaces values below `lod` by fill-in values.

## References

Helsel, D.R. (1990). Less Than Obvious - Statistical Treatment of Data
Below the Detection Limit. Environmental Science & Technology, 24(12),
1766-1774.

## Examples

``` r
# Example dataset
set.seed(113)
values <- c(0.05, 0.1, 0.2, NA, 0.4, 0.5)
lod <- rep(0.1, length(values))

# Apply fill_in function
imputed_values <- fill_in(values, lod)
print(imputed_values)
#> [1] 0.07635231 0.10000000 0.20000000         NA 0.40000000 0.50000000
```
