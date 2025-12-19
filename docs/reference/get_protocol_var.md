# Get Protocol Variables for Standardisation

This function identifies which among the candidate protocol variables
are significantly associated with the exposure to be standardised. The
selection is based on ANOVA p-values (\< 0.2). Model outputs are
exported to files for transparency and quality control.

## Usage

``` r
get_protocol_var(
  data,
  var_to_std,
  protocol_vars,
  covariates,
  folder,
  group,
  export_check_model,
  p_thresh
)
```

## Arguments

- data:

  A data.frame containing the exposure, protocol variables, and optional
  covariates.

- var_to_std:

  A character string giving the name of the variable to be standardised.

- protocol_vars:

  A character vector of names of potential protocol variables.

- covariates:

  A character vector of additional covariates to include in the model
  (not corrected for).

- folder:

  A character string giving the directory where outputs (model summary,
  ANOVA table, optional model check) will be saved.

- group:

  A character string or vector used for naming the output files
  (typically the exposure variable name).

- export_check_model:

  Logical. If `TRUE`, exports model diagnostic plots from
  [`performance::check_model()`](https://easystats.github.io/performance/reference/check_model.html).
  Default is `FALSE`.

- p_thresh:

  Numeric. The p-value threshold used to select protocol variables for
  adjustment. Default is `0.05`.

## Value

A character vector containing the names of the protocol variables with
p-value \< 0.2, to be used for standardisation.

## See also

[`Anova`](https://rdrr.io/pkg/car/man/Anova.html),
[`tidy`](https://generics.r-lib.org/reference/tidy.html),
[`check_model`](https://easystats.github.io/performance/reference/check_model.html),
[`export`](http://gesistsa.github.io/rio/reference/export.md)
