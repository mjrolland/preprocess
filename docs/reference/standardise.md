# Apply Standardisation on Protocol Variables

This function standardises the exposure data on selected protocol
variables (either categorical or continuous) as per the SEPAGES pipeline
guide. The function is designed to be called within a `mutate` call on
grouped data, where each group represents a different exposure.

## Usage

``` r
standardise(
  data = dplyr::pick(everything()),
  var_to_std,
  protocol_vars,
  covariates = NULL,
  folder,
  group = dplyr::cur_group(),
  force = NULL,
  export_check_model = FALSE,
  p_thresh = 0.05
)
```

## Arguments

- data:

  A data.frame in tidy format containing the data to standardise, i.e.,
  one row per ID and per exposure.

- var_to_std:

  A character string representing the variable to standardise. This
  variable should be normally distributed.

- protocol_vars:

  A character vector of potential protocol variables on which to
  standardise.

- covariates:

  A character vector of names of model covariates.

- folder:

  A character string representing the folder where regression outputs
  will be saved.

- group:

  A character string or vector representing the exposure group. Used for
  naming output files and tracking grouping structure. Defaults to the
  current grouping (via
  [`dplyr::cur_group()`](https://dplyr.tidyverse.org/reference/context.html)).

- force:

  A character vector of protocol variables to forcibly include in the
  model regardless of statistical significance. Default is `NULL`.

- export_check_model:

  Logical. If `TRUE`, exports model diagnostic plots from
  [`performance::check_model()`](https://easystats.github.io/performance/reference/check_model.html).
  Default is `FALSE`.

- p_thresh:

  Numeric. The p-value threshold used to select protocol variables for
  adjustment. Default is `0.05`.

## Value

A numeric vector of corrected values after standardisation. Length
matches the number of rows in `data`; values for which residuals cannot
be computed are returned as `NA`.

## Details

The function selects protocol variables with p-value \< 0.2 in ANOVA to
adjust for protocol effects. A linear model is fitted including selected
protocol variables and covariates. Residual-based corrections are
computed and applied to the exposure variable. If no protocol variable
is selected, the original variable is returned unmodified. Forced
variables are always included.

## See also

[`mutate`](https://dplyr.tidyverse.org/reference/mutate.html),
[`lm`](https://rdrr.io/r/stats/lm.html),
[`predict`](https://rdrr.io/r/stats/predict.html),
[`check_model`](https://easystats.github.io/performance/reference/check_model.html),
[`Anova`](https://rdrr.io/pkg/car/man/Anova.html)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- data.frame(
  exposure = rnorm(100, mean = 5, sd = 2),
  protocol_var1 = sample(c("A", "B", "C"), 100, replace = TRUE),
  protocol_var2 = runif(100, 0, 1),
  age = rnorm(100, mean = 35, sd = 10)
)

data |>
  dplyr::mutate(
    exposure_std = standardise(
      var_to_std = "exposure",
      protocol_vars = c("protocol_var1", "protocol_var2"),
      covariates = "age",
      folder = "outputs/",
      export_check_model = TRUE
    )
  )
} # }
```
