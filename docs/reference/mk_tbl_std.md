# Create a Standardized Table of Significant Variables from AOV Results

This function scans a directory for AOV result files (in `.xlsx` format,
ending with `_aov.xlsx`), extracts protocol variables with a p-value \<
0.2, and creates a wide-format summary table. Each row corresponds to an
exposure and each column to a protocol variable.

## Usage

``` r
mk_tbl_std(path, lst_prot_vars)
```

## Arguments

- path:

  Character string. Path to the folder containing AOV result files in
  `.xlsx` format.

- lst_prot_vars:

  Character vector. List of protocol variables to check for
  significance.

## Value

A tibble in wide format with one row per exposure and one column per
protocol variable. The value `"X"` indicates that the variable was
significantly associated (p \< 0.2) with the exposure. Empty cells mean
non-significant or missing data.

## Examples

``` r
if (FALSE) { # \dontrun{
path_to_results <- "outputs/aov_results/"
prot_vars <- c("batch", "season", "storage_time")
mk_tbl_std(path_to_results, prot_vars)
} # }
```
