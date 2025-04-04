# helper functions
# M. Rolland
# 20-03-2025

#' Create a Standardized Table of Significant Variables from AOV Results
#'
#' This function scans a directory for AOV result files (in `.xlsx` format, ending with `_aov.xlsx`),
#' extracts protocol variables with a p-value < 0.2, and creates a wide-format summary table.
#' Each row corresponds to an exposure and each column to a protocol variable.
#'
#' @param path Character string. Path to the folder containing AOV result files in `.xlsx` format.
#' @param lst_prot_vars Character vector. List of protocol variables to check for significance.
#'
#' @return A tibble in wide format with one row per exposure and one column per protocol variable.
#'   The value `"X"` indicates that the variable was significantly associated (p < 0.2) with the exposure.
#'   Empty cells mean non-significant or missing data.
#'
#' @import dplyr tidyr stringr rio
#' @export
#'
#' @examples
#' \dontrun{
#' path_to_results <- "outputs/aov_results/"
#' prot_vars <- c("batch", "season", "storage_time")
#' mk_tbl_std(path_to_results, prot_vars)
#' }

mk_tbl_std <- function(path, lst_prot_vars){
  # List all files in the folder
  files <- list.files(path, full.names = TRUE, pattern = "_aov\\.xlsx$")  # Adjust pattern if needed

  lst_exp <- str_remove_all(files, path) |>
    str_remove("_aov.xlsx") |>
    str_remove_all("/")

  # Initialize an empty list to store results
  results_df <- expand_grid(
    exposure = lst_exp,
    variable = lst_prot_vars
  ) |>
    mutate(value = "")

  # Loop through each file
  for (file in files) {
    # Read the AOV results
    df <- import(file)

    # Extract the exposure name from the filename
    exp <- file |> str_remove(path) |> str_remove("_aov.xlsx") |> str_remove_all("/")

    # Filter variables with p < 0.2 (excluding "Residuals")
    prot_vars <- df |>
      filter(term %in% lst_prot_vars) |>
      filter(p.value < 0.2) |>
      pull(term)

    results_df$value[results_df$exposure == exp & results_df$variable %in% prot_vars] <- "X"
  }

  results_df <- results_df  |>
    pivot_wider(
      names_from = variable,
      values_from = value,
      id_cols = exposure
    )

  return(results_df)
}
