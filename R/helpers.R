# helper functions
# M. Rolland
# 20-03-2025

#' Create a Standardized Table of Significant Variables from AOV Results
#'
#' This function scans a directory for files ending in `_aov.csv`, extracts
#' variables with a p-value below 0.2, and creates a summary table indicating
#' which variables were found to be significant.
#'
#' @param path Character string. Path to the folder containing AOV result files.
#' @param lst_prot_vars Character vector. List of variable names to check in the AOV results.
#'
#' @return A tibble with exposures as rows and variables as columns, where `"X"`
#' indicates a significant association (p < 0.2).
#'
#' @import dplyr tidyr stringr rio
#' @export
#'
#' @examples
#' \dontrun{
#' # Define the path to the folder with AOV results
#' path_to_results <- "path/to/aov/files/"
#'
#' # Define the list of variables to check
#' prot_vars <- c("var1", "var2", "var3")
#'
#' # Generate the summary table
#' mk_tbl_std(path_to_results, prot_vars)
#' }
mk_tbl_std <- function(path, lst_prot_vars){
  # List all files in the folder
  files <- list.files(path, full.names = TRUE, pattern = "_aov\\.csv$")  # Adjust pattern if needed

  lst_exp <- str_remove_all(files, path) |>
    str_remove("_aov.csv")

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
    exp <- file |> str_remove(path) |> str_remove("_aov.csv") |> str_remove("/")

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
