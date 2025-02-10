# Functions for the preprocessing of SEPAGES data
# 25/09/23
# M. Rolland

#' Fill-in Data Below Limit Of Detection (LOD)
#'
#' This function imputes values below the Limit of Detection (LOD) using the
#' fill-in method by Helsel (1990). It assumes a normal distribution of data.
#'
#' @param var_to_fill A numeric vector of values to be imputed.
#' @param lod A numeric vector of the limit of detection (LOD), the same length as `var_to_fill`.
#' @param loq A numeric vector of the limit of quantification (LOQ), the same length as `var_to_fill`. Defaults to NULL.
#'
#' @details
#' The function identifies censored values (below LOD or LOQ), estimates distribution parameters
#' using \code{NADA::cenros}, and then replaces censored values with random samples from a
#' truncated normal distribution.
#'
#' @return A numeric vector with imputed values replacing those below LOD.
#'
#' @references
#' Helsel, D.R. (1990). Less Than Obvious - Statistical Treatment of
#' Data Below the Detection Limit. Environmental Science & Technology, 24(12), 1766-1774.
#'
#' @export
#' @importFrom NADA cenros mean sd
#' @importFrom msm rtnorm
fill_in <- function(var_to_fill, lod, loq = NULL) {
  # Determine the threshold for imputation (either LOD or LOQ)
  max_lim <- ifelse(!is.null(loq), loq, lod)

  # Identify values below the threshold (censored values)
  censored <- var_to_fill < max_lim
  censored[is.na(censored)] <- FALSE  # Set NA to FALSE to avoid indexing issues

  # Count number of values to be imputed
  n_impute <- sum(censored)

  # Compute distribution parameters using regression on order statistics (ROS)
  stats_ros <- NADA::cenros(var_to_fill, censored, forwardT = NULL)  # No log transform
  dist_mean <- NADA::mean(stats_ros)
  dist_sd   <- NADA::sd(stats_ros)

  # Handle cases where values between LOD and LOQ should also be imputed
  if (!is.null(loq)) {
    n_impute_loq <- sum(dplyr::between(var_to_fill, lod, loq), na.rm = TRUE)
    censored_loq <- dplyr::between(var_to_fill, lod, loq)
    censored_loq[is.na(censored_loq)] <- FALSE

    # Generate fill-in values within the LOQ range
    fill_in_values_loq <- msm::rtnorm(
      n = n_impute_loq, mean = dist_mean,
      sd = dist_sd, lower = lod, upper = loq
    )

    # Exclude LOQ-imputed values from further LOD imputations
    censored <- ifelse(censored_loq, FALSE, censored)
    n_impute <- n_impute - n_impute_loq
  }

  # Generate fill-in values below LOD using truncated normal distribution
  fill_in_values <- msm::rtnorm(
    n = n_impute,
    mean = dist_mean,
    sd = dist_sd,
    lower = rep(-Inf, n_impute),
    upper = lod
  )

  # Replace values below LOD with imputed values
  vec_filled_in <- var_to_fill
  vec_filled_in[censored] <- fill_in_values
  if (!is.null(loq)) vec_filled_in[censored_loq] <- fill_in_values_loq

  return(vec_filled_in)
}

#' Standardize Data on Protocol Variables
#'
#' Standardizes exposure data based on selected protocol variables.
#' This function is designed for use within a `mutate` call in grouped data.
#'
#' @param data A data frame in tidy format containing the data to standardize.
#' @param var_to_std A character string specifying the variable to standardize.
#' @param protocol_vars A character vector of potential protocol variables for standardization.
#' @param covariates A character vector of model covariates. Defaults to NULL.
#' @param folder A character string specifying the folder where model outputs will be saved.
#' @param group A character string indicating the grouping variable for exposure. Defaults to `dplyr::cur_group()`.
#'
#' @details
#' The function selects relevant protocol variables (p < 0.2), builds a model,
#' extracts residuals, and standardizes values accordingly.
#'
#' @return A numeric vector of standardized values.
#'
#' @export
#' @importFrom dplyr mutate cur_group left_join select everything pull
#' @importFrom stringr str_c
#' @importFrom broom tidy
standardise <- function(data = dplyr::pick(everything()), var_to_std,
                        protocol_vars, covariates = NULL, folder,
                        group = dplyr::cur_group()) {

  # Identify protocol variables that influence the variable to standardize
  final_std_vars <- get_protocol_var(
    data,
    var_to_std,
    protocol_vars,
    covariates,
    folder,
    group
  )

  # Construct linear model formula for standardization
  form <- stringr::str_c(
    var_to_std, "~", paste(final_std_vars, collapse = "+"),
    if (!is.null(covariates)) { "+" }, paste(covariates, collapse = "+")
  ) |>
    as.formula()

  # Fit linear model
  lm1 <- lm(form, data = data)

  # Extract model coefficients (betas)
  betas <- broom::tidy(lm1)

  # Apply correction based on protocol variables
  for (prot_var in final_std_vars) {
    class_i <- class(data[[prot_var]])

    if (class_i == "numeric") {
      # Numeric variable: apply median-centered correction
      beta_i <- betas |>
        dplyr::filter(variable == prot_var) |>
        dplyr::pull(estimate)
      data <- data |>
        dplyr::mutate(
          !!stringr::str_c("correct_", prot_var) := beta_i * (.data[[prot_var]] - median(data[[prot_var]], na.rm = TRUE))
        )

    } else if (class_i == "factor") {
      # Factor variable: apply categorical correction
      beta_i <- betas |>
        dplyr::filter(variable == prot_var) |>
        dplyr::select(estimate, level) |>
        dplyr::rename(!!prot_var := level, !!stringr::str_c("correct_", prot_var) := estimate)

      data <- data |>
        dplyr::left_join(beta_i, by = prot_var)

    } else {
      stop(stringr::str_c(
        "Protocol variables must be continuous or factor, not: ", class_i,
        " (variable ", prot_var, ")"
      )
      )
    }
  }

  # Compute standardized values by removing protocol-based corrections
  data <- data |>
    dplyr::mutate(
      val_std = .data[[var_to_std]] - rowSums(dplyr::pick(dplyr::starts_with("correct")), na.rm = TRUE)
    )

  return(data$val_std)
}
