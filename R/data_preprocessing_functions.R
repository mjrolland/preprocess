# Functions for the preprocessing of SEPAGES data
# 25/09/23
# M. Rolland

#' Fill-in Data Below Limit Of Detection (LOD)
#'
#' This function fills in data below LOD using the fill-in method by Helsel
#' (1990). It is designed to work with data that should have a normal
#' distribution.
#'
#' @param var_to_fill A numeric vector of values to fill.
#' @param lod A numeric vector of the limit of detection (LOD), should be of the
#'   same length as `var_to_fill`.
#' @param loq A numeric vector of the limit of quantification (LOQ), should be of the
#'   same length as `var_to_fill`.
#' @details The function flags values that are to be imputed and computes
#'   distribution parameters using the \code{NADA::cenros}. It then computes
#'   fill-in values by performing a random sample between 0 and \code{lod} for a
#'   normal distribution with the previously computed parameters. Finally, the
#'   function replaces values below \code{lod} by fill-in values.
#'
#' @return A numeric vector with the original values and where those are below
#'   LOD replaced by fill-in values below the limit of detection.
#'
#' @references Helsel, D.R. (1990). Less Than Obvious - Statistical Treatment of
#'   Data Below the Detection Limit. Environmental Science & Technology, 24(12),
#'   1766-1774.
#'
#' @export
#'
#' @importFrom NADA cenros mean sd
#' @importFrom msm rtnorm
#'
#' @examples
#' # Example dataset
#' set.seed(113)
#' values <- c(0.05, 0.1, 0.2, NA, 0.4, 0.5)
#' lod <- rep(0.1, length(values))
#'
#' # Apply fill_in function
#' imputed_values <- fill_in(values, lod)
#' print(imputed_values)
fill_in <- function(var_to_fill, lod, loq = NULL){
  # Define maximum censoring limit (LOQ if provided, otherwise LOD)
  max_lim <- ifelse(!is.null(loq), loq, lod)

  # Identify censored observations (below maximum limit)
  censored <- var_to_fill < max_lim

  # Estimate distribution parameters using censored regression on order statistics (ROS)
  stats_ros <- NADA::cenros(var_to_fill, censored, forwardT = NULL)
  dist_mean <- NADA::mean(stats_ros)
  dist_sd <- NADA::sd(stats_ros)

  # Set censored NA values to FALSE (no imputation needed)
  censored[is.na(censored)] <- FALSE

  # Count number of censored observations below LOD
  n_impute <- sum(censored)

  # If LOQ provided, handle additional censoring between LOD and LOQ
  if (!is.null(loq)) {
    # Identify values between LOD and LOQ
    censored_loq <- var_to_fill >= lod & var_to_fill < loq
    censored_loq[is.na(censored_loq)] <- FALSE
    n_impute_loq <- sum(censored_loq)

    # Generate fill-in values between LOD and LOQ from truncated normal distribution
    fill_in_values_loq <- msm::rtnorm(
      n = n_impute_loq,
      mean = dist_mean,
      sd = dist_sd,
      lower = lod,
      upper = loq
    )

    # Update censoring status to exclude those already imputed between LOD and LOQ
    censored <- ifelse(censored_loq, FALSE, censored)
    n_impute <- sum(censored)
  }

  # Generate fill-in values below LOD from truncated normal distribution
  fill_in_values <- msm::rtnorm(
    n = n_impute,
    mean = dist_mean,
    sd = dist_sd,
    lower = rep(-Inf, n_impute),
    upper = lod
  )

  # Replace censored values with generated fill-in values
  vec_filled_in <- var_to_fill
  vec_filled_in[censored] <- fill_in_values

  # If LOQ was provided, replace censored values between LOD and LOQ
  if (!is.null(loq)) {
    vec_filled_in[censored_loq] <- fill_in_values_loq
  }

  return(vec_filled_in)
}


#' Apply Standardisation on Protocol Variables
#'
#' This function standardises the exposure data on selected protocol variables
#' (either categorical or continuous) as per the SEPAGES pipeline guide.
#' The function is designed to be called within a `mutate` call on grouped data,
#' where each group represents a different exposure.
#'
#' @param data A data.frame in tidy format containing the data to standardise,
#'   i.e., one row per ID and per exposure.
#' @param var_to_std A character string representing the variable to standardise.
#'   This variable should be normally distributed.
#' @param protocol_vars A character vector of potential protocol variables on
#'   which to standardise.
#' @param covariates A character vector of names of model covariates.
#' @param folder A character string representing the folder where regression
#'   outputs will be saved.
#' @param group A character string or vector representing the exposure group.
#'   Used for naming output files and tracking grouping structure.
#'   Defaults to the current grouping (via `dplyr::cur_group()`).
#' @param force A character vector of protocol variables to forcibly include in
#'   the model regardless of statistical significance. Default is `NULL`.
#' @param export_check_model Logical. If `TRUE`, exports model diagnostic plots
#'   from `performance::check_model()`. Default is `FALSE`.
#'
#' @details
#' The function selects protocol variables with p-value < 0.2 in ANOVA to adjust
#' for protocol effects. A linear model is fitted including selected protocol
#' variables and covariates. Residual-based corrections are computed and applied
#' to the exposure variable. If no protocol variable is selected, the original
#' variable is returned unmodified. Forced variables are always included.
#'
#' @return
#' A numeric vector of corrected values after standardisation. Length matches
#' the number of rows in `data`; values for which residuals cannot be computed
#' are returned as `NA`.
#'
#' @seealso
#' \code{\link[dplyr]{mutate}}, \code{\link[stats]{lm}}, \code{\link[stats]{predict}},
#' \code{\link[performance]{check_model}}, \code{\link[car]{Anova}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   exposure = rnorm(100, mean = 5, sd = 2),
#'   protocol_var1 = sample(c("A", "B", "C"), 100, replace = TRUE),
#'   protocol_var2 = runif(100, 0, 1),
#'   age = rnorm(100, mean = 35, sd = 10)
#' )
#'
#' data |>
#'   dplyr::mutate(
#'     exposure_std = standardise(
#'       var_to_std = "exposure",
#'       protocol_vars = c("protocol_var1", "protocol_var2"),
#'       covariates = "age",
#'       folder = "outputs/",
#'       export_check_model = TRUE
#'     )
#'   )
#' }

standardise <- function(
    data = dplyr::pick(everything()),
    var_to_std,
    protocol_vars,
    covariates = NULL,
    folder,
    group = dplyr::cur_group(),
    force = NULL,
    export_check_model = FALSE
){
  # Select protocol variables for standardisation (p < 0.2)
  final_std_vars <- get_protocol_var(data, var_to_std, protocol_vars, covariates, folder, group, export_check_model)
  if(!is.null(force)){final_std_vars = c(final_std_vars, force)}

  # --> EXIT EARLY IF NOTHING TO ADJUST FOR
  if (length(final_std_vars) == 0) {
    return(data[[var_to_std]])
  }

  # Construct model formula with final protocol variables (p < 0.2)
  rhs <- c(final_std_vars, covariates) |> purrr::compact() |> paste(collapse = " + ")
  form <- as.formula(paste(var_to_std, "~", rhs))

  # Fit model
  lm1 <- lm(form, data = data)

  # Get betas
  betas <- lm1 |>
    broom::tidy() |>
    tidycat::tidy_categorical(m = lm1)

  # add all correction factor
  for(prot_var in final_std_vars){

    # get class (numeric or categorical)
    class_i <- class(data[[prot_var]])

    if(class_i == "numeric"){

      # get beta
      beta_i <- betas |>
        dplyr::filter(variable == prot_var) |>
        dplyr::pull(estimate)

      # prepare correction factor
      data <- data |>
        mutate(
          !!stringr::str_c("correct_", prot_var) := beta_i * (.data[[prot_var]] - median(data[[prot_var]], na.rm = TRUE))
        )

    }else if(class_i == "factor"){

      # get betas and prepare correction factor (==beta)
      beta_i <- betas |>
        dplyr::filter(variable == prot_var) |>
        dplyr::select(estimate, level) |>
        dplyr::rename(
          !!prot_var := level,
          !!stringr::str_c("correct_", prot_var) := estimate
        )

      # add betas to data frame
      data <- data |>
        dplyr::left_join(beta_i, by = prot_var)

    }else{

      stop(
        stringr::str_c(
          "Protocol variables need to be coded as continuous or factor, not: ",
          class_i, "(variable ", prot_var, ")"
        )
      )

    }
  }

  # standardise
  data <- data |>
    dplyr::mutate(
      val_std = .data[[var_to_std]] - rowSums(dplyr::pick(starts_with("correct")), na.rm = TRUE)
    )

  # Return the vector of standardised values
  return(data$val_std)
}

#' Get Protocol Variables for Standardisation
#'
#' This function identifies which among the candidate protocol variables are
#' significantly associated with the exposure to be standardised. The selection
#' is based on ANOVA p-values (< 0.2). Model outputs are exported to files for
#' transparency and quality control.
#'
#' @param data A data.frame containing the exposure, protocol variables, and
#'   optional covariates.
#' @param var_to_std A character string giving the name of the variable to be
#'   standardised.
#' @param protocol_vars A character vector of names of potential protocol
#'   variables.
#' @param covariates A character vector of additional covariates to include in
#'   the model (not corrected for).
#' @param folder A character string giving the directory where outputs (model
#'   summary, ANOVA table, optional model check) will be saved.
#' @param group A character string or vector used for naming the output files
#'   (typically the exposure variable name).
#' @param export_check_model Logical. If `TRUE`, exports model diagnostic plots
#'   from `performance::check_model()`. Default is `FALSE`.
#'
#' @return A character vector containing the names of the protocol variables
#'   with p-value < 0.2, to be used for standardisation.
#'
#' @seealso \code{\link[car]{Anova}}, \code{\link[broom]{tidy}},
#' \code{\link[performance]{check_model}}, \code{\link[rio]{export}}
#'
#' @export

get_protocol_var <- function(data, var_to_std, protocol_vars, covariates, folder, group, export_check_model){
  # Construct linear model formula
  model_formula <- stringr::str_c(
    var_to_std,
    "~",
    paste(protocol_vars, collapse = "+"),
    if(!is.null(covariates)){"+"},
    paste(covariates, collapse = "+")
  )

  # Fit the full linear model
  lm_full <- lm(as.formula(model_formula), data = data)

  # Extract betas from the model and export to CSV
  betas <- broom::tidy(lm_full, conf.int = TRUE)

  # export
  filename <- stringr::str_c(paste(group, collapse = "_"), ".xlsx")
  rio::export(betas, file.path(folder, filename))
  if(export_check_model){
    filename <- stringr::str_c(paste(group, collapse = "_"), ".png")
    p <- invisible(plot(performance::check_model(lm_full)))
    ggsave(p, filename = file.path(folder, filename), height = 10, width = 8)
  }

  # Perform ANOVA and get p-values
  aov_output <- car::Anova(lm_full)

  # export anova output
  filename <- stringr::str_c(paste(group, collapse = "_"), "_aov.xlsx")
  rio::export(broom::tidy(aov_output), file.path(folder, filename))

  # Identify and return protocol variables with p < 0.2
  final_std_vars <- aov_output |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "term") |>
    dplyr::rename(p = "Pr(>F)") |>
    dplyr::filter(term %in% protocol_vars & p < 0.2) |>
    dplyr::pull(term)

  return(final_std_vars)
}
