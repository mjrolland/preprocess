% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{mk_tbl_std}
\alias{mk_tbl_std}
\title{Create a Standardized Table of Significant Variables from AOV Results}
\usage{
mk_tbl_std(path, lst_prot_vars)
}
\arguments{
\item{path}{Character string. Path to the folder containing AOV result files.}

\item{lst_prot_vars}{Character vector. List of variable names to check in the AOV results.}
}
\value{
A tibble with exposures as rows and variables as columns, where \code{"X"}
indicates a significant association (p < 0.2).
}
\description{
This function scans a directory for files ending in \verb{_aov.csv}, extracts
variables with a p-value below 0.2, and creates a summary table indicating
which variables were found to be significant.
}
\examples{
\dontrun{
# Define the path to the folder with AOV results
path_to_results <- "path/to/aov/files/"

# Define the list of variables to check
prot_vars <- c("var1", "var2", "var3")

# Generate the summary table
mk_tbl_std(path_to_results, prot_vars)
}
}
