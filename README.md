# preprocess

An R package for preprocessing environmental exposure data, with a focus on handling left-censored measurements and standardising for protocol-related variability.

## Installation

You can install the development version from GitHub using `devtools`:

```r
devtools::install_github("mjrolland/preprocess")
```

## Quick Start

```r
library(preprocess)

# Impute values below the limit of detection
imputed <- fill_in(var_to_fill = x, lod = lod)

# Standardise exposure values on protocol variables
standardised <- standardise(
  data = your_data,
  var_to_std = "x",
  protocol_vars = c("batch", "storage_time"),
  covariates = c("age", "season"),
  folder = "model_outputs"
)
```

## Core Functions

- `fill_in()`: Imputes values below the limit of detection (LOD) using censored regression on order statistics.
- `standardise()`: Corrects for protocol effects using linear models with variable selection.
- `mk_tbl_std()`: Summarises protocol variables associated with exposures across models.

## Documentation

The package includes detailed articles to guide your use:

- [Preprocessing Environmental Data](https://mjrolland.github.io/preprocess/articles/preprocess-intro.html)  
  *Walkthrough of a typical preprocessing pipeline using simulated data.*

- [Preprocessing Methodology](https://bookdown.org/mj_rolland/sepages_pipeline_doc/)  
  *Overview of the theoretical framework and main steps involved.*

You can also access the vignette from R:

```r
vignette("preprocess-intro", package = "preprocess")
```

## Contributing

Feedback, issues and contributions are welcome. Please open an issue or submit a pull request.

## License

This package is licensed under the MIT License.
