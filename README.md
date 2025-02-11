# preprocess

An R package for preprocessing environmental datasets.

## Installation

You can install this package from GitHub using `devtools`:

```r
devtools::install_github("mjrolland/preprocess")
```

## Usage

```r
library(preprocess)
```

## Functions

- `fill_in()`: Handles values below the limit of detection.
- `standardise()`: Standardizes exposure data based on protocol variables.

## Documentation

For a practical **step-by-step guide** on using this package, refer to the **vignette** or the **accompanying article**:

- Vignette: `vignette("preprocess-intro", package = "preprocess")`
- Article: [Preprocessing Environmental Data](https://mjrolland.github.io/preprocess/articles/preprocess-intro.html)

For a **theoretical background** on the preprocessing methodology, see the detailed documentation:

[https://bookdown.org/mj_rolland/sepages_pipeline_doc/](https://bookdown.org/mj_rolland/sepages_pipeline_doc/)

## Contributing

Contributions are welcome! Feel free to submit issues or pull requests.

## License

This package is licensed under the MIT License.
