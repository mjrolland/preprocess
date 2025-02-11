# preprocess

An R package for preprocessing environmental datasets.

## Installation

You can install this package from GitHub using `devtools`:

```
devtools::install_github("mjrolland/preprocess")
```

## Usage

```
library(preprocess)
```

## Functions

- `fill_in()`: Handles values below the limit of detection.
- `standardise()`: Standardizes exposure data based on protocol variables.

## Example

### Filling in Values Below LOD

```

# Example data
data <- c(1, 2, 0.5, 3, 0.2, 0.8)
lod <- c(NA, NA, 1, NA, 1, 1)

# Apply fill-in method
result <- fill_in(data, lod)
print(result)

```

### Standardizing Exposure Data

```

library(preprocess)
library(dplyr)

# Simulated dataset
df <- data.frame(
  id = 1:100,
  exp1 = rnorm(100, mean = 3, sd = 1),
  batch = factor(sample(1:5, 100, replace = TRUE)),
  prot1 = rnorm(100, mean = 0, sd = 0.3)
)

# Standardize the data
df_standardized <- df |>
  group_by(batch) |>
  mutate(
    exp1_std = standardise(., var_to_std = "exp1", protocol_vars = "prot1")
  )

head(df_standardized)

```

## Contributing

Contributions are welcome! Feel free to submit issues or pull requests.

## License

This package is licensed under the MIT License.
