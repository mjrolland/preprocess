# Simulating data to test fill-in and standardise function
# M. Rolland
# 20/09/2023

# Setting random seed for reproducibility
set.seed(123)

# Set dataset parameters
n <- 1000
exp_means <- c(exp1 = 3, exp2 = 2, exp3 = 0.5)
batch_effect <- c(0, 1, -1, 0, 2)
season_effect <- c(2, 3, -1, 1)

# Simulate data
df <- data.frame(id = str_c("id", 1:n)) |>                         # init
  mutate(
    batch = factor(sample(1:5, n, replace = TRUE)),                 # create random batch var
    exp1  = rlnorm(n, meanlog = log(exp_means[1]), sdlog = 0.5),    # simulate lognormally distributed data for exp1
    exp2  = rlnorm(n, meanlog = log(exp_means[2]), sdlog = 0.5),    # simulate lognormally distributed data for exp2
    exp3  = rlnorm(n, meanlog = log(exp_means[3]), sdlog = 0.5),    # simulate lognormally distributed data for exp3
    exp1  = exp1 + batch_effect[as.numeric(as.character(batch))],   # add batch effect to exp1
    prot1 = rnorm(n, mean = 0, sd = 0.3),                           # simulate protocol variable 1
    exp2  = exp2 + prot1,                                           # add protocol variable effect to exp2
    prot2 = rnorm(n, mean = 50, sd = 10)                            # simulate protocol variable 2
  ) |>
  dplyr::arrange(batch)

# Adding season that is associated with batch
df$season <- NA
df$season[df$batch %in% c(1, 2)] <- sample(c("Summer", "Spring", "Autumn", "Winter"), size = sum(df$batch %in% c(1, 2)), prob = c(0.25, 0.25, 0.25, 0.25), replace = TRUE)
df$season[df$batch %in% c(3, 4, 5)] <- sample(c("Summer", "Spring", "Autumn", "Winter"), size = sum(df$batch %in% c(3, 4, 5)), prob = c(0.1, 0.6, 0.2, 0.1), replace = TRUE)
df$season <- as.factor(df$season)

# Add season effect + other covariate
df <- df |>
  mutate(
    exp4 = exp1 + season_effect[as.numeric(season)],
    cov1 = rnorm(n, mean = 0, sd = 0.3)
  )

# Simulate LODs for a given % of data < LOD (5%, 15% and 35%)
lods <- c(
  lod1 = round(as.numeric(quantile(df$exp1, 0.05)), 1),
  lod2 = round(as.numeric(quantile(df$exp2, 0.15)), 1),
  lod3 = round(as.numeric(quantile(df$exp3, 0.35)), 1),
  lod4 = round(as.numeric(quantile(df$exp4, 0.05)), 1)
)

# Generate lists of ids to generate missing data
missing1 <- sample(1:1000, 40)
missing2 <- sample(1:1000, 10)
missing3 <- sample(1:1000, 10)

# Set "ND" and NA data
df <- df |>
  mutate(
    # set NA data for each exp
    exp1 = ifelse(dplyr::row_number() %in% missing1, NA, exp1),
    exp1 = ifelse(dplyr::row_number() %in% missing2, NA, exp1),
    exp2 = ifelse(dplyr::row_number() %in% missing1, NA, exp2),
    exp2 = ifelse(dplyr::row_number() %in% missing3, NA, exp2),
    exp3 = ifelse(dplyr::row_number() %in% missing1, NA, exp3),
    exp3 = ifelse(dplyr::row_number() %in% missing2, NA, exp3),
    # set "ND" (non detect) data for each exp
    exp1 = ifelse(exp1 < lods[1], "<LOD", as.character(exp1)),
    exp2 = ifelse(exp2 < lods[2], "<LOD", as.character(exp2)),
    exp3 = ifelse(exp3 < lods[3], "<LOD", as.character(exp3)),
    exp4 = ifelse(exp4 < lods[4], "<LOD", as.character(exp4))
  )

# Export simulated data
# readr::write_csv(df, "data/simulated_data.csv")

# Export LODs
lod_data <- data.frame(
  exp = c("exp1", "exp2", "exp3", "exp4"),
  lod = lods
)
# readr::write_csv(lod_data, "data/simulated_lod_data.csv")
