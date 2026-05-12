library(tidyverse)
library(here)

# 1. load data ----

data <- read_rds(here("data/intermediate/l50_raw.rds"))

species_info <- read_rds(here("data/metadata/species_info.rds")) |>
  mutate(species_clean = str_replace_all(species, " ", "_")) |>
  select(species_clean, aphia_id,
         l_50_fb_median = l_50_median,
         l_50_fb_min    = l_50_min,
         l_50_fb_max    = l_50_max)

# 2. clean estimates ----

data_clean <- data |>
  # extract sample size from nested data before any filtering
  mutate(n = map_int(data, nrow)) |>
  left_join(species_info, by = "species_clean") |>
  filter(
    l50_est  > 0,                                    # L50 must be positive
    l50_se   < 10,                                   # SE biologically reasonable
    l50_lower > 0,                                   # lower CI must be positive
    !is.na(l50_se),                                  # drop rows where SE failed
    (l50_upper - l50_lower) <= 0.5 * l50_est,        # CI width at most 50% of L50
    l50_est  > 0.7 * l_50_fb_min,                    # not suspiciously smaller than FishBase min
    l50_est  < 1.3 * l_50_fb_max                     # not suspiciously larger than FishBase max
  ) |>
  select(-data, -l_50_fb_min, -l_50_fb_median, -l_50_fb_max)

# 3. save cleaned estimates ----

write_rds(data_clean, here("data/final/l50_clean.rds"))

# 4. prepare csv for data repository ----

data_export <- data_clean |>
  ungroup() |>
  mutate(species = str_replace_all(species_clean, "_", " ")) |>
  select(-species_clean) |>
  rename(
    spatial_scale   = model,
    ices_stock      = species_stock,
    time_period     = period,
    sample_size     = n,
    l25_cm          = l25_est,
    l25_se_cm       = l25_se,
    l25_ci_lower_cm = l25_lower,
    l25_ci_upper_cm = l25_upper,
    l50_cm          = l50_est,
    l50_se_cm       = l50_se,
    l50_ci_lower_cm = l50_lower,
    l50_ci_upper_cm = l50_upper,
    l75_cm          = l75_est,
    l75_se_cm       = l75_se,
    l75_ci_lower_cm = l75_lower,
    l75_ci_upper_cm = l75_upper,
    model_auc       = auc,
    model_tjur_r2   = tjur_r2
  ) |>
  select(
    aphia_id, species, sex, spatial_scale, ices_stock, ices_area, time_period,
    sample_size,
    l25_cm, l25_se_cm, l25_ci_lower_cm, l25_ci_upper_cm,
    l50_cm, l50_se_cm, l50_ci_lower_cm, l50_ci_upper_cm,
    l75_cm, l75_se_cm, l75_ci_lower_cm, l75_ci_upper_cm,
    model_auc, model_tjur_r2
  ) |>
  # round all numeric columns to 2 decimal places
  mutate(across(where(is.numeric) & !sample_size, ~ round(.x, 2)))

# 5. save ----

write_csv(data_export, here("data/final/l50_estimates.csv"))
