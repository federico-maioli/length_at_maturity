# libraries
library(tidyverse)
library(here)

# load data
data <- read_rds(here('data/intermediate/l50_raw.rds'))

species_info <- read_rds(here('data/metadata/species_info.rds')) |>
  mutate(species_clean = str_replace_all(species, " ", "_")) %>%
  select(species_clean, aphia_id, l_50_fb_median = l_50_median, l_50_fb_min = l_50_min, l_50_fb_max = l_50_max)

# clean data
data_clean <- data %>%
  # extract sample size from nested data before any filtering
  mutate(n = map_int(data, nrow)) %>%
  # join by scientific name
  left_join(
    species_info,
    by = 'species_clean'
  ) %>%
  filter(
    # remove models that failed or have infinite/extreme SEs
    l50_est > 0,               # L50 must be positive
    l50_se < 10,           # SE should be biologically reasonable (say < 10cm)
    l50_lower > 0,         # l50_lower should be > 0
    !is.na(l50_se),         # drop anything where SE couldn't calculate
    # CI width (upper - lower) cannot exceed 50% of the L50 estimate
    (l50_upper - l50_lower) <= 0.5 * l50_est,
    # ensure l50 is within a realistic biological value
    #   l50 < l_inf * 1.2,                       # cannot mature larger than max size
    l50_est > (0.7 * l_50_fb_min),                   # not suspiciously smaller than min fishbase
    l50_est < (1.3 * l_50_fb_max)                    # not suspiciously larger than max fishbase
  ) %>%
  # cleanup columns
  select(-data, -l_50_fb_min, -l_50_fb_median, -l_50_fb_min, -l_50_fb_max)


write_rds(data_clean, here("data/final/l50_clean.rds"))

# report Data 


# prepare CSV for data repository -----------------------------------------

data_export <- data_clean %>% ungroup() |> 
  mutate(species = str_replace_all(species_clean, "_", " ")) %>%
  select(-species_clean) %>%
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
  ) %>%
  select(
    aphia_id, species, sex, spatial_scale, ices_stock, ices_area, time_period,
    sample_size,
    l25_cm, l25_se_cm, l25_ci_lower_cm, l25_ci_upper_cm,
    l50_cm, l50_se_cm, l50_ci_lower_cm, l50_ci_upper_cm,
    l75_cm, l75_se_cm, l75_ci_lower_cm, l75_ci_upper_cm,
    model_auc, model_tjur_r2
  )

write_csv(data_export, here("data/final/l50_estimates.csv"))
