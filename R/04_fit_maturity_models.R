# libraries ---------------------------------------------------------------
library(tidyverse)
library(here)
library(pROC)
library(performance)

set.seed(2)

# data -----------------------------------------------------
data <- readRDS(here('data/intermediate/maturity_clean.rds')) %>%
  dplyr::select(
    survey, quarter, country, haul_id, year, month, lon, lat,
    stat_rec, species, valid_aphia, ices_area, stock, species_clean,
    species_stock, sex, lngt_cm, mature
  ) %>%
  mutate(period = paste0((year %/% 5) * 5, "-", (year %/% 5) * 5 + 4))

# prepare data for different scales ------------------------------------

## global ------------------------------------------------------------------
data_global_comb <- data %>%
  group_by(species_clean) %>%
  filter(
    sum(mature == 0) >= 15, 
    sum(mature == 1) >= 15
  ) %>%
  nest() %>%
  mutate(model = "Global", sex = "Combined")

data_global_sex <- data %>% 
  filter(sex %in% c('F', 'M')) %>%
  group_by(species_clean) %>%
  # check that both sexes meet the threshold within the species
  filter(
    sum(sex == "M" & mature == 0) >= 15,
    sum(sex == "M" & mature == 1) >= 15,
    sum(sex == "F" & mature == 0) >= 15,
    sum(sex == "F" & mature == 1) >= 15 
  ) %>%
  # regroup to create separate nested dataframes for each sex
  group_by(species_clean, sex) %>%
  nest() %>%
  mutate(model = "Global")


## stock -------------------------------------------------------------------
data_stock_comb <- data %>% 
  filter(str_detect(species_stock, "\\.27\\.")) %>%
  group_by(species_clean, species_stock) %>% 
  filter(
    sum(mature == 0) >= 15, 
    sum(mature == 1) >= 15
  ) %>%
  nest() %>%
  mutate(model = "Stock", sex = "Combined")

data_stock_sex <- data %>% 
  filter(sex %in% c('F', 'M'), str_detect(species_stock, "\\.27\\.")) %>%
  group_by(species_stock) %>%
  filter(
    sum(sex == "M" & mature == 0) >= 15, 
    sum(sex == "M" & mature == 1) >= 15,
    sum(sex == "F" & mature == 0) >= 15, 
    sum(sex == "F" & mature == 1) >= 15
  ) %>%
  group_by(species_clean, species_stock, sex) %>% 
  nest() %>%
  mutate(model = "Stock")

## area --------------------------------------------------------------------
data_area_comb <- data %>%
  group_by(species_clean, ices_area) %>%
  filter(
    sum(mature == 0) >= 15,
    sum(mature == 1) >= 15
  ) %>%
  nest() %>%
  mutate(model = "Area", sex = "Combined")

data_area_sex <- data %>%
  filter(sex %in% c('F', 'M')) %>%
  group_by(species_clean, ices_area) %>%
  filter(
    sum(sex == "M" & mature == 0) >= 15,
    sum(sex == "M" & mature == 1) >= 15,
    sum(sex == "F" & mature == 0) >= 15,
    sum(sex == "F" & mature == 1) >= 15
  ) %>%
  group_by(species_clean, ices_area, sex) %>%
  nest() %>%
  mutate(model = "Area")

## period ------------------------------------------------------------------
data_global_comb_period <- data %>%
  group_by(species_clean, period) %>%
  filter(
    sum(mature == 0) >= 15,
    sum(mature == 1) >= 15
  ) %>%
  nest() %>%
  mutate(model = "Global", sex = "Combined")

data_global_sex_period <- data %>%
  filter(sex %in% c('F', 'M')) %>%
  group_by(species_clean, period) %>%
  filter(
    sum(sex == "M" & mature == 0) >= 15,
    sum(sex == "M" & mature == 1) >= 15,
    sum(sex == "F" & mature == 0) >= 15,
    sum(sex == "F" & mature == 1) >= 15
  ) %>%
  group_by(species_clean, period, sex) %>%
  nest() %>%
  mutate(model = "Global")

data_stock_comb_period <- data %>%
  filter(str_detect(species_stock, "\\.27\\.")) %>%
  group_by(species_clean, species_stock, period) %>%
  filter(
    sum(mature == 0) >= 15,
    sum(mature == 1) >= 15
  ) %>%
  nest() %>%
  mutate(model = "Stock", sex = "Combined")

data_stock_sex_period <- data %>%
  filter(sex %in% c('F', 'M'), str_detect(species_stock, "\\.27\\.")) %>%
  group_by(species_stock, period) %>%
  filter(
    sum(sex == "M" & mature == 0) >= 15,
    sum(sex == "M" & mature == 1) >= 15,
    sum(sex == "F" & mature == 0) >= 15,
    sum(sex == "F" & mature == 1) >= 15
  ) %>%
  group_by(species_clean, species_stock, period, sex) %>%
  nest() %>%
  mutate(model = "Stock")

data_area_comb_period <- data %>%
  group_by(species_clean, ices_area, period) %>%
  filter(
    sum(mature == 0) >= 15,
    sum(mature == 1) >= 15
  ) %>%
  nest() %>%
  mutate(model = "Area", sex = "Combined")

data_area_sex_period <- data %>%
  filter(sex %in% c('F', 'M')) %>%
  group_by(species_clean, ices_area, period) %>%
  filter(
    sum(sex == "M" & mature == 0) >= 15,
    sum(sex == "M" & mature == 1) >= 15,
    sum(sex == "F" & mature == 0) >= 15,
    sum(sex == "F" & mature == 1) >= 15
  ) %>%
  group_by(species_clean, ices_area, period, sex) %>%
  nest() %>%
  mutate(model = "Area")


# fit the models ----------------------------------------------------------

# function for length at probability
length_at_probability <- function(cf, p = 0.5) {
  (qlogis(p) - cf[1]) / cf[2]
}

# helper to build an all-NA result table for a given probability vector
.na_lp_result <- function(p, est = NA_real_) {
  tibble::tibble(
    label   = paste0("l", round(p * 100)),
    est     = est,
    se      = NA_real_,
    lower   = NA_real_,
    upper   = NA_real_,
    auc     = NA_real_,
    tjur_r2 = NA_real_
  )
}

# fit a binomial glm on resampled rows; returns coefficients or c(NA, NA)
.boot_one <- function(sub_df) {
  idx <- sample.int(nrow(sub_df), replace = TRUE)
  fit <- tryCatch(
    suppressWarnings(
      glm(mature ~ lngt_cm, data = sub_df[idx, ], family = binomial)
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || !isTRUE(fit$converged)) return(c(NA_real_, NA_real_))
  unname(coef(fit))
}

# bootstrap CIs around L25 / L50 / L75 (or any p)
get_lp_ci <- function(sub_df, p = c(0.25, 0.5, 0.75), n_boot = 100) { # do 1k for final run
  
  # check maturity variation
  if (length(unique(sub_df$mature)) < 2 ||
      sum(sub_df$mature == 0, na.rm = TRUE) < 15 ||
      sum(sub_df$mature == 1, na.rm = TRUE) < 15) {
    return(.na_lp_result(p))
  }
  
  # fit model on full data
  mod <- tryCatch(
    suppressWarnings(glm(mature ~ lngt_cm, data = sub_df, family = binomial)),
    error = function(e) NULL
  )
  if (is.null(mod) || !isTRUE(mod$converged)) return(.na_lp_result(p))
  
  cf <- coef(mod)

  # model validation metrics
  pred    <- predict(mod, type = "response")
  roc_obj <- pROC::roc(sub_df$mature, pred, quiet = TRUE)
  auc_val     <- as.numeric(pROC::auc(roc_obj))
  tjur_r2_val <- as.numeric(performance::r2_tjur(mod))

  # manual bootstrap: resample rows, refit, collect coefficients.
  boot_coefs <- replicate(n_boot, .boot_one(sub_df))   # 2 x n_boot matrix

  # compute results per probability
  purrr::map_dfr(p, function(prob) {
    label  <- paste0("l", round(prob * 100))
    lp_est <- length_at_probability(cf, prob)

    bLp <- apply(boot_coefs, 2, function(b) {
      if (anyNA(b) || !is.finite(b[1]) || !is.finite(b[2]) || abs(b[2]) < 1e-6) {
        return(NA_real_)
      }
      length_at_probability(b, prob)
    })
    bLp <- bLp[is.finite(bLp)]

    if (length(bLp) < 10) {
      return(tibble::tibble(
        label = label, est = lp_est,
        se = NA_real_, lower = NA_real_, upper = NA_real_,
        auc = auc_val, tjur_r2 = tjur_r2_val
      ))
    }

    ci <- quantile(bLp, c(0.025, 0.975), na.rm = TRUE)
    tibble::tibble(
      label   = label,
      est     = lp_est,
      se      = sd(bLp),
      lower   = unname(ci[1]),
      upper   = unname(ci[2]),
      auc     = auc_val,
      tjur_r2 = tjur_r2_val
    )
  })
}

maturity_data <- bind_rows(
  data_global_comb        %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Global combined")) %>% unnest(res, keep_empty = TRUE),
  data_global_sex         %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Global by sex")) %>% unnest(res, keep_empty = TRUE),
  data_stock_comb         %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Stock combined")) %>% unnest(res, keep_empty = TRUE),
  data_stock_sex          %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Stock by sex")) %>% unnest(res, keep_empty = TRUE),
  data_area_comb          %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Area combined")) %>% unnest(res, keep_empty = TRUE),
  data_area_sex           %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Area by sex")) %>% unnest(res, keep_empty = TRUE),
  data_global_comb_period %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Global combined period")) %>% unnest(res, keep_empty = TRUE),
  data_global_sex_period  %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Global by sex period")) %>% unnest(res, keep_empty = TRUE),
  data_stock_comb_period  %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Stock combined period")) %>% unnest(res, keep_empty = TRUE),
  data_stock_sex_period   %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Stock by sex period")) %>% unnest(res, keep_empty = TRUE),
  data_area_comb_period   %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Area combined period")) %>% unnest(res, keep_empty = TRUE),
  data_area_sex_period    %>% mutate(res = map(data, ~ get_lp_ci(.x), .progress = "Area by sex period")) %>% unnest(res, keep_empty = TRUE)
) %>%
  tidyr::pivot_wider(
    names_from = label,
    values_from = c(est, se, lower, upper),
    names_glue = "{label}_{.value}"
  )


# save data  --------------------------------------------------------------
write_rds(maturity_data, here('data/intermediate/l50_raw.rds'))
