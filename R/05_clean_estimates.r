library(tidyverse)
library(here)
# library(icesdata)  # only needed for the ICES stock-assessment comparison (commented out below)

# 01 Load data ----

data <- read_rds(here("data/intermediate/l50_raw.rds"))

# species aphia id for the export (FishBase l_50_median no longer loaded; comparison removed)
species_info <- read_rds(here("data/metadata/species_info.rds")) |>
  mutate(species_clean = str_replace_all(species, " ", "_")) |>
  select(species_clean, aphia_id)

# 02 Clean estimates ----
# sanity check: keep only statistically sound estimates

data_clean <- data |>
  left_join(species_info, by = "species_clean") |>
  filter(
    l50_est   > 0,                              # L50 must be positive
    l25_est   > 0,                              # L25 must be positive
    l50_lower > 0,                              # lower CI must be positive
    !is.na(l50_se),                             # drop rows where the SE failed
    l50_se    < 10,                             # SE biologically reasonable
    (l50_upper - l50_lower) <= 0.5 * l50_est,   # CI no wider than 50% of L50
    l25_est < l50_est, l50_est < l75_est,       # increasing ogive (positive slope)
    auc >= 0.5                                  # model discriminates no worse than chance
  )

# report how many estimates the sanity filter drops
n_total <- nrow(data)
n_na    <- sum(is.na(data$l50_est))
n_kept  <- nrow(data_clean)
cat(sprintf("Estimates: %d total, %d unfittable (NA), %d dropped by sanity filters, %d retained\n",
            n_total, n_na, n_total - n_na - n_kept, n_kept))

# 03 Save cleaned estimates ----

write_rds(data_clean, here("data/final/l50_clean.rds"))

# 04 Prepare csv for data repository ----

data_export <- data_clean |>
  ungroup() |>
  mutate(
    species = str_replace_all(species_clean, "_", " "),
    # ICES stock code without the species prefix (Gadus_morhua_cod.27.21 -> cod.27.21)
    ices_stock = str_remove(species_stock, paste0(species_clean, "_"))
  ) |>
  select(-species_clean, -species_stock) |>
  rename(
    spatial_scale = model,
    time_period = period,
    sample_size = n,
    l25_cm = l25_est,
    l25_se_cm = l25_se,
    l25_ci_lower_cm = l25_lower,
    l25_ci_upper_cm = l25_upper,
    l50_cm = l50_est,
    l50_se_cm = l50_se,
    l50_ci_lower_cm = l50_lower,
    l50_ci_upper_cm = l50_upper,
    l75_cm = l75_est,
    l75_se_cm = l75_se,
    l75_ci_lower_cm = l75_lower,
    l75_ci_upper_cm = l75_upper,
    model_auc = auc,
    model_tjur_r2 = tjur_r2
  ) |>
  select(
    aphia_id, species, sex, spatial_scale, ices_stock, ices_area, time_period,
    sample_size,
    l25_cm, l25_se_cm, l25_ci_lower_cm, l25_ci_upper_cm,
    l50_cm, l50_se_cm, l50_ci_lower_cm, l50_ci_upper_cm,
    l75_cm, l75_se_cm, l75_ci_lower_cm, l75_ci_upper_cm,
    model_auc, model_tjur_r2
  ) |>
  mutate(across(where(is.numeric) & !sample_size, ~ round(.x, 2)))

write_csv(data_export, here("data/final/l50_estimates.csv"))

# number of species in the final estimates
cat(sprintf("Species in the final estimates: %d\n", n_distinct(data_clean$species_clean)))

# 05 Validate: fit quality by aggregation level ----
# AUC and Tjur R2 by aggregation level are plotted in 09_plot_performance.R

# 06 Global L50 vs FishBase (not used - commented out) ----
# our global (all-year, combined-sex) L50 per species, against the FishBase median Lm
#
# our_global <- data_clean |>
#   ungroup() |>
#   filter(model == "Global", sex == "Combined", is.na(period)) |>
#   select(species_clean, l50_est)
#
# global_compare <- our_global |>
#   inner_join(species_info |> select(species_clean, l_50_fb_median), by = "species_clean") |>
#   filter(!is.na(l_50_fb_median))
#
# r_global <- cor(global_compare$l50_est, global_compare$l_50_fb_median)
#
# ggplot(global_compare, aes(x = l_50_fb_median, y = l50_est)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
#   geom_point() +
#   geom_text(aes(label = str_replace_all(species_clean, "_", " ")),
#             size = 2.3, vjust = -0.6, check_overlap = TRUE) +
#   coord_equal() +
#   labs(
#     x = "FishBase median length at maturity (cm)",
#     y = "Our global L50 (combined sexes, cm)",
#     title = "Our global L50 vs FishBase",
#     subtitle = sprintf("n = %d species | r = %.2f", nrow(global_compare), r_global)
#   ) +
#   theme_light()
#
# ggsave(here("outputs/supp/check_l50_vs_fishbase.png"), width = 16, height = 15,
#        units = "cm", dpi = 300)

# 07 Stock L50 vs ICES size at maturity (not used - commented out) ----
#
# # ICES size at maturity per stock, from age at maturity and von Bertalanffy growth
# out_list <- list()
# for (stock_name in names(icesdata)) {
#   stk <- icesdata[[stock_name]]
#   mat <- mat50(stk)
#   fl  <- attributes(stk)$fishlife
#   if (is.null(fl)) next
#   linf <- unname(fl["linf"])
#   k    <- unname(fl["k"])
#   out_list[[stock_name]] <- data.frame(
#     stock         = stock_name,
#     species       = attributes(stk)$species,
#     year          = as.integer(names(mat)),
#     age_maturity  = as.numeric(mat),
#     linf          = linf,
#     k             = k,
#     size_maturity = linf * (1 - exp(-k * as.numeric(mat)))
#   )
# }
# out <- bind_rows(out_list)
#
# # mean ICES size at maturity per stock over years > 2000
# out_mean <- out |>
#   filter(year > 2000) |>
#   summarise(size_maturity_ices = mean(size_maturity, na.rm = TRUE), .by = stock)
#
# # our stock-level L50 (combined sexes, all-years); recover the ICES stock code from species_stock
# our_l50_stock <- data_clean |>
#   ungroup() |>
#   filter(model == "Stock", is.na(period), sex == "Combined",
#          str_detect(species_stock, "\\.27\\.")) |>
#   mutate(stock = str_sub(species_stock, nchar(species_clean) + 2)) |>
#   select(species_clean, species_stock, stock, l50_est)
#
# stock_compare <- our_l50_stock |>
#   inner_join(out_mean, by = "stock")
#
# r_stock <- cor(stock_compare$l50_est, stock_compare$size_maturity_ices)
#
# ggplot(stock_compare, aes(x = size_maturity_ices, y = l50_est)) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
#   geom_point() +
#   geom_text(aes(label = stock), size = 2.3, vjust = -0.6, check_overlap = TRUE) +
#   coord_equal() +
#   labs(
#     x = "ICES size at maturity (mean of years > 2000, cm)",
#     y = "Our stock L50 (combined sexes, cm)",
#     title = "Our stock L50 vs ICES size at maturity",
#     subtitle = sprintf("n = %d stocks | r = %.2f", nrow(stock_compare), r_stock)
#   ) +
#   theme_light()
#
# ggsave(here("outputs/supp/check_l50_vs_ices.png"), width = 16, height = 15,
#        units = "cm", dpi = 300)
