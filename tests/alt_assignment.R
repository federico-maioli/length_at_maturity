library(sf)
library(purrr)
library(tidyverse)
library(here)
library(DATRASextra)

# 1. load data ----

data <- readRDS(here("data/raw/raw_datras.rds")) |>
  DATRASextra::correct_species()

species_info <- read_rds(here("data/processed/species_info.rds")) |>
  mutate(aphia_id = as.numeric(aphia_id))

ca <- data[["CA"]] |>
  select(Survey, Quarter, Country, Ship, HaulNo, haul.id, Year,
         AreaType, AreaCode, Sex, MaturityScale, Maturity,
         Species, Valid_Aphia, LngtCm) |>
  mutate(
    Year = as.integer(as.character(Year)),
    MaturityScale = as.character(MaturityScale) |> na_if(""),
    Maturity = as.character(Maturity) |> na_if("")
  ) |>
  droplevels() |>
  left_join(species_info |> select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) |>
  mutate(Species = coalesce(species, Species)) |>
  select(-species) |>
  filter(!is.na(Species), class %in% c("Teleostei", "Elasmobranchii"))

hh <- data[["HH"]] |>
  select(Survey, Quarter, Year, Month, Country, Ship, lon, lat, StatRec, haul.id) |>
  mutate(Year = as.integer(as.character(Year))) |>
  droplevels()

m_data <- ca |>
  left_join(hh) |>
  janitor::clean_names() |>
  filter(!is.na(maturity)) |>
  mutate(maturity = as.character(maturity)) |>
  filter(
    year >= 2000,
    species %in% c("Gadus morhua", "Clupea harengus", "Platichthys flesus", "Merlangius merlangus")
  )

# 2. ices area + spawning season filter ----

sf_use_s2(FALSE)

ices_sf <- read_sf(
  here("data/extra/ICES_StatRec_mapto_ICES_Areas/StatRec_map_Areas_Full_20170124.shp")
) |>
  group_by(Area_27) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

m_data <- m_data |>
  filter(!is.na(lon), !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_join(ices_sf |> select(ices_area = Area_27), left = TRUE) |>
  st_drop_geometry()

spawn_lookup <- read_rds(here("data/processed/spawning_lookup.rds"))

m_data <- m_data |>
  left_join(
    spawn_lookup |> select(species, ices_area, spawn_months, spawn_months_ext),
    by = c("species", "ices_area")
  ) |>
  mutate(in_spawn = map2_lgl(month, spawn_months_ext, ~ !is.null(.y) && .x %in% .y)) |>
  filter(map_lgl(spawn_months, ~ is.null(.x) || length(.x) == 0) | in_spawn) # nolint

# 3. classify code types ----

# detect I/M binomial scale within each country-survey-year
im_scale <- m_data |>
  group_by(country, survey, year) |>
  summarise(has_M = any(maturity == "M", na.rm = TRUE), .groups = "drop")

m_data <- m_data |>
  left_join(im_scale, by = c("country", "survey", "year")) |>
  mutate(
    code_type = case_when(
      grepl("^R1_", maturity) ~ "R1",
      grepl("^R2_", maturity) ~ "R2",
      grepl("^6[0-9]$", maturity) ~ "M6",
      maturity %in% c("A","B","Ba","Bb","C","Ca","Cb","D","Da","Db","E","F") ~ "SMSF",
      maturity %in% c("1","2","3","4","5","6") ~ "Numeric",
      maturity %in% c("I","M") & has_M ~ "IM",
      grepl("^[IVX]+$", maturity) & !has_M ~ "Roman",
      TRUE ~ "Unknown"
    )
  ) |>
  group_by(country, code_type) |>
  mutate(
    has_high_roman = any(maturity %in% c("VII","VIII","IX","X")),
    has_VI_only = any(maturity == "VI") & !has_high_roman,
    # subdivide Roman into full (up to X) vs truncated (up to VI) vs other
    roman_subtype = case_when(
      code_type != "Roman" ~ NA_character_,
      has_high_roman ~ "Roman_full",
      has_VI_only ~ "Roman_trunc",
      TRUE ~ "Roman_other"
    )
  ) |>
  ungroup() |>
  mutate(scale = coalesce(roman_subtype, code_type))

# 3b. scale overview by country and year ----
# one row per species x country x scale; years packed as a list column

scale_overview <- m_data |>
  group_by(species, country, scale) |>
  summarise(
    years = list(sort(unique(year))),
    year_min = min(year),
    year_max = max(year),
    n_obs = n(),
    .groups = "drop"
  ) |>
  arrange(species, country, scale)

# 4. maturity assignment function ----

# parameterise the ambiguous boundary decisions for each scale
assign_mature <- function(df,
                          smsf_Ba = 0,
                          roman_full_II = 0,
                          roman_trunc_II = 1,
                          numeric_2 = 1) {
  df |>
    mutate(
      mature = case_when(
        scale == "R1" & maturity == "R1_1" ~ 0,
        scale == "R1" & maturity %in% c("R1_2","R1_3","R1_4") ~ 1,
        scale == "R2" & maturity == "R2_2" ~ 0,
        scale == "R2" & maturity %in% c("R2_4","R2_6","R2_8") ~ 1,
        scale == "M6" & maturity == "61" ~ 0,
        scale == "M6" & maturity %in% c("62","63","64","65") ~ 1,
        scale == "M6" & maturity == "66" ~ NA_real_,
        scale == "SMSF" & maturity == "A" ~ 0,
        scale == "SMSF" & maturity == "Ba" ~ smsf_Ba,
        scale == "SMSF" & maturity %in% c("B","Bb","C","Ca","Cb","D","Da","Db","E") ~ 1,
        scale == "SMSF" & maturity == "F" ~ NA_real_,
        scale == "Roman_full" & maturity == "I" ~ 0,
        scale == "Roman_full" & maturity == "II" ~ roman_full_II,
        scale == "Roman_full" & maturity %in% c("III","IV","V","VI","VII","VIII") ~ 1,
        scale == "Roman_full" & maturity %in% c("IX","X") ~ NA_real_,
        scale == "Roman_trunc" & maturity == "I" ~ 0,
        scale == "Roman_trunc" & maturity == "II" ~ roman_trunc_II,
        scale == "Roman_trunc" & maturity %in% c("III","IV","V") ~ 1,
        scale == "Roman_trunc" & maturity == "VI" ~ NA_real_,
        scale == "Roman_other" & maturity == "I" ~ 0,
        scale == "Roman_other" & maturity %in% c("II","III","IV","V") ~ 1,
        scale == "Numeric" & maturity == "1" ~ 0,
        scale == "Numeric" & maturity == "2" ~ numeric_2,
        scale == "Numeric" & maturity %in% c("3","4","5") ~ 1,
        scale == "Numeric" & maturity == "6" ~ NA_real_,
        scale == "IM" & maturity == "I" ~ 0,
        scale == "IM" & maturity == "M" ~ 1,
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(mature))
}

# 5. scenarios ----

scenarios <- tribble(
  ~scenario,           ~label,                      ~smsf_Ba, ~roman_full_II, ~roman_trunc_II, ~numeric_2,
  "current",           "current",                   0,        0,              1,               1,
  "SMSF_Ba_mat",       "SMSF Ba = mature",          1,        0,              1,               1,
  "RomanFull_II_mat",  "Roman full II = mature",    0,        1,              1,               1,
  "RomanTrunc_II_imm", "Roman trunc II = immature", 0,        0,              0,               1,
  "Numeric_2_imm",     "Numeric 2 = immature",      0,        0,              1,               0,
)

# for each uncertain scale, only the scenarios that change its L50
scale_scenarios <- tribble(
  ~scale,        ~relevant,
  "SMSF",        c("current", "SMSF_Ba_mat"),
  "Roman_full",  c("current", "RomanFull_II_mat"),
  "Roman_trunc", c("current", "RomanTrunc_II_imm"),
  "Numeric",     c("current", "Numeric_2_imm"),
)

# display names, palette, and boundary descriptions used across plots
scale_display <- c(
  "SMSF" = "SMSF",
  "Roman_full" = "Roman (7+ stages)",
  "Roman_trunc" = "Roman (max VI)",
  "Numeric" = "Numeric (6 stages)"
)

scale_pal <- c(
  "SMSF" = "#E69F00",
  "Roman (7+ stages)" = "#56B4E9",
  "Roman (max VI)" = "#009E73",
  "Numeric (6 stages)" = "#CC79A7",
  "M6 (reference)" = "#444444"
)

boundary_labels <- tribble(
  ~scale,        ~scenario,           ~boundary,
  "SMSF",        "current",           "Ba = immature",
  "SMSF",        "SMSF_Ba_mat",       "Ba = mature",
  "Roman_full",  "current",           "II = immature",
  "Roman_full",  "RomanFull_II_mat",  "II = mature",
  "Roman_trunc", "current",           "II = mature",
  "Roman_trunc", "RomanTrunc_II_imm", "II = immature",
  "Numeric",     "current",           "2 = mature",
  "Numeric",     "Numeric_2_imm",     "2 = immature",
) |>
  mutate(assignment_type = if_else(scenario == "current", "current", "alternative"))

# 6. fit l50 helper ----

fit_l50 <- function(df) {
  if (nrow(df) < 20 || sum(df$mature == 0) < 5 || sum(df$mature == 1) < 5)
    return(NULL)
  tryCatch({
    fit <- glm(mature ~ lngt_cm, data = df, family = binomial)
    if (!fit$converged) return(NULL)
    b <- coef(fit)
    tibble(l50 = as.numeric(-b[1] / b[2]),
           n = nrow(df),
           n_imm = sum(df$mature == 0),
           n_mat = sum(df$mature == 1))
  }, error = function(e) NULL)
}

# 7. specific comparison cases ----
# M6 always included as the unambiguous reference

uncertain_scales <- c("SMSF", "Roman_full", "Roman_trunc", "Numeric")

cases <- tribble(
  ~species,               ~country, ~scales,
  "Clupea harengus",      "DE",     c("M6", "Numeric", "SMSF"),
  "Clupea harengus",      "SE",     c("M6", "Numeric"),
  "Clupea harengus",      "DK",     c("M6", "Numeric", "SMSF"),
  "Gadus morhua",         "DE",     c("M6", "Numeric", "Roman_full"),
  "Gadus morhua",         "SE",     c("M6", "Numeric", "Roman_full", "SMSF"),
  "Gadus morhua",         "DK",     c("M6", "Numeric", "SMSF"),
  "Gadus morhua",         "FR",     c("M6", "Numeric", "SMSF"),
  "Merlangius merlangus", "DE",     c("M6", "Numeric", "SMSF"),
  "Merlangius merlangus", "DK",     c("M6", "Numeric", "SMSF"),
)

# 7b. find best ices area per case ----
# best area = most observations among areas where >= 2 required scales co-occur

best_area <- map_dfr(seq_len(nrow(cases)), function(i) {
  m_data |>
    filter(
      species == cases$species[i],
      country == cases$country[i],
      scale %in% cases$scales[[i]]
    ) |>
    group_by(ices_area) |>
    summarise(
      n_scales = n_distinct(scale),
      n_obs = n(),
      .groups = "drop"
    ) |>
    filter(n_scales >= 2, !is.na(ices_area)) |>
    slice_max(n_obs, n = 1, with_ties = FALSE) |>
    mutate(species = cases$species[i], country = cases$country[i])
})

cat("\nBest ICES area per case:\n")
print(best_area |> select(species, country, ices_area, n_scales, n_obs))

# 7c. filter data to best overlapping areas ----

m_overlap <- map_dfr(seq_len(nrow(best_area)), function(i) {
  req_scales <- cases |>
    filter(species == best_area$species[i], country == best_area$country[i]) |>
    pull(scales) |> .[[1]]
  m_data |>
    filter(
      species == best_area$species[i],
      country == best_area$country[i],
      ices_area == best_area$ices_area[i],
      scale %in% req_scales
    )
})

# 8. run scenarios: pooled and per year ----

run_l50 <- function(df, grp_vars) {
  map_dfr(seq_len(nrow(scenarios)), function(i) {
    sc <- scenarios[i, ]
    assign_mature(df,
                  smsf_Ba = sc$smsf_Ba,
                  roman_full_II = sc$roman_full_II,
                  roman_trunc_II = sc$roman_trunc_II,
                  numeric_2 = sc$numeric_2) |>
      filter(!is.na(lngt_cm)) |>
      group_by(across(all_of(grp_vars))) |>
      group_modify(~ {
        r <- fit_l50(.x)
        if (is.null(r)) tibble() else r
      }) |>
      ungroup() |>
      mutate(scenario = sc$scenario, label = sc$label)
  })
}

cat("Running pooled scenarios...\n")
l50_pooled <- run_l50(m_overlap, c("species", "country", "ices_area", "scale"))

cat("Running per-year scenarios...\n")
l50_yearly <- run_l50(m_overlap, c("species", "country", "ices_area", "scale", "year"))

# 9. m6 reference and deviations ----

ref_m6 <- l50_pooled |>
  filter(scale == "M6", scenario == "current") |>
  select(species, country, ices_area, ref_l50 = l50)

uncert_pooled <- l50_pooled |>
  filter(scale %in% uncertain_scales) |>
  inner_join(
    scale_scenarios |> unnest(relevant) |> rename(scenario = relevant),
    by = c("scale", "scenario")
  ) |>
  left_join(ref_m6, by = c("species", "country", "ices_area")) |>
  filter(!is.na(ref_l50)) |>
  mutate(dev = l50 - ref_l50)

# 10. output directory and shared aesthetics ----

out_dir <- here("output/04_check_maturity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sp_abbr <- c(
  "Gadus morhua" = "Cod",
  "Clupea harengus" = "Herring",
  "Merlangius merlangus" = "Whiting"
)

# x-axis order for plot A: group by scale, immature variant first within each
x_order <- c(
  "Ba = immature", "Ba = mature",
  "II = immature", "II = mature",
  "2 = mature", "2 = immature"
)

boundary_linetypes <- c(
  "M6 reference" = "solid",
  "Ba = immature" = "solid",
  "Ba = mature" = "dashed",
  "II = immature" = "solid",
  "II = mature" = "dashed",
  "2 = mature" = "solid",
  "2 = immature" = "dashed"
)
boundary_shapes <- c(
  "M6 reference" = 16,
  "Ba = immature" = 16,
  "Ba = mature" = 17,
  "II = immature" = 16,
  "II = mature" = 17,
  "2 = mature" = 16,
  "2 = immature" = 17
)

# 11. plot a: l50 vs m6, colour by scale ----

plot_data_a <- uncert_pooled |>
  left_join(boundary_labels, by = c("scale", "scenario")) |>
  mutate(
    scale_disp = scale_display[scale],
    boundary = factor(boundary, levels = x_order),
    case = paste0(sp_abbr[species], " (", country, ")\n", ices_area)
  )

ggplot(plot_data_a, aes(x = boundary, y = dev, color = scale_disp)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_point(alpha = 0.9, size = 2.5) +
  scale_color_manual(values = scale_pal, name = "Scale") +
  # scale_shape_manual(
  #   values = c("current" = 16, "alternative" = 17),
  #   labels = c("current"     = "current assignment",
  #              "alternative" = "alternative assignment"),
  #   name = NULL
  # ) +
  # scale_size_continuous(range = c(4, 9), guide = "none") +
  scale_x_discrete(drop = TRUE) +
  facet_wrap(~ case, ncol = 3, scales = "free") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    # strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  labs(x = NULL, y = "L50 - L50 (M6)  [cm]")

ggsave(here("output/supp/L50_vs_M6.png"), width = 7, height = 9, dpi = 450)
cat("Saved L50_vs_M6.png\n")

# 12. plot b: l50 per year, colour by scale, linetype by assignment ----

m6_yearly <- l50_yearly |>
  filter(scale == "M6", scenario == "current") |>
  mutate(
    scale_disp = "M6 (reference)",
    assignment_type = "current",
    boundary = "M6 reference"
  )

uncert_yearly <- l50_yearly |>
  filter(scale %in% uncertain_scales) |>
  inner_join(
    scale_scenarios |> unnest(relevant) |> rename(scenario = relevant),
    by = c("scale", "scenario")
  ) |>
  left_join(boundary_labels, by = c("scale", "scenario")) |>
  mutate(scale_disp = scale_display[scale])

yearly_plot_data <- bind_rows(m6_yearly, uncert_yearly) |>
  mutate(
    case = paste0(sp_abbr[species], " (", country, ")\n", ices_area),
    legend_group = factor(
      paste(scale_disp, boundary, sep = " - "),
      levels = c(
        "M6 (reference) - M6 reference",
        "Numeric (6 stages) - 2 = immature",
        "Numeric (6 stages) - 2 = mature",
        "SMSF - Ba = immature",
        "SMSF - Ba = mature",
        "Roman (7+ stages) - II = immature",
        "Roman (7+ stages) - II = mature"
      )
    )
  )

ggplot(
  yearly_plot_data |> filter(!is.na(l50)),
  aes(x = year, y = l50, color = legend_group, linetype = legend_group,
      group = interaction(scale, scenario))
) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(
    name = NULL,
    values = c(
      "M6 (reference) - M6 reference" = "black",
      "Numeric (6 stages) - 2 = immature" = "#E78AC3",
      "Numeric (6 stages) - 2 = mature" = "#E78AC3",
      "SMSF - Ba = immature" = "#66C2A5",
      "SMSF - Ba = mature" = "#66C2A5",
      "Roman (7+ stages) - II = immature" = "#8DA0CB",
      "Roman (7+ stages) - II = mature" = "#8DA0CB"
    )
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Numeric (6 stages) - 2 = immature" = "dotted",
      "Numeric (6 stages) - 2 = mature" = "solid",
      "SMSF - Ba = immature" = "dotted",
      "SMSF - Ba = mature" = "solid",
      "Roman (7+ stages) - II = immature" = "dotted",
      "Roman (7+ stages) - II = mature" = "solid",
      "M6 (reference) - M6 reference" = "solid"
    )
  ) +
  facet_wrap(~ case, scales = "free_y", ncol = 3) +
  labs(
    # title = "L50 per year under alternative maturity assignments",
    x = "Year",
    y = expression(L[50] ~ "(cm)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(size = 9, face = "bold", lineheight = 0.9),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.spacing.y = unit(0, "cm"),
    legend.key.height = unit(0.35, "cm"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  guides(
    color = guide_legend(ncol = 1, byrow = TRUE),
    linetype = guide_legend(ncol = 1, byrow = TRUE)
  )

ggsave(here("output/supp/L50_comparison_by_year.png"), width = 10, height = 9, dpi = 450)

# 13. console summary ----

cat("\n=== Best scenario per species x country x area x scale (smallest |dev| from M6) ===\n")
uncert_pooled |>
  mutate(abs_dev = abs(dev)) |>
  group_by(species, country, ices_area, scale) |>
  slice_min(abs_dev, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(species, country, ices_area, scale, scenario, l50, ref_l50, dev) |>
  arrange(species, country, ices_area, scale) |>
  print(n = 50)
