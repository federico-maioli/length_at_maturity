library(tidyverse)
library(here)
library(sf)
library(DATRASextra)
library(knitr)
library(kableExtra)

# Diagnostics supporting the Baltic (BITS) maturity harmonisation. For BITS records coded on
# the Roman-numeral national keys we test how the ambiguous stage II should be treated, by
# comparing L50 under three rules (stage II excluded, I+II immature, II mature) against the
# unambiguous reference keys (BITS numeric, M6, SMSF). Three figures:
#  (1) check_bits_roman_sensitivity: two well-sampled cod cases, L50 per year under each rule
#      against the reference keys used in other years (same ICES area).
#  (2) check_bits_rult_immature: Lithuania and Russia, which record almost no stage I, so
#      stages I and II together are immature; checked against their reference keys.
#  (3) check_bits_pooled_by_key: one pooled L50 (delta-method CI) per species and rule, across
#      several ICES subdivisions, plus a numeric summary of the deviation from the reference.

# 01 Load data ----

raw <- readRDS(here("data/raw/raw_datras.rds")) |> DATRASextra::correct_species()

species_info <- read_rds(here("data/metadata/species_info.rds")) |>
  mutate(aphia_id = as.numeric(aphia_id))

# 02 Build BITS maturity records ----

ca <- raw[["CA"]] |>
  select(Survey, Quarter, Country, Ship, HaulNo, haul.id, Year,
         Sex, Maturity, Species, Valid_Aphia, LngtCm) |>
  mutate(
    Year = as.integer(as.character(Year)),
    Maturity = na_if(as.character(Maturity), "")
  ) |>
  droplevels() |>
  left_join(species_info |> select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) |>
  mutate(Species = coalesce(species, Species)) |>
  select(-species) |>
  filter(!is.na(Species), class %in% c("Teleostei", "Elasmobranchii"))

# haul positions, so each record can be assigned an ICES area
hh <- raw[["HH"]] |>
  select(Survey, Quarter, Year, Country, Ship, lon, lat, haul.id) |>
  mutate(Year = as.integer(as.character(Year))) |>
  droplevels()

m_data <- ca |>
  left_join(hh) |>
  janitor::clean_names() |>
  filter(survey == "BITS", !is.na(maturity), year >= 2000,
         !country %in% c("RU", "LT")) |>
  mutate(maturity = as.character(maturity))

# 03 Classify coding scheme ----

# detect the I/M binomial scale within each country-year (not used in BITS, but keeps
# the Roman detection identical to the main pipeline)
im_scale <- m_data |>
  group_by(country, year) |>
  summarise(has_M = any(maturity == "M", na.rm = TRUE), .groups = "drop")

m_data <- m_data |>
  left_join(im_scale, by = c("country", "year"))

# SE reported its 2000-2005 BITS maturity in Roman numerals, but this is the 6-stage
# numeric scale; recode to numeric before classifying, as in the main pipeline
m_data <- m_data |>
  mutate(maturity = if_else(
    country == "SE" & year >= 2000 & year <= 2005 & grepl("^[IVX]+$", maturity),
    recode(maturity, "I" = "1", "II" = "2", "III" = "3", "IV" = "4", "V" = "5", "VI" = "6"),
    maturity
  ))

m_data <- m_data |>
  mutate(
    code_type = case_when(
      grepl("^6[0-9]$", maturity) ~ "M6",
      maturity %in% c("A","B","Ba","Bb","C","Ca","Cb","D","Da","Db","E","F") ~ "SMSF",
      maturity %in% c("1","2","3","4","5","6") ~ "Numeric",
      maturity %in% c("I","M") & has_M ~ "IM",
      grepl("^[IVX]+$", maturity) & !has_M ~ "Roman",
      TRUE ~ "Unknown"
    )
  )

# 04 Assign an ICES area to every record ----

ices_sf <- read_sf(here("data/metadata/ices_areas/ICES_Statrec_mapto_ICES_Areas.shp")) |>
  group_by(Area_27) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

m_data <- m_data |>
  filter(!is.na(lon), !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_join(ices_sf |> select(ices_area = Area_27), left = TRUE) |>
  st_drop_geometry()

# 05 Shared helpers ----

# standard assignment for the non-Roman scales; series named after the scale
assign_std <- function(df) {
  df |>
    filter(code_type != "Roman") |>
    mutate(
      mature = case_when(
        code_type == "M6" & maturity == "61" ~ 0,
        code_type == "M6" & maturity %in% c("62","63","64","65") ~ 1,
        code_type == "M6" & maturity == "66" ~ NA_real_,
        code_type == "Numeric" & maturity == "1" ~ 0,
        code_type == "Numeric" & maturity %in% c("2","3","4","5") ~ 1,
        code_type == "Numeric" & maturity == "6" ~ NA_real_,
        code_type == "SMSF" & maturity %in% c("A","Ba") ~ 0,
        code_type == "SMSF" & maturity %in% c("B","Bb","C","Ca","Cb","D","Da","Db","E") ~ 1,
        code_type == "SMSF" & maturity == "F" ~ NA_real_,
        TRUE ~ NA_real_
      ),
      series = code_type
    )
}

# Roman assignment: stage I immature, III-VIII mature, IX/X dropped; stage II set by
# ii_value (NA drops it), which is the decision being tested
assign_roman <- function(df, ii_value, series_label) {
  df |>
    filter(code_type == "Roman") |>
    mutate(
      mature = case_when(
        maturity == "I" ~ 0,
        maturity == "II" ~ ii_value,
        maturity %in% c("III","IV","V","VI","VII","VIII") ~ 1,
        maturity %in% c("IX","X") ~ NA_real_,
        TRUE ~ NA_real_
      ),
      series = series_label
    )
}

# logistic L50 with a delta-method 95% CI, requiring at least 15 immature and 15 mature
fit_l50 <- function(df) {
  if (sum(df$mature == 0) < 15 || sum(df$mature == 1) < 15) return(tibble())
  fit <- tryCatch(suppressWarnings(glm(mature ~ lngt_cm, data = df, family = binomial)),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(tibble())
  b <- coef(fit)
  if (b[2] <= 0) return(tibble())
  l50  <- as.numeric(-b[1] / b[2])
  grad <- -c(1, l50) / b[2]                       # d(L50)/d(coef), delta method
  se   <- sqrt(as.numeric(t(grad) %*% vcov(fit) %*% grad))
  tibble(l50 = l50, se = se,
         lower = l50 - 1.96 * se, upper = l50 + 1.96 * se, n = nrow(df))
}

# split a per-year series into segments of consecutive years, for line-breaking at gaps
add_segments <- function(df) {
  df |>
    arrange(panel, grp, year) |>
    group_by(panel, grp) |>
    mutate(seg = cumsum(c(TRUE, diff(year) != 1))) |>
    ungroup()
}

# common colour scheme used across every figure (reference + three Roman rules)
key_levels <- c("Reference key", "Roman: II excluded",
                "Roman: I+II immature", "Roman: II mature")
key_pal <- c(
  "Reference key" = "#009E73",
  "Roman: II excluded" = "#000000",
  "Roman: I+II immature" = "#0072B2",
  "Roman: II mature" = "#D55E00"
)

# where the reference key is shown by scale, the scale is distinguished by point shape
scale_levels <- c("BITS (numeric)", "M6", "SMSF")
scale_shapes <- c("BITS (numeric)" = 16, "M6" = 17, "SMSF" = 15)

# 06 Sensitivity: stage II decision for two cod cases ----

cases <- tribble(
  ~country, ~species,       ~ices_area,
  "SE",     "Gadus morhua", "3.d.25",
  "DE",     "Gadus morhua", "3.d.24"
)

m_sub <- m_data |>
  semi_join(cases, by = c("country", "species", "ices_area"))

# colour = the four coding rules (reference pooled into one green series); the reference
# points additionally carry a shape for which scale they came from
dat_series <- bind_rows(
  assign_std(m_sub) |>
    mutate(series = "Reference key",
           scale_lab = recode(code_type, Numeric = "BITS (numeric)")),
  assign_roman(m_sub, NA_real_, "Roman: II excluded")  |> mutate(scale_lab = NA_character_),
  assign_roman(m_sub, 0,        "Roman: I+II immature") |> mutate(scale_lab = NA_character_),
  assign_roman(m_sub, 1,        "Roman: II mature")     |> mutate(scale_lab = NA_character_)
) |>
  filter(!is.na(mature), !is.na(lngt_cm))

l50_year <- dat_series |>
  group_by(country, species, ices_area, series, scale_lab, year) |>
  group_modify(~ fit_l50(.x)) |>
  ungroup() |>
  filter(l50 > 0, l50 < 120)

l50_plot <- l50_year |>
  mutate(
    series = factor(series, levels = key_levels),
    scale_lab = factor(scale_lab, levels = scale_levels),
    # plotmath label so the species name renders in italics (parsed by label_parsed)
    panel = paste0("atop('", country, " — '*italic('", species, "'), '", ices_area, "')"),
    grp = as.character(series)
  ) |>
  add_segments()

line_segs <- l50_plot |>
  group_by(panel, grp, seg) |>
  filter(n() > 1) |>
  ungroup()

p_sens <- ggplot(l50_plot, aes(x = year, y = l50, colour = series)) +
  geom_line(data = line_segs, aes(group = interaction(series, seg)), linewidth = 0.6) +
  geom_point(data = function(d) filter(d, series == "Reference key"),
             aes(shape = scale_lab), size = 1.8) +
  geom_point(data = function(d) filter(d, series != "Reference key"),
             shape = 16, size = 1.6) +
  scale_colour_manual(values = key_pal, name = NULL, drop = TRUE) +
  scale_shape_manual(values = scale_shapes, name = "Reference scale", drop = TRUE) +
  facet_wrap(~ panel, scales = "free", ncol = 2, labeller = label_parsed) +
  labs(x = "Year", y = expression(L[50] ~ "(cm)")) +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

ggsave(here("outputs/supp/check_bits_roman_sensitivity.png"), p_sens,
       width = 20, height = 12, units = "cm", dpi = 300)

# 07 Lithuania and Russia: I+II immature vs their reference keys ----
# RU and LT (excluded above) record almost no stage I, so stages I and II together stand
# in for the immature class. Dropping stage II or calling it mature leaves essentially no
# immature fish and cannot be fitted, so the only workable Roman rule is I+II immature.
# Here we check it reproduces the L50 from the same country's reference keys (M6, Numeric,
# SMSF) in the same species and ICES area.

rult <- ca |>
  filter(Country %in% c("RU", "LT"), Survey == "BITS") |>
  left_join(hh) |>
  janitor::clean_names() |>
  filter(!is.na(maturity), year >= 2000) |>
  mutate(maturity = as.character(maturity))

im_rult <- rult |>
  group_by(country, year) |>
  summarise(has_M = any(maturity == "M", na.rm = TRUE), .groups = "drop")

rult <- rult |>
  left_join(im_rult, by = c("country", "year")) |>
  mutate(
    code_type = case_when(
      grepl("^6[0-9]$", maturity) ~ "M6",
      maturity %in% c("A","B","Ba","Bb","C","Ca","Cb","D","Da","Db","E","F") ~ "SMSF",
      maturity %in% c("1","2","3","4","5","6") ~ "Numeric",
      maturity %in% c("I","M") & has_M ~ "IM",
      grepl("^[IVX]+$", maturity) & !has_M ~ "Roman",
      TRUE ~ "Unknown"
    )
  ) |>
  filter(!is.na(lon), !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_join(ices_sf |> select(ices_area = Area_27), left = TRUE) |>
  st_drop_geometry()

# per country, the species x area best covered by both Roman and a reference scale
rult_cases <- rult |>
  filter(!is.na(ices_area)) |>
  group_by(country, species, ices_area) |>
  summarise(
    n_roman = sum(code_type == "Roman"),
    n_other = sum(code_type %in% c("M6", "Numeric", "SMSF")),
    n_obs = n(),
    .groups = "drop"
  ) |>
  filter(n_roman > 0, n_other > 0) |>
  group_by(country) |>
  slice_max(n_obs, n = 2, with_ties = FALSE) |>
  ungroup() |>
  # RU flounder has too few Roman years to be informative; keep RU cod, LT cod, LT flounder
  filter(!(country == "RU" & species == "Platichthys flesus"))

cat("RU/LT cases (Roman + reference scale):\n")
print(rult_cases)

rult_sub <- rult |>
  semi_join(rult_cases, by = c("country", "species", "ices_area"))

# keep the specific reference key on each point (shape), pooling the black line;
# the numeric scale is the BITS numeric key
rult_series <- bind_rows(
  assign_std(rult_sub) |>
    mutate(series = "Reference key",
           scale_lab = recode(code_type, Numeric = "BITS (numeric)")),
  assign_roman(rult_sub, 0, "Roman: I+II immature") |>
    mutate(scale_lab = "Roman (I+II immature)")
) |>
  filter(!is.na(mature), !is.na(lngt_cm))

rult_l50 <- rult_series |>
  group_by(country, species, ices_area, series, scale_lab, year) |>
  group_modify(~ fit_l50(.x)) |>
  ungroup() |>
  filter(l50 > 0, l50 < 120)

rult_levels <- c("Reference key", "Roman: I+II immature")
rult_pal    <- key_pal[rult_levels]   # green reference, blue I+II immature (shared palette)

# reference keys shown by point shape (scale_levels/scale_shapes from the helpers); the Roman
# series is drawn as a fixed diamond so it stays out of the shape legend
rult_plot <- rult_l50 |>
  mutate(
    series = factor(series, levels = rult_levels),
    scale_lab = factor(scale_lab, levels = c(scale_levels, "Roman (I+II immature)")),
    # plotmath label so the species name renders in italics (parsed by label_parsed)
    panel = paste0("atop('", country, " — '*italic('", species, "'), '", ices_area, "')"),
    grp = as.character(series)
  ) |>
  add_segments()

rult_line_segs <- rult_plot |>
  group_by(panel, grp, seg) |>
  filter(n() > 1) |>
  ungroup()

p_rult <- ggplot(rult_plot, aes(x = year, y = l50, colour = series)) +
  geom_line(data = rult_line_segs, aes(group = interaction(series, seg)), linewidth = 0.6) +
  geom_point(data = function(d) filter(d, series == "Reference key"),
             aes(shape = scale_lab), size = 2) +
  geom_point(data = function(d) filter(d, series != "Reference key"),
             shape = 18, size = 2.3) +
  scale_colour_manual(values = rult_pal, name = NULL) +
  scale_shape_manual(values = scale_shapes, name = "Reference scale", drop = TRUE) +
  facet_wrap(~ panel, scales = "free", ncol = 2, labeller = label_parsed) +
  labs(x = "Year", y = expression(L[50] ~ "(cm)")) +
  theme_light(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

n_rult <- n_distinct(rult_plot$panel)
ggsave(here("outputs/supp/check_bits_rult_immature.png"), p_rult,
       width = 20, height = 2 + 6 * ceiling(n_rult / 2), units = "cm", dpi = 300)

# 08 One pooled estimate per key and species, across several ICES areas ----
# collapse to a single L50 (with a delta-method CI) per species and rule, in each ICES
# subdivision that carries enough species; restrict the reference records to each species'
# Roman coverage window so the pooled estimates are time-comparable. Records are pooled
# across countries within an area (RU and LT excluded, as they have their own rule).

# species x area carrying enough Roman and reference records
area_species <- m_data |>
  filter(!is.na(ices_area)) |>
  group_by(ices_area, species) |>
  summarise(
    n_roman = sum(code_type == "Roman"),
    n_other = sum(code_type %in% c("M6", "Numeric", "SMSF")),
    .groups = "drop"
  ) |>
  filter(n_roman >= 200, n_other >= 200)

# the subdivisions with the most such species
target_areas <- area_species |>
  count(ices_area, name = "n_species") |>
  filter(n_species >= 3) |>
  arrange(desc(n_species)) |>
  slice_head(n = 4) |>
  pull(ices_area)

cat("Areas used for the pooled comparison:\n")
print(target_areas)

multi_dat <- m_data |>
  filter(ices_area %in% target_areas) |>
  semi_join(area_species, by = c("ices_area", "species"))

# Roman coverage window per area x species, to time-match the reference records
roman_window <- multi_dat |>
  filter(code_type == "Roman") |>
  group_by(ices_area, species) |>
  summarise(y_min = min(year), y_max = max(year), .groups = "drop")

multi_pooled_dat <- multi_dat |>
  inner_join(roman_window, by = c("ices_area", "species")) |>
  filter(code_type == "Roman" | (year >= y_min & year <= y_max))

pooled_l50 <- bind_rows(
  assign_std(multi_pooled_dat)              |> mutate(series = "Reference key"),
  assign_roman(multi_pooled_dat, NA_real_, "Roman: II excluded"),
  assign_roman(multi_pooled_dat, 0,        "Roman: I+II immature"),
  assign_roman(multi_pooled_dat, 1,        "Roman: II mature")
) |>
  filter(!is.na(mature), !is.na(lngt_cm)) |>
  group_by(ices_area, species, series) |>
  group_modify(~ fit_l50(.x)) |>
  ungroup() |>
  filter(l50 > 0, l50 < 120)

# keep area x species with both a reference and a II-excluded estimate
pooled_keep <- pooled_l50 |>
  filter(series %in% c("Reference key", "Roman: II excluded")) |>
  count(ices_area, species) |>
  filter(n == 2) |>
  select(ices_area, species)

# 08b Quantitative summary: deviation of each Roman rule from the reference key ----
dev_summary <- pooled_l50 |>
  semi_join(pooled_keep, by = c("ices_area", "species")) |>
  select(ices_area, species, series, l50) |>
  pivot_wider(names_from = series, values_from = l50) |>
  transmute(
    ices_area, species,
    `Roman: II excluded`   = `Roman: II excluded`   - `Reference key`,
    `Roman: I+II immature` = `Roman: I+II immature` - `Reference key`,
    `Roman: II mature`     = `Roman: II mature`     - `Reference key`
  ) |>
  pivot_longer(-c(ices_area, species), names_to = "rule", values_to = "dev") |>
  filter(!is.na(dev)) |>
  group_by(rule) |>
  summarise(
    n_cases = n(),
    median_dev = round(median(dev), 1),
    median_abs_dev = round(median(abs(dev)), 1),
    max_abs_dev = round(max(abs(dev)), 1),
    .groups = "drop"
  ) |>
  arrange(median_abs_dev)

cat("\nDeviation from the reference-key L50 (cm), pooled per area x species:\n")
print(dev_summary)

# write a booktabs LaTeX table for the supplement
dev_tex <- dev_summary |>
  mutate(rule = recode(rule,
    "Roman: II excluded" = "Stage II excluded",
    "Roman: I+II immature" = "Stages I and II immature",
    "Roman: II mature" = "Stage II mature")) |>
  rename(
    `Stage-II rule` = rule,
    `\\textit{n}` = n_cases,
    `Median dev.\\ (cm)` = median_dev,
    `Median $|\\text{dev.}|$ (cm)` = median_abs_dev,
    `Max $|\\text{dev.}|$ (cm)` = max_abs_dev
  )

dev_caption <- paste0(
  "Deviation of each Roman-numeral stage-II rule from the reference-key $L_{50}$ ",
  "(BITS numeric, M6, SMSF), pooled per species and ICES subdivision. A positive ",
  "deviation means the rule overestimates $L_{50}$ relative to the reference key. ",
  "\\textit{n} is the number of species-by-subdivision comparisons."
)

dev_kable <- kable(
  dev_tex, format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",
  align = "lrrrr", caption = dev_caption, label = "roman_deviation"
) |>
  kable_styling(latex_options = "hold_position")

writeLines(as.character(dev_kable), here("outputs/supp/bits_roman_deviation_summary.tex"))

# 08c Figure: one pooled estimate per key and species ----

# order species by the reference L50 within each area (reorder-within via a keyed factor)
order_lvls <- pooled_l50 |>
  semi_join(pooled_keep, by = c("ices_area", "species")) |>
  filter(series == "Reference key") |>
  arrange(ices_area, l50) |>
  mutate(sp_area = paste(ices_area, species, sep = "___")) |>
  pull(sp_area)

pooled_plot <- pooled_l50 |>
  semi_join(pooled_keep, by = c("ices_area", "species")) |>
  mutate(
    series = factor(series, levels = key_levels),
    sp_area = factor(paste(ices_area, species, sep = "___"), levels = order_lvls)
  )

p_pooled <- ggplot(pooled_plot, aes(x = l50, y = sp_area, colour = series)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = 0.4, linewidth = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.7)) +
  scale_colour_manual(values = key_pal, name = NULL) +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  facet_wrap(~ ices_area, scales = "free", ncol = 2) +
  labs(x = expression(L[50] ~ "(cm)"), y = NULL) +
  theme_light(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "italic"),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

n_pooled_area <- n_distinct(pooled_plot$ices_area)
ggsave(here("outputs/supp/check_bits_pooled_by_key.png"), p_pooled,
       width = 22, height = 4 + 8 * ceiling(n_pooled_area / 2), units = "cm", dpi = 300)
