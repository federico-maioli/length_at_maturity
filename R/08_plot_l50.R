library(tidyverse)
library(here)
library(rnaturalearth)
library(ggrepel)
library(paletteer)
library(rphylopic)
library(sf)
library(patchwork)
library(ggstats)
library(geomtextpath)

# 01 Load data ----

raw <- read_rds(here("data/final/l50_clean.rds")) |>
  rename(l50 = l50_est, l25 = l25_est, l75 = l75_est)
data_raw <- read_rds(here("data/intermediate/maturity_clean.rds"))

# shared colour palette
col_mid <- "#756BB1"
col_dark <- "#3F007D"
col_light <- "#BCBDDC"

# 02 Panel a: l50 by stock ----

coeff <- raw |>
  filter(model == "Stock", sex == "Combined", is.na(period)) |>
  mutate(
    parts = strsplit(species_stock, "_"),
    genus = sapply(parts, `[`, 1),
    sp = sapply(parts, `[`, 2),
    stock = sapply(parts, function(x) if (length(x) > 2) paste(x[3:length(x)], collapse = "_") else NA_character_),
    species_name = paste(genus, sp),
    # build italic expression label, appending stock code when present
    label = ifelse(
      is.na(stock),
      paste0("italic('", species_name, "')"),
      paste0("italic('", species_name, "')~' [", stock, "]'")
    )
  ) |>
  group_by(species_name) |>
  mutate(mean_l50 = mean(l50, na.rm = TRUE)) |>
  ungroup() |>
  arrange(desc(mean_l50), desc(l50)) |>
  mutate(label = factor(label, levels = rev(unique(label))))

pA <- ggplot(coeff, aes(x = l50, y = label)) +
  ggstats::geom_stripped_rows(odd = "white", even = "grey90", alpha = 0.3) +
  geom_errorbarh(aes(xmin = l50_lower, xmax = l50_upper),
                 color = col_mid, height = 0, linewidth = 0.5) +
  geom_hline(
    data = coeff |> group_by(species_name) |> summarise(y = min(as.numeric(label))),
    aes(yintercept = y - 0.5), color = "grey70", linewidth = 0.3
  ) +
  geom_point(color = col_mid, fill = "white", shape = 21, size = 1.4, stroke = 0.8) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(L[50]~"(cm)"), y = NULL) +
  theme_light(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(size = 8, lineheight = 0.9),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.minor.x = element_blank()
  )

# 03 Spatial setup ----

ices_sf <- read_sf(here("data/metadata/ices_areas/ICES_Statrec_mapto_ICES_Areas.shp"))
coast <- ne_countries(scale = "medium", returnclass = "sf")

# pre-computed lookup of species x ICES area from the cleaned maturity data
area_lookup <- readRDS(here("data/intermediate/maturity_clean.rds")) |>
  dplyr::select(species, ices_area, species_stock) |>
  distinct()

# 04 Panel b: cod stock map ----

cod_sf <- coeff |>
  select(-ices_area) |>
  filter(species_name == "Gadus morhua") |>
  left_join(area_lookup |> filter(species == "Gadus morhua") |>
              distinct(ices_area, species_stock), by = "species_stock") |>
  left_join(ices_sf, by = c("ices_area" = "Area_27")) |>
  filter(!is.na(l50)) |>
  group_by(species_stock, stock, l50) |>
  summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_as_sf() |>
  st_transform(3035)

centroids <- st_point_on_surface(cod_sf)
bbox <- st_bbox(cod_sf)
coast_3035 <- st_transform(coast, 3035)
buf <- 300000

cod_global_l50 <- raw |>
  filter(species_clean == "Gadus_morhua", model == "Global",
         sex == "Combined", is.na(period)) |>
  pull(l50)

p_cod <- ggplot() +
  geom_sf(data = cod_sf, aes(fill = l50), color = NA) +
  geom_sf(data = coast_3035, fill = "grey90", color = NA, linewidth = 0.05) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(geometry = geometry, label = stock),
    stat = "sf_coordinates",
    size = 3.2, force = 100, box.padding = 0.9,
    min.segment.length = 0.1, seed = 1,
    color = "black", bg.color = "white", bg.r = 0.2
  ) +
  paletteer::scale_fill_paletteer_c("ggthemes::Purple", direction = 1,
                                    name = expression(L[50]~"(cm)"),
                                    breaks = scales::breaks_pretty(n = 5)) +
  coord_sf(crs = 3035,
           xlim = c(bbox$xmin, bbox$xmax),
           ylim = c(bbox$ymin - buf, bbox$ymax + buf),
           expand = FALSE) +
  labs(title = expression(bolditalic("Gadus morhua")), x = NULL, y = NULL) +
  guides(fill = guide_colorbar(
    direction = "horizontal", title.position = "left",
    barwidth = 6, barheight = 0.3, ticks.colour = "white",
    title.theme = element_text(size = 8, vjust = 1)
  )) +
  theme_light(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = c(0.99, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7)
  )

# 05 Panel c: cod.27.46a7d20 time series ----

period_levels <- c("2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2024")

cod_ts <- raw |>
  filter(species_clean == "Gadus_morhua", model == "Stock",
         sex == "Combined", !is.na(period),
         grepl("cod.27.46a7d20", species_stock)) |>
  mutate(period = factor(period, levels = period_levels)) |>
  filter(!is.na(period))

pC <- ggplot(cod_ts, aes(x = period, y = l50, group = 1)) +
  geom_ribbon(aes(ymin = l50_lower, ymax = l50_upper),
              fill = col_mid, alpha = 0.2, color = NA) +
  geom_line(color = col_mid, linewidth = 0.6) +
  geom_point(color = col_mid, size = 1.6, shape = 21, fill = "white", stroke = 0.8) +
  labs(x = NULL, y = expression(L[50]~"(cm)"), title = "cod.27.46a7d20") +
  theme_light(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
  )

# 06 Panel d: maturity ogive by sex for cod.27.46a7d20 ----

sex_colors <- c("F" = col_dark, "M" = col_light)

# length grid for the fitted curves (starts at 0 so the ogive spans the full axis)
x_seq <- seq(0, 100, length.out = 300)

# raw individual observations for the cod.27.46a7d20 stock, by sex
cod21_obs <- data_raw |>
  filter(grepl("cod.27.46a7d20", species_stock), sex %in% c("M", "F"),
         !is.na(mature), !is.na(lngt_cm))

# refit the logistic maturity model per sex on the raw observations
cod21_models <- cod21_obs |>
  nest(data = -sex) |>
  mutate(fit = map(data, ~ glm(mature ~ lngt_cm, data = .x, family = binomial)))

# predicted ogive per sex over the length grid
cod21_fits <- cod21_models |>
  mutate(pred = map(fit, ~ tibble(
    lngt_cm = x_seq,
    p = predict(.x, tibble(lngt_cm = x_seq), type = "response")
  ))) |>
  select(sex, pred) |>
  unnest(pred)

# L50 per sex from the refitted models (for the dashed reference lines)
l50_lines <- cod21_models |>
  mutate(l50 = map_dbl(fit, ~ -coef(.x)[1] / coef(.x)[2])) |>
  select(sex, l50) |>
  mutate(label = sprintf("%.1f", l50))

# bin raw observations to overlay proportions as points
cod21_binned <- cod21_obs |>
  rename(model_sex = sex) |>
  mutate(lngt_bin = floor(lngt_cm / 3) * 3 + 1.5) |>
  group_by(model_sex, lngt_bin) |>
  summarise(n = n(), prop = mean(mature, na.rm = TRUE), .groups = "drop") |>
  filter(n >= 5)

pD <- ggplot() +
  geom_line(data = cod21_fits,
            aes(x = lngt_cm, y = p, color = sex), linewidth = 0.7) +
  geom_point(data = cod21_binned,
             aes(x = lngt_bin, y = prop, color = model_sex, size = n), alpha = 0.7) +
  geom_textvline(data = l50_lines,
                 aes(xintercept = l50, color = sex, label = label),
                 linetype = "dashed", linewidth = 0.4,
                 hjust = 0.08, vjust = 0.5, size = 3,
                 show.legend = FALSE) +
  scale_color_manual(values = sex_colors,
                     labels = c("F" = "Female", "M" = "Male"), name = NULL) +
  scale_size_continuous(range = c(0.8, 3.5), guide = "none") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Total length (cm)", y = "Proportion mature", title = "cod.27.46a7d20") +
  theme_light(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = c(0.99, 0.03),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
  )

# 07 Combine and save ----

pA + (p_cod / pC / pD + plot_layout(heights = c(2, 0.8, 2))) +
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1, 1.3)) &
  theme(plot.tag = element_text(size = 10))

ggsave(
  here("outputs/main/data_examples.png"),
  width = 24, height = 24, units = "cm", dpi = 600
)
