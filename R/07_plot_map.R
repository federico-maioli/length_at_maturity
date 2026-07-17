library(here)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

# 1. load and prepare haul-level data ----

# keep only variables needed for mapping
data <- readRDS(here("data/intermediate/maturity_clean.rds")) |>
  select(
    survey, quarter, country, haul_id, year, month,
    lon, lat, stat_rec, species, valid_aphia, ices_area,
    stock, species_clean, species_stock, lngt_cm, mature
  )

# one row per haul; n_obs = sampled fish records, n_species = distinct species
haul_sf <- data |>
  group_by(haul_id, lon, lat) |>
  summarise(
    n_obs = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(3035)

# 2. coastline background ----

world <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035)

# 3. main map ----

p1 <- ggplot() +
  geom_sf(data = world, fill = "gray90", linewidth = 0) +
  # point size = observations per haul, fill = species richness
  geom_sf(
    data = haul_sf,
    aes(size = n_obs, fill = n_species),
    shape = 21, color = "white", stroke = 0.2, alpha = 0.5
  ) +
  scale_size_continuous(
    range = c(0.5, 4),
    name = "Samples\nper haul",
    breaks = c(500, 100, 50, 10),
    guide = guide_legend(override.aes = list(shape = 19, color = "black"))
  ) +
  scale_fill_viridis_c(option = "C", name = "Number of\nspecies sampled") +
  coord_sf(
    crs = st_crs(3035),
    xlim = c(2.3e6, 5.4e6),
    ylim = c(1.6e6, 4.7e6)
  ) +
  theme_light() +
  theme(
    legend.position = c(0.80, 0.3),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.1),
    panel.grid.minor = element_blank()
  )

# 4. time-series inset ----

# total sampled records per year
ts_df <- data |>
  group_by(year) |>
  summarise(n_obs = n(), .groups = "drop")

p2 <- ggplot(ts_df, aes(x = year, y = n_obs)) +
  geom_col(fill = scales::viridis_pal(option = "C")(1), width = 0.8) +
  labs(x = "Year", y = "Samples") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), panel.background = element_blank())

# 5. combine and save ----

final_plot <- p1 +
  inset_element(p2, left = 0, bottom = 0.8, right = 0.45, top = 0.99)

ggsave(
  here("output/map.png"),
  plot = final_plot,
  bg = "white",
  width = 18,
  height = 13,
  units = "cm",
  dpi = 600
)
