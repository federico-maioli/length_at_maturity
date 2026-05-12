# Libraries ----------------------------------------------------------------

library(here)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

# Load and prepare haul-level data -----------------------------------------

# Read the maturity dataset and keep only the variables needed here
data <- readRDS(here("data/intermediate/maturity_clean.rds")) %>%
  select(
    survey, quarter, country, haul_id, year, month,
    lon, lat, stat_rec, species, valid_aphia, ices_area,
    stock, species_clean, species_stock, lngt_cm, mature
  )

# Collapse observations to one row per haul location
# n_obs = total number of sampled fish records in the haul
# n_species = number of distinct species observed in the haul
haul_sf <- data %>%
  group_by(haul_id, lon, lat) %>%
  summarise(
    n_obs = n(),
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3035)

# Load coastline / land background -----------------------------------------

# Natural Earth world polygons
world <- ne_countries(scale = "medium", returnclass = "sf")

# Reproject to the same CRS as the haul data
world_laea <- world %>%
  st_transform(3035)

# Main map -----------------------------------------------------------------

p1 <- ggplot() +
  
  # Background land polygons
  geom_sf(
    data = world_laea,
    fill = "gray90",
    linewidth = 0
  ) +
  
  # Haul locations
  # point size shows number of observations per haul
  # fill color shows number of species sampled in that haul
  geom_sf(
    data = haul_sf,
    aes(size = n_obs, fill = n_species),
    shape = 21,
    color = "white",
    stroke = 0.2,
    alpha = 0.5
  ) +
  
  # Size legend for haul sample size
  scale_size_continuous(
    range = c(0.5, 4),
    name = "Samples\nper haul",
    breaks = c(500, 100, 50, 10),
    guide = guide_legend(
      override.aes = list(shape = 19, color = "black")
    )
  ) +
  
  # Fill legend for species richness per haul
  scale_fill_viridis_c(
    option = "C",
    name = "Number of\nspecies sampled"
  ) +
  
  # Map extent and projection
  coord_sf(
    crs = st_crs(3035),
    xlim = c(2.3e6, 5.4e6),
    ylim = c(1.6e6, 4.7e6)
  ) +
  
  # Theme and legend styling
  theme_light() +
  theme(
    legend.position = c(0.80, 0.3),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.1
    ),
    panel.grid.minor = element_blank()
  )

# Time-series inset --------------------------------------------------------

# Total number of sampled records per year
ts_df <- data %>%
  group_by(year) %>%
  summarise(n_obs = n(), .groups = "drop")

# Small inset barplot
p2 <- ggplot(ts_df, aes(x = year, y = n_obs)) +
  geom_col(
    fill = scales::viridis_pal(option = "C")(1),
    width = 0.8
  ) +
  labs(
    x = "Year",
    y = "Samples"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# Combine map and inset ----------------------------------------------------

final_plot <- p1 +
  inset_element(
    p2,
    left = 0,
    bottom = 0.8,
    right = 0.45,
    top = 0.99
  )

# Save output --------------------------------------------------------------

ggsave(
  here("output/map.png"),
  plot = final_plot,
  bg = "white",
  width = 18,
  height = 13,
  units = "cm",
  dpi = 600
)

