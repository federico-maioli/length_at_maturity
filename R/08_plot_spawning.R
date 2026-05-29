library(tidyverse)
library(here)

# 1. load data ----

spawn_lookup <- read_rds(here("data/metadata/spawning_lookup.rds"))
l50_clean <- read_rds(here("data/final/l50_clean.rds"))

# species with at least one validated L50 estimate
target_species <- l50_clean |>
  distinct(species_clean) |>
  pull(species_clean) |>
  str_replace_all("_", " ")

# 2. build species x month matrix ----

# keep only direct sources (gofish, wkmat, gofish+wkmat)
spawn_long <- spawn_lookup |>
  filter(
    source %in% c("gofish", "wkmat", "gofish+wkmat"),
    species %in% target_species,
    lengths(spawn_months) > 0
  ) |>
  select(species, ices_area, spawn_months) |>
  unnest_longer(spawn_months, values_to = "month")

# total distinct ICES areas per species (denominator for proportion)
species_n_areas <- spawn_long |>
  group_by(species) |>
  summarise(n_areas = n_distinct(ices_area), .groups = "drop")

# proportion of ICES areas where each species spawns in each month
spawn_summary <- spawn_long |>
  group_by(species, month) |>
  summarise(n_spawn = n_distinct(ices_area), .groups = "drop") |>
  complete(species, month = 1:12, fill = list(n_spawn = 0L)) |>
  left_join(species_n_areas, by = "species") |>
  mutate(
    prop = n_spawn / n_areas,
    species = fct_rev(factor(species)),   # alphabetical, A at top
    month = factor(month, levels = 1:12, labels = month.abb)
  )

# 3. plot ----

# purple palette matching 06_plot_l50 (white = no spawning, capped at medium-dark)
purple_pal <- c("white", "#BCBDDC", "#756BB1", "#54278F")

ggplot(spawn_summary, aes(x = month, y = species, fill = prop)) +
  geom_tile(color = "grey92", linewidth = 0.25) +
  scale_fill_gradientn(
    colors = purple_pal,
    values = c(0, 0.01, 0.5, 1),   # white only at exactly 0
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = "Proportion of \nICES areas"
  ) +
  scale_x_discrete(expand = c(0, 0), position = "top") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9, face = "italic", hjust = 1),
    axis.text.x = element_text(size = 10),
    panel.grid = element_blank(),
    legend.key.height = unit(0.9, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_colorbar(ticks.colour = "white"))

# 4. save ----

ggsave(
  here("outputs/supp/spawning_heatmap.png"),
  width = 18,
  height = 15,
  units = "cm",
  dpi = 450,
  limitsize = FALSE
)
