library(tidyverse)
library(here)
library(ggrepel)

# Internal sanity check: compare our global (all-year, combined-sex) L50 per species with
# the median length at maturity reported in FishBase. This is an interspecific comparison,
# so it is dominated by the wide range of body sizes across species; it is meant only to
# confirm our estimates are the right order of magnitude, not as a formal validation.

# 01 Load data ----

clean <- read_rds(here("data/final/l50_clean.rds"))

# FishBase median (and range) length at maturity per species
fb <- read_rds(here("data/metadata/species_info.rds")) |>
  mutate(species_clean = str_replace_all(species, " ", "_")) |>
  select(species_clean, l_50_median, l_50_min, l_50_max)

# 02 Our global estimates ----

our_global <- clean |>
  ungroup() |>
  filter(model == "Global", sex == "Combined", is.na(period)) |>
  select(species_clean, l50_est)

# 03 Join our estimates with FishBase ----

compare <- our_global |>
  inner_join(fb, by = "species_clean") |>
  filter(!is.na(l_50_median)) |>
  mutate(
    species = str_replace_all(species_clean, "_", " "),
    # label only the well-separated species (FishBase L50 >= 20 cm) plus any clear
    # off-diagonal outliers; the dense small-species cluster is left unlabelled
    lab = if_else(l_50_median >= 20 | abs(l50_est - l_50_median) > 8, species, NA_character_)
  )

r_pearson  <- cor(compare$l50_est, compare$l_50_median)
r_spearman <- cor(compare$l50_est, compare$l_50_median, method = "spearman")

cat(sprintf("Global L50 vs FishBase: n = %d species, Pearson r = %.2f, Spearman = %.2f\n",
            nrow(compare), r_pearson, r_spearman))

# 04 Correlation plot ----

lim <- c(0, max(compare$l50_est, compare$l_50_median, na.rm = TRUE) * 1.05)

p_fb <- ggplot(compare, aes(x = l_50_median, y = l50_est)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(shape = 21, fill = "#756BB1", colour = "white", size = 2.4, stroke = 0.6) +
  ggrepel::geom_text_repel(aes(label = lab), size = 2.4, fontface = "italic",
                           colour = "grey30", na.rm = TRUE, seed = 1,
                           max.overlaps = Inf, box.padding = 0.5, force = 5,
                           force_pull = 0.2, min.segment.length = 0,
                           segment.colour = "grey75", segment.size = 0.2) +
  # Pearson r and sample size annotated in the empty upper-left corner
  annotate("text", x = lim[1] + 0.03 * diff(lim), y = lim[2] * 0.98,
           hjust = 0, vjust = 1, size = 4.2,
           label = sprintf("italic(r) == %.2f", r_pearson), parse = TRUE) +
  annotate("text", x = lim[1] + 0.03 * diff(lim), y = lim[2] * 0.91,
           hjust = 0, vjust = 1, size = 3.4,
           label = sprintf("italic(n) == %d~species", nrow(compare)), parse = TRUE) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(
    x = expression("FishBase median " * italic(L)[50] * " (cm)"),
    y = expression("Our global " * italic(L)[50] * " (cm)")
  ) +
  theme_light(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(here("outputs/main/check_l50_vs_fishbase.png"), p_fb,
       width = 16, height = 16, units = "cm", dpi = 300)
