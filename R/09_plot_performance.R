library(tidyverse)
library(here)

# Model performance of the retained maturity ogives: discrimination (AUC) and fit
# (Tjur R2) by spatial aggregation level. Reads the cleaned estimates written by
# 05_clean_estimates.R and produces the supplementary performance figure.

# 01 Load data ----

data_clean <- read_rds(here("data/final/l50_clean.rds"))

# 02 Fit quality by aggregation level ----

fit_quality <- data_clean |>
  ungroup() |>
  select(model, auc, tjur_r2) |>
  pivot_longer(c(auc, tjur_r2), names_to = "metric", values_to = "value") |>
  mutate(
    level = factor(model, levels = c("Global", "Stock", "Area")),
    metric = recode(metric, auc = "AUC", tjur_r2 = "Tjur~italic(R)^2")
  )

# median fit quality per level (for the text)
fit_quality_summary <- data_clean |>
  ungroup() |>
  mutate(level = factor(model, levels = c("Global", "Stock", "Area"))) |>
  group_by(level) |>
  summarise(
    n = n(),
    median_auc = median(auc, na.rm = TRUE),
    median_tjur_r2 = median(tjur_r2, na.rm = TRUE),
    .groups = "drop"
  )
print(fit_quality_summary)

# 03 Plot ----

p_quality <- ggplot(fit_quality, aes(x = level, y = value)) +
  geom_violin(fill = "#BCBDDC", colour = NA, alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.4, fill = "white") +
  facet_wrap(~ metric, scales = "free_y", labeller = label_parsed) +
  labs(x = "Aggregation level", y = NULL) +
  theme_light() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = "black", size = 11))

ggsave(here("outputs/supp/model_performance.png"), p_quality,
       width = 18, height = 9, units = "cm", dpi = 600)
