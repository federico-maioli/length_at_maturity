library(tidyverse)
library(here)
library(rfishbase)
library(worrms)
library(taxize)

# safely wrap the WoRMS classification call so failures return NULL
safe_wm <- safely(wm_classification)

# 01 Load species from datras ----

data <- readRDS(here("data/raw/raw_datras.rds")) |> DATRASextra::correct_species()

data <- data[["HL"]]

sp <- data |> distinct(Species, Valid_Aphia)

sp <- sp |>
  janitor::clean_names() |>
  # keep only names with exactly two words
  filter(str_count(species, "\\S+") == 2) |>
  # genus capitalised, species lowercase
  filter(str_detect(species, "^[A-Z][a-z]+ [a-z]+$")) |>
  droplevels()

# 02 Get accepted names and taxonomy ----

# unique aphia ids from the list
out <- id2name(sp$valid_aphia, db = "worms")

species_lookup <- out |>
  map_df(~ as.data.frame(.x), .id = "input_id") |>
  select(valid_aphia = id, accepted_name = name) |>
  distinct()

class_results <- map(as.numeric(species_lookup$valid_aphia), safe_wm)

# name the list by aphia id so results can be tracked
names(class_results) <- species_lookup$valid_aphia

taxonomy_wide <- class_results |>
  map(~ .x$result) |>
  compact() |>
  bind_rows(.id = "aphia_id") |>
  filter(rank %in% c("Class", "Order", "Family", "Genus", "Species")) |>
  select(aphia_id, rank, scientificname) |>
  pivot_wider(names_from = rank, values_from = scientificname)

# 03 Keep bony fish and sharks ----

sp <- taxonomy_wide |>
  filter(Class %in% c("Teleostei", "Elasmobranchii")) |>
  janitor::clean_names()

# 04 Get traits from fishbase ----

my_sp_names <- sp$species

# maximum length (Lmax) from the species table
df_max <- species(my_sp_names) |>
  select(Species, max_l = Length)

# length at maturity (L50) from the maturity table, summarised across studies
df_l50 <- maturity(my_sp_names) |>
  select(Species, Lm) |>
  filter(!is.na(Lm)) |>
  group_by(Species) |>
  summarize(l_50_median = median(Lm, na.rm = TRUE),
            l_50_min = min(Lm, na.rm = TRUE),
            l_50_max = max(Lm, na.rm = TRUE))

# 05 Join and save ----

sp_final <- sp |>
  left_join(df_max, by = c("species" = "Species")) |>
  left_join(df_l50, by = c("species" = "Species"))

write_rds(sp_final, here("data/metadata/species_info.rds"))
