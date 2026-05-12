# libraries
library(tidyverse)
library(here)
library(purrr)
library(rfishbase)
library(worrms)
library(taxize)

data <- readRDS(here('data/raw/raw_datras.rds')) %>% DATRASextra::correct_species()

data <- data[['HL']]

sp <- data |> distinct(Species,Valid_Aphia)

sp <- sp %>% janitor::clean_names() |> 
  # keep only names with exactly two words
  filter(str_count(species, "\\S+") == 2) %>%
  # proper formatting (Genus capitalized, species lowercase)
  filter(str_detect(species, "^[A-Z][a-z]+ [a-z]+$")) |> droplevels()

# get the unique AphiaIDs from your list
out <- id2name(sp$valid_aphia, db = "worms")

species_lookup <- out %>%
  # This map handles the list structure safely
  map_df(~ as.data.frame(.x), .id = "input_id") %>%
  select(
    valid_aphia = id,
    accepted_name = name
  ) |> distinct()

# Map the function over your IDs
safe_wm <- safely(wm_classification)

class_results <- map(as.numeric(species_lookup$valid_aphia), safe_wm)

# set the names of the list to the AphiaIDs so we can track them
names(class_results) <- species_lookup$valid_aphia

taxonomy_wide <- class_results %>%
  map(~ .x$result) %>%
  compact() %>% 
  bind_rows(.id = "aphia_id") %>%
  filter(rank %in% c("Class", "Order", "Family", "Genus", "Species")) %>%
  select(aphia_id, rank, scientificname) %>%
  pivot_wider(names_from = rank, values_from = scientificname)

# filter for bony fish and sharks

sp <- taxonomy_wide |> filter(Class %in% c('Teleostei','Elasmobranchii')) |> janitor::clean_names()

# now get l inf etc form rfishbase
my_sp_names <- sp$species

# Pull Lmax (from the species table)
# FishBase uses 'Length' as the common maximum length
df_max <- species(my_sp_names) %>%
  select(Species, max_l = Length)

# Pull Linf (from the popgrowth table)
# We take the median Linf (Loo) to avoid outliers from a single bad study
df_linf <- popgrowth(my_sp_names) %>%
  select(Species, Loo) %>%
  filter(!is.na(Loo)) %>%
  group_by(Species) %>%
  summarize(l_inf_median = median(Loo, na.rm = TRUE),l_inf_min = min(Loo, na.rm = TRUE),l_inf_max = max(Loo, na.rm = TRUE))

# Pull L50 (from the maturity table)
# FishBase column 'lm' is the length at first maturity (L50)
df_l50 <- maturity(my_sp_names) %>%
  select(Species, Lm) %>%
  filter(!is.na(Lm)) %>%
  group_by(Species) %>%
  summarize(l_50_median = median(Lm, na.rm = TRUE), l_50_min = min(Lm, na.rm = TRUE),l_50_max = max(Lm, na.rm = TRUE))

# join everything back to my original taxonomy table
sp_final <- sp %>%
  left_join(df_max, by = c("species" = "Species")) %>%
  left_join(df_linf, by = c("species" = "Species")) %>%
  left_join(df_l50, by = c("species" = "Species"))

# flag species where linf is less or more tha 95 % lmax
write_rds(sp_final, here('data/metadata/species_info.rds'))
