
library(sf)
library(purrr)
library(tidyverse)
library(here)
library(DATRASextra)
library(readxl)

sf_use_s2(FALSE)

# ── 1. Species present in DATRAS CA ──────────────────────────────────────────

raw_data <- readRDS(here("data/raw/raw_datras.rds")) %>%
  DATRASextra::correct_species()

species_info <- read_rds(here("data/metadata/species_info.rds")) %>%
  mutate(aphia_id = as.numeric(aphia_id))

# corrected scientific names actually observed in CA
dataset_species <- raw_data[["CA"]] %>%
  as_tibble() %>%
  select(Valid_Aphia, Species) %>%
  left_join(species_info %>% select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) %>%
  mutate(species_clean = coalesce(species, as.character(Species))) %>%
  filter(!is.na(species_clean),
         class %in% c("Teleostei", "Elasmobranchii")) %>%
  distinct(species_clean) %>%
  pull(species_clean)

# ── 2. ICES statistical areas ─────────────────────────────────────────────────

ices_sf <- read_sf(
  here("data/metadata/ices_areas/StatRec_map_Areas_Full_20170124.shp")
) %>%
  group_by(Area_27) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# ── 3. Species × ICES area combinations in the dataset ───────────────────────
# Derive from raw haul positions so the lookup is pre-computed for exactly
# those combinations that appear in 03_clean_maturity.R after the spatial join.

hh_sf <- raw_data[["HH"]] %>%
  as_tibble() %>%
  select(haul.id, lon, lat) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

hh_ices <- st_join(hh_sf, ices_sf %>% select(ices_area = Area_27), left = TRUE) %>%
  st_drop_geometry() %>%
  distinct(haul.id, ices_area) %>%
  filter(!is.na(ices_area))

ca_haul_species <- raw_data[["CA"]] %>%
  as_tibble() %>%
  select(haul.id, Valid_Aphia, Species) %>%
  left_join(species_info %>% select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) %>%
  mutate(species = coalesce(species, as.character(Species))) %>%
  filter(!is.na(species), class %in% c("Teleostei", "Elasmobranchii")) %>%
  distinct(haul.id, species)

species_area_combos <- ca_haul_species %>%
  inner_join(hh_ices, by = "haul.id") %>%
  distinct(species, ices_area)

# ── 4. GO-FISH: spawning months per species × ICES area ──────────────────────

gofish_raw <- read_sf(here("data/metadata/go_fish_spawning/GO-FISH-hs.shp")) %>%
  filter(species %in% dataset_species)

gofish_months <- gofish_raw %>%
  rowwise() %>%
  mutate(spawn_months = list(as.integer(which(c_across(starts_with("X")) > 0)))) %>%
  ungroup() %>%
  select(species, spawn_months, geometry)

spawn_gofish <- st_intersection(gofish_months,
                                ices_sf %>% select(ices_area = Area_27)) %>%
  st_drop_geometry() %>%
  filter(!is.na(ices_area), lengths(spawn_months) > 0) %>%
  group_by(species, ices_area) %>%
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    .groups      = "drop"
  ) %>%
  mutate(spawn_months = map(spawn_months, as.integer),
         source       = "gofish")


# ── 5. WKMAT07: spawning months per species × ICES area ──────────────────────

wkmat07_raw <- read_excel(here("data/metadata/spawning_season_wkmat07.xlsx"))

spawn_wkmat <- wkmat07_raw %>%
  rowwise() %>%
  mutate(spawn_months = list(as.integer(which(c_across(Jan:Dec) == 1)))) %>%
  ungroup() %>%
  filter(lengths(spawn_months) > 0) %>%
  select(Scientific_name, Area, spawn_months) %>%
  separate_rows(Area, sep = ",") %>%
  mutate(
    Area         = trimws(Area),
    n_months     = lengths(spawn_months)
  ) %>%
  # if a species × area pair appears more than once, keep the narrower window
  group_by(Scientific_name, Area) %>%
  slice_min(n_months, with_ties = FALSE) %>%
  ungroup() %>%
  select(-n_months) %>%
  rename(species = Scientific_name, ices_area = Area) %>%
  mutate(
    spawn_months = map(spawn_months, as.integer),
    source       = "wkmat"
  ) %>%
  filter(species %in% dataset_species)


# ── 6. Merge sources ─────────────────────────────────────────────────────────

spawn_direct <- bind_rows(spawn_gofish, spawn_wkmat) %>%
  filter(ices_area != "all") %>%
  group_by(species, ices_area) %>%
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    # will be "gofish", "wkmat", or "gofish+wkmat"
    source       = paste(sort(unique(source)), collapse = "+"),
    .groups      = "drop"
  )

# ── 7. Propagate months to all ancestor ICES areas ───────────────────────────

get_area_parents <- function(area) {
  parts <- strsplit(area, "\\.")[[1]]
  # returns the area itself AND all its ancestors, from least to most specific
  sapply(seq_along(parts), function(i) paste(parts[1:i], collapse = "."))
}

# Roll up child area data to every ancestor level (including the area itself)
spawn_all <- spawn_direct %>%
  mutate(ancestors = map(ices_area, get_area_parents)) %>%
  unnest(ancestors) %>%
  group_by(species, ices_area = ancestors) %>%
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    source       = paste(sort(unique(source)), collapse = "+"),
    .groups      = "drop"
  )

# Region-level fallback (just the leading integer of the ICES code)
get_region <- function(area) as.numeric(sub("\\..*", "", area))

spawn_region <- spawn_all %>%
  mutate(region = get_region(ices_area)) %>%
  group_by(species, region) %>%
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    source       = paste(sort(unique(source)), collapse = "+"),
    .groups      = "drop"
  )

# ── 8. Resolver ───────────────────────────────────────────────────────────────

get_area_hierarchy <- function(area) {
  # from most specific to least specific (reverse of get_area_parents)
  parts <- strsplit(area, "\\.")[[1]]
  sapply(length(parts):1, function(i) paste(parts[1:i], collapse = "."))
}

resolve_spawn <- function(sp, area, spawn_direct, spawn_all, spawn_region) {

  hierarchy <- get_area_hierarchy(area)

  for (a in hierarchy) {

    hit_all    <- filter(spawn_all,    species == sp, ices_area == a)
    if (nrow(hit_all) == 0) next

    hit_direct <- filter(spawn_direct, species == sp, ices_area == a)
    has_direct <- nrow(hit_direct) > 0

    src <- if (a == area) {
      # Found at the exact requested area
      if (has_direct) hit_direct$source[1]
      else            paste0("child_of:", a, "|", hit_all$source[1])
    } else {
      # Found at a parent area
      if (has_direct) paste0("parent_area:", a, "|", hit_direct$source[1])
      else            paste0("child_of_parent:", a, "|", hit_all$source[1])
    }

    return(list(
      spawn_months = hit_all$spawn_months[[1]],
      source       = src,
      source_area  = a
    ))
  }

  # Last resort: borrow from the nearest numbered ICES region
  target_region   <- get_region(area)
  regions_to_try  <- c(target_region,
                       target_region - 1, target_region + 1,
                       target_region - 2, target_region + 2)
  regions_to_try  <- regions_to_try[regions_to_try > 0]

  region_hit <- spawn_region %>%
    filter(species == sp, region %in% regions_to_try) %>%
    mutate(dist = abs(region - target_region)) %>%
    slice_min(dist, n = 1, with_ties = FALSE)

  if (nrow(region_hit) == 0)
    return(list(spawn_months  = NULL,
                source        = NA_character_,
                source_area   = NA_character_))

  list(
    spawn_months = sort(unique(unlist(region_hit$spawn_months))),
    source       = paste0("borrowed_region:", region_hit$region[1],
                          "|", region_hit$source[1]),
    source_area  = as.character(region_hit$region[1])
  )
}

# ── 9. Build the lookup ───────────────────────────────────────────────────────

spawn_lookup <- species_area_combos %>%
  mutate(
    resolved = map2(
      species, ices_area,
      resolve_spawn,
      spawn_direct = spawn_direct,
      spawn_all    = spawn_all,
      spawn_region = spawn_region
    ),
    spawn_months = map(resolved, "spawn_months"),
    source       = map_chr(resolved, ~ .x$source  %||% NA_character_),
    source_area  = map_chr(resolved, ~ .x$source_area %||% NA_character_)
  ) %>%
  select(-resolved)

# ── 10. Extend spawning months by 3-month pre-spawn window ───────────────────
# include the 3 calendar months immediately before first
# spawning month (wraps around December → January).

add_pre_spawn <- function(months) {
  if (is.null(months) || length(months) == 0) return(months)
  pre <- ((months - 4) %% 12) + 1
  sort(unique(c(months, pre)))
}

spawn_lookup <- spawn_lookup %>%
  mutate(spawn_months_ext = map(spawn_months, add_pre_spawn))

# ── 11. Summary ───────────────────────────────────────────────────────────────

source_summary <- spawn_lookup %>%
  mutate(
    source_type = case_when(
      source %in% c("gofish", "wkmat", "gofish+wkmat") ~ source,
      grepl("^child_of_parent", source)                ~ "child_of_parent_area",
      grepl("^child_of:", source)                      ~ "child_areas_rollup",
      grepl("^parent_area:", source)                   ~ "parent_area",
      grepl("^borrowed_region:", source)               ~ "borrowed_region",
      is.na(source)                                    ~ "no_data",
      TRUE                                             ~ "other"
    )
  ) %>%
  count(source_type, name = "n_entries")

# ── 12. Save ─────────────────────────────────────────────────────────────────

write_rds(spawn_lookup, here("data/metadata/spawning_lookup.rds"))
