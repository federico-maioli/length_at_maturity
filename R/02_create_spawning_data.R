
library(sf)
library(purrr)
library(tidyverse)
library(here)
library(DATRASextra)
library(readxl)
library(rfishbase)

sf_use_s2(FALSE)

# 1. species present in datras ca ----

raw_data <- readRDS(here("data/raw/raw_datras.rds")) |>
  DATRASextra::correct_species()

species_info <- read_rds(here("data/metadata/species_info.rds")) |>
  mutate(aphia_id = as.numeric(aphia_id))

# fish and elasmobranch species observed in CA with valid names
dataset_species <- raw_data[["CA"]] |>
  as_tibble() |>
  select(Valid_Aphia, Species) |>
  left_join(species_info |> select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) |>
  mutate(species_clean = coalesce(species, as.character(Species))) |>
  filter(!is.na(species_clean),
         class %in% c("Teleostei", "Elasmobranchii")) |>
  distinct(species_clean) |>
  pull(species_clean)

# 2. ices statistical areas ----

ices_sf <- read_sf(
  here("data/metadata/ices_areas/StatRec_map_Areas_Full_20170124.shp")
) |>
  group_by(Area_27) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

# 3. species x ices area combinations in the dataset ----
# pre-compute the exact combinations that appear after the spatial join in 03_clean_maturity.R

hh_sf <- raw_data[["HH"]] |>
  as_tibble() |>
  select(haul.id, lon, lat) |>
  filter(!is.na(lon), !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# spatial join to assign each haul to an ICES area
hh_ices <- st_join(hh_sf, ices_sf |> select(ices_area = Area_27), left = TRUE) |>
  st_drop_geometry() |>
  distinct(haul.id, ices_area) |>
  filter(!is.na(ices_area))

ca_haul_species <- raw_data[["CA"]] |>
  as_tibble() |>
  select(haul.id, Valid_Aphia, Species) |>
  left_join(species_info |> select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) |>
  mutate(species = coalesce(species, as.character(Species))) |>
  filter(!is.na(species), class %in% c("Teleostei", "Elasmobranchii")) |>
  distinct(haul.id, species)

species_area_combos <- ca_haul_species |>
  inner_join(hh_ices, by = "haul.id") |>
  distinct(species, ices_area)

# 4. go-fish: spawning months per species x ices area ----

gofish_raw <- read_sf(here("data/metadata/go_fish_spawning/GO-FISH-hs.shp")) |>
  filter(species %in% dataset_species)

# convert binary month columns to a list of active month indices
gofish_months <- gofish_raw |>
  rowwise() |>
  mutate(spawn_months = list(as.integer(which(c_across(starts_with("X")) > 0)))) |>
  ungroup() |>
  select(species, spawn_months, geometry)

spawn_gofish <- st_intersection(gofish_months,
                                ices_sf |> select(ices_area = Area_27)) |>
  st_drop_geometry() |>
  filter(!is.na(ices_area), lengths(spawn_months) > 0) |>
  group_by(species, ices_area) |>
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    .groups      = "drop"
  ) |>
  mutate(spawn_months = map(spawn_months, as.integer),
         source       = "gofish")

# 5. wkmat07: spawning months per species x ices area ----

wkmat07_raw <- read_excel(here("data/metadata/spawning_season_wkmat07.xlsx"))

spawn_wkmat <- wkmat07_raw |>
  rowwise() |>
  mutate(spawn_months = list(as.integer(which(c_across(Jan:Dec) == 1)))) |>
  ungroup() |>
  filter(lengths(spawn_months) > 0) |>
  select(Scientific_name, Area, spawn_months) |>
  separate_rows(Area, sep = ",") |>
  mutate(
    Area     = trimws(Area),
    n_months = lengths(spawn_months)
  ) |>
  # keep narrower spawning window when a species x area pair appears more than once
  group_by(Scientific_name, Area) |>
  slice_min(n_months, with_ties = FALSE) |>
  ungroup() |>
  select(-n_months) |>
  rename(species = Scientific_name, ices_area = Area) |>
  mutate(
    spawn_months = map(spawn_months, as.integer),
    source       = "wkmat"
  ) |>
  filter(species %in% dataset_species)

# 5b. fishbase: spawning months mapped to ICES regions ----
# FishBase has obviously no ICES key! Spawningarea and SpawningGround are free text. Two tiers:
#  - specific seas (North Sea, Baltic, Bay of Biscay, ...) map to ONE ICES region
#    (source "fishbase").
#  - generic basin strings ("Northeast Atlantic", "Eastern Atlantic", "Europe", ...) can
#    not be pinned to a region. They are assigned to *all* dataset regions but kept only
#    as a *last resort* where every other source is empty (source "fishbase_basin").
# Non-ICES / western / tropical localities are dropped entirely.

month_cols <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# keyword (regex, lower-case) -> ICES region; evaluated top to bottom, first match retained
region_keywords <- tribble(
  ~pattern,                                                                                  ~region,
  "barents",                                                                                  1,
  "norwegian sea|lofoten|vester|vestfjord",                                                   2,
  "baltic|bornholm|bothnia|gotland|belt sea|skagerrak|kattegat|the sound|oresund|øresund", 3,
  "north sea|german bight|southern bight|dogger",                                             4,
  "iceland|faroe",                                                                            5,
  "rockall|hebrides|west.*scotland|scotland.*west|minch",                                     6,
  "celtic sea|irish sea|english channel|bristol channel|porcupine|west.*ireland",             7,
  "bay of biscay|golfe de gascogne|gascogne|cantabrian|armorican",                            8,
  "portug|iberian|galicia|gulf of cadiz|cadiz",                                               9,
  "azores",                                                                                  10
)

classify_region <- function(x) {
  x <- str_to_lower(x %||% "")
  for (i in seq_len(nrow(region_keywords))) {
    if (str_detect(x, region_keywords$pattern[i])) return(region_keywords$region[i])
  }
  NA_real_
}

# generic NE-Atlantic / Europe basin strings (kept) vs western / tropical / other (dropped)
generic_incl <- "northeast.*atlantic|north-east.*atlantic|eastern.*atlantic|east atlantic|ne atlantic|europe|uk|france"
generic_excl <- "west|tropical|south|america|gulf of maine|carib|brazil|argentin|africa|indian|pacific"

# collapse the 0/1 month flag columns into a list of active month indices
months_from <- function(df) {
  df |>
    rowwise() |>
    mutate(spawn_months = list(as.integer(which(c_across(all_of(month_cols)) == 1)))) |>
    ungroup() |>
    select(-all_of(month_cols))
}

fb_records <- tryCatch(
  rfishbase::spawning(species_list = dataset_species) |>
    as_tibble() |>
    mutate(
      loc    = coalesce(Spawningarea, SpawningGround),
      low    = str_to_lower(coalesce(loc, "")),
      region = map_dbl(loc, classify_region),
      across(all_of(month_cols), ~ as.integer(!is.na(.) & . != 0))
    ) |>
    select(species = Species, loc, low, region, all_of(month_cols)),
  error = function(e) {
    warning("FishBase query failed: ", conditionMessage(e))
    tibble(species = character(), loc = character(), low = character(),
           region = numeric())
  }
)

# specific tier: union of months per species x ICES region
spawn_fishbase <- fb_records |>
  filter(!is.na(region)) |>
  group_by(species, region) |>
  summarise(across(all_of(month_cols), max), .groups = "drop") |>
  months_from() |>
  filter(lengths(spawn_months) > 0) |>
  mutate(source = "fishbase") |>
  select(species, region, spawn_months, source)

# generic tier: union of months per species from basin-scale NE-Atlantic / Europe strings
spawn_fishbase_generic <- fb_records |>
  filter(is.na(region), str_detect(low, generic_incl), !str_detect(low, generic_excl)) |>
  group_by(species) |>
  summarise(across(all_of(month_cols), max), .groups = "drop") |>
  months_from() |>
  filter(lengths(spawn_months) > 0) |>
  mutate(source = "fishbase_basin") |>
  select(species, spawn_months, source)

# 6. merge sources ----

# merge a vector of (possibly already compound) source labels into a single string of
# unique, sorted atomic tokens, e.g. c("gofish", "gofish+wkmat") -> "gofish+wkmat"
combine_sources <- function(x) {
  paste(sort(unique(unlist(strsplit(x, "\\+")))), collapse = "+")
}

spawn_direct <- bind_rows(spawn_gofish, spawn_wkmat) |>
  filter(ices_area != "all") |>
  group_by(species, ices_area) |>
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    # contributing sources, e.g. "gofish", "wkmat", or "gofish+wkmat"
    source       = combine_sources(source),
    .groups      = "drop"
  )

# 7. propagate months to all ancestor ices areas ----

get_area_parents <- function(area) {
  parts <- strsplit(area, "\\.")[[1]]
  # returns the area itself and all its ancestors, from least to most specific
  sapply(seq_along(parts), function(i) paste(parts[1:i], collapse = "."))
}

# roll up child area data to every parent level (including the area itself)
spawn_all <- spawn_direct |>
  mutate(ancestors = map(ices_area, get_area_parents)) |>
  unnest(ancestors) |>
  group_by(species, ices_area = ancestors) |>
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    source       = combine_sources(source),
    .groups      = "drop"
  )

# region-level fallback using the leading integer of the ICES code
get_region <- function(area) as.numeric(sub("\\..*", "", area))

# the region tier feeds the resolver's region fallback at two priority levels:
#  - "good"  : gofish/wkmat rolled up to region, plus region-specific FishBase where the
#              survey sources are absent. Searched first, at the target region and its
#              neighbours (nearest first).
#  - "basin" : generic basin-scale FishBase, broadcast to every dataset region but kept
#              only where no "good" data exists. Used only as a final fallback, so a
#              basin-wide window never outranks real data from a neighbouring region.

# gofish/wkmat rolled up to region
spawn_region_survey <- spawn_all |>
  mutate(region = get_region(ices_area)) |>
  select(species, region, spawn_months, source)

# region-specific FishBase, only where the survey sources have nothing for that region
spawn_fishbase_keep <- spawn_fishbase |>
  anti_join(distinct(spawn_region_survey, species, region),
            by = c("species", "region"))

spawn_region_good <- bind_rows(spawn_region_survey, spawn_fishbase_keep) |>
  group_by(species, region) |>
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    source       = combine_sources(source),
    .groups      = "drop"
  )

# generic FishBase basin window: broadcast to every dataset region, kept only where no
# "good" region data exists for that species x region
dataset_regions <- sort(unique(get_region(species_area_combos$ices_area)))

spawn_region_basin <- spawn_fishbase_generic |>
  cross_join(tibble(region = dataset_regions)) |>
  anti_join(distinct(spawn_region_good, species, region),
            by = c("species", "region")) |>
  group_by(species, region) |>
  summarise(
    spawn_months = list(sort(unique(unlist(spawn_months)))),
    source       = combine_sources(source),
    .groups      = "drop"
  )

# 8. resolver ----

# search from most specific to least specific area
get_area_hierarchy <- function(area) {
  parts <- strsplit(area, "\\.")[[1]]
  sapply(length(parts):1, function(i) paste(parts[1:i], collapse = "."))
}

# resolution priority (most to least specific / preferred):
#   exact_area > parent_area > region > neighbour_region > basin
# `source` records the contributing data set(s); `match_type` records how the value was
# resolved; `source_area` records the ICES area or region it was taken from.
no_match <- list(spawn_months = NULL, source = NA_character_,
                 match_type = NA_character_, source_area = NA_character_)

resolve_spawn <- function(sp, area, spawn_direct, spawn_all,
                          spawn_region_good, spawn_region_basin) {

  # 1. walk the ICES area hierarchy, most specific first (survey sources only)
  for (a in get_area_hierarchy(area)) {

    hit_all <- filter(spawn_all, species == sp, ices_area == a)
    if (nrow(hit_all) == 0) next

    hit_direct <- filter(spawn_direct, species == sp, ices_area == a)
    has_direct <- nrow(hit_direct) > 0

    # "_subareas" = value comes from finer sub-areas rolled up to this area
    match_type <- paste0(
      if (a == area) "exact_area" else "parent_area",
      if (has_direct) "" else "_subareas"
    )

    return(list(
      spawn_months = hit_all$spawn_months[[1]],
      source       = if (has_direct) hit_direct$source[1] else hit_all$source[1],
      match_type   = match_type,
      source_area  = a
    ))
  }

  # 2. region fallback: target region and its neighbours (nearest first), using the
  #    "good" region data (rolled-up survey data + region-specific FishBase)
  target_region  <- get_region(area)
  regions_to_try <- c(target_region,
                      target_region - 1, target_region + 1,
                      target_region - 2, target_region + 2)
  regions_to_try <- regions_to_try[regions_to_try > 0]

  good_hit <- spawn_region_good |>
    filter(species == sp, region %in% regions_to_try) |>
    mutate(dist = abs(region - target_region)) |>
    slice_min(dist, n = 1, with_ties = FALSE)

  if (nrow(good_hit) == 1)
    return(list(
      spawn_months = good_hit$spawn_months[[1]],
      source       = good_hit$source[1],
      match_type   = if (good_hit$dist == 0) "region" else "neighbour_region",
      source_area  = as.character(good_hit$region[1])
    ))

  # 3. final fallback: basin-scale FishBase window for the target region
  basin_hit <- filter(spawn_region_basin, species == sp, region == target_region)

  if (nrow(basin_hit) == 1)
    return(list(
      spawn_months = basin_hit$spawn_months[[1]],
      source       = basin_hit$source[1],
      match_type   = "basin",
      source_area  = as.character(basin_hit$region[1])
    ))

  no_match
}

# 9. build the lookup ----

spawn_lookup <- species_area_combos |>
  mutate(
    resolved = map2(
      species, ices_area,
      resolve_spawn,
      spawn_direct       = spawn_direct,
      spawn_all          = spawn_all,
      spawn_region_good  = spawn_region_good,
      spawn_region_basin = spawn_region_basin
    ),
    spawn_months = map(resolved, "spawn_months"),
    source       = map_chr(resolved, ~ .x$source      %||% NA_character_),
    match_type   = map_chr(resolved, ~ .x$match_type  %||% NA_character_),
    source_area  = map_chr(resolved, ~ .x$source_area %||% NA_character_)
  ) |>
  select(-resolved) |>
  # keep only combinations for which spawning months were resolved
  filter(lengths(spawn_months) > 0)

# 10. manual additions ----
# species missing from all automated sources, filled in by hand. Each entry is applied to
# every ICES area where the species occurs in the survey (species_area_combos).

manual_spawn <- tibble(
  species      = "Chelidonichthys lucerna",
  spawn_months = list(c(1L, 2L)),
  source       = "fishbase",
  match_type   = "basin",
  source_area  = "Mediterranean"
) |>
  inner_join(distinct(species_area_combos, species, ices_area), by = "species") |>
  # do not overwrite anything already resolved for that species x area
  anti_join(distinct(spawn_lookup, species, ices_area),
            by = c("species", "ices_area"))

spawn_lookup <- bind_rows(spawn_lookup, manual_spawn)

# 11. extend spawning months by 3-month pre-spawn window ----
# for each spawning month, include the three calendar months immediately before it
# (wraps across the year boundary)

add_pre_spawn <- function(months) {
  if (is.null(months) || length(months) == 0) return(months)
  pre <- unlist(lapply(1:3, function(k) ((months - 1 - k) %% 12) + 1))
  sort(unique(c(months, pre)))
}

spawn_lookup <- spawn_lookup |>
  mutate(spawn_months_ext = map(spawn_months, add_pre_spawn))

# 12. summary ----

# how many species x area combinations were resolved at each level of the hierarchy
match_type_levels <- c("exact_area", "exact_area_subareas",
                       "parent_area", "parent_area_subareas",
                       "region", "neighbour_region", "basin")

source_summary <- spawn_lookup |>
  mutate(match_type = factor(match_type, levels = match_type_levels)) |>
  count(match_type, source, name = "n_entries") |>
  arrange(match_type, desc(n_entries))

# 13. species with no spawning information ----
# a species "has information" if at least one of its species x area combos
# resolved to a non-empty set of spawning months (source not NA)

species_with_info <- spawn_lookup |>
  filter(!is.na(source), lengths(spawn_months) > 0) |>
  distinct(species) |>
  pull(species)

# species we care about (dataset_species) that never resolved any spawning months
species_no_info <- setdiff(dataset_species, species_with_info)

# 14. save ----

write_rds(spawn_lookup, here("data/metadata/spawning_lookup.rds"))
