# libraries ---------------------------------------------------------------
library(sf)
library(purrr)
library(tidyverse)
library(here)
library(DATRASextra)
library(rfishbase)
library(readxl)

# read in
data <- readRDS(here('data/raw/raw_datras.rds')) %>% DATRASextra::correct_species()

# clean CA ----------------------------------------------------------------

# get CA and HH and select only columns we care about
ca <- data[['CA']]

# get rid of empty values etc.
ca <- ca %>% select(Survey,Quarter,Country,Ship,HaulNo,haul.id,Year,AreaType,AreaCode,Sex,MaturityScale,Maturity,Species,Valid_Aphia,LngtCm)  %>%  mutate(
  Year = as.integer(as.character(Year)),
  MaturityScale = as.character(MaturityScale),
  MaturityScale = ifelse(MaturityScale == "", NA, MaturityScale),
  Maturity = as.character(Maturity),
  Maturity = ifelse(Maturity == "", NA, Maturity)
) %>% droplevels()

# fix species names
species_info <- read_rds(here('data/metadata/species_info.rds')) |> mutate(aphia_id=as.numeric(aphia_id))

ca <- ca %>%
  left_join(
    species_info %>%
      select(aphia_id, species, class),
    by = c("Valid_Aphia" = "aphia_id")
  ) %>%
  mutate(
    Species = if_else(
      !is.na(species),
      species,
      Species
    )
  ) %>%
  select(-species)

# remove non fish and NA
ca <- ca %>% filter(!is.na(Species), class %in% c('Teleostei','Elasmobranchii'))

# merge haul level data
hh <- data[['HH']] %>% select(Survey,Quarter,Year, Month,Country,Ship,lon,lat, StatRec,haul.id) %>% mutate(Year = as.integer(as.character(Year))) %>% droplevels()

# left_join and tidy up names
m_data <- ca %>% left_join(hh) %>% janitor::clean_names() %>% filter(
  !(is.na(species) & is.na(valid_aphia))
) %>% filter(!is.na(maturity))

# clean maturity codes and assign binary maturity
m_data <- m_data %>%
  # ensure character
  mutate(
    maturity = as.character(maturity)
  )

# fix maturity scales -----------------------------------------------------

# detect if I/M binomial scale exists within survey-country-year
im_scale <- m_data %>%
  group_by(country, survey, year) %>%
  summarise(
    has_M = any(maturity == "M", na.rm = TRUE),
    .groups = "drop"
  )

# join this information
m_data <- m_data %>%
  left_join(im_scale, by = c("country","survey","year"))

# detect code type
m_data <- m_data %>%
  mutate(
    code_type = case_when(

      grepl("^R1_", maturity) ~ "R1",
      grepl("^R2_", maturity) ~ "R2",

      grepl("^6[0-9]$", maturity) ~ "M6",

      maturity %in% c("A","B","Ba","Bb","C","Ca","Cb","D","Da","Db","E","F") ~ "SMSF",

      maturity %in% c("1","2","3","4","5","6") ~ "Numeric",

      maturity %in% c("I","M") & has_M ~ "IM",

      grepl("^[IVX]+$", maturity) & !has_M ~ "Roman",

      TRUE ~ "Unknown"
    )
  )

m_data = m_data |> filter(year>=2000)

# assign binary maturity
m_data <- m_data %>%
  mutate(
    maturity_chr = as.character(maturity)
  ) %>%
  group_by(country, code_type) %>%
  mutate(

    # detect structure
    has_high_roman = any(maturity_chr %in% c("VII","VIII","IX","X")),
    has_VI_only    = any(maturity_chr == "VI") & !has_high_roman,

    has_numeric6   = any(maturity_chr == "6")

  ) %>%
  ungroup() %>%
  mutate(
    mature = case_when(

      # ---------------- R1 ----------------
      code_type == "R1" & maturity_chr == "R1_1" ~ 0,
      code_type == "R1" & maturity_chr %in% c("R1_2","R1_3","R1_4") ~ 1,

      # ---------------- R2 ----------------
      code_type == "R2" & maturity_chr == "R2_2" ~ 0,
      code_type == "R2" & maturity_chr %in% c("R2_4","R2_6","R2_8") ~ 1,

      # ---------------- M6 ----------------
      code_type == "M6" & maturity_chr == "61" ~ 0,
      code_type == "M6" & maturity_chr %in% c("62","63","64","65") ~ 1,
      code_type == "M6" & maturity_chr == "66" ~ NA_real_,

      # ---------------- SMSF ----------------
      code_type == "SMSF" & maturity_chr %in% c("A","Ba") ~ 0,
      code_type == "SMSF" & maturity_chr %in% c("B","Bb","C","Ca","Cb","D","Da","Db","E") ~ 1,
      code_type == "SMSF" & maturity_chr == "F" ~ NA_real_,

      # ---------------- Roman (FULL scale) ----------------
      code_type == "Roman" & has_high_roman & maturity_chr %in% c("I") ~ 0,
      code_type == "Roman" & has_high_roman & maturity_chr %in% c("II","III","IV","V","VI","VII","VIII") ~ 1,
      code_type == "Roman" & has_high_roman & maturity_chr %in% c("IX","X") ~ NA_real_,

      # ---------------- Roman (TRUNCATED scale) ----------------
      code_type == "Roman" & has_VI_only & maturity_chr == "I" ~ 0,
      code_type == "Roman" & has_VI_only & maturity_chr %in% c("II","III","IV","V") ~ 1,
      code_type == "Roman" & has_VI_only & maturity_chr == "VI" ~ NA_real_,

      # ---------------- Numeric ----------------
      code_type == "Numeric" & maturity_chr == "1" ~ 0,
      code_type == "Numeric" & maturity_chr %in% c("2","3","4","5") ~ 1,
      code_type == "Numeric" & maturity_chr == "6" ~ NA_real_,

      # ---------------- I/M ----------------
      code_type == "IM" & maturity_chr == "I" ~ 0,
      code_type == "IM" & maturity_chr == "M" ~ 1,

      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(mature))

# get all the contrasts
# maturity_overview <- m_data %>%
#   distinct(year, country, code_type, maturity) %>%
#   group_by(year, country, code_type) %>%
#   summarise(
#     maturity_values = paste(sort(unique(maturity)), collapse = ", "),
#     .groups = "drop"
#   )
# maturity_overview <- m_data %>%
#   distinct(year, country, code_type, maturity) %>%
#   group_by(country, code_type) %>%
#   summarise(
#     start_year = min(year),
#     end_year   = max(year),
#     maturity_values = paste(sort(unique(maturity)), collapse = ", "),
#     .groups = "drop"
#   )


# add area information ----------------------------------------------------

ices_sf <- read_sf(here('data/metadata/ices_areas/StatRec_map_Areas_Full_20170124.shp'))

ices_sf <- ices_sf %>%
  group_by(Area_27) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# ggplot(sf) +
#   geom_sf(aes(fill = Area_27)) +
#   theme_minimal()

m_data_sf <- m_data %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# spatial join
m_joined <- st_join(
  m_data_sf,
  ices_sf %>% select(Area_27),
  left = TRUE
)

m_joined <- m_joined %>%
  mutate(
    ices_area = Area_27
  ) %>%
  select( -Area_27)

m_data <- m_joined %>%
  st_drop_geometry()

head(m_data)

#  add stock information --------------------------------------------------

stock_table <- read_csv(here('data/metadata/species_area.csv'))

m_data <- m_data %>% left_join(stock_table %>% select(species,ices_area,stock))

# filter to spawning / pre-spawning period --------------------------------
# lookup built by create_spawning_data.R (GO-FISH + WKMAT07, with ICES area
# hierarchy inference and 3-month pre-spawn extension)

spawn_lookup <- read_rds(here('data/metadata/spawning_lookup.rds'))

m_data <- m_data %>%
  left_join(
    spawn_lookup %>% select(species, ices_area, spawn_months, spawn_months_ext,
                            source, source_area),
    by = c("species", "ices_area")
  ) %>%
  mutate(
    spawning_season = map2_lgl(month, spawn_months_ext,
                               ~ !is.null(.y) && .x %in% .y)
  ) %>%
  filter(
    # keep rows with no spawning info (cannot filter them)
    map_lgl(spawn_months, ~ is.null(.x) || length(.x) == 0) |
      spawning_season == TRUE
  )


# check for length outliers ----------------------------------------------

m_data <- m_data %>%
  left_join(
    species_info %>%
      select(species, max_l, l_inf_max)
  )

m_data <- m_data %>%
  mutate(
    length_flag =
      lngt_cm > 1.4 * max_l |
      lngt_cm < 1
  )

m_data <- m_data |>
  filter(grepl("^[A-Z][a-z]+ [a-z]+$", species))

# length_check <- m_data %>%
#   group_by(species) %>%
#   summarise(
#     n = n(),
#     min_length = min(lngt_cm, na.rm = TRUE),
#     max_length = max(lngt_cm, na.rm = TRUE),
#     max_allowed = max(maxL, na.rm = TRUE),
#     prop_flagged = mean(length_flag, na.rm = TRUE),
#     .groups = "drop"
#   )

m_data <- m_data |> filter(length_flag == FALSE)

# create species stock identifier
m_data <- m_data |>
  mutate(
    species_clean = str_replace_all(species, " ", "_"),
    species_stock = case_when(
      !is.na(stock) ~ paste0(species_clean, "_", stock),
      TRUE ~ species_clean
    )
  )


# set min number of observations ------------------------------------------


# species_counts <- m_data %>%
#   count(species, name = "n_obs") |> arrange(n_obs)
#
# species_counts

m_data_filtered <- m_data %>%
  group_by(species_stock) %>%
  filter(n() >= 50,
         sum(mature == 0, na.rm = TRUE) >= 15,
         sum(mature == 1, na.rm = TRUE) >= 15) %>%
  ungroup()


write_rds(m_data_filtered, here('data/intermediate/maturity_clean.rds'))
#unique(m_data_filtered$species)
