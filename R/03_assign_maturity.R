library(tidyverse)
library(here)
library(sf)
library(DATRASextra)
library(rfishbase)
library(readxl)
library(knitr)
library(kableExtra)

# 01 Load data ----

raw <- readRDS(here("data/raw/raw_datras.rds")) |> DATRASextra::correct_species()

species_info <- read_rds(here("data/metadata/species_info.rds")) |>
  mutate(aphia_id = as.numeric(aphia_id))

# 02 Clean ca ----

# select relevant columns and normalise types
ca <- raw[["CA"]] |>
  select(Survey, Quarter, Country, Ship, HaulNo, haul.id, Year, AreaType, AreaCode,
         Sex, MaturityScale, Maturity, Species, Valid_Aphia, LngtCm) |>
  mutate(
    Year = as.integer(as.character(Year)),
    MaturityScale = as.character(MaturityScale),
    MaturityScale = ifelse(MaturityScale == "", NA, MaturityScale),
    Maturity = as.character(Maturity),
    Maturity = ifelse(Maturity == "", NA, Maturity)
  ) |>
  droplevels()

# replace raw species labels with validated names from species_info
ca <- ca |>
  left_join(species_info |> select(aphia_id, species, class),
            by = c("Valid_Aphia" = "aphia_id")) |>
  mutate(Species = if_else(!is.na(species), species, Species)) |>
  select(-species)

# keep only fish and elasmobranchs with a valid species name
ca <- ca |> filter(!is.na(Species), class %in% c("Teleostei", "Elasmobranchii"))

# merge haul-level metadata
hh <- raw[["HH"]] |>
  select(Survey, Quarter, Year, Month, Country, Ship, lon, lat, StatRec, haul.id) |>
  mutate(Year = as.integer(as.character(Year))) |>
  droplevels()

m_data <- ca |>
  left_join(hh) |>
  janitor::clean_names() |>
  filter(!(is.na(species) & is.na(valid_aphia))) |>
  filter(!is.na(maturity))

m_data <- m_data |> mutate(maturity = as.character(maturity))

# 03 Filter to study years ----

m_data <- m_data |> filter(year >= 2000)

# 04 Taxonomic validation ----
# keep only valid binomial species names (drops codes and partial / genus-only names)

m_data <- m_data |> filter(grepl("^[A-Z][a-z]+ [a-z]+$", species))

# 05 Remove length outliers ----

m_data <- m_data |>
  left_join(species_info |> select(species, max_l)) |>
  # flag a missing or sub-1 cm length as invalid, and flag as too large only when an Lmax is
  # available; species without an Lmax are never flagged on the upper bound
  mutate(length_flag = is.na(lngt_cm) | lngt_cm < 1 |
           (!is.na(max_l) & lngt_cm > 1.2 * max_l))

# how many individuals the length filter removes, and how many lack an Lmax reference
cat(sprintf(
  "Length filter: %d of %d records removed (%.2f%%); %d have no Lmax (not flagged as too large)\n",
  sum(m_data$length_flag), nrow(m_data),
  100 * mean(m_data$length_flag),
  sum(is.na(m_data$max_l))
))

m_data <- m_data |> filter(length_flag == FALSE)

# 06 Fix maturity scales ----

# detect whether the I/M binomial scale is used within each country-survey-year
im_scale <- m_data |>
  group_by(country, survey, year) |>
  summarise(has_M = any(maturity == "M", na.rm = TRUE), .groups = "drop")

m_data <- m_data |>
  left_join(im_scale, by = c("country", "survey", "year"))

# SE reported its 2000-2005 BITS maturity in Roman numerals, but this is the same 6-stage
# scale it otherwise reports numerically. Recode these to the numeric codes before the coding
# scheme is classified, so they are handled by the Numeric rule (I->1, II->2, ..., VI->6).
m_data <- m_data |>
  mutate(maturity = if_else(
    country == "SE" & survey == "BITS" & year >= 2000 & year <= 2005 & grepl("^[IVX]+$", maturity),
    dplyr::recode(maturity, "I" = "1", "II" = "2", "III" = "3", "IV" = "4", "V" = "5", "VI" = "6"),
    maturity
  ))

# classify each record's maturity coding scheme
m_data <- m_data |>
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

# 07 Assign binary maturity ----

m_data <- m_data |>
  mutate(maturity_chr = as.character(maturity)) |>
  mutate(
    mature = case_when(

      # R1
      code_type == "R1" & maturity_chr == "R1_1" ~ 0,
      code_type == "R1" & maturity_chr %in% c("R1_2","R1_3","R1_4") ~ 1,

      # R2
      code_type == "R2" & maturity_chr == "R2_2" ~ 0,
      code_type == "R2" & maturity_chr %in% c("R2_4","R2_6","R2_8") ~ 1,

      # M6
      code_type == "M6" & maturity_chr == "61" ~ 0,
      code_type == "M6" & maturity_chr %in% c("62","63","64","65") ~ 1,
      code_type == "M6" & maturity_chr == "66" ~ NA_real_,

      # SMSF
      code_type == "SMSF" & maturity_chr %in% c("A","Ba") ~ 0,
      code_type == "SMSF" & maturity_chr %in% c("B","Bb","C","Ca","Cb","D","Da","Db","E") ~ 1,
      code_type == "SMSF" & maturity_chr == "F" ~ NA_real_,

      # Roman (Baltic national keys). Stage I is virgin (immature) and stages III and above are
      # developing, spawning, or spent (mature). Stage II is ambiguous -- it labels both young
      # developing fish and resting post-spawning adults -- so we drop it. Exceptions: (i) Russia
      # and Lithuania record almost no stage I, so for them stage II stands in for the immature
      # class (I and II immature); (ii) in the SE-SOUND survey stage II is kept as mature (stage I
      # immature, all later stages mature); (iii) for herring in the Polish and Swedish BITS,
      # stage II is not ambiguous and is a genuine immature stage, so I and II are both immature.
      # IX/X are rare terminal/abnormal codes and are dropped. (SE's 2000-2005 BITS Roman was
      # recoded to Numeric above, so it is not handled here.)
      code_type == "Roman" & survey == "BITS" & species == "Clupea harengus" &
        country %in% c("PL","SE") & maturity_chr %in% c("I","II") ~ 0,
      code_type == "Roman" & country %in% c("RU","LT") & maturity_chr %in% c("I","II") ~ 0,
      code_type == "Roman" & survey == "SE-SOUND" & maturity_chr == "II" ~ 1,
      code_type == "Roman" & !(country %in% c("RU","LT")) & maturity_chr == "I" ~ 0,
      code_type == "Roman" & !(country %in% c("RU","LT")) & maturity_chr == "II" ~ NA_real_,
      code_type == "Roman" & maturity_chr %in% c("III","IV","V","VI","VII","VIII") ~ 1,
      code_type == "Roman" & maturity_chr %in% c("IX","X") ~ NA_real_,

      # Numeric
      code_type == "Numeric" & maturity_chr == "1" ~ 0,
      code_type == "Numeric" & maturity_chr %in% c("2","3","4","5") ~ 1,
      code_type == "Numeric" & maturity_chr == "6" ~ NA_real_,

      # I/M binomial
      code_type == "IM" & maturity_chr == "I" ~ 0,
      code_type == "IM" & maturity_chr == "M" ~ 1,

      TRUE ~ NA_real_
    )
  ) |>
  filter(!is.na(mature))

# 08 Summarise maturity coding ----
# collapsed overview of the coding schemes: per survey, country and code type, the years and
# codes used, the number of years and observations, and which codes were classed immature and
# mature. Saved as a csv and as a latex table for the supplement.

# order codes within a scale rather than alphabetically
code_order <- c("1", "2", "3", "4", "5", "6",
                "61", "62", "63", "64", "65", "66",
                "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                "A", "Ba", "B", "Bb", "C", "Ca", "Cb", "D", "Da", "Db", "E", "F", "M",
                "R1_1", "R1_2", "R1_3", "R1_4", "R2_2", "R2_4", "R2_6", "R2_8")

# collapse a set of codes into a single ordered, comma-separated string
order_codes <- function(x) {
  u <- unique(x)
  paste(u[order(match(u, code_order))], collapse = ", ")
}

# list years, merging consecutive runs into ranges (e.g. 2016, 2019--2021)
collapse_years <- function(years) {
  y <- sort(unique(years))
  runs <- cumsum(c(TRUE, diff(y) != 1))
  parts <- tapply(y, runs, function(g) if (length(g) == 1) as.character(g) else paste0(min(g), "--", max(g)))
  paste(parts, collapse = ", ")
}

maturity_summary <- m_data |>
  group_by(survey, country, code_type) |>
  summarise(
    years = collapse_years(year),
    codes = order_codes(maturity_chr),
    n_years = n_distinct(year),
    n_obs = n(),
    immature_codes = order_codes(maturity_chr[mature == 0]),
    mature_codes = order_codes(maturity_chr[mature == 1]),
    .groups = "drop"
  ) |>
  arrange(survey, country, code_type)

# share of records assigned from the Roman scale, of BITS and of all surveys
pct_roman_bits <- 100 * mean(m_data$code_type[m_data$survey == "BITS"] == "Roman")
pct_roman_all <- 100 * mean(m_data$code_type == "Roman")
# cat(sprintf("Roman-coded records: %.1f%% of BITS, %.1f%% of all surveys\n",
#             pct_roman_bits, pct_roman_all))

write_csv(maturity_summary, here("data/metadata/maturity_scale_summary.csv"))

# latex table for the supplement (escape underscores in the RIVO codes)
maturity_tex <- maturity_summary |>
  mutate(across(c(codes, immature_codes, mature_codes),
                ~ str_replace_all(.x, "_", "\\\\_"))) |>
  rename(
    Survey = survey,
    Country = country,
    Scale = code_type,
    Years = years,
    Codes = codes,
    `N years` = n_years,
    `N obs` = n_obs,
    Immature = immature_codes,
    Mature = mature_codes
  )

# country-code legend for the caption
country_names <- c(
  BE = "Belgium", DE = "Germany", DK = "Denmark", EE = "Estonia", ES = "Spain",
  FR = "France", GB = "United Kingdom", `GB-ENG` = "England", `GB-NIR` = "Northern Ireland",
  `GB-SCT` = "Scotland", IE = "Ireland", LT = "Lithuania", LV = "Latvia",
  NL = "Netherlands", NO = "Norway", PL = "Poland", PT = "Portugal",
  RU = "Russia", SE = "Sweden"
)
codes_present <- sort(unique(as.character(maturity_summary$country)))
legend_names <- coalesce(unname(country_names[codes_present]), codes_present)
country_legend <- paste(paste0(codes_present, ", ", legend_names), collapse = "; ")
tab_caption <- paste0(
  "Maturity coding schemes in the DATRAS data, by survey, country and scale, showing the ",
  "codes classified as immature and mature. Country codes: ", country_legend, "."
)

# small font, wrapping widths on the text columns (years, codes, immature, mature) and reduced
# inter-column padding so all columns fit the page width
maturity_kable <- kable(
  maturity_tex, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
  linesep = "", caption = tab_caption, label = "maturity_scales"
) |>
  kable_styling(latex_options = c("repeat_header", "striped"), font_size = 8,
                stripe_color = "gray!10") |>
  column_spec(4, width = "2.3cm") |>
  column_spec(5, width = "2.4cm") |>
  column_spec(8, width = "1.2cm") |>
  column_spec(9, width = "2.4cm")

# tighten the column padding, then save
maturity_kable <- sub("\\begin{longtable}",
                      "\\setlength{\\tabcolsep}{3pt}\n\\begin{longtable}",
                      as.character(maturity_kable), fixed = TRUE)
writeLines(maturity_kable, here("data/metadata/maturity_scale_summary.tex"))

# 09 Add area information ----

ices_sf <- read_sf(here("data/metadata/ices_areas/ICES_Statrec_mapto_ICES_Areas.shp")) |>
  group_by(Area_27) |>
  summarise(geometry = st_union(geometry), .groups = "drop")

m_data_sf <- m_data |>
  filter(!is.na(lon), !is.na(lat)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# spatial join to assign an ICES area to each record
m_data <- st_join(m_data_sf, ices_sf |> select(Area_27), left = TRUE) |>
  mutate(ices_area = Area_27) |>
  select(-Area_27) |>
  st_drop_geometry()

# 10 Add stock information ----

stock_table <- read_csv(here("data/metadata/species_area.csv"))

m_data <- m_data |>
  left_join(stock_table |> select(species, ices_area, stock))

# 11 Filter to spawning period ----
# lookup built by 02_create_spawning_data.R (GO-FISH + WKMAT07, with ICES area
# hierarchy inference and 3-month pre-spawn extension)

spawn_lookup <- read_rds(here("data/metadata/spawning_lookup.rds"))

m_data <- m_data |>
  left_join(
    spawn_lookup |> select(species, ices_area, spawn_months, spawn_months_ext,
                           source, source_area),
    by = c("species", "ices_area")
  ) |>
  mutate(
    spawning_season = map2_lgl(month, spawn_months_ext, ~ !is.null(.y) && .x %in% .y)
  ) |>
  filter(
    # retain records with no spawning info rather than dropping them
    map_lgl(spawn_months, ~ is.null(.x) || length(.x) == 0) |
      spawning_season == TRUE
  )

# 12 Save ----

# create species-stock identifier
m_data <- m_data |>
  mutate(
    species_clean = str_replace_all(species, " ", "_"),
    species_stock = case_when(
      !is.na(stock) ~ paste0(species_clean, "_", stock),
      TRUE ~ species_clean
    )
  )

write_rds(m_data, here("data/intermediate/maturity_clean.rds"))

# final dataset after taxonomic validation, length-outlier removal and spawning-window filtering
cat(sprintf(
  "Final dataset: %d individuals, %d species, %d hauls\n",
  nrow(m_data), n_distinct(m_data$species), n_distinct(m_data$haul_id)
))
