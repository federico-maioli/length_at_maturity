library(DATRASextra) # remotes::install_github("tokami/DATRASextra")
library(tidyverse)
library(here)
library(sf)
library(jsonlite)

# 01 Create project folders ----

dirs <- c(
  here("data", "raw", "surveys"),
  here("data", "intermediate"),
  here("data", "final"),
  here("data", "metadata", "ices_areas"),
  here("data", "metadata", "go_fish_spawning"),
  here("outputs", "main"),
  here("outputs", "supp")
)

walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# 02 Select surveys ----

surveys <- list_surveys() |>
  pull(survey) |>
  setdiff(c("BTS-GSA17", "Can-Mar", "IS-IDPS", "NS-IDPS"))

# 03 Download DATRAS survey data ----

dir <- here("data", "raw", "surveys")

# this takes quite a while
download_datras(surveys = surveys, dir = dir)

data <- read_datras(file.path(dir, surveys))

write_rds(data, here("data", "raw", "raw_datras.rds"))

# 04 Download ICES area polygons ----

url <- paste0(
  "https://gis.ices.dk/gis/rest/services/",
  "Mapping_layers/ICES_Statrec_mapto_ICES_Areas/",
  "MapServer/0/query?",
  "where=1%3D1&outFields=*&f=geojson"
)

ices <- st_read(url)

st_write(
  ices,
  here("data", "metadata", "ices_areas", "ICES_Statrec_mapto_ICES_Areas.shp"),
  delete_layer = TRUE
)

# 05 Download spawning data (GO-FISH) ----
# https://zenodo.org/records/11098993

record <- fromJSON("https://zenodo.org/api/records/11098993")
files <- record$files

for (i in seq_len(nrow(files))) {
  download.file(
    url = files$links$self[i],
    destfile = here("data", "metadata", "go_fish_spawning", files$key[i]),
    mode = "wb"
  )
}
