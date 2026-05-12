# libraries
library(DATRASextra) #remotes::install_github("tokami/DATRASextra")
library(tidyverse)
library(here)
library(sf)
library(jsonlite)

# select surveys
surveys <- list_surveys() %>% pull(survey) %>% 
  setdiff(c("BTS-GSA17",
              "Can-Mar",
              "IS-IDPS",
              "NS-IDPS",
              "IS-IDPS"))

# create directory
dir <- here('data/raw/surveys/')

# download
download_datras(surveys = surveys, dir = dir) # this takes quite a while!

# read in
data <- read_datras(file.path(dir,surveys))

# save them
saveRDS(data, here('data/raw/raw_datras.rds'))

# download ices data
url <- paste0(
  "https://gis.ices.dk/gis/rest/services/",
  "Mapping_layers/ICES_Statrec_mapto_ICES_Areas/",
  "MapServer/0/query?",
  "where=1%3D1&outFields=*&f=geojson"
)

ices <- st_read(url)

# save as shapefile
st_write(
  ices, here("data/metadata/ices_areas/ICES_Statrec_mapto_ICES_Areas.shp"),
  delete_layer = TRUE
)

# download fish spawning
#https://zenodo.org/records/11098993
dir.create(
  here::here("data", "metadata", "go_fish_spawning"),
  recursive = TRUE,
  showWarnings = FALSE
)

record <- fromJSON("https://zenodo.org/api/records/11098993")

# Extract file table
files <- record$files

# Download all files
for(i in seq_len(nrow(files))) {
  download.file(
    url = files$links$self[i],
    destfile = here("data", "metadata","go_fish_spawning", files$key[i]),
    mode = "wb"
  )
}
