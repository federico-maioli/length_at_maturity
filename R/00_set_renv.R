# 00 Set up renv ----
# One-time setup of a reproducible, project-local package library. Run this once before
# the rest of the pipeline to create renv.lock. Collaborators cloning the repository can
# instead run renv::restore() to install the exact package versions recorded in renv.lock.

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# initialise a project-local library (first run only; skipped if renv.lock already exists)
if (!file.exists("renv.lock")) renv::init(bare = TRUE)

# CRAN packages used across the pipeline
cran_pkgs <- c(
  "tidyverse", "here", "sf", "janitor", "jsonlite", "rfishbase", "worrms",
  "taxize", "readxl", "pROC", "performance", "knitr", "kableExtra",
  "rnaturalearth", "patchwork", "ggrepel", "paletteer", "geomtextpath",
  "ggstats", "DiagrammeR", "DiagrammeRsvg", "rsvg"
)

renv::install(cran_pkgs)
renv::install("tokami/DATRASextra")

# record the exact versions in renv.lock
renv::snapshot()
