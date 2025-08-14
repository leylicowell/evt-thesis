# The functions found in src/R/ are from Engelke S, Gnecco N, Rottger F (2025)
# https://github.com/nicolagnecco/extremeSCM
# This setup was taken from Engelke S, Gnecco N, Rottger F (2025)
# to load required functions

# Define required packages
cran_packs <- c(
  "pcalg", "dHSIC", "ranger",
  "tidyverse",
  "graph", "MASS",
  "furrr", "future", "progressr",
  "mvtnorm",
  "glue", "ggpubr", "latex2exp", 
  "VGAM", "skewsamp", "rngtools",
  "ggExtra", "igraph", "geosphere", "egg",
  "gamlss.dist", "ggpubr"
)

# Automatically load all installed packages
invisible(lapply(cran_packs, library, character.only = TRUE))
library(graphicalExtremes)  # Manually load GitHub package

# Load local functions from "R" directory
purrr::walk(list.files(here::here("src","R"), full.names = TRUE), source)

# Confirmation message
message("Packages and functions loaded successfully!")
