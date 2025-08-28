# The functions found in src/R/ are from Engelke S, Gnecco N, Rottger F (2025)
# This setup was taken from Engelke S, Gnecco N, Rottger F (2025)
# https://github.com/nicolagnecco/extremeSCM 
# to load required functions

# Set CRAN mirror (pak handles this, but keeping for safety)
chooseCRANmirror(ind = 1)

# Install {pak} if not already installed
if (!"pak" %in% rownames(installed.packages())) {
  install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")
}

# Load pak
library(pak)

# Define required packages
cran_packs <- c(
  "pcalg", "dHSIC", "ranger",
  "tidyverse", "here",
  "graph", "MASS",
  "furrr", "future", "progressr",
  "mvtnorm",
  "glue", "ggpubr", "latex2exp", 
  "VGAM", "skewsamp", "rngtools",
  "ggExtra", "igraph", "geosphere", "egg", "gamlss.dist",
  "ggpubr"
)

# Install missing packages (CRAN + GitHub)
pak::pkg_install(c(cran_packs, 
                   "sebastian-engelke/graphicalExtremes",
                   "nicolagnecco/causalXtreme"))

# Setup completed
message("Setup complete! Now run `source('main/load_packages.R')` in your scripts.")
