# R script to initialize renv and install dependencies
# Save this script as setup_renv.R in your project root directory.

# --- Helper function to install and load renv ---
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
library(renv)

# --- Initialize renv for the project ---
# This will create an renv.lock file and a project-local library.
# Run this from your project's root directory in an R console.
# If you have an existing .Rprofile that loads renv, it might do this automatically.
# renv::init()
# If init() was already run or you prefer manual setup:
# renv::activate()

# --- Set Bioconductor version for R 4.2.x ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install(version = "3.16", update = FALSE, ask = FALSE)

# --- Install CRAN packages ---
# Based on your sessionInfo and script comments
# renv will automatically find packages used in your scripts if you run renv::snapshot()
# However, to be explicit and ensure all are captured from sessionInfo:

cran_packages <- c(
  "yaml", "WikidataQueryServiceR", "wesanderson", "webchem", "vegan", "lattice", "permute",
  "lubridate", "forcats", "stringr", "purrr", "tidyr", "tibble", "tidyverse", "rockchalk",
  "rfPermute", "readr", "pls", "plotly", "orca", "janitor", "iheatmapr", "htmltools",
  "ggrepel", "ggh4x", "ggplot2", "DT", "dplyr", "digest", "crosstalk", "websocket",
  "minqa", "colorspace", "snakecase", "xml2", "codetools", "splines", "knitr", "itertools",
  "jsonlite", "ratelimitr", "nloptr", "cluster", "missForest", "data.tree", "httr",
  "assertthat", "Matrix", "fastmap", "lazyeval", "cli", "later", "gtable", "glue",
  "reshape2", "ggthemes", "doRNG", "Rcpp", "carData", "vctrs", "nlme", "iterators",
  "xfun", "fastcluster", "ps", "openxlsx", "lme4", "rvest", "timechange", "lifecycle",
  "rngtools", "MASS", "scales", "promises", "hms", "kutils", 
  "RColorBrewer", "gridExtra", "stringi", "foreach", "randomForest", "boot", "zip",
  "rlang", "pkgconfig", "matrixStats", "bitops", "evaluate", "htmlwidgets", "processx",
  "cowplot", "tidyselect", "plyr", "magrittr", "R6", "generics", "pillar", "foreign",
  "withr", "mgcv", "sp", "RCurl", "tzdb", "data.table", "xtable", "munsell", "viridisLite"
)

print("Installing CRAN packages...")
renv::install(cran_packages)

# --- Install Bioconductor packages ---
bioc_packages <- c(
  "structToolbox", "pmp", "XVector", "GenomicRanges", "impute", "ontologyIndex", "Biobase",
  "GenomeInfoDbData", "zlibbioc", "S4Vectors", "BiocGenerics", "GenomeInfoDb", "IRanges",
  "DelayedArray", "SummarizedExperiment", "pcaMethods", "ropls", "BiocFileCache",
  "MatrixGenerics" # Moved MatrixGenerics here
)

print("Installing Bioconductor packages...")
# BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)
# Using renv::install for Bioconductor packages is preferred when using renv
renv::install(paste0("bioc::", bioc_packages))


# --- Install GitHub packages ---
github_packages <- c(
  "mapp-metabolomics-unit/MAPPstructToolbox", # As per your script and sessionInfo
  "KarstensLab/microshades",                 # As per your script and sessionInfo
  "jcheng5/d3scatter",                     # From script comments
  "mikemc/speedyseq"
  )

print("Installing GitHub packages...")
renv::install(github_packages)

# --- Create a lockfile ---
# This captures the state of your project's library.
print("Taking renv snapshot...")
renv::snapshot()

print("renv setup script finished. Check renv.lock for captured dependencies.")
print("You can now share your project, including the renv.lock file, .Rprofile, and the renv folder (optional but can speed up restore for others).")

