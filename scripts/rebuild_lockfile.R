#!/usr/bin/env Rscript

snapshot_repo <- Sys.getenv(
  "BIOSTAT_TOOLBOX_CRAN_REPO",
  unset = "https://packagemanager.posit.co/cran/2025-04-10"
)
bioc_version <- "3.16"

cran_packages <- c(
  "crosstalk",
  "digest",
  "dplyr",
  "DT",
  "forcats",
  "ggh4x",
  "ggplot2",
  "ggrepel",
  "htmltools",
  "htmlwidgets",
  "iheatmapr",
  "janitor",
  "jsonlite",
  "optparse",
  "plotly",
  "plotrix",
  "pls",
  "purrr",
  "readr",
  "rfPermute",
  "rlang",
  "rockchalk",
  "stringr",
  "svglite",
  "tibble",
  "tidyr",
  "vegan",
  "emmeans",
  "mapdata",
  "maps",
  "viridisLite",
  "webchem",
  "wesanderson",
  "WikidataQueryServiceR",
  "yaml"
)

bioc_packages <- c("pmp")

github_packages <- c(
  "KarstensLab/microshades",
  "mapp-metabolomics-unit/MAPPstructToolbox@e78fb8f645d75d5d060fbdf0282b1e6df02d46d3"
)

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
if (length(script_path)) {
  script_path <- normalizePath(script_path[1])
} else {
  script_path <- normalizePath(file.path(getwd(), "scripts", "rebuild_lockfile.R"), mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."))

args <- commandArgs(trailingOnly = TRUE)
clean <- "--clean" %in% args

project_library <- file.path(repo_root, "renv", "library")
lockfile_path <- file.path(repo_root, "renv.lock")

Sys.setenv(RENV_CONFIG_PAK_ENABLED = "FALSE")
Sys.setenv(R_LIBS_USER = "")
options(
  renv.config.pak.enabled = FALSE,
  repos = c(CRAN = snapshot_repo)
)

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = getOption("repos"))
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = getOption("repos"))
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = getOption("repos"))
}

if (clean) {
  unlink(project_library, recursive = TRUE, force = TRUE)
  unlink(lockfile_path, force = TRUE)
}

renv::activate(project = repo_root)
renv::consent(provided = TRUE)
renv::settings$use.cache(FALSE, project = repo_root)
renv::settings$bioconductor.version(bioc_version, project = repo_root)

project_lib <- renv::paths$library(project = repo_root)
dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_lib, .Library))

install.packages(
  cran_packages,
  lib = project_lib,
  repos = getOption("repos"),
  dependencies = c("Depends", "Imports", "LinkingTo")
)

BiocManager::install(
  bioc_packages,
  ask = FALSE,
  update = FALSE
)

for (pkg in github_packages) {
  remotes::install_github(
    pkg,
    upgrade = "never",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

renv::snapshot(project = repo_root, prompt = FALSE, force = TRUE)

cat("Lockfile rebuild complete.\n")
cat(sprintf("Repository: %s\n", repo_root))
cat(sprintf("CRAN snapshot: %s\n", snapshot_repo))
cat(sprintf("Project library: %s\n", project_lib))
