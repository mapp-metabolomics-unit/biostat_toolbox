#!/usr/bin/env Rscript
# Rebuilds renv.lock from scratch. Run as maintainer when adding/updating packages.
# Usage: Rscript scripts/rebuild_lockfile.R [--clean]
#
# To add a package: edit packages.yaml only, then re-run this script.

get_script_path <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
  if (length(script_path)) return(normalizePath(script_path[1]))
  normalizePath(file.path(getwd(), "scripts", "rebuild_lockfile.R"), mustWork = FALSE)
}

repo_root    <- normalizePath(file.path(dirname(get_script_path()), ".."))
packages_yml <- file.path(repo_root, "packages.yaml")
lockfile_path <- file.path(repo_root, "renv.lock")

if (!file.exists(packages_yml)) {
  stop("packages.yaml not found at: ", packages_yml, call. = FALSE)
}

# Bootstrap minimal deps before renv is active
for (pkg in c("renv", "remotes", "BiocManager", "yaml")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

cfg <- yaml::read_yaml(packages_yml)

snapshot_repo <- Sys.getenv("BIOSTAT_TOOLBOX_CRAN_REPO", unset = cfg$cran_snapshot)
bioc_version  <- cfg$bioc_version
cran_packages <- cfg$cran
bioc_packages <- cfg$bioc
github_packages <- cfg$github

Sys.setenv(RENV_CONFIG_PAK_ENABLED = "FALSE")
Sys.setenv(R_LIBS_USER = "")
options(
  renv.config.pak.enabled = FALSE,
  repos = c(CRAN = snapshot_repo)
)

args  <- commandArgs(trailingOnly = TRUE)
clean <- "--clean" %in% args

project_library <- file.path(repo_root, "renv", "library")

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
  lib          = project_lib,
  repos        = getOption("repos"),
  dependencies = c("Depends", "Imports", "LinkingTo")
)

BiocManager::install(bioc_packages, ask = FALSE, update = FALSE)

for (pkg in github_packages) {
  remotes::install_github(
    pkg,
    upgrade      = "never",
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

renv::snapshot(project = repo_root, prompt = FALSE, force = TRUE)

cat("Lockfile rebuild complete.\n")
cat(sprintf("Repository:   %s\n", repo_root))
cat(sprintf("CRAN snapshot: %s\n", snapshot_repo))
cat(sprintf("Project lib:  %s\n", project_lib))
