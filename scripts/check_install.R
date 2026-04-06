#!/usr/bin/env Rscript
# Checks that all explicitly required packages (from packages.yaml) are installed.
# Run: Rscript scripts/check_install.R

get_script_path <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
  if (length(script_path)) return(normalizePath(script_path[1]))
  normalizePath(file.path(getwd(), "scripts", "check_install.R"), mustWork = FALSE)
}

repo_root    <- normalizePath(file.path(dirname(get_script_path()), ".."))
packages_yml <- file.path(repo_root, "packages.yaml")
issues   <- character()
warnings <- character()

# Prepend renv library so installed packages are visible
if (requireNamespace("renv", quietly = TRUE)) {
  project_library <- renv::paths$library(project = repo_root)
  if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
}

if (!file.exists(packages_yml)) {
  stop("packages.yaml not found. Ensure you are running from the repository root.", call. = FALSE)
}

cfg <- yaml::read_yaml(packages_yml)

# R version check
if (!startsWith(as.character(getRversion()), cfg$r_version)) {
  issues <- c(issues, sprintf(
    "R %s is required; found %s.", cfg$r_version, getRversion()
  ))
}

# Build list of required top-level package names
pkg_names <- function(github_specs) {
  # "user/repo@sha" or "user/repo" -> "repo"
  repos <- sub("@.*$", "", github_specs)       # strip @sha
  basename(repos)                               # strip user/
}

required_packages <- c(
  cfg$cran,
  cfg$bioc,
  pkg_names(cfg$github)
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  issues <- c(issues, sprintf(
    "Missing R packages: %s.\nRun `Rscript install.R` from the repository root.",
    paste(missing_packages, collapse = ", ")
  ))
}

# Pandoc check
if (!nzchar(Sys.which("pandoc"))) {
  issues <- c(issues, "Pandoc not found on PATH. Install pandoc before generating HTML outputs.")
}

# Config file check
required_paths <- c(
  file.path(repo_root, "params", "params.yaml"),
  file.path(repo_root, "params", "params_user.yaml")
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  warnings <- c(warnings, sprintf(
    "Missing config files: %s. Copy the templates in `params/` before running the toolbox.",
    paste(basename(missing_paths), collapse = ", ")
  ))
}

cat("biostat_toolbox install check\n")
cat(sprintf("Repository: %s\n", repo_root))

if (length(warnings)) {
  cat("\nWarnings:\n")
  for (item in warnings) cat(sprintf("  - %s\n", item))
}

if (length(issues)) {
  cat("\nFailed checks:\n")
  for (item in issues) cat(sprintf("  - %s\n", item))
  quit(status = 1)
}

cat("\nAll required checks passed.\n")
