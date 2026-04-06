#!/usr/bin/env Rscript
# Quick sanity check: R version, renv library, pandoc, and param files.
# Package installation is handled entirely by renv — run `Rscript install.R` if needed.

required_r_version <- "4.2"

get_script_path <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
  if (length(script_path)) return(normalizePath(script_path[1]))
  normalizePath(file.path(getwd(), "scripts", "check_install.R"), mustWork = FALSE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."))
issues   <- character()
warnings <- character()

# Prepend renv library
if (requireNamespace("renv", quietly = TRUE)) {
  lib <- renv::paths$library(project = repo_root)
  if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
}

# R version
if (!startsWith(as.character(getRversion()), required_r_version)) {
  issues <- c(issues, sprintf("R %s is required; found %s.", required_r_version, getRversion()))
}

# renv library exists
lockfile <- file.path(repo_root, "renv.lock")
if (!file.exists(lockfile)) {
  issues <- c(issues, "renv.lock not found. Clone the repository cleanly and re-run.")
} else if (!dir.exists(renv::paths$library(project = repo_root))) {
  issues <- c(issues, "renv library not found. Run `Rscript install.R` to provision packages.")
}

# Pandoc
if (!nzchar(Sys.which("pandoc"))) {
  issues <- c(issues, "Pandoc not found on PATH. Install pandoc before generating HTML outputs.")
}

# Param files
required_paths <- c(
  file.path(repo_root, "params", "params.yaml"),
  file.path(repo_root, "params", "params_user.yaml")
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  warnings <- c(warnings, sprintf(
    "Missing config files: %s. Copy the templates in params/ before running.",
    paste(basename(missing_paths), collapse = ", ")
  ))
}

cat("biostat_toolbox install check\n")
cat(sprintf("Repository: %s\n", repo_root))

if (length(warnings)) {
  cat("\nWarnings:\n")
  for (w in warnings) cat(sprintf("  - %s\n", w))
}

if (length(issues)) {
  cat("\nFailed checks:\n")
  for (i in issues) cat(sprintf("  - %s\n", i))
  quit(status = 1)
}

cat("\nAll required checks passed.\n")
