#!/usr/bin/env Rscript

# Rebuild the locked MAPPstructToolbox revision after correcting its one stale
# reference to the legacy structToolbox namespace. Keeping both packages
# installed causes their duplicate PLSDA S4 classes to collide at runtime.

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1]) || !nzchar(x[1])) y else x
}

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
script_path <- script_path[1] %||% file.path(getwd(), "scripts", "install_patched_mappstructtoolbox.R")
repo_root <- normalizePath(file.path(dirname(normalizePath(script_path)), ".."))

lockfile_path <- file.path(repo_root, "renv.lock")
lockfile <- renv::lockfile_read(lockfile_path)
record <- lockfile$Packages$MAPPstructToolbox

if (is.null(record) || is.null(record$RemoteSha)) {
  stop("MAPPstructToolbox GitHub revision is missing from renv.lock.")
}

sha <- record$RemoteSha
archive_name <- sprintf("MAPPstructToolbox_%s.tar.gz", sha)
archive_path <- file.path(
  path.expand("~/.cache/R/renv/source/github/MAPPstructToolbox"),
  archive_name
)

temp_dir <- tempfile("patched-mappstructtoolbox-")
dir.create(temp_dir, recursive = TRUE)
on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

if (!file.exists(archive_path)) {
  archive_path <- file.path(temp_dir, archive_name)
  archive_url <- sprintf(
    "https://github.com/mapp-metabolomics-unit/MAPPstructToolbox/archive/%s.tar.gz",
    sha
  )
  download.file(archive_url, archive_path, mode = "wb", quiet = FALSE)
}

utils::untar(archive_path, exdir = temp_dir)
source_candidates <- list.files(
  temp_dir,
  pattern = "^oplsr_class[.]R$",
  recursive = TRUE,
  full.names = TRUE
)
source_candidates <- source_candidates[basename(dirname(source_candidates)) == "R"]

if (length(source_candidates) != 1L) {
  stop(sprintf("Expected one oplsr_class.R in the source archive; found %d.", length(source_candidates)))
}

oplsr_path <- source_candidates[[1]]
package_dir <- dirname(dirname(oplsr_path))
source_lines <- readLines(oplsr_path, warn = FALSE)
patched_lines <- gsub(
  "structToolbox:::ents$factor_name",
  "ents$factor_name",
  source_lines,
  fixed = TRUE
)

if (identical(source_lines, patched_lines)) {
  stop("The expected stale structToolbox reference was not found; refusing to install an unverified source tree.")
}
writeLines(patched_lines, oplsr_path, useBytes = TRUE)

library_path <- renv::paths$library(project = repo_root)
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
install_libraries <- unique(c(library_path, .libPaths()))
install_env <- sprintf("R_LIBS=%s", paste(install_libraries, collapse = .Platform$path.sep))
install_args <- c(
  "CMD",
  "INSTALL",
  paste0("--library=", shQuote(library_path)),
  shQuote(package_dir)
)

status <- system2(file.path(R.home("bin"), "R"), install_args, env = install_env)
if (!identical(status, 0L)) {
  stop(sprintf("Patched MAPPstructToolbox installation failed with status %s.", status))
}

.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(MAPPstructToolbox))
class_package <- getClassDef("PLSDA")@package
if (!identical(class_package, "MAPPstructToolbox")) {
  stop(sprintf("PLSDA resolved to %s instead of MAPPstructToolbox.", class_package))
}

message(sprintf("Installed patched MAPPstructToolbox %s; PLSDA resolves correctly.", record$Version))
