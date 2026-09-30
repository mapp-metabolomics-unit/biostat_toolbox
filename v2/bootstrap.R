#!/usr/bin/env Rscript

if (getRversion() < "4.6.0") stop("The replacement pipeline requires R >= 4.6.")
project <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project, "renv.lock")) || !file.exists(file.path(project, "bootstrap.R"))) {
  stop("Run this script from the repository's v2 directory.")
}
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
lock_version <- renv::lockfile_read(file.path(project, "renv.lock"))$R$Version
if (as.character(getRversion()) != lock_version) {
  stop("This lockfile requires R ", lock_version, "; current R is ", as.character(getRversion()), ".")
}
renv::restore(project = project, prompt = FALSE)
required <- c("yaml", "digest", "jsonlite", "ggplot2", "shiny", "pls", "plotly")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("The pinned environment is missing: ", paste(missing, collapse = ", "))
message("Replacement pipeline environment ready: ", project)
