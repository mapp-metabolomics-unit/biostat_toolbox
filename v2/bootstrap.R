#!/usr/bin/env Rscript

if (getRversion() < "4.6.0") {
  stop("V2 requires R >= 4.6. Install it side by side with `rig add 4.6`, then run `rig run -r 4.6 -f bootstrap.R`. This does not change the default R 4.2 used by legacy runs.")
}
project <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project, "bootstrap.R"))) {
  stop("Run bootstrap.R from the repository's v2 directory.")
}
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")
if (!file.exists(file.path(project, "renv", "activate.R"))) {
  renv::init(project = project, bare = TRUE, restart = FALSE)
}
# renv::init() can prepend an activation line when the repository already
# provides one. Keep exactly one line so activation cannot lock against itself.
writeLines('source("renv/activate.R")', file.path(project, ".Rprofile"))
renv::install(
  c("bioc::structToolbox", "yaml", "digest", "jsonlite", "ggplot2", "shiny", "callr"),
  project = project
)
renv::snapshot(project = project, prompt = FALSE)
message("V2 environment created at: ", project)
