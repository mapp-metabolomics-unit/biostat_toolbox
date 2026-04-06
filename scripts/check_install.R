#!/usr/bin/env Rscript

required_r_version <- "4.2"
required_packages <- c(
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
  "MAPPstructToolbox",
  "microshades",
  "optparse",
  "plotly",
  "plotrix",
  "pls",
  "pmp",
  "purrr",
  "readr",
  "rfPermute",
  "rlang",
  "rockchalk",
  "stringr",
  "structToolbox",
  "svglite",
  "tibble",
  "tidyr",
  "vegan",
  "webchem",
  "wesanderson",
  "WikidataQueryServiceR",
  "yaml"
)

get_script_path <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
  if (length(script_path)) {
    return(normalizePath(script_path[1]))
  }
  normalizePath(file.path(getwd(), "scripts", "check_install.R"), mustWork = FALSE)
}

repo_root <- normalizePath(file.path(dirname(get_script_path()), ".."))
issues <- character()
warnings <- character()

if (requireNamespace("renv", quietly = TRUE)) {
  project_library <- renv::paths$library(project = repo_root)
  if (dir.exists(project_library)) {
    .libPaths(c(project_library, .libPaths()))
  }
}

if (!startsWith(as.character(getRversion()), required_r_version)) {
  issues <- c(
    issues,
    sprintf("R %s is required; found %s.", required_r_version, getRversion())
  )
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  issues <- c(
    issues,
    sprintf(
      "Missing R packages: %s. Run `Rscript install.R` from the repository root.",
      paste(missing_packages, collapse = ", ")
    )
  )
}

if (!nzchar(Sys.which("pandoc"))) {
  issues <- c(
    issues,
    "Pandoc was not found on PATH. Install pandoc before generating HTML outputs."
  )
}

required_paths <- c(
  file.path(repo_root, "params", "params.yaml"),
  file.path(repo_root, "params", "params_user.yaml")
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  warnings <- c(
    warnings,
    sprintf(
      "Missing config files: %s. Copy the templates in `params/` before running the toolbox.",
      paste(basename(missing_paths), collapse = ", ")
    )
  )
}

cat("biostat_toolbox install check\n")
cat(sprintf("Repository: %s\n", repo_root))

if (length(warnings)) {
  cat("\nWarnings:\n")
  for (item in warnings) {
    cat(sprintf("- %s\n", item))
  }
}

if (length(issues)) {
  cat("\nFailed checks:\n")
  for (item in issues) {
    cat(sprintf("- %s\n", item))
  }
  quit(status = 1)
}

cat("\nAll required checks passed.\n")
