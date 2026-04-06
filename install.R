#!/usr/bin/env Rscript

snapshot_repo <- Sys.getenv(
  "BIOSTAT_TOOLBOX_CRAN_REPO",
  unset = "https://packagemanager.posit.co/cran/2025-04-10"
)

Sys.setenv(RENV_CONFIG_PAK_ENABLED = "FALSE")
options(
  renv.config.pak.enabled = FALSE,
  repos = c(CRAN = snapshot_repo)
)

repo_root <- normalizePath(getwd())
lockfile_path <- file.path(repo_root, "renv.lock")

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = getOption("repos"))
}

if (!file.exists(lockfile_path)) {
  stop(
    paste(
      "renv.lock is missing.",
      "Maintainers should rebuild it with",
      "`Rscript scripts/rebuild_lockfile.R --clean` from the repository root."
    ),
    call. = FALSE
  )
}

renv::restore(project = repo_root, prompt = FALSE)
