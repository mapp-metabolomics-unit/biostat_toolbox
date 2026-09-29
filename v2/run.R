#!/usr/bin/env Rscript
args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else file.path(getwd(), "run.R")
v2_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(v2_dir, ".."), mustWork = TRUE)
Sys.setenv(MAPP_STATS_REPO_ROOT = repo_root)
source(file.path(repo_root, "src", "mapp_stats.R"))
