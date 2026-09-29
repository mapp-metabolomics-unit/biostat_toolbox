#!/usr/bin/env Rscript
args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else file.path(getwd(), "app.R")
v2_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
shiny::runApp(file.path(v2_dir, "..", "app"), launch.browser = TRUE)
