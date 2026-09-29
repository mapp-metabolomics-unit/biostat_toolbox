#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else file.path(getwd(), "src", "mapp_stats.R")
repo_root_override <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
repo_root <- if (nzchar(repo_root_override)) {
  normalizePath(repo_root_override, mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
}
source(file.path(repo_root, "src", "v2", "utils.R"))
source(file.path(repo_root, "src", "v2", "runner.R"))
source_v2_modules(repo_root)

usage <- function() {
  cat(paste(
    "MAPP statistics pipeline V2",
    "",
    "Usage:",
    "  Rscript src/mapp_stats.R validate --dataset DATASET.yaml --recipe RECIPE.yaml",
    "  Rscript src/mapp_stats.R plan     --dataset DATASET.yaml --recipe RECIPE.yaml",
    "  Rscript src/mapp_stats.R run      --dataset DATASET.yaml --recipe RECIPE.yaml",
    sep = "\n"
  ), "\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] %in% c("-h", "--help", "help")) { usage(); quit(status = 0) }
command <- args[1]
flag <- function(name) {
  index <- match(name, args)
  if (is.na(index) || index == length(args)) stop("Missing required option: ", name, call. = FALSE)
  args[index + 1]
}
dataset_yaml <- flag("--dataset")
recipe_yaml <- flag("--recipe")
config <- resolve_v2_config(dataset_yaml, recipe_yaml, repo_root)

if (command == "validate") {
  checked <- validate_v2_run(config)
  cat(jsonlite::toJSON(checked$validation, auto_unbox = TRUE, pretty = TRUE), "\n")
} else if (command == "plan") {
  checked <- validate_v2_run(config)
  preview <- preprocess_dataset(checked$dataset, config$recipe)
  preprocessing_summary <- list(
    raw_features = ncol(checked$dataset$raw),
    blank_filter_removed = if (nrow(preview$diagnostics$blank)) sum(preview$diagnostics$blank$removed_as_blank) else 0,
    qc_rsd_removed = if (nrow(preview$diagnostics$qc)) sum(preview$diagnostics$qc$removed_by_qc_rsd) else 0,
    retained_features = ncol(preview$matrices$scaled),
    analysis_samples = nrow(preview$matrices$scaled)
  )
  cat(yaml::as.yaml(list(output = file.path(config$paths$output_root, checked$identity$run_hash), hashes = checked$identity[c("dataset_hash", "preprocess_hash", "run_hash")], validation = checked$validation, preprocessing_preview = preprocessing_summary)))
} else if (command == "run") {
  execute_v2_run(config)
} else {
  usage()
  stop("Unknown command: ", command, call. = FALSE)
}
