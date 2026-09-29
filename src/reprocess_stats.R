#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(digest)
  library(dplyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
if (length(script_path)) {
  script_path <- normalizePath(script_path[1])
} else {
  script_path <- normalizePath(file.path(getwd(), "src", "reprocess_stats.R"), mustWork = FALSE)
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

source(file.path(script_dir, "helpers.r"), local = TRUE)

option_list <- list(
  make_option(c("-s", "--stats-dir"), default = NULL, help = "Stats root containing hash subdirectories with params.yaml"),
  make_option(c("-o", "--output-root"), default = NULL, help = "New stats output root"),
  make_option(c("-u", "--params-user"), default = NULL, help = "params_user.yaml recursively merged over each archived params_user.yaml"),
  make_option(c("--defaults-yaml"), default = file.path(repo_root, "params", "params.yaml"), help = "YAML used for missing plot-default sections such as ordination and npc_summed_intensity"),
  make_option(c("--override-yaml"), default = NULL, help = "YAML file recursively merged into each archived params.yaml"),
  make_option(c("--include"), default = NULL, help = "Comma-separated original hashes to process"),
  make_option(c("--exclude"), default = NULL, help = "Comma-separated original hashes to skip"),
  make_option(c("--dry-run"), action = "store_true", default = FALSE, help = "Show planned runs without writing outputs"),
  make_option(c("--overwrite"), action = "store_true", default = FALSE, help = "Rerun even when the predicted output already exists"),
  make_option(c("--stop-on-error"), action = "store_true", default = FALSE, help = "Stop at the first failed run")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

normalize_option_name <- function(opt, underscore_name, hyphen_name) {
  if (is.null(opt[[underscore_name]]) && !is.null(opt[[hyphen_name]])) {
    opt[[underscore_name]] <- opt[[hyphen_name]]
  }
  opt
}

opt <- normalize_option_name(opt, "stats_dir", "stats-dir")
opt <- normalize_option_name(opt, "output_root", "output-root")
opt <- normalize_option_name(opt, "params_user", "params-user")
opt <- normalize_option_name(opt, "defaults_yaml", "defaults-yaml")
opt <- normalize_option_name(opt, "override_yaml", "override-yaml")
opt <- normalize_option_name(opt, "dry_run", "dry-run")
opt <- normalize_option_name(opt, "stop_on_error", "stop-on-error")

has_value <- function(value) {
  !is.null(value) && length(value) && !is.na(value[1]) && nzchar(trimws(as.character(value[1])))
}

resolve_path <- function(path_value, fallback_dir = getwd(), must_work = FALSE) {
  if (!has_value(path_value)) {
    return(NULL)
  }
  path_value <- trimws(as.character(path_value[1]))
  if (grepl("^/", path_value)) {
    return(normalizePath(path_value, mustWork = must_work))
  }
  if (file.exists(path_value)) {
    return(normalizePath(path_value, mustWork = must_work))
  }
  normalizePath(file.path(fallback_dir, path_value), mustWork = must_work)
}

project_library_path <- function(project_root) {
  r_version <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  candidate <- file.path(project_root, "renv", "library", paste0("R-", r_version), R.version$platform)
  if (dir.exists(candidate)) {
    return(normalizePath(candidate, mustWork = TRUE))
  }
  NULL
}

parse_hash_list <- function(value) {
  if (!has_value(value)) {
    return(character())
  }
  values <- unlist(strsplit(as.character(value), ",", fixed = TRUE))
  values <- trimws(values)
  values[nzchar(values)]
}

deep_merge <- function(base, override) {
  if (is.null(override) || (is.list(override) && !length(override))) {
    return(base)
  }
  if (is.null(base) || !is.list(base) || !is.list(override) || is.null(names(override))) {
    return(override)
  }
  for (name in names(override)) {
    base[[name]] <- deep_merge(base[[name]], override[[name]])
  }
  base
}

default_params_user <- function(params) {
  list(
    paths = list(
      docs = params$paths$docs %||% "",
      output = ""
    ),
    operating_system = list(
      system = params$operating_system$system %||% "unix",
      pandoc = params$operating_system$pandoc %||% "pandoc"
    )
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x) || is.na(x[1])) {
    return(y)
  }
  x
}

prepare_hash_params <- function(params, params_user) {
  params$paths$docs <- params_user$paths$docs
  params$paths$output <- params_user$paths$output
  params$operating_system$system <- params_user$operating_system$system
  params$operating_system$pandoc <- params_user$operating_system$pandoc
  params
}

merge_color_palette <- function(params, override) {
  override_keys <- override$colors$all$key
  override_values <- override$colors$all$value
  current_keys <- params$colors$all$key
  current_values <- params$colors$all$value

  if (is.null(override_keys) || is.null(override_values) || is.null(current_keys) || !length(current_keys)) {
    return(list(params = params, override = override))
  }
  if (length(override_keys) != length(override_values)) {
    stop("Override colors.all.key and colors.all.value must have the same length.")
  }
  if (is.null(current_values) || length(current_values) != length(current_keys)) {
    return(list(params = params, override = override))
  }

  palette_lookup <- stats::setNames(as.character(override_values), as.character(override_keys))
  matched_keys <- as.character(current_keys) %in% names(palette_lookup)
  current_values[matched_keys] <- unname(palette_lookup[as.character(current_keys)[matched_keys]])
  params$colors$all$value <- current_values

  override$colors$all$key <- NULL
  override$colors$all$value <- NULL
  if (!length(override$colors$all)) {
    override$colors$all <- NULL
  }
  if (!length(override$colors)) {
    override$colors <- NULL
  }

  list(params = params, override = override)
}

merge_params_override <- function(params, override) {
  color_merge <- merge_color_palette(params, override)
  params <- color_merge$params
  override <- color_merge$override
  if (!is.null(override$npc_summed_intensity)) {
    params$npc_summed_intensity <- NULL
  }
  deep_merge(params, override)
}

split_override_params <- function(override_params) {
  if (is.null(override_params)) {
    override_params <- list()
  }
  by_target <- override_params$by_target
  override_params$by_target <- NULL
  list(global = override_params, by_target = by_target %||% list())
}

apply_plot_defaults <- function(params, defaults) {
  if (is.null(defaults)) {
    return(params)
  }
  for (section_name in c("ordination", "npc_summed_intensity")) {
    if (!is.null(defaults[[section_name]])) {
      params[[section_name]] <- deep_merge(defaults[[section_name]], params[[section_name]])
    }
  }
  params
}

apply_target_override <- function(params, override_parts) {
  params <- merge_params_override(params, override_parts$global)
  target_name <- params$target$sample_metadata_header
  if (!is.null(target_name) && target_name %in% names(override_parts$by_target)) {
    params <- merge_params_override(params, override_parts$by_target[[target_name]])
  }
  params
}

write_manifest <- function(manifest, manifest_path) {
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  write.table(manifest, manifest_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

copy_yaml_if_available <- function(source_path, destination_path) {
  if (is.null(source_path) || !file.exists(source_path)) {
    return(invisible(FALSE))
  }
  dir.create(dirname(destination_path), recursive = TRUE, showWarnings = FALSE)
  file.copy(source_path, destination_path, overwrite = TRUE)
}

write_reprocess_inputs <- function(destination_dir, params, params_user, source_params, source_params_user, params_user_override_path, defaults_path, override_path) {
  inputs_dir <- file.path(destination_dir, "_reprocess_inputs")
  dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)

  yaml::write_yaml(params, file.path(inputs_dir, "merged_params.yaml"))
  yaml::write_yaml(params_user, file.path(inputs_dir, "merged_params_user.yaml"))
  copy_yaml_if_available(source_params, file.path(inputs_dir, "source_params.yaml"))
  copy_yaml_if_available(source_params_user, file.path(inputs_dir, "source_params_user.yaml"))
  copy_yaml_if_available(params_user_override_path, file.path(inputs_dir, "params_user_override.yaml"))
  copy_yaml_if_available(defaults_path, file.path(inputs_dir, "defaults.yaml"))
  copy_yaml_if_available(override_path, file.path(inputs_dir, "override.yaml"))

  invisible(inputs_dir)
}

if (!has_value(opt$stats_dir)) {
  stop("--stats-dir is required.")
}
if (!has_value(opt$output_root)) {
  stop("--output-root is required.")
}

stats_dir <- resolve_path(opt$stats_dir, must_work = TRUE)
output_root <- resolve_path(opt$output_root, must_work = FALSE)
params_user_override_path <- resolve_path(opt$params_user, must_work = TRUE)
defaults_path <- resolve_path(opt$defaults_yaml, must_work = FALSE)
override_path <- resolve_path(opt$override_yaml, must_work = TRUE)

params_user_override <- list()
if (!is.null(params_user_override_path)) {
  params_user_override <- yaml.load_file(params_user_override_path)
  if (is.null(params_user_override)) {
    params_user_override <- list()
  }
}

default_params <- list()
if (!is.null(defaults_path) && file.exists(defaults_path)) {
  default_params <- yaml.load_file(defaults_path)
  if (is.null(default_params)) {
    default_params <- list()
  }
}

override_params <- list()
if (!is.null(override_path)) {
  override_params <- yaml.load_file(override_path)
  if (is.null(override_params)) {
    override_params <- list()
  }
}
override_parts <- split_override_params(override_params)

include_hashes <- parse_hash_list(opt$include)
exclude_hashes <- parse_hash_list(opt$exclude)

candidate_dirs <- list.dirs(stats_dir, full.names = TRUE, recursive = FALSE)
candidate_dirs <- candidate_dirs[file.exists(file.path(candidate_dirs, "params.yaml"))]
if (length(include_hashes)) {
  candidate_dirs <- candidate_dirs[basename(candidate_dirs) %in% include_hashes]
}
if (length(exclude_hashes)) {
  candidate_dirs <- candidate_dirs[!basename(candidate_dirs) %in% exclude_hashes]
}
candidate_dirs <- sort(candidate_dirs)

if (!length(candidate_dirs)) {
  stop("No stats result folders with params.yaml matched the selection.")
}

work_root <- file.path(output_root, "_reprocess_work")
log_root <- file.path(output_root, "_reprocess_logs")
manifest_path <- file.path(output_root, "reprocess_manifest.tsv")
biostat_script <- file.path(script_dir, "biostat_toolbox.r")
renv_library <- project_library_path(repo_root)
child_env <- character()
if (!is.null(renv_library)) {
  child_libraries <- c(.libPaths(), renv_library)
  child_libraries <- child_libraries[nzchar(child_libraries)]
  child_env <- sprintf("R_LIBS=%s", paste(child_libraries, collapse = .Platform$path.sep))
  message(sprintf("Adding project R library as child-run fallback: %s", renv_library))
}

if (!isTRUE(opt$dry_run)) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  batch_inputs_dir <- file.path(output_root, "_reprocess_inputs")
  dir.create(batch_inputs_dir, recursive = TRUE, showWarnings = FALSE)
  copy_yaml_if_available(params_user_override_path, file.path(batch_inputs_dir, "params_user_override.yaml"))
  copy_yaml_if_available(defaults_path, file.path(batch_inputs_dir, "defaults.yaml"))
  copy_yaml_if_available(override_path, file.path(batch_inputs_dir, "override.yaml"))
}

manifest <- data.frame(
  original_hash = character(),
  new_hash = character(),
  description = character(),
  source_params = character(),
  output_dir = character(),
  status = character(),
  started_at = character(),
  ended_at = character(),
  log_path = character(),
  stringsAsFactors = FALSE
)

message(sprintf("Discovered %d stats run(s).", length(candidate_dirs)))

for (result_dir in candidate_dirs) {
  original_hash <- basename(result_dir)
  source_params <- file.path(result_dir, "params.yaml")
  source_params_user <- file.path(result_dir, "params_user.yaml")
  started_at <- as.character(Sys.time())
  log_path <- file.path(log_root, paste0(original_hash, ".log"))

  params <- yaml.load_file(source_params)
  params_user <- if (file.exists(source_params_user)) {
    yaml.load_file(source_params_user)
  } else {
    default_params_user(params)
  }
  params_user <- deep_merge(params_user, params_user_override)
  params <- apply_plot_defaults(params, default_params)
  params <- apply_target_override(params, override_parts)
  params_user$paths$output <- output_root

  if (is.null(params_user$paths$docs) || !nzchar(params_user$paths$docs)) {
    stop(sprintf("Missing paths.docs for %s; provide it in params_user.yaml or the archived run.", original_hash))
  }
  if (is.null(params_user$operating_system$system) || !nzchar(params_user$operating_system$system)) {
    params_user$operating_system$system <- "unix"
  }
  if (is.null(params_user$operating_system$pandoc) || !nzchar(params_user$operating_system$pandoc)) {
    params_user$operating_system$pandoc <- "pandoc"
  }

  hash_params <- prepare_hash_params(params, params_user)
  hash_row <- convert_yaml_to_single_row_df_with_hash(hash_params)
  new_hash <- hash_row$hash
  description <- hash_row$description
  output_dir <- file.path(output_root, new_hash)

  if (isTRUE(opt$dry_run)) {
    message(sprintf("[dry-run] %s -> %s", original_hash, new_hash))
    manifest <- bind_rows(manifest, data.frame(
      original_hash = original_hash,
      new_hash = new_hash,
      description = description,
      source_params = source_params,
      output_dir = output_dir,
      status = "dry-run",
      started_at = started_at,
      ended_at = as.character(Sys.time()),
      log_path = log_path,
      stringsAsFactors = FALSE
    ))
    next
  }

  # session_info.txt is written at the end of a successful stats run. DE.rds is
  # written much earlier and can therefore exist in an incomplete output folder.
  if (dir.exists(output_dir) && file.exists(file.path(output_dir, "session_info.txt")) && !isTRUE(opt$overwrite)) {
    message(sprintf("[skip] %s -> %s already exists", original_hash, new_hash))
    write_reprocess_inputs(output_dir, params, params_user, source_params, source_params_user, params_user_override_path, defaults_path, override_path)
    manifest <- bind_rows(manifest, data.frame(
      original_hash = original_hash,
      new_hash = new_hash,
      description = description,
      source_params = source_params,
      output_dir = output_dir,
      status = "skipped-existing",
      started_at = started_at,
      ended_at = as.character(Sys.time()),
      log_path = log_path,
      stringsAsFactors = FALSE
    ))
    write_manifest(manifest, manifest_path)
    next
  }

  work_dir <- file.path(work_root, original_hash)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(params, file.path(work_dir, "params.yaml"))
  yaml::write_yaml(params_user, file.path(work_dir, "params_user.yaml"))

  message(sprintf("[run] %s -> %s", original_hash, new_hash))
  command_args <- c(
    "--vanilla",
    biostat_script,
    "--params", file.path(work_dir, "params.yaml"),
    "--params-user", file.path(work_dir, "params_user.yaml")
  )
  status <- system2("Rscript", command_args, stdout = log_path, stderr = log_path, env = child_env)
  run_status <- if (identical(status, 0L)) "success" else sprintf("failed:%s", status)
  write_reprocess_inputs(output_dir, params, params_user, source_params, source_params_user, params_user_override_path, defaults_path, override_path)

  manifest <- bind_rows(manifest, data.frame(
    original_hash = original_hash,
    new_hash = new_hash,
    description = description,
    source_params = source_params,
    output_dir = output_dir,
    status = run_status,
    started_at = started_at,
    ended_at = as.character(Sys.time()),
    log_path = log_path,
    stringsAsFactors = FALSE
  ))
  write_manifest(manifest, manifest_path)

  fatal_signal <- !identical(status, 0L) && !is.na(status) && status >= 128
  if (!identical(status, 0L) && (isTRUE(opt$stop_on_error) || fatal_signal)) {
    stop(sprintf("Reprocessing failed for %s. See log: %s", original_hash, log_path))
  }
}

if (!isTRUE(opt$dry_run)) {
  write_manifest(manifest, manifest_path)
  message(sprintf("Manifest written to %s", manifest_path))
}
