#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3L && identical(args[[1L]], "--")) args <- args[-1L]
if (length(args) != 2L || !identical(args[[1L]], "--batch") || !nzchar(args[[2L]])) {
  stop("Usage: Rscript scripts/publish_stats_dashboard.R --batch <batch directory>", call. = FALSE)
}

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg[[1L]]) else {
  file_option <- match("-f", full_args)
  if (is.na(file_option) || file_option == length(full_args)) stop("Run this script from its file path", call. = FALSE)
  full_args[[file_option + 1L]]
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(script_dir, "..", "src", "v2", "utils.R"), local = TRUE)
source(file.path(script_dir, "..", "src", "v2", "runner.R"), local = TRUE)
source_dir <- file.path(script_dir, "..", "dashboard")
batch_dir <- normalizePath(args[[2L]], mustWork = TRUE)
if (!dir.exists(batch_dir)) stop("Batch must be a directory: ", batch_dir, call. = FALSE)
output_root <- file.path(batch_dir, "results", "stats_v2")
if (!dir.exists(output_root)) stop("Missing V2 results directory: ", output_root, call. = FALSE)
dashboard_dir <- file.path(output_root, "dashboard")
canonical_dataset_id <- paste(basename(dirname(batch_dir)), basename(batch_dir), sep = "/")

hex_hash <- "^[0-9a-f]{64}$"
is_symlink <- function(path) {
  link <- Sys.readlink(path)
  !is.na(link) && nzchar(link)
}

required_string <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("Invalid ", label, call. = FALSE)
  }
  value
}

valid_description <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || nchar(value, type = "chars") > 160L ||
      !identical(value, trimws(value)) || grepl("[[:cntrl:]]", value)) {
    stop("Invalid ", label, " (expected one trimmed line of at most 160 characters)", call. = FALSE)
  }
  value
}

derived_description <- function(manifest) {
  recipe <- manifest$effective_recipe
  counts <- manifest$validation$group_counts
  groups <- names(counts)
  if (is.null(groups) || !length(groups) || anyNA(groups) || any(!nzchar(groups))) {
    stop("Missing validated group counts for legacy run", call. = FALSE)
  }
  group_names <- paste(groups, collapse = ", ")
  if (nchar(group_names, type = "chars") > 64L) {
    group_names <- paste(length(groups), "groups")
  }
  sample_count <- manifest$validation$role_counts$sample
  if (!is.numeric(sample_count) || length(sample_count) != 1L || is.na(sample_count) ||
      sample_count < 1 || sample_count != as.integer(sample_count)) {
    stop("Invalid validated sample count for legacy run", call. = FALSE)
  }
  contrasts <- recipe$design$contrasts
  comparison <- if (identical(contrasts, "all")) {
    "all contrasts"
  } else if (is.list(contrasts) && length(contrasts)) {
    pairs <- vapply(contrasts, function(pair) {
      paste(required_string(pair$numerator, "contrast numerator"), "vs",
            required_string(pair$denominator, "contrast denominator"))
    }, character(1))
    paste(pairs, collapse = ", ")
  } else {
    "no contrasts"
  }
  if (nchar(comparison, type = "chars") > 80L && is.list(contrasts)) {
    comparison <- paste(length(contrasts), "planned contrasts")
  }
  analyses <- c(pca = "PCA", pcoa = "PCoA", omnibus = "omnibus",
                differential = "differential", plsda = "PLS-DA")
  enabled <- names(analyses)[vapply(names(analyses), function(key)
    isTRUE(recipe$analyses[[key]]$enabled), logical(1))]
  details <- paste0(recipe$design$group, ": ", group_names, " (", sample_count,
                    " samples); ", comparison, "; ",
                    if (length(enabled)) paste(analyses[enabled], collapse = ", ") else "no enabled analyses")
  if (nchar(details, type = "chars") > 160L) {
    details <- paste0(recipe$design$group, ": ", group_names, " (", sample_count,
                      " samples); ", comparison)
  }
  if (nchar(details, type = "chars") > 160L) {
    details <- paste0(length(groups), " groups (", sample_count, " samples); ", comparison)
  }
  valid_description(details, "derived description")
}

scientific_fields <- function(recipe) {
  recipe$presentation <- NULL
  recipe$output_root <- NULL
  recipe
}

# Report changed leaf paths, but keep positional arrays (e.g. contrast pairs) together.
changed_paths <- function(before, after, prefix) {
  if (identical(before, after)) return(character())
  if (is.list(before) && is.list(after) &&
      !is.null(names(before)) && !is.null(names(after)) &&
      !anyDuplicated(names(before)) && !anyDuplicated(names(after))) {
    keys <- sort(union(names(before), names(after)))
    return(unlist(lapply(keys, function(key) {
      path <- paste(prefix, key, sep = ".")
      if (!(key %in% names(before)) || !(key %in% names(after))) return(path)
      changed_paths(before[[key]], after[[key]], path)
    }), use.names = FALSE))
  }
  prefix
}

changes_since <- function(previous, current) {
  if (is.null(previous)) return("First published run for dataset")
  parts <- c(
    changed_paths(previous$recipe, current$recipe, "recipe"),
    changed_paths(previous$inputs, current$inputs, "inputs"),
    changed_paths(previous$mapping, current$mapping, "dataset"))
  old_env <- previous$environment
  new_env <- current$environment
  code <- changed_paths(old_env$code_sha256, new_env$code_sha256, "")
  if (length(code)) {
    # The checksum keys name source files, not hashes; display just their filenames.
    parts <- c(parts, paste0("code: ", paste(basename(sub("^\\.", "", code)), collapse = ", ")))
  }
  other_env <- changed_paths(old_env[names(old_env) != "code_sha256"],
                             new_env[names(new_env) != "code_sha256"], "environment")
  parts <- c(parts, other_env)
  if (!length(parts)) "No scientific, input, dataset mapping, or environment differences recorded" else
    paste(parts, collapse = "; ")
}

safe_ledger_text <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value) || nchar(value, type = "bytes") > 8192L ||
      !identical(value, trimws(value)) || grepl("[[:cntrl:]]", value)) {
    stop("Invalid ", label, " in completed-job index", call. = FALSE)
  }
  value
}

read_run <- function(run_dir, hash) {
  marker_path <- file.path(run_dir, "COMPLETE")
  manifest_path <- file.path(run_dir, "manifest.yaml")
  if (is_symlink(marker_path) || is_symlink(manifest_path)) {
    stop("Symlink in completed run: ", run_dir, call. = FALSE)
  }
  # Staging directories have non-hash names; a manifest in a final hash directory
  # without its seal is a damaged publication, not an ignorable staging run.
  if (!file.exists(marker_path)) {
    if (file.exists(manifest_path)) {
      stop("Manifest without COMPLETE marker: ", run_dir, call. = FALSE)
    }
    return(NULL)
  }
  if (!file.exists(manifest_path)) stop("Missing manifest for completed run: ", run_dir, call. = FALSE)
  marker <- tryCatch(yaml::read_yaml(marker_path), error = function(e) NULL)
  manifest <- tryCatch(yaml::read_yaml(manifest_path), error = function(e) NULL)
  if (!is.list(marker) || !identical(marker$run_hash, hash) ||
      !is.character(marker$content_sha256) || length(marker$content_sha256) != 1L ||
      is.na(marker$content_sha256) || !grepl(hex_hash, marker$content_sha256)) {
    stop("Invalid COMPLETE marker: ", run_dir, call. = FALSE)
  }
  if (!is.list(manifest) || !identical(manifest$schema_version, "1.0.0") ||
      !identical(manifest$status, "complete") || !is.list(manifest$hashes) ||
      !identical(manifest$hashes$run_hash, hash)) {
    stop("Unsupported or invalid completed manifest: ", run_dir, call. = FALSE)
  }

  created_at <- required_string(manifest$created_at, paste0(hash, " created_at"))
  if (!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} UTC$", created_at)) {
    stop("Invalid created_at: ", run_dir, call. = FALSE)
  }
  dataset_id <- required_string(manifest$dataset_id, paste0(hash, " dataset_id"))
  design_group <- required_string(manifest$effective_recipe$design$group, paste0(hash, " design.group"))
  intensity_measure <- required_string(manifest$effective_dataset$intensity_measure,
                                       paste0(hash, " intensity_measure"))
  transformation <- required_string(manifest$effective_recipe$preprocessing$transformation$method,
                                    paste0(hash, " transformation"))
  description <- if ("description" %in% names(manifest)) {
    valid_description(manifest$description, paste0(hash, " description"))
  } else {
    derived_description(manifest)
  }
  if (!completed_run(run_dir, hash)) stop("Completed run is missing its COMPLETE seal: ", run_dir, call. = FALSE)

  outputs <- manifest$outputs
  paths <- names(outputs)
  if (!is.list(outputs) || !length(outputs) || is.null(paths) || anyNA(paths) ||
      anyDuplicated(paths) || any(!nzchar(paths))) {
    stop("Invalid output map: ", run_dir, call. = FALSE)
  }
  for (path in paths) {
    # Output keys are browser-relative paths, never absolute paths or traversal.
    if (!grepl("^(tables|plots|objects)/[^/]+(/[^/]+)*$", path) ||
        grepl("(^|/)\\.{1,2}(/|$)|\\\\", path)) {
      stop("Unsafe output path in ", run_dir, ": ", path, call. = FALSE)
    }
    checksum <- outputs[[path]]
    if (!is.character(checksum) || length(checksum) != 1L || is.na(checksum) ||
        !grepl(hex_hash, checksum)) {
      stop("Invalid output checksum in ", run_dir, ": ", path, call. = FALSE)
    }
    components <- strsplit(path, "/", fixed = TRUE)[[1L]]
    relative <- run_dir
    for (component in components) {
      relative <- file.path(relative, component)
      if (is_symlink(relative)) stop("Symlink in completed output: ", relative, call. = FALSE)
    }
    if (!file.exists(relative) || dir.exists(relative)) {
      stop("Missing completed output: ", relative, call. = FALSE)
    }
  }

  recipe <- manifest$effective_recipe
  source <- NULL
  if (identical(dataset_id, canonical_dataset_id) &&
      identical(manifest$effective_dataset$id, canonical_dataset_id) &&
      is.null(recipe$sample_plan) && is.null(manifest$effective_dataset$sample_plan) &&
      !is.null(outputs[["tables/sample_metadata.tsv"]]) &&
      is.list(recipe$roles) && is.list(recipe$preprocessing) && is.list(recipe$analyses) &&
      !is.null(recipe$seed) && !is.null(recipe$design$group)) {
    # The COMPLETE seal was checked above for every published run.
    source <- list(run_hash = hash, dataset_id = dataset_id,
                   group_column = design_group, roles = recipe$roles,
                   defaults = list(seed = recipe$seed, preprocessing = recipe$preprocessing,
                                   analyses = recipe$analyses))
  }
  mapping <- manifest$effective_dataset
  mapping[c("batch_dir", "metadata", "quantification", "annotations", "input_sha256")] <- NULL
  mapping$feature_id_column <- mapping$feature_id_column %||% "row ID"
  mapping$annotation_sources <- names(Filter(Negate(is.null), manifest$effective_dataset$annotations %||% list()))
  recipe <- scientific_fields(recipe)
  list(hash = hash, dataset_id = dataset_id, created_at = created_at,
       description = description, recipe_sha256 = sha256_object(recipe),
       comparison = list(recipe = recipe, inputs = manifest$inputs, mapping = mapping,
                         environment = manifest$environment),
       schema_version = manifest$schema_version, outputs = outputs,
       design_group = design_group, intensity_measure = intensity_measure,
       transformation = transformation,
       preprocessing_filters = list(
         blank_filter = recipe$preprocessing$blank_filter,
         qc_rsd_filter = recipe$preprocessing$qc_rsd_filter,
         missing_values = recipe$preprocessing$missing_values,
         normalization = recipe$preprocessing$normalization),
       builder_candidate = source)
}

entries <- list.files(output_root, all.files = FALSE, no.. = TRUE, full.names = FALSE)
hashes <- entries[grepl(hex_hash, entries)]
runs <- list()
for (hash in hashes) {
  run_dir <- file.path(output_root, hash)
  if (is_symlink(run_dir)) stop("Symlinked run directory: ", run_dir, call. = FALSE)
  if (!dir.exists(run_dir)) stop("Run hash is not a directory: ", run_dir, call. = FALSE)
  run <- read_run(run_dir, hash)
  if (!is.null(run)) runs[[length(runs) + 1L]] <- run
}
if (length(runs)) {
  dates <- vapply(runs, `[[`, character(1), "created_at")
  ids <- vapply(runs, `[[`, character(1), "hash")
  runs <- runs[order(dates, ids, decreasing = TRUE)]
}
if (length(runs)) {
  latest_by_dataset <- new.env(parent = emptyenv())
  # Runs are newest first (hash breaks second-resolution timestamp ties).
  for (i in rev(seq_along(runs))) {
    dataset <- runs[[i]]$dataset_id
    previous <- latest_by_dataset[[dataset]]
    runs[[i]]$changes_vs_previous_dataset_run <-
      safe_ledger_text(changes_since(previous, runs[[i]]$comparison),
                       "changes_vs_previous_dataset_run")
    latest_by_dataset[[dataset]] <- runs[[i]]$comparison
  }
}


params_log <- file.path(output_root, "params_log.tsv")
if (is_symlink(params_log)) stop("Symlinked completed-job index: ", params_log, call. = FALSE)
local({
  log_rows <- data.frame(
    run_hash = vapply(runs, `[[`, character(1), "hash"),
    created_at = vapply(runs, `[[`, character(1), "created_at"),
    dataset_id = vapply(runs, `[[`, character(1), "dataset_id"),
    description = vapply(runs, `[[`, character(1), "description"),
    recipe_sha256 = vapply(runs, `[[`, character(1), "recipe_sha256"),
    changes_vs_previous_dataset_run = vapply(runs, `[[`, character(1),
                                                  "changes_vs_previous_dataset_run"))
  if (any(grepl("[[:cntrl:]]", log_rows$dataset_id))) {
    stop("Invalid dataset ID in completed-job index", call. = FALSE)
  }
  temporary <- tempfile(pattern = ".params-log-", tmpdir = output_root)
  on.exit(unlink(temporary), add = TRUE)
  write.table(log_rows, temporary, sep = "\t", row.names = FALSE,
              col.names = TRUE, quote = FALSE, na = "")
  if (!file.rename(temporary, params_log)) {
    stop("Could not atomically publish params_log.tsv", call. = FALSE)
  }
})

assets <- c("index.html", "dashboard.js", "dashboard.css")
asset_sources <- file.path(source_dir, assets)
plotly_source <- system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package = "plotly")
if (any(!file.exists(asset_sources)) || !nzchar(plotly_source) || !file.exists(plotly_source)) {
  stop("Missing dashboard source asset or installed plotly.min.js", call. = FALSE)
}
if (is_symlink(dashboard_dir)) stop("Symlinked dashboard directory: ", dashboard_dir, call. = FALSE)
if (!dir.exists(dashboard_dir) && !dir.create(dashboard_dir, recursive = TRUE)) {
  stop("Could not create dashboard directory: ", dashboard_dir, call. = FALSE)
}
if (!all(file.copy(c(asset_sources, plotly_source), dashboard_dir, overwrite = TRUE))) {
  stop("Could not copy dashboard assets to: ", dashboard_dir, call. = FALSE)
}
plotly_destination <- file.path(dashboard_dir, "plotly-latest.min.js")
if (!file.rename(plotly_destination, file.path(dashboard_dir, "plotly.min.js"))) {
  stop("Could not name vendored plotly.min.js", call. = FALSE)
}

local({
  candidates <- which(vapply(runs, function(run) !is.null(run$builder_candidate), logical(1)))
  builder_source <- if (length(candidates)) runs[[candidates[[1L]]]]$builder_candidate else NULL
  runs <- lapply(runs, function(run) {
    run$builder_candidate <- NULL
    run$comparison <- NULL
    run
  })
  catalogue <- list(schema_version = 1L, runs = runs)
  if (!is.null(builder_source)) catalogue$builder_source <- builder_source
  temporary <- tempfile(pattern = ".runs-", tmpdir = dashboard_dir)
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(catalogue, temporary, auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(temporary, file.path(dashboard_dir, "runs.json"))) {
    stop("Could not atomically publish runs.json", call. = FALSE)
  }
})
message("Published ", length(runs), " sealed runs to ", dashboard_dir)
