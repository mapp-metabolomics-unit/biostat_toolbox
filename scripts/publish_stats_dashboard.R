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
source(file.path(script_dir, "..", "src", "biostats", "utils.R"), local = TRUE)
source(file.path(script_dir, "..", "src", "biostats", "runner.R"), local = TRUE)
source_dir <- file.path(script_dir, "..", "dashboard")
batch_dir <- normalizePath(args[[2L]], mustWork = TRUE)
if (!dir.exists(batch_dir)) stop("Batch must be a directory: ", batch_dir, call. = FALSE)
output_root <- file.path(batch_dir, "results", "stats")
if (!dir.exists(output_root)) stop("Missing Biostats results directory: ", output_root, call. = FALSE)
dashboard_dir <- file.path(output_root, "dashboard")

# GitHub Pages serves the repository's docs/ tree at the repository root.
# Never construct a published URL from an arbitrary local directory name.
project_dir <- dirname(batch_dir)
docs_dir <- dirname(project_dir)
if (!identical(basename(docs_dir), "docs") ||
    !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", basename(dirname(docs_dir))) ||
    !grepl("^mapp_project_[A-Za-z0-9_-]+$", basename(project_dir)) ||
    !grepl("^mapp_batch_[A-Za-z0-9_-]+$", basename(batch_dir))) {
  stop("Batch must be under <group-repo>/docs/mapp_project_*/mapp_batch_* with URL-safe names",
       call. = FALSE)
}
pages_root <- paste0("https://mapp-metabolomics-unit.github.io/",
                     basename(dirname(docs_dir)), "/", basename(project_dir),
                     "/", basename(batch_dir), "/results/stats/dashboard/")

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

validated_count <- function(value, label, minimum = 0) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum || value != floor(value)) {
    stop("Invalid validated ", label, call. = FALSE)
  }
  value
}

# Sealed manifest facts may still include user-entered text. Escape it as text,
# rather than allowing it to supply Markdown links or HTML.
markdown_text <- function(value, label) {
  value <- valid_description(value, label)
  if (grepl("https?://|www\\.|mailto:", value, ignore.case = TRUE)) {
    stop("URL in ", label, " cannot be published as Markdown text", call. = FALSE)
  }
  characters <- strsplit(value, "", fixed = TRUE)[[1L]]
  escaped <- characters %in% c("\\", "`", "*", "_", "{", "}", "[", "]",
                                "(", ")", "#", "+", "-", ".", "!", "|", ">", "<", "~")
  paste0(ifelse(escaped, "\\", ""), characters, collapse = "")
}

readme_facts <- function(manifest) {
  validation <- manifest$validation
  roles <- validation$role_counts
  groups <- validation$group_counts
  if (!is.list(roles) || !setequal(names(roles), c("sample", "qc", "blank")) ||
      anyDuplicated(names(roles)) ||
      !is.list(groups) || !length(groups) || is.null(names(groups)) ||
      anyNA(names(groups)) || anyDuplicated(names(groups)) || any(!nzchar(names(groups)))) {
    stop("Invalid validated role or group counts", call. = FALSE)
  }
  biological <- validated_count(roles$sample, "biological sample count", 1)
  qc <- validated_count(roles$qc, "QC count")
  blank <- validated_count(roles$blank, "blank count")
  injections <- validated_count(validation$dimensions$samples, "injection count", 1)
  features <- validated_count(validation$dimensions$features, "measured feature count", 1)
  group_sizes <- vapply(groups, validated_count, numeric(1), label = "group count", minimum = 1)
  if (sum(group_sizes) > biological || biological + qc + blank > injections) {
    stop("Inconsistent validated injection and group counts", call. = FALSE)
  }
  analysis_labels <- c(pca = "PCA", pcoa = "PCoA", omnibus = "omnibus",
                       differential = "differential", plsda = "PLS-DA")
  enabled <- names(analysis_labels)[vapply(names(analysis_labels), function(name)
    isTRUE(manifest$effective_recipe$analyses[[name]]$enabled), logical(1))]
  list(
    counts = paste0("Validated injections: ", injections, "; biological samples: ", biological,
                    "; QC injections: ", qc, "; blank injections: ", blank,
                    "; measured features: ", features, "."),
    groups = paste0(markdown_text(manifest$effective_recipe$design$group %||%
                                    manifest$effective_recipe$display_group %||% "sample_type",
                                  "group column"),
                    ": ", paste0(vapply(names(groups), markdown_text, character(1),
                                         label = "group name"), " (", group_sizes, ")",
                                  collapse = ", ")),
    analyses = if (length(enabled)) paste(analysis_labels[enabled], collapse = ", ") else
      "none enabled",
    selected = sum(group_sizes)
  )
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

study_description <- function(manifest) {
  count <- validated_count(manifest$validation$role_counts$sample, "biological sample count", 1)
  valid_description(paste0("Shared study (", count, " biological samples); ",
                           manifest$effective_recipe$display_group %||% "sample_type",
                           " display groups"), "study description")
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

reject_tree_symlinks <- function(directory) {
  pending <- directory
  while (length(pending)) {
    current <- pending[[length(pending)]]
    pending <- pending[-length(pending)]
    children <- list.files(current, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    for (child in children) {
      if (is_symlink(child)) stop("Symlink in completed run: ", child, call. = FALSE)
      if (dir.exists(child)) pending <- c(pending, child)
    }
  }
}

read_run <- function(run_dir, hash, child = FALSE) {
  marker_path <- file.path(run_dir, "COMPLETE")
  manifest_path <- file.path(run_dir, "manifest.yaml")
  if (is_symlink(marker_path) || is_symlink(manifest_path)) {
    stop("Symlink in completed run: ", run_dir, call. = FALSE)
  }
  # Staging directories have non-hash names; a manifest in a final hash directory
  # without its seal is a damaged publication, not an ignorable staging run.
  if (!file.exists(marker_path)) {
    if (child || file.exists(manifest_path)) {
      stop("Missing COMPLETE marker in final run directory: ", run_dir, call. = FALSE)
    }
    return(NULL)
  }
  if (!file.exists(manifest_path)) stop("Missing manifest for completed run: ", run_dir, call. = FALSE)
  reject_tree_symlinks(run_dir)
  if (!completed_run(run_dir, hash)) stop("Completed run is missing its COMPLETE seal: ", run_dir, call. = FALSE)
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
  kind <- if (is.null(manifest$kind)) "legacy" else
    required_string(manifest$kind, paste0(hash, " kind"))
  if (!is.null(manifest$kind) && identical(kind, "legacy"))
    stop("Invalid manifest kind: ", run_dir, call. = FALSE)
  if (!kind %in% c("legacy", "study", "comparison") ||
      (child && !identical(kind, "comparison")) ||
      (!child && identical(kind, "comparison")) ||
      (!identical(kind, "legacy") && !identical(manifest$effective_recipe$kind, kind))) {
    stop("Invalid manifest kind for run location: ", run_dir, call. = FALSE)
  }
  if (child && (!identical(manifest$effective_recipe$study_hash, manifest$study_hash) ||
                !is.character(manifest$study_hash) || length(manifest$study_hash) != 1L ||
                is.na(manifest$study_hash) || !grepl(hex_hash, manifest$study_hash))) {
    stop("Invalid comparison study_hash: ", run_dir, call. = FALSE)
  }
  for (field in c("dataset_hash", "preprocess_hash")) {
    if (!is.character(manifest$hashes[[field]]) || length(manifest$hashes[[field]]) != 1L ||
        is.na(manifest$hashes[[field]]) || !grepl(hex_hash, manifest$hashes[[field]])) {
      stop("Invalid ", field, " in completed manifest: ", run_dir, call. = FALSE)
    }
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
    if (identical(kind, "study")) study_description(manifest) else derived_description(manifest)
  }

  outputs <- manifest$outputs
  paths <- names(outputs)
  if (!is.list(outputs) || !length(outputs) || is.null(paths) || anyNA(paths) ||
      anyDuplicated(paths) || any(!nzchar(paths))) {
    stop("Invalid output map: ", run_dir, call. = FALSE)
  }
  for (path in paths) {
    # Output keys are browser-relative paths, never absolute paths or traversal.
    if (!grepl("^(tables|plots|objects)/[A-Za-z0-9][A-Za-z0-9._/-]*$", path) ||
        grepl("(^|/)\\.{1,2}(/|$)|//|\\\\", path)) {
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
    if (!identical(sha256_file(relative), checksum)) {
      stop("Completed output checksum mismatch: ", relative, call. = FALSE)
    }
  }

  mapping <- manifest$effective_dataset
  mapping[c("batch_dir", "metadata", "quantification", "annotations", "input_sha256")] <- NULL
  mapping$feature_id_column <- mapping$feature_id_column %||% "row ID"
  mapping$annotation_sources <- names(Filter(Negate(is.null), manifest$effective_dataset$annotations %||% list()))
  recipe <- scientific_fields(manifest$effective_recipe)
  list(hash = hash, kind = if (identical(kind, "legacy")) NULL else kind,
       study_hash = if (identical(kind, "study")) hash else if (child) manifest$study_hash else NULL,
       hashes = manifest$hashes, effective_dataset = manifest$effective_dataset,
       preprocessing = manifest$effective_recipe$preprocessing,
       dataset_id = dataset_id, created_at = created_at,
       description = description, recipe_sha256 = sha256_object(recipe),
       readme = c(readme_facts(manifest), list(
         description = if ("description" %in% names(manifest))
           markdown_text(manifest$description, "run description") else NULL)),
       comparison = list(recipe = recipe, inputs = manifest$inputs, mapping = mapping,
                         environment = manifest$environment),
       schema_version = manifest$schema_version, outputs = outputs,
       design_group = design_group, intensity_measure = intensity_measure,
       transformation = transformation,
       preprocessing_filters = list(
         blank_filter = recipe$preprocessing$blank_filter,
         qc_rsd_filter = recipe$preprocessing$qc_rsd_filter,
         missing_values = recipe$preprocessing$missing_values,
         normalization = recipe$preprocessing$normalization))
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

comparison_root <- file.path(output_root, "comparisons")
if (is_symlink(comparison_root) || (file.exists(comparison_root) && !dir.exists(comparison_root))) {
  stop("Invalid comparisons directory: ", comparison_root, call. = FALSE)
}
if (dir.exists(comparison_root)) {
  parent_hashes <- list.files(comparison_root, all.files = FALSE, no.. = TRUE)
  for (study_hash in parent_hashes) {
    if (!grepl(hex_hash, study_hash)) {
      stop("Invalid study directory in comparisons: ", study_hash, call. = FALSE)
    }
    parent_dir <- file.path(comparison_root, study_hash)
    if (is_symlink(parent_dir) || !dir.exists(parent_dir)) {
      stop("Invalid comparison parent directory: ", parent_dir, call. = FALSE)
    }
    parents <- Filter(function(run) identical(run$hash, study_hash) &&
                        identical(run$kind, "study"), runs)
    if (length(parents) != 1L) {
      stop("Comparisons require a sealed published study: ", parent_dir, call. = FALSE)
    }
    parent <- parents[[1L]]
    if (!"tables/sample_metadata.tsv" %in% names(parent$outputs)) {
      stop("Study is missing shared sample metadata: ", parent_dir, call. = FALSE)
    }
    children <- list.files(parent_dir, all.files = FALSE, no.. = TRUE)
    for (child_hash in children) {
      if (!grepl(hex_hash, child_hash)) {
        stop("Invalid comparison directory: ", child_hash, call. = FALSE)
      }
      child_dir <- file.path(parent_dir, child_hash)
      if (is_symlink(child_dir) || !dir.exists(child_dir)) {
        stop("Invalid comparison directory: ", child_dir, call. = FALSE)
      }
      child <- read_run(child_dir, child_hash, child = TRUE)
      if (!identical(child$study_hash, study_hash) ||
          !identical(child$dataset_id, parent$dataset_id) ||
          !identical(child$hashes$dataset_hash, parent$hashes$dataset_hash) ||
          !identical(child$hashes$preprocess_hash, parent$hashes$preprocess_hash) ||
          !identical(child$effective_dataset, parent$effective_dataset) ||
          !identical(child$preprocessing, parent$preprocessing)) {
        stop("Comparison does not match its sealed study: ", child_dir, call. = FALSE)
      }
      child_paths <- names(child$outputs)
      if (any(child_paths %in% names(parent$outputs)) ||
          !all(grepl("^(tables/(comparison_sample_metadata|omnibus|differential|plsda_[A-Za-z_]+)\\.tsv|plots/volcano_[^/]+\\.(png|pdf)|objects/[^/]+)$",
                     child_paths)) ||
          !"tables/comparison_sample_metadata.tsv" %in% child_paths) {
        stop("Invalid child-owned outputs: ", child_dir, call. = FALSE)
      }
      child$sources <- as.list(stats::setNames(
        rep(paste("comparisons", study_hash, child_hash, sep = "/"), length(child_paths)), child_paths))
      child$outputs <- c(parent$outputs, child$outputs)
      runs[[length(runs) + 1L]] <- child
    }
  }
}
if (anyDuplicated(vapply(runs, `[[`, character(1), "hash"))) {
  stop("Duplicate run hash across studies and comparisons", call. = FALSE)
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


params_log <- file.path(output_root, "biostats_params_log.tsv")
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
    stop("Could not atomically publish biostats_params_log.tsv", call. = FALSE)
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
retired_assets <- file.path(dashboard_dir,
                            c("recipe-builder.html", "recipe-builder.js", "recipe-builder.css"))
for (path in retired_assets) {
  if (is_symlink(path) || dir.exists(path)) stop("Unexpected retired dashboard asset: ", path, call. = FALSE)
  if (file.exists(path) && !file.remove(path)) {
    stop("Could not remove retired dashboard asset: ", path, call. = FALSE)
  }
}

local({
  runs <- lapply(runs, function(run) {
    for (field in c("comparison", "readme", "effective_dataset", "preprocessing"))
      run[[field]] <- NULL
    if (is.null(run$kind)) {
      run$kind <- NULL
      run$hashes <- NULL
    }
    if (is.null(run$study_hash)) run$study_hash <- NULL
    run
  })
  catalogue <- list(schema_version = 1L, runs = runs)
  temporary <- tempfile(pattern = ".runs-", tmpdir = dashboard_dir)
  on.exit(unlink(temporary), add = TRUE)
  jsonlite::write_json(catalogue, temporary, auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(temporary, file.path(dashboard_dir, "runs.json"))) {
    stop("Could not atomically publish runs.json", call. = FALSE)
  }
})

readme_path <- file.path(output_root, "README.md")
if (is_symlink(readme_path)) stop("Symlinked Biostats README: ", readme_path, call. = FALSE)
dashboard_published <- file.exists(file.path(dashboard_dir, "index.html")) &&
  file.exists(file.path(dashboard_dir, "runs.json"))
readme_lines <- c("# Biostats results", "",
                  if (length(runs)) paste(length(runs), "completed sealed",
                                          if (length(runs) == 1L) "run." else "runs.") else
                    "No completed runs are available.", "")
if (dashboard_published) {
  readme_lines <- c(readme_lines, paste0("[Open the Biostats dashboard](", pages_root, ")"), "")
}
if (length(runs)) {
  readme_lines <- c(readme_lines, "## Available statistics", "")
  studies <- Filter(function(run) identical(run$kind, "study"), runs)
  standalone <- Filter(function(run) is.null(run$kind), runs)
  append_facts <- function(lines, run) {
    facts <- run$readme
    lines <- c(lines, paste0("- ", facts$counts),
               paste0("- Biological groups (", facts$groups, ")."),
               if (identical(run$kind, "comparison"))
                 paste0("- Selected biological samples: ", facts$selected, "."),
               paste0("- Enabled analyses: ", facts$analyses, "."))
    if (!is.null(facts$description))
      lines <- c(lines, paste0("- Description: ", facts$description))
    c(lines, "")
  }
  for (run in c(studies, standalone)) {
    readme_lines <- c(readme_lines,
                      paste0("### [", if (identical(run$kind, "study")) "Study " else "Run ",
                             run$hash, "](", pages_root, "?run=", run$hash, ")"), "")
    readme_lines <- append_facts(readme_lines, run)
    if (identical(run$kind, "study")) {
      children <- Filter(function(candidate) identical(candidate$kind, "comparison") &&
                           identical(candidate$study_hash, run$hash), runs)
      if (length(children)) readme_lines <- c(readme_lines, "#### Comparisons", "")
      for (child in children) {
        readme_lines <- c(readme_lines,
                          paste0("- [Comparison ", child$hash, "](", pages_root,
                                 "?run=", child$hash, ") — ",
                                 markdown_text(child$description, "comparison description")), "")
        readme_lines <- append_facts(readme_lines, child)
      }
    }
  }
}
local({
  temporary <- tempfile(pattern = ".biostats-readme-", tmpdir = output_root)
  on.exit(unlink(temporary), add = TRUE)
  writeLines(readme_lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, readme_path)) {
    stop("Could not atomically publish Biostats README.md", call. = FALSE)
  }
})
message("Published ", length(runs), " sealed runs to ", dashboard_dir)
