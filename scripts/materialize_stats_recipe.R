#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg) else {
  file_option <- match("-f", args_full)
  if (is.na(file_option) || file_option == length(args_full)) stop("Run this script from its file path", call. = FALSE)
  args_full[[file_option + 1L]]
}
repo_root <- normalizePath(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), ".."), mustWork = TRUE)
source(file.path(repo_root, "src", "v2", "utils.R"))
source(file.path(repo_root, "src", "v2", "runner.R"))
source_v2_modules(repo_root)

exact_keys <- function(value, keys, label) {
  if (!is.list(value) || is.null(names(value)) || anyNA(names(value)) ||
      anyDuplicated(names(value)) || !setequal(names(value), keys) || length(value) != length(keys))
    stop(label, " requires exactly: ", paste(keys, collapse = ", "), call. = FALSE)
}
text_value <- function(value, label) {
  config_text(value, label)
  if (!identical(value, trimws(value))) stop(label, " must have no surrounding whitespace", call. = FALSE)
  value
}
number_value <- function(value, label, minimum = 0, maximum = Inf, integer = FALSE) {
  config_number(value, label, minimum, maximum, integer)
  value
}
boolean_value <- function(value, label) { config_boolean(value, label); value }
hex_hash <- "^[0-9a-f]{64}$"
check_hash <- function(hash, label) {
  if (!is.character(hash) || length(hash) != 1L || is.na(hash) || !grepl(hex_hash, hash))
    stop("Invalid ", label, call. = FALSE)
  hash
}
contained <- function(path, directory) {
  path <- normalizePath(path, mustWork = TRUE)
  directory <- normalizePath(directory, mustWork = TRUE)
  if (!startsWith(path, paste0(directory, "/")) || dir.exists(path))
    stop("Input must be a file inside the batch: ", path, call. = FALSE)
  path
}
source_input <- function(path, old_batch, batch, label) {
  text_value(path, label)
  if (grepl("^/", path)) {
    if (!startsWith(path, paste0(old_batch, "/")))
      stop(label, " was not inside the original batch", call. = FALSE)
    path <- substring(path, nchar(old_batch) + 2L)
  }
  if (grepl("\\\\|(^|/)\\.{1,2}(/|$)|^/", path)) stop("Unsafe batch-relative ", label, call. = FALSE)
  contained(file.path(batch, path), batch)
}
validate_plan <- function(plan) {
  exact_keys(plan, c("schema_version", "source_run_hash", "slug", "group_column", "groups",
                     "seed", "preprocessing", "analyses", "contrasts",
                     intersect("description", names(plan))), "plan")
  if (!is.numeric(plan$schema_version) || length(plan$schema_version) != 1L ||
      is.na(plan$schema_version) || plan$schema_version != 1) stop("Unsupported plan schema_version", call. = FALSE)
  check_hash(plan$source_run_hash, "source_run_hash")
  if (!grepl("^[a-z][a-z0-9_-]{2,50}$", text_value(plan$slug, "slug")))
    stop("slug must be 3–51 lowercase filesystem-safe characters starting with a letter", call. = FALSE)
  if ("description" %in% names(plan)) {
    description <- plan$description
    text_value(description, "description")
    if (nchar(description, type = "chars") > 160L || grepl("[[:cntrl:]]", description))
      stop("description must be one line of at most 160 characters without controls", call. = FALSE)
  }
  text_value(plan$group_column, "group_column")
  if (!identical(plan$contrasts, "all")) stop("Only contrasts: 'all' is supported", call. = FALSE)
  number_value(plan$seed, "seed", 0, .Machine$integer.max, TRUE)
  if (!is.list(plan$groups) || !is.null(names(plan$groups)) || !length(plan$groups))
    stop("groups must be a nonempty array", call. = FALSE)
  for (entry in plan$groups) {
    exact_keys(entry, c("source", "target"), "group mapping")
    config_text(entry$source, "group source")
    target <- text_value(entry$target, "group target")
    if (!grepl("^[A-Za-z][A-Za-z0-9_-]*$", target) || identical(target, "NA"))
      stop("group target must be a safe non-NA label", call. = FALSE)
  }
  sources <- vapply(plan$groups, `[[`, character(1), "source")
  if (anyDuplicated(sources)) stop("group source values must be unique", call. = FALSE)
  sections <- list(blank_filter = c("enabled", "minimum_sample_to_blank_ratio", "minimum_blank_prevalence"),
                   qc_rsd_filter = c("enabled", "threshold_percent"), missing_values = "method",
                   normalization = "method", transformation = "method", scaling = "method")
  exact_keys(plan$preprocessing, names(sections), "preprocessing")
  for (section in names(sections)) exact_keys(plan$preprocessing[[section]], sections[[section]], paste0("preprocessing.", section))
  for (section in c("blank_filter", "qc_rsd_filter"))
    boolean_value(plan$preprocessing[[section]]$enabled, paste0(section, ".enabled"))
  number_value(plan$preprocessing$blank_filter$minimum_sample_to_blank_ratio, "blank filter ratio")
  number_value(plan$preprocessing$blank_filter$minimum_blank_prevalence, "blank filter prevalence", 0, 1)
  number_value(plan$preprocessing$qc_rsd_filter$threshold_percent, "QC RSD threshold")
  methods <- list(missing_values = c("none", "half_minimum"), normalization = c("none", "total_sum", "median", "pqn"),
                  transformation = c("none", "log2", "log10", "log1p"), scaling = c("none", "pareto", "autoscale"))
  for (section in names(methods)) config_choice(plan$preprocessing[[section]]$method, methods[[section]], paste0(section, ".method"))
  options <- list(pca = c("enabled", "components"), pcoa = c("enabled", "components", "stage", "distance"),
                  omnibus = "enabled", differential = "enabled",
                  plsda = c("enabled", "components", "folds", "repeats", "permutations"))
  exact_keys(plan$analyses, names(options), "analyses")
  for (name in names(options)) {
    exact_keys(plan$analyses[[name]], options[[name]], paste0("analyses.", name))
    boolean_value(plan$analyses[[name]]$enabled, paste0(name, ".enabled"))
    if ("components" %in% options[[name]]) number_value(plan$analyses[[name]]$components, paste0(name, ".components"), 1, Inf, TRUE)
  }
  config_choice(plan$analyses$pcoa$stage, c("imputed", "normalized", "transformed", "scaled"), "PCoA stage")
  config_choice(plan$analyses$pcoa$distance, c("bray", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), "PCoA distance")
  for (field in c("folds", "repeats", "permutations"))
    number_value(plan$analyses$plsda[[field]], paste0("PLS-DA ", field),
                 if (field == "permutations") 19 else 2, Inf, TRUE)
  invisible(plan)
}

materialize_stats_recipe <- function(batch, plan_path) {
  batch <- normalizePath(batch, mustWork = TRUE)
  if (!dir.exists(batch)) stop("Batch must be a directory", call. = FALSE)
  plan <- jsonlite::read_json(plan_path, simplifyVector = FALSE)
  validate_plan(plan)
  hash <- plan$source_run_hash
  run_dir <- file.path(batch, "results", "stats_v2", hash)
  if (!dir.exists(run_dir) || nzchar(Sys.readlink(run_dir))) stop("Missing sealed source run: ", hash, call. = FALSE)
  if (!completed_run(run_dir, hash)) stop("Source run has no COMPLETE seal: ", hash, call. = FALSE)
  manifest <- read_yaml_file(file.path(run_dir, "manifest.yaml"))
  canonical_id <- paste(basename(dirname(batch)), basename(batch), sep = "/")
  if (!is.list(manifest) || !identical(manifest$schema_version, mapp_v2_schema_version) ||
      !identical(manifest$status, "complete") || !identical(manifest$hashes$run_hash, hash) ||
      !identical(manifest$dataset_id, canonical_id) || !identical(manifest$effective_dataset$id, canonical_id) ||
      is.null(manifest$outputs[["tables/sample_metadata.tsv"]]))
    stop("Source must be a sealed canonical full-batch run with sample metadata", call. = FALSE)
  source_recipe <- manifest$effective_recipe
  source_dataset <- manifest$effective_dataset
  if (!is.null(source_recipe$sample_plan) || !is.null(source_dataset$sample_plan))
    stop("Source must be a canonical full-batch run without a sample plan", call. = FALSE)
  source_roles <- source_recipe$roles %||% dataset_roles(source_dataset)
  old_batch <- text_value(source_dataset$batch_dir, "source batch_dir")
  if (!grepl("^/", old_batch) || grepl("/\\.{1,2}(/|$)", old_batch)) stop("Invalid sealed batch_dir", call. = FALSE)
  paths <- list(metadata = source_input(source_dataset$metadata, old_batch, batch, "metadata"),
                quantification = source_input(source_dataset$quantification, old_batch, batch, "quantification"))
  annotation_fields <- Filter(Negate(is.null), source_dataset$annotations %||% list())
  if (!is.list(annotation_fields) || (length(annotation_fields) &&
      (is.null(names(annotation_fields)) || anyNA(names(annotation_fields)) ||
       anyDuplicated(names(annotation_fields)) || any(!nzchar(names(annotation_fields))))))
    stop("Invalid source annotations", call. = FALSE)
  paths$annotations <- lapply(names(annotation_fields), function(name)
    source_input(annotation_fields[[name]], old_batch, batch, paste0("annotation ", name)))
  names(paths$annotations) <- names(annotation_fields)
  inputs <- manifest$inputs
  if (!is.list(inputs) || !is.list(inputs$annotations) ||
      !setequal(names(inputs$annotations), names(paths$annotations))) stop("Missing original input hashes", call. = FALSE)
  for (field in c("metadata", "quantification")) {
    expected <- check_hash(inputs[[field]], paste0(field, " SHA256"))
    if (!identical(sha256_file(paths[[field]]), expected)) stop("Original ", field, " SHA256 mismatch", call. = FALSE)
  }
  for (name in names(paths$annotations)) {
    expected <- check_hash(inputs$annotations[[name]], paste0(name, " SHA256"))
    if (!identical(sha256_file(paths$annotations[[name]]), expected)) stop("Original annotation SHA256 mismatch: ", name, call. = FALSE)
  }
  sealed_metadata <- file.path(run_dir, "tables", "sample_metadata.tsv")
  if (!identical(sha256_file(sealed_metadata), check_hash(manifest$outputs[["tables/sample_metadata.tsv"]], "sealed metadata SHA256")))
    stop("Sealed sample metadata checksum mismatch", call. = FALSE)

  original <- utils::read.delim(paths$metadata, check.names = FALSE, stringsAsFactors = FALSE,
                                na.strings = c("", "NA"))
  mapping <- source_dataset$metadata_columns %||% list()
  filename_column <- mapping$filename %||% "filename"
  role_column <- if (identical(source_roles$column, "sample_type"))
    mapping$sample_type %||% "sample_type" else source_roles$column
  role_sample <- source_roles$levels$sample
  role_qc <- source_roles$levels$qc
  role_blank <- source_roles$levels$blank
  if (!all(c(filename_column, role_column, plan$group_column) %in% names(original)) ||
      "analysis_group" %in% names(original) || anyDuplicated(names(original)))
    stop("Missing or conflicting original metadata columns", call. = FALSE)
  filenames <- original[[filename_column]]
  if (anyNA(filenames) || any(!nzchar(filenames)) || anyDuplicated(filenames))
    stop("Original metadata filenames must be unique and nonempty", call. = FALSE)
  sealed <- utils::read.delim(sealed_metadata, check.names = FALSE, stringsAsFactors = FALSE,
                              na.strings = c("", "NA"))
  if (!"filename" %in% names(sealed) || !setequal(as.character(sealed$filename), as.character(filenames)))
    stop("Sealed metadata does not cover the full original batch", call. = FALSE)
  roles <- original[[role_column]]
  biological <- !is.na(roles) & roles == role_sample
  other <- !is.na(roles) & roles %in% c(role_qc, role_blank)
  if (any(!biological & !other)) stop("Unknown source sample roles", call. = FALSE)
  groups <- as.character(original[[plan$group_column]][biological])
  if (anyNA(groups) || any(!nzchar(trimws(groups)))) stop("Missing biological source groups", call. = FALSE)
  sources <- vapply(plan$groups, `[[`, character(1), "source")
  targets <- vapply(plan$groups, `[[`, character(1), "target")
  if (any(!sources %in% groups)) stop("Plan contains unknown biological group source", call. = FALSE)
  selected <- match(as.character(original[[plan$group_column]]), sources)
  keep <- other | (biological & !is.na(selected))
  metadata <- original[keep, , drop = FALSE]
  metadata$analysis_group <- ifelse(other[keep], NA_character_, targets[selected[keep]])
  group_counts <- table(metadata$analysis_group[biological[keep]])
  if (length(group_counts) < 2L || any(group_counts < 3L)) stop("At least two target groups with three samples each are required", call. = FALSE)
  if (isTRUE(plan$analyses$plsda$enabled) &&
      (any(group_counts < 4L) || any(group_counts < plan$analyses$plsda$folds)))
    stop("PLS-DA requires four samples per target group and folds <= smallest group", call. = FALSE)

  quant <- utils::read.csv(paths$quantification, check.names = FALSE, nrows = 0)
  if (anyDuplicated(names(quant))) stop("Duplicate original quantification column names", call. = FALSE)
  intensity <- grep(" Peak (height|area)$", names(quant))
  sample_names <- sub(" Peak (height|area)$", "", names(quant)[intensity])
  if (any(!sample_names %in% filenames))
    stop("Original quantification contains unknown sample columns", call. = FALSE)
  measure_columns <- grep(paste0(" Peak ", source_dataset$intensity_measure, "$"), names(quant), value = TRUE)
  measure_names <- sub(paste0(" Peak ", source_dataset$intensity_measure, "$"), "", measure_columns)
  if (anyDuplicated(measure_names) || !setequal(measure_names, filenames))
    stop("Original metadata and MZmine filenames do not match exactly", call. = FALSE)
  parent <- file.path(batch, "recipes")
  final <- file.path(parent, plan$slug)
  final_link <- Sys.readlink(final)
  if (file.exists(final) || dir.exists(final) || (!is.na(final_link) && nzchar(final_link)))
    stop("Refusing to overwrite existing materialization: ", final, call. = FALSE)
  ancestors <- parent
  links <- Sys.readlink(ancestors)
  if (any(!is.na(links) & nzchar(links)))
    stop("Refusing symlinked materialization directory", call. = FALSE)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) stop("Cannot create plans directory", call. = FALSE)
  if (nzchar(Sys.readlink(parent)) ||
      !startsWith(normalizePath(parent, mustWork = TRUE), paste0(batch, "/")))
    stop("Materialization directory must be inside the batch", call. = FALSE)
  staging <- tempfile(pattern = ".materialize-", tmpdir = parent)
  if (!dir.create(staging)) stop("Cannot create temporary materialization", call. = FALSE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
  generated_dataset <- list(id = paste(canonical_id, plan$slug, sep = "/"), batch_dir = "../..",
                            metadata = substring(paths$metadata, nchar(batch) + 2L),
                            quantification = substring(paths$quantification, nchar(batch) + 2L))
  for (field in c("metadata_columns", "feature_id_column", "intensity_measure"))
    if (!is.null(source_dataset[[field]])) generated_dataset[[field]] <- source_dataset[[field]]
  if (length(paths$annotations)) generated_dataset$annotations <- lapply(paths$annotations, function(path)
    substring(path, nchar(batch) + 2L))
  generated_dataset$input_sha256 <- inputs[c("metadata", "quantification")]
  # The metadata alias is written once, in metadata_columns; roles refer to its normalized name.
  if (identical(source_roles$column, mapping$sample_type))
    source_roles$column <- "sample_type"
  if (!identical(source_roles, default_v2_roles()))
    generated_dataset$roles <- source_roles
  recipe <- list(seed = plan$seed,
                 sample_plan = list(group_column = plan$group_column, groups = plan$groups),
                 preprocessing = plan$preprocessing, design = list(group = "analysis_group", contrasts = "all"),
                 analyses = plan$analyses)
  if ("description" %in% names(plan)) recipe$presentation <- list(description = plan$description)
  write_yaml_file(generated_dataset, file.path(staging, "dataset.yaml"))
  write_yaml_file(recipe, file.path(staging, "recipe.yaml"))
  staged <- resolve_v2_config(file.path(staging, "dataset.yaml"), file.path(staging, "recipe.yaml"), repo_root)
  report <- validate_dataset(read_mzmine_dataset(staged$dataset, staged$recipe$sample_plan), staged$recipe)
  if (!isTRUE(report$valid))
    stop("Dataset validation failed: ", paste(report$issues, collapse = " "), call. = FALSE)
  atomic_publish(staging, final)
  cat("Materialized ", final, "\n", sep = "")
  for (command in c("validate", "plan", "run")) {
    cat("rig run -r 4.6.1 -f ", shQuote(file.path(repo_root, "v2", "run.R")),
        " -- ", command, " --dataset ", shQuote(file.path(final, "dataset.yaml")),
        " --recipe ", shQuote(file.path(final, "recipe.yaml")), "\n", sep = "")
  }
  invisible(final)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && identical(args[[1L]], "--")) args <- args[-1L]
  if (length(args) != 4L || !identical(args[c(1L, 3L)], c("--batch", "--plan")))
    stop("Usage: Rscript scripts/materialize_stats_recipe.R --batch <batch directory> --plan <downloaded plan.json>", call. = FALSE)
  materialize_stats_recipe(args[[2L]], args[[4L]])
}
