source_biostats_modules <- function(repo_root) {
  module_dir <- file.path(repo_root, "src", "biostats")
  for (module in c("utils.R", "config.R", "dataset.R", "preprocess.R", "analysis.R", "supervised.R", "annotations.R", "export.R"))
    source(file.path(module_dir, module), local = .GlobalEnv)
}

study_environment <- function(config) {
  lock_path <- file.path(config$repo_root, "biostats", "renv.lock")
  if (!file.exists(lock_path)) stop("Pinned Biostats environment not found: ", lock_path, call. = FALSE)
  list(
    schema_version = mapp_biostats_schema_version,
    code_sha256 = biostats_code_checksums(config$repo_root),
    renv_lock_sha256 = sha256_file(lock_path),
    packages = package_versions(c("ropls", "yaml", "digest", "ggplot2", "jsonlite")),
    r_version = as.character(getRversion())
  )
}

sealed_study <- function(config, dataset_hash, environment) {
  hash <- config$recipe$study_hash
  directory <- file.path(config$paths$output_root, hash)
  if (!dir.exists(directory) || nzchar(Sys.readlink(directory)) ||
      !completed_run(directory, hash))
    stop("No sealed study for comparison: ", hash, call. = FALSE)
  manifest_path <- file.path(directory, "manifest.yaml")
  objects <- file.path(directory, "objects", c("dataset.rds", "processed.rds"))
  if (any(!file.exists(c(manifest_path, objects))) ||
      any(nzchar(Sys.readlink(c(manifest_path, objects)))))
    stop("Sealed study lacks safe saved inputs: ", hash, call. = FALSE)
  manifest <- read_yaml_file(manifest_path)
  if (!identical(manifest$kind, "study") ||
      !identical(manifest$effective_recipe$kind, "study") ||
      !identical(manifest$hashes$run_hash, hash) ||
      !identical(manifest$hashes$dataset_hash, dataset_hash) ||
      !is.character(manifest$hashes$preprocess_hash) ||
      length(manifest$hashes$preprocess_hash) != 1L ||
      is.na(manifest$hashes$preprocess_hash) ||
      !grepl("^[0-9a-f]{64}$", manifest$hashes$preprocess_hash) ||
      !identical(manifest$environment[names(manifest$environment) != "git_commit"], environment))
    stop("Comparison dataset or code/environment differs from its sealed study; create a new study.", call. = FALSE)
  list(directory = directory, manifest = manifest)
}

build_run_identity <- function(config) {
  annotation_paths <- Filter(Negate(is.null), config$dataset$annotations %||% list())
  checksums <- list(
    metadata = sha256_file(config$dataset$metadata),
    quantification = sha256_file(config$dataset$quantification),
    annotations = lapply(annotation_paths, sha256_file)
  )
  verify_input_sha256(config$dataset, checksums)
  # Dataset input locations and checksum guards do not affect scientific identity.
  mapping <- config$dataset[setdiff(names(config$dataset),
                                    c("batch_dir", "metadata", "quantification", "annotations", "input_sha256"))]
  mapping$feature_id_column <- mapping$feature_id_column %||% "row ID"
  mapping$roles <- dataset_roles(config$dataset)
  mapping$annotation_sources <- names(annotation_paths)
  environment <- study_environment(config)
  dataset_hash <- sha256_object(list(schema = mapp_biostats_schema_version, mapping = mapping, inputs = checksums))
  scientific <- scientific_recipe(config$recipe)
  if (identical(scientific$kind, "comparison")) {
    parent <- sealed_study(config, dataset_hash, environment)
    preprocess_hash <- parent$manifest$hashes$preprocess_hash
    run_hash <- sha256_object(list(study_hash = scientific$study_hash,
      sample_plan = scientific$sample_plan, design = scientific$design,
      analyses = scientific$analyses, seed = scientific$seed %||% 1,
      environment = environment))
    return(list(dataset_hash = dataset_hash, preprocess_hash = preprocess_hash,
      run_hash = run_hash, inputs = checksums, environment = environment, study = parent))
  }
  preprocess_hash <- sha256_object(list(dataset_hash = dataset_hash, preprocessing = scientific$preprocessing,
                                        roles = mapping$roles, environment = environment))
  run_hash <- sha256_object(list(preprocess_hash = preprocess_hash, kind = "study",
                                analysis = scientific$analyses, display_group = scientific$design$group,
                                seed = scientific$seed %||% 1, environment = environment))
  list(dataset_hash = dataset_hash, preprocess_hash = preprocess_hash, run_hash = run_hash,
       inputs = checksums, environment = environment)
}

validate_biostats_run <- function(config) {
  validate_recipe_config(config$recipe)
  if (isTRUE(config$recipe$analyses$plsda$enabled) ||
      isTRUE(config$recipe$analyses$pca$enabled)) assert_packages("ropls")
  identity <- build_run_identity(config)
  recipe <- config$recipe
  if (identical(recipe$kind, "comparison")) {
    parent <- identity$study
    dataset <- readRDS(file.path(parent$directory, "objects", "dataset.rds"))
    if (!is.null(recipe$sample_plan))
      dataset <- subset_biostats_dataset(dataset, recipe$sample_plan)
    recipe$preprocessing <- parent$manifest$effective_recipe$preprocessing
    validate_recipe_config(recipe)
  } else {
    dataset <- read_mzmine_dataset(config$dataset)
  }
  report <- validate_dataset(dataset, recipe)
  list(dataset = dataset, validation = report, identity = identity, effective_recipe = recipe)
}

run_contents <- function(directory) {
  files <- sort(list.files(directory, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE))
  files <- files[!dir.exists(files) & files != file.path(directory, "COMPLETE")]
  names(files) <- substring(files, nchar(directory) + 2L)
  as.list(vapply(files, sha256_file, character(1)))
}

completed_run <- function(directory, run_hash) {
  marker_path <- file.path(directory, "COMPLETE")
  if (!file.exists(marker_path)) return(FALSE)
  marker <- tryCatch(read_yaml_file(marker_path), error = function(e) NULL)
  if (!is.list(marker) || !identical(marker$run_hash, run_hash) ||
      !identical(marker$content_sha256, sha256_object(run_contents(directory))))
    stop("Completed run is corrupt or has an incompatible COMPLETE marker: ", directory, call. = FALSE)
  TRUE
}

execute_biostats_run <- function(config, force = FALSE) {
  commit <- committed_biostats_code(config$repo_root)
  started <- proc.time()[["elapsed"]]
  elapsed <- function() {
    seconds <- floor(proc.time()[["elapsed"]] - started)
    sprintf("%dm%02ds", seconds %/% 60, seconds %% 60)
  }
  say <- function(stage) message("[Biostats ", elapsed(), "] ", stage)
  say("Validating metadata, quantification, annotations, recipe and input checksums")
  validation <- validate_biostats_run(config)
  identity <- validation$identity
  output_root <- config$paths$output_root
  kind <- config$recipe$kind
  final_dir <- if (identical(kind, "comparison"))
    file.path(output_root, "comparisons", config$recipe$study_hash, identity$run_hash) else
    file.path(output_root, identity$run_hash)
  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(final_dir, ".lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) stop("Another process holds this run lock: ", lock_dir, call. = FALSE)
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  write_yaml_file(list(run_hash = identity$run_hash, started_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                       pid = Sys.getpid()), file.path(lock_dir, "owner.yaml"))
  if (file.exists(final_dir)) {
    if (completed_run(final_dir, identity$run_hash)) {
      if (isTRUE(force)) stop("Refusing to overwrite a completed immutable run: ", final_dir, call. = FALSE)
      message("Completed run already exists: ", final_dir)
      return(invisible(final_dir))
    }
    stop("An incomplete run directory already exists; inspect or move it before rerunning: ", final_dir, call. = FALSE)
  }
  staging_dir <- tempfile(pattern = paste0(".", substr(identity$run_hash, 1, 12), ".staging-"),
                          tmpdir = dirname(final_dir))
  if (!dir.create(staging_dir)) stop("Could not create staging directory: ", staging_dir, call. = FALSE)
  success <- FALSE
  on.exit(if (!success && dir.exists(staging_dir)) message("Failed staging directory retained for diagnosis: ", staging_dir), add = TRUE)

  terminal <- isatty(stderr())
  last_progress <- -Inf
  next_percent <- 5L
  last_width <- 0L
  open_line <- FALSE
  on.exit(if (open_line) cat("\n", file = stderr()), add = TRUE)
  report <- function(phase, completed, total) {
    now <- proc.time()[["elapsed"]]
    percent <- floor(100 * completed / total)
    if (terminal) {
      if (completed > 0 && completed < total && now - last_progress < 0.25) return(invisible(NULL))
      filled <- floor(20 * completed / total)
      line <- sprintf("[Biostats %s] PLS-DA [%-20s] %3d%% (%.0f/%.0f fits) %s",
                      elapsed(), paste0(strrep("=", filled), strrep(" ", 20L - filled)),
                      percent, completed, total, phase)
      cat("\r", line, strrep(" ", max(0L, last_width - nchar(line))),
          if (completed == total) "\n" else "", sep = "", file = stderr())
      last_width <<- nchar(line)
      open_line <<- completed < total
      flush(stderr())
    } else if (completed == 0 || completed == total || percent >= next_percent || now - last_progress >= 10) {
      message("[Biostats ", elapsed(), "] PLS-DA ", completed, "/", total,
              " fits (", percent, "%) — ", phase)
      next_percent <<- (percent %/% 5L + 1L) * 5L
    } else return(invisible(NULL))
    last_progress <<- now
    invisible(NULL)
  }
  say(paste0("Starting new ", kind, " ", substr(identity$run_hash, 1L, 12L)))
  set.seed(as.integer(config$recipe$seed %||% 1))
  recipe <- validation$effective_recipe
  if (identical(kind, "study")) {
    say(sprintf("Preprocessing %d injections and %d measured features",
                nrow(validation$dataset$raw), ncol(validation$dataset$raw)))
    processed <- preprocess_dataset(validation$dataset, recipe)
    say(sprintf("Retained %d features across %d biological samples",
                ncol(processed$matrices$scaled), nrow(processed$matrices$scaled)))
  } else {
    say(paste0("Reusing sealed study ", substr(config$recipe$study_hash, 1L, 12L),
               " and its common preprocessing and ordinations"))
    study_processed <- readRDS(file.path(identity$study$directory, "objects", "processed.rds"))
    samples <- rownames(validation$dataset$sample_metadata)[
      role_rows(validation$dataset$sample_metadata, recipe$roles, "sample")]
    transformed <- study_processed$matrices$transformed[samples, , drop = FALSE]
    processed <- list(matrices = list(transformed = transformed),
                      sample_metadata = validation$dataset$sample_metadata[samples, , drop = FALSE],
                      methods = study_processed$methods)
  }
  analyses <- run_analyses(processed, recipe, function(stage) say(paste("Running", stage)))
  if (isTRUE(recipe$analyses$plsda$enabled))
    say("Running PLS-DA with fold-local preprocessing and label permutations")
  supervised <- run_supervised(validation$dataset, recipe, progress = report)
  if (!is.null(supervised)) {
    analyses$plsda <- supervised
    if (identical(supervised$status, "withheld")) say(paste("PLS-DA withheld:", supervised$reason))
  }
  manifest <- list(
    schema_version = mapp_biostats_schema_version,
    status = "complete",
    kind = kind,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    dataset_id = config$dataset$id,
    hashes = identity[c("dataset_hash", "preprocess_hash", "run_hash")],
    inputs = identity$inputs,
    environment = c(identity$environment, list(git_commit = commit)),
    validation = validation$validation,
    effective_dataset = if (identical(kind, "comparison"))
      identity$study$manifest$effective_dataset else config$dataset,
    effective_recipe = scientific_recipe(recipe)
  )
  if (identical(kind, "comparison")) manifest$study_hash <- config$recipe$study_hash
  if (!is.null(config$recipe$presentation))
    manifest$description <- config$recipe$presentation$description
  say("Writing saved tables, objects, plots and output checksums")
  if (identical(kind, "comparison"))
    export_comparison_results(validation$dataset, analyses, manifest, staging_dir)
  else
    export_results(validation$dataset, processed, analyses, manifest, staging_dir)
  committed_biostats_code(config$repo_root, identity$environment$code_sha256, commit)
  say("Sealing and publishing the completed run")
  write_yaml_file(list(run_hash = identity$run_hash, content_sha256 = sha256_object(run_contents(staging_dir)),
                       completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(staging_dir, "COMPLETE"))
  atomic_publish(staging_dir, final_dir)
  success <- TRUE
  message("Completed Biostats run: ", final_dir)
  invisible(final_dir)
}
