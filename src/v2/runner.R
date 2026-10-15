source_v2_modules <- function(repo_root) {
  module_dir <- file.path(repo_root, "src", "v2")
  for (module in c("utils.R", "config.R", "dataset.R", "preprocess.R", "analysis.R", "supervised.R", "annotations.R", "export.R"))
    source(file.path(module_dir, module), local = .GlobalEnv)
}

build_run_identity <- function(config) {
  annotation_paths <- Filter(Negate(is.null), config$dataset$annotations %||% list())
  checksums <- list(
    metadata = sha256_file(config$dataset$metadata),
    quantification = sha256_file(config$dataset$quantification),
    annotations = lapply(annotation_paths, sha256_file)
  )
  verify_input_sha256(config$dataset, checksums)
  # Input locations and expected-checksum guards are not analysis choices; actual file content is.
  mapping <- config$dataset[setdiff(names(config$dataset),
                                    c("batch_dir", "metadata", "quantification", "annotations", "input_sha256"))]
  mapping$feature_id_column <- mapping$feature_id_column %||% "row ID"
  mapping$roles <- dataset_roles(config$dataset)
  mapping$annotation_sources <- names(annotation_paths)
  lock_path <- file.path(config$repo_root, "v2", "renv.lock")
  if (!file.exists(lock_path)) stop("Pinned V2 environment not found: ", lock_path, call. = FALSE)
  environment <- list(
    schema_version = mapp_v2_schema_version,
    code_sha256 = v2_code_checksums(config$repo_root),
    renv_lock_sha256 = sha256_file(lock_path),
    packages = package_versions(c("ropls", "yaml", "digest", "ggplot2", "jsonlite")),
    r_version = as.character(getRversion())
  )
  dataset_hash <- sha256_object(list(schema = mapp_v2_schema_version, mapping = mapping, inputs = checksums))
  scientific <- scientific_recipe(config$recipe)
  preprocess_hash <- sha256_object(list(dataset_hash = dataset_hash, preprocessing = config$recipe$preprocessing,
                                        roles = mapping$roles, sample_plan = scientific$sample_plan,
                                        environment = environment))
  run_hash <- sha256_object(list(preprocess_hash = preprocess_hash, analysis = scientific$analyses,
                                   design = scientific$design, seed = scientific$seed %||% 1, environment = environment))
  list(dataset_hash = dataset_hash, preprocess_hash = preprocess_hash, run_hash = run_hash,
       inputs = checksums, environment = environment)
}

validate_v2_run <- function(config) {
  validate_recipe_config(config$recipe)
  if (isTRUE(config$recipe$analyses$plsda$enabled)) assert_packages("ropls")
  dataset <- read_mzmine_dataset(config$dataset, config$recipe$sample_plan)
  report <- validate_dataset(dataset, config$recipe)
  list(dataset = dataset, validation = report, identity = build_run_identity(config))
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

execute_v2_run <- function(config, force = FALSE) {
  started <- proc.time()[["elapsed"]]
  elapsed <- function() {
    seconds <- floor(proc.time()[["elapsed"]] - started)
    sprintf("%dm%02ds", seconds %/% 60, seconds %% 60)
  }
  say <- function(stage) message("[V2 ", elapsed(), "] ", stage)
  say("Validating metadata, quantification, annotations, recipe and input checksums")
  validation <- validate_v2_run(config)
  identity <- validation$identity
  output_root <- config$paths$output_root
  final_dir <- file.path(output_root, identity$run_hash)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
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
  staging_dir <- tempfile(pattern = paste0(".", substr(identity$run_hash, 1, 12), ".staging-"), tmpdir = output_root)
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
      line <- sprintf("[V2 %s] PLS-DA [%-20s] %3d%% (%.0f/%.0f fits) %s",
                      elapsed(), paste0(strrep("=", filled), strrep(" ", 20L - filled)),
                      percent, completed, total, phase)
      cat("\r", line, strrep(" ", max(0L, last_width - nchar(line))),
          if (completed == total) "\n" else "", sep = "", file = stderr())
      last_width <<- nchar(line)
      open_line <<- completed < total
      flush(stderr())
    } else if (completed == 0 || completed == total || percent >= next_percent || now - last_progress >= 10) {
      message("[V2 ", elapsed(), "] PLS-DA ", completed, "/", total,
              " fits (", percent, "%) — ", phase)
      next_percent <<- (percent %/% 5L + 1L) * 5L
    } else return(invisible(NULL))
    last_progress <<- now
    invisible(NULL)
  }
  say(paste0("Starting new run ", substr(identity$run_hash, 1L, 12L)))
  set.seed(as.integer(config$recipe$seed %||% 1))
  say(sprintf("Preprocessing %d injections and %d measured features",
              nrow(validation$dataset$raw), ncol(validation$dataset$raw)))
  processed <- preprocess_dataset(validation$dataset, config$recipe)
  say(sprintf("Retained %d features across %d biological samples",
              ncol(processed$matrices$scaled), nrow(processed$matrices$scaled)))
  analyses <- run_analyses(processed, config$recipe, function(stage) say(paste("Running", stage)))
  if (isTRUE(config$recipe$analyses$plsda$enabled)) say("Running PLS-DA with fold-local preprocessing and label permutations")
  supervised <- run_supervised(validation$dataset, config$recipe, progress = report)
  if (!is.null(supervised)) {
    analyses$plsda <- supervised
    if (identical(supervised$status, "withheld")) say(paste("PLS-DA withheld:", supervised$reason))
  }
  manifest <- list(
    schema_version = mapp_v2_schema_version,
    status = "complete",
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    dataset_id = config$dataset$id,
    hashes = identity[c("dataset_hash", "preprocess_hash", "run_hash")],
    inputs = identity$inputs,
    environment = identity$environment,
    validation = validation$validation,
    effective_dataset = config$dataset,
    effective_recipe = scientific_recipe(config$recipe)
  )
  if (!is.null(config$recipe$presentation))
    manifest$description <- config$recipe$presentation$description
  say("Writing saved tables, objects, plots and output checksums")
  export_results(validation$dataset, processed, analyses, manifest, staging_dir)
  say("Sealing and publishing the completed run")
  write_yaml_file(list(run_hash = identity$run_hash, content_sha256 = sha256_object(run_contents(staging_dir)),
                       completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(staging_dir, "COMPLETE"))
  atomic_publish(staging_dir, final_dir)
  success <- TRUE
  message("Completed V2 run: ", final_dir)
  invisible(final_dir)
}
