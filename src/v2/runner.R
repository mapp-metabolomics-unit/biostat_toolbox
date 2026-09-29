source_v2_modules <- function(repo_root) {
  module_dir <- file.path(repo_root, "src", "v2")
  for (module in c("utils.R", "config.R", "dataset.R", "preprocess.R", "analysis.R", "export.R")) source(file.path(module_dir, module), local = .GlobalEnv)
}

build_run_identity <- function(config) {
  inputs <- c(metadata = config$dataset$metadata, quantification = config$dataset$quantification)
  annotation_paths <- unlist(config$dataset$annotations %||% list(), use.names = TRUE)
  annotation_paths <- annotation_paths[!vapply(annotation_paths, is.null, logical(1))]
  inputs <- c(inputs, annotation_paths)
  checksums <- as.list(vapply(inputs, sha256_file, character(1)))
  v2_lock <- file.path(config$repo_root, "v2", "renv.lock")
  lock_path <- if (file.exists(v2_lock)) v2_lock else file.path(config$repo_root, "renv.lock")
  environment <- list(
    schema_version = mapp_v2_schema_version,
    git = git_fingerprint(config$repo_root),
    code_sha256 = v2_code_checksums(config$repo_root),
    renv_lock = sub(paste0("^", normalizePath(config$repo_root), "/?"), "", normalizePath(lock_path, mustWork = FALSE)),
    renv_lock_sha256 = if (file.exists(lock_path)) sha256_file(lock_path) else NA_character_,
    packages = package_versions(c("structToolbox", "struct", "yaml", "digest")),
    r_version = as.character(getRversion())
  )
  dataset_hash <- sha256_object(list(schema = mapp_v2_schema_version, dataset_id = config$dataset$id, inputs = checksums))
  preprocess_hash <- sha256_object(list(dataset_hash = dataset_hash, preprocessing = config$recipe$preprocessing, roles = config$recipe$roles, environment = environment))
  run_hash <- sha256_object(list(preprocess_hash = preprocess_hash, analysis = scientific_recipe(config$recipe)$analyses, design = config$recipe$design, seed = config$recipe$seed %||% 1, environment = environment))
  list(dataset_hash = dataset_hash, preprocess_hash = preprocess_hash, run_hash = run_hash, inputs = checksums, environment = environment)
}

validate_v2_run <- function(config) {
  dataset <- read_mzmine_dataset(config$dataset)
  report <- validate_dataset(dataset, config$recipe)
  list(dataset = dataset, validation = report, identity = build_run_identity(config))
}

execute_v2_run <- function(config, force = FALSE) {
  validation <- validate_v2_run(config)
  identity <- validation$identity
  output_root <- config$paths$output_root
  final_dir <- file.path(output_root, identity$run_hash)
  if (file.exists(file.path(final_dir, "COMPLETE")) && !force) {
    message("Completed run already exists: ", final_dir)
    return(invisible(final_dir))
  }
  if (dir.exists(final_dir)) stop("An incomplete run directory already exists; inspect or move it before rerunning: ", final_dir, call. = FALSE)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(final_dir, ".lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) stop("Another process holds this run lock: ", lock_dir, call. = FALSE)
  on.exit(if (dir.exists(lock_dir)) unlink(lock_dir, recursive = TRUE), add = TRUE)
  staging_dir <- tempfile(pattern = paste0(".", substr(identity$run_hash, 1, 12), ".staging-"), tmpdir = output_root)
  dir.create(staging_dir, recursive = TRUE)
  success <- FALSE
  on.exit(if (!success && dir.exists(staging_dir)) message("Failed staging directory retained for diagnosis: ", staging_dir), add = TRUE)

  set.seed(as.integer(config$recipe$seed %||% 1))
  processed <- preprocess_dataset(validation$dataset, config$recipe)
  analyses <- run_analyses(processed, config$recipe)
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
  export_results(validation$dataset, processed, analyses, manifest, staging_dir)
  file.create(file.path(staging_dir, "COMPLETE"))
  atomic_publish(staging_dir, final_dir)
  success <- TRUE
  message("Completed V2 run: ", final_dir)
  invisible(final_dir)
}
