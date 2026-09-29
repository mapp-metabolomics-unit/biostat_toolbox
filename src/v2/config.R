resolve_v2_config <- function(dataset_yaml, recipe_yaml, repo_root) {
  dataset_path <- normalize_existing_path(dataset_yaml)
  recipe_path <- normalize_existing_path(recipe_yaml)
  dataset <- read_yaml_file(dataset_path)
  recipe <- read_yaml_file(recipe_path)
  dataset_base <- dirname(dataset_path)

  required_dataset <- c("id", "batch_dir", "metadata", "quantification")
  missing_dataset <- setdiff(required_dataset, names(dataset))
  if (length(missing_dataset)) stop("Dataset YAML lacks: ", paste(missing_dataset, collapse = ", "), call. = FALSE)

  batch_dir <- normalize_existing_path(dataset$batch_dir, dataset_base)
  resolve_input <- function(value) normalize_existing_path(value, batch_dir)
  dataset$batch_dir <- batch_dir
  dataset$metadata <- resolve_input(dataset$metadata)
  dataset$quantification <- resolve_input(dataset$quantification)
  if (!is.null(dataset$annotations)) {
    dataset$annotations <- lapply(dataset$annotations, function(value) {
      if (is.null(value) || !nzchar(value)) return(NULL)
      resolve_input(value)
    })
  }

  output_root <- recipe$output_root %||% file.path(batch_dir, "results", "stats_v2")
  if (!grepl("^/", output_root)) output_root <- file.path(batch_dir, output_root)

  list(
    dataset = dataset,
    recipe = recipe,
    paths = list(dataset_yaml = dataset_path, recipe_yaml = recipe_path, output_root = normalizePath(output_root, mustWork = FALSE)),
    repo_root = repo_root
  )
}

scientific_recipe <- function(recipe) {
  recipe$output_root <- NULL
  recipe$presentation <- NULL
  recipe
}

