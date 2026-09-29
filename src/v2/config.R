config_fields <- function(value, allowed, context) {
  if (!is.list(value) || (length(value) && (is.null(names(value)) || anyNA(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value)))))
    stop(context, " must be a mapping with unique keys.", call. = FALSE)
  unknown <- setdiff(names(value), allowed)
  if (length(unknown)) stop("Unsupported ", context, " key(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

config_text <- function(value, context) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value)))
    stop(context, " must be a non-empty string.", call. = FALSE)
  invisible(value)
}

config_choice <- function(value, allowed, context) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !value %in% allowed)
    stop("Unsupported ", context, "; expected one of: ", paste(allowed, collapse = ", "), call. = FALSE)
}

config_boolean <- function(value, context) {
  if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop(context, " must be true or false.", call. = FALSE)
}

config_number <- function(value, context, minimum = 0, maximum = Inf, integer = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < minimum || value > maximum ||
      (integer && value != floor(value)))
    stop(context, " must be ", if (integer) "an integer" else "a number", " in [", minimum, ", ", maximum, "].", call. = FALSE)
}

validate_recipe_config <- function(recipe) {
  config_fields(recipe, c("seed", "output_root", "presentation", "roles", "preprocessing", "design", "analyses"), "recipe")
  if (!is.null(recipe$seed)) config_number(recipe$seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  if (!is.null(recipe$output_root)) config_text(recipe$output_root, "output_root")

  roles <- recipe$roles %||% list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank"))
  config_fields(roles, c("column", "levels"), "roles")
  config_text(roles$column, "roles.column")
  config_fields(roles$levels, c("sample", "qc", "blank"), "roles.levels")
  for (role in c("sample", "qc", "blank")) config_text(roles$levels[[role]], paste0("roles.levels.", role))
  if (anyDuplicated(unlist(roles$levels, use.names = FALSE))) stop("Role levels must be distinct.", call. = FALSE)

  preprocessing <- recipe$preprocessing %||% list()
  config_fields(preprocessing, c("blank_filter", "qc_rsd_filter", "missing_values", "normalization", "transformation", "scaling"), "preprocessing")
  sections <- list(blank_filter = c("enabled", "minimum_sample_to_blank_ratio", "minimum_blank_prevalence"),
                   qc_rsd_filter = c("enabled", "threshold_percent"), missing_values = "method",
                   normalization = "method", transformation = "method", scaling = "method")
  for (section in names(sections)) config_fields(preprocessing[[section]] %||% list(), sections[[section]], paste0("preprocessing.", section))
  for (section in c("blank_filter", "qc_rsd_filter")) {
    enabled <- preprocessing[[section]]$enabled
    if (!is.null(enabled)) config_boolean(enabled, paste0("preprocessing.", section, ".enabled"))
  }
  blank <- preprocessing$blank_filter
  qc <- preprocessing$qc_rsd_filter
  if (!is.null(blank$minimum_sample_to_blank_ratio)) config_number(blank$minimum_sample_to_blank_ratio, "blank filter ratio", 0)
  if (!is.null(blank$minimum_blank_prevalence)) config_number(blank$minimum_blank_prevalence, "blank filter prevalence", 0, 1)
  if (!is.null(qc$threshold_percent)) config_number(qc$threshold_percent, "QC RSD threshold", 0)
  methods <- list(missing_values = c("none", "half_minimum"), normalization = c("none", "total_sum", "median", "pqn"),
                  transformation = c("none", "log2", "log10", "log1p"), scaling = c("none", "pareto", "autoscale"))
  for (section in names(methods)) {
    method <- preprocessing[[section]]$method
    if (!is.null(method)) config_choice(tolower(method), methods[[section]], paste0("preprocessing.", section, ".method"))
  }

  design <- recipe$design
  config_fields(design, c("group", "block", "block_verified", "contrasts"), "design")
  config_text(design$group, "design.group")
  if (!is.null(design$block)) {
    config_text(design$block, "design.block")
    if (identical(design$block, design$group)) stop("The group and block columns must differ.", call. = FALSE)
    if (!identical(design$block_verified, TRUE)) stop("design.block requires explicit block_verified: true after experimental verification.", call. = FALSE)
  } else if (!is.null(design$block_verified)) {
    stop("design.block_verified requires a design.block column.", call. = FALSE)
  }
  contrasts <- design$contrasts %||% list()
  if (!is.list(contrasts) || (length(contrasts) && !is.null(names(contrasts)))) stop("design.contrasts must be a sequence.", call. = FALSE)
  contrast_names <- character()
  for (contrast in contrasts) {
    config_fields(contrast, c("name", "numerator", "denominator"), "contrast")
    config_text(contrast$numerator, "contrast.numerator")
    config_text(contrast$denominator, "contrast.denominator")
    if (identical(contrast$numerator, contrast$denominator)) stop("Contrast levels must differ.", call. = FALSE)
    name <- contrast$name %||% paste(contrast$numerator, "vs", contrast$denominator, sep = "_")
    config_text(name, "contrast.name")
    contrast_names <- c(contrast_names, name)
  }
  if (anyDuplicated(contrast_names) || anyDuplicated(vapply(contrast_names, safe_name, character(1))))
    stop("Contrast names (including filesystem-safe names) must be unique.", call. = FALSE)

  analyses <- recipe$analyses %||% list()
  config_fields(analyses, c("pca", "pcoa", "omnibus", "differential", "plsda"), "analyses")
  options <- list(pca = c("enabled", "components"), pcoa = c("enabled", "components", "stage", "distance"),
                  omnibus = "enabled", differential = "enabled",
                  plsda = c("enabled", "components", "folds", "repeats", "permutations"))
  for (analysis in names(options)) {
    value <- analyses[[analysis]] %||% list()
    config_fields(value, options[[analysis]], paste0("analyses.", analysis))
    if (!is.null(value$enabled)) config_boolean(value$enabled, paste0("analyses.", analysis, ".enabled"))
    if (!is.null(value$components)) config_number(value$components, paste0("analyses.", analysis, ".components"), 1, integer = TRUE)
  }
  pcoa <- analyses$pcoa
  if (!is.null(pcoa$stage)) config_choice(pcoa$stage, c("imputed", "normalized", "transformed", "scaled"), "PCoA stage")
  if (!is.null(pcoa$distance)) config_choice(pcoa$distance, c("bray", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), "PCoA distance")
  if (isTRUE(pcoa$enabled) && identical(pcoa$distance %||% "bray", "bray") &&
      !(pcoa$stage %||% "normalized") %in% c("imputed", "normalized"))
    stop("Bray-Curtis PCoA requires the imputed or normalized nonnegative stage.", call. = FALSE)
  if (isTRUE(analyses$differential$enabled) && !length(contrasts))
    stop("Differential analysis requires explicit planned design.contrasts.", call. = FALSE)
  plsda <- analyses$plsda
  if (!is.null(plsda$folds)) config_number(plsda$folds, "PLS-DA folds", 2, integer = TRUE)
  if (!is.null(plsda$repeats)) config_number(plsda$repeats, "PLS-DA repeats", 2, integer = TRUE)
  if (!is.null(plsda$permutations)) config_number(plsda$permutations, "PLS-DA permutations", 19, integer = TRUE)
  if (isTRUE(plsda$enabled) && !is.null(design$block))
    stop("PLS-DA does not model verified blocks; disable PLS-DA or remove the block design.", call. = FALSE)
  invisible(recipe)
}

resolve_v2_config <- function(dataset_yaml, recipe_yaml, repo_root) {
  dataset_path <- normalize_existing_path(dataset_yaml)
  recipe_path <- normalize_existing_path(recipe_yaml)
  dataset <- read_yaml_file(dataset_path)
  recipe <- read_yaml_file(recipe_path)
  config_fields(dataset, c("id", "batch_dir", "metadata", "quantification", "feature_id_column", "intensity_measure", "annotations"), "dataset")
  for (field in c("id", "batch_dir", "metadata", "quantification")) config_text(dataset[[field]], paste0("dataset.", field))
  if (!is.null(dataset$feature_id_column)) config_text(dataset$feature_id_column, "dataset.feature_id_column")
  if (!is.null(dataset$intensity_measure)) config_choice(dataset$intensity_measure, c("height", "area"), "dataset.intensity_measure")
  annotations <- dataset$annotations %||% list()
  if (!is.list(annotations) || (length(annotations) && (is.null(names(annotations)) || anyNA(names(annotations)) ||
      any(!nzchar(names(annotations))) || anyDuplicated(names(annotations)))))
    stop("dataset.annotations must map unique source names to paths.", call. = FALSE)
  for (name in names(annotations)) if (!is.null(annotations[[name]])) config_text(annotations[[name]], paste0("dataset.annotations.", name))
  validate_recipe_config(recipe)

  dataset_base <- dirname(dataset_path)
  batch_dir <- normalize_existing_path(dataset$batch_dir, dataset_base)
  resolve_input <- function(value) normalize_existing_path(value, batch_dir)
  dataset$batch_dir <- batch_dir
  dataset$metadata <- resolve_input(dataset$metadata)
  dataset$quantification <- resolve_input(dataset$quantification)
  if (length(annotations)) dataset$annotations <- lapply(annotations, function(value) if (is.null(value)) NULL else resolve_input(value))

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
