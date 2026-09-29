read_mzmine_dataset <- function(dataset_config) {
  metadata <- utils::read.delim(dataset_config$metadata, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  quant <- utils::read.csv(dataset_config$quantification, check.names = FALSE, stringsAsFactors = FALSE)

  required_metadata <- c("filename", "sample_id", "sample_type")
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata)) stop("Metadata lacks required column(s): ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  if (anyDuplicated(metadata$filename)) stop("Metadata filename values must be unique.", call. = FALSE)

  intensity_columns <- grep(" Peak height$", names(quant), value = TRUE)
  if (!length(intensity_columns)) stop("No MZmine 'Peak height' columns were found.", call. = FALSE)
  sample_names <- sub(" Peak height$", "", intensity_columns)
  if (anyDuplicated(sample_names)) stop("Quantification sample names must be unique.", call. = FALSE)

  feature_column <- dataset_config$feature_id_column %||% "row ID"
  if (!feature_column %in% names(quant)) stop("Feature ID column not found: ", feature_column, call. = FALSE)
  feature_ids <- as.character(quant[[feature_column]])
  if (anyNA(feature_ids) || any(!nzchar(feature_ids)) || anyDuplicated(feature_ids)) stop("Feature IDs must be non-empty and unique.", call. = FALSE)

  missing_metadata_rows <- setdiff(sample_names, metadata$filename)
  metadata_only_rows <- setdiff(metadata$filename, sample_names)
  if (length(missing_metadata_rows)) stop("Quantified samples missing from metadata: ", paste(missing_metadata_rows, collapse = ", "), call. = FALSE)
  if (length(metadata_only_rows)) message("Metadata-only samples will be excluded: ", paste(metadata_only_rows, collapse = ", "))

  metadata <- metadata[match(sample_names, metadata$filename), , drop = FALSE]
  rownames(metadata) <- metadata$filename
  x <- t(as.matrix(data.frame(lapply(quant[intensity_columns], as.numeric), check.names = FALSE)))
  rownames(x) <- sample_names
  colnames(x) <- feature_ids
  storage.mode(x) <- "double"
  if (any(x < 0, na.rm = TRUE)) stop("Negative raw intensities are not supported.", call. = FALSE)

  variable_indices <- which(!names(quant) %in% intensity_columns & !is.na(names(quant)) & nzchar(names(quant)))
  variable_metadata <- quant[variable_indices]
  variable_metadata$feature_id <- feature_ids
  rownames(variable_metadata) <- feature_ids

  annotation_tables <- list()
  annotation_config <- dataset_config$annotations %||% list()
  for (annotation_name in names(annotation_config)) {
    annotation_path <- annotation_config[[annotation_name]]
    if (is.null(annotation_path)) next
    annotation_tables[[annotation_name]] <- utils::read.delim(annotation_path, check.names = FALSE, stringsAsFactors = FALSE)
  }

  list(
    raw = x,
    sample_metadata = metadata,
    variable_metadata = variable_metadata,
    annotations = annotation_tables,
    source = dataset_config
  )
}

validate_dataset <- function(dataset, recipe = list()) {
  sm <- dataset$sample_metadata
  x <- dataset$raw
  role_column <- recipe$roles$column %||% "sample_type"
  if (!role_column %in% names(sm)) stop("Sample-role column not found: ", role_column, call. = FALSE)
  roles <- recipe$roles$levels %||% list(sample = "sample", qc = "QC", blank = "blank")
  missing_roles <- setdiff(c("sample", "qc", "blank"), names(roles))
  if (length(missing_roles)) stop("Role mapping lacks: ", paste(missing_roles, collapse = ", "), call. = FALSE)
  counts <- vapply(roles, function(level) sum(sm[[role_column]] == level, na.rm = TRUE), integer(1))
  issues <- list()
  if (!counts[["sample"]]) issues <- c(issues, "No biological samples were found.")
  if (isTRUE(recipe$preprocessing$blank_filter$enabled) && !counts[["blank"]]) issues <- c(issues, "Blank filtering is enabled but no blanks were found.")
  if (isTRUE(recipe$preprocessing$qc_rsd_filter$enabled) && counts[["qc"]] < 3) issues <- c(issues, "QC RSD filtering needs at least three QCs.")

  group_column <- recipe$design$group %||% NULL
  group_counts <- integer()
  if (!is.null(group_column)) {
    if (!group_column %in% names(sm)) issues <- c(issues, paste("Design group column not found:", group_column)) else {
      sample_rows <- sm[[role_column]] == roles$sample
      group_counts <- table(sm[[group_column]][sample_rows], useNA = "ifany")
      if (any(group_counts < 3)) issues <- c(issues, "At least one group has fewer than three samples.")
    }
  }

  report <- list(
    valid = !length(issues),
    dimensions = list(samples = nrow(x), features = ncol(x)),
    role_counts = as.list(counts),
    group_counts = as.list(group_counts),
    zero_fraction = mean(x == 0, na.rm = TRUE),
    missing_fraction = mean(is.na(x)),
    issues = unname(issues)
  )
  if (!report$valid) stop("Dataset validation failed: ", paste(unlist(issues), collapse = " "), call. = FALSE)
  report
}
