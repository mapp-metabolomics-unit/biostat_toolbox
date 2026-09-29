read_mzmine_dataset <- function(dataset_config) {
  required_metadata <- c("filename", "sample_id", "sample_type")
  metadata_columns <- names(utils::read.delim(dataset_config$metadata, check.names = FALSE, nrows = 0))
  feature_column <- dataset_config$feature_id_column %||% "row ID"
  quant_columns <- names(utils::read.csv(dataset_config$quantification, check.names = FALSE, nrows = 0))
  missing_metadata <- setdiff(required_metadata, metadata_columns)
  if (length(missing_metadata)) stop("Metadata lacks required column(s): ", paste(missing_metadata, collapse = ", "), call. = FALSE)
  if (!feature_column %in% quant_columns) stop("Feature ID column not found: ", feature_column, call. = FALSE)
  if (anyDuplicated(metadata_columns) || anyDuplicated(quant_columns)) stop("Input column names must be unique.", call. = FALSE)
  metadata <- utils::read.delim(dataset_config$metadata, check.names = FALSE, stringsAsFactors = FALSE,
                                na.strings = c("", "NA"), colClasses = setNames(rep("character", 3L), required_metadata))
  quant <- utils::read.csv(dataset_config$quantification, check.names = FALSE, stringsAsFactors = FALSE,
                           na.strings = c("", "NA"), colClasses = setNames("character", feature_column))
  for (field in required_metadata) {
    values <- metadata[[field]]
    if (anyNA(values) || any(!nzchar(trimws(as.character(values)))))
      stop("Metadata ", field, " values must be non-empty.", call. = FALSE)
  }
  if (anyDuplicated(metadata$filename) || anyDuplicated(metadata$sample_id))
    stop("Metadata filename and sample_id values must each be unique.", call. = FALSE)

  intensity_pattern <- " Peak (height|area)$"
  all_intensity_columns <- grep(intensity_pattern, names(quant), value = TRUE)
  if (!length(all_intensity_columns)) stop("No MZmine 'Peak height' or 'Peak area' columns were found.", call. = FALSE)
  available_measures <- unique(sub("^.* Peak (height|area)$", "\\1", all_intensity_columns))
  measure <- dataset_config$intensity_measure
  if (is.null(measure)) {
    if (length(available_measures) != 1L)
      stop("Both MZmine Peak height and Peak area are present; set dataset intensity_measure to 'height' or 'area'.", call. = FALSE)
    measure <- available_measures
  }
  if (!is.character(measure) || length(measure) != 1L || is.na(measure) || !measure %in% c("height", "area"))
    stop("intensity_measure must be 'height' or 'area'.", call. = FALSE)
  if (!measure %in% available_measures) stop("No MZmine Peak ", measure, " columns were found.", call. = FALSE)
  intensity_columns <- grep(paste0(" Peak ", measure, "$"), names(quant), value = TRUE)
  sample_names <- sub(paste0(" Peak ", measure, "$"), "", intensity_columns)
  if (any(!nzchar(trimws(sample_names))) || anyDuplicated(sample_names))
    stop("Quantification sample names must be non-empty and unique.", call. = FALSE)

  feature_ids <- as.character(quant[[feature_column]])
  if (anyNA(feature_ids) || any(!nzchar(trimws(feature_ids))) || anyDuplicated(feature_ids))
    stop("Feature IDs must be non-empty and unique.", call. = FALSE)
  if (!length(feature_ids)) stop("Quantification contains no features.", call. = FALSE)

  missing_metadata_rows <- setdiff(sample_names, metadata$filename)
  metadata_only_rows <- setdiff(metadata$filename, sample_names)
  if (length(missing_metadata_rows) || length(metadata_only_rows)) {
    stop("Quantification and metadata filenames must match exactly. Quantification only: ",
         paste(missing_metadata_rows, collapse = ", "), "; metadata only: ",
         paste(metadata_only_rows, collapse = ", "), call. = FALSE)
  }
  metadata <- metadata[match(sample_names, metadata$filename), , drop = FALSE]
  rownames(metadata) <- metadata$filename
  intensities <- lapply(intensity_columns, function(column) {
    values <- quant[[column]]
    if (!is.numeric(values) && !is.character(values))
      stop("Non-numeric intensity in column: ", column, call. = FALSE)
    parsed <- suppressWarnings(as.numeric(values))
    if (any(!is.na(values) & (is.na(parsed) | !is.finite(parsed))))
      stop("Non-numeric or non-finite intensity in column: ", column, call. = FALSE)
    parsed
  })
  x <- t(do.call(cbind, intensities))
  rownames(x) <- sample_names
  colnames(x) <- feature_ids
  storage.mode(x) <- "double"
  if (any(x < 0, na.rm = TRUE)) stop("Negative raw intensities are not supported.", call. = FALSE)

  variable_metadata <- quant[!names(quant) %in% all_intensity_columns]
  variable_metadata$feature_id <- feature_ids
  rownames(variable_metadata) <- feature_ids

  annotation_tables <- list()
  for (annotation_name in names(dataset_config$annotations %||% list())) {
    annotation_path <- dataset_config$annotations[[annotation_name]]
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
  roles <- recipe$roles %||% list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank"))
  role_column <- roles$column
  if (!role_column %in% names(sm)) stop("Sample-role column not found: ", role_column, call. = FALSE)
  counts <- vapply(roles$levels, function(level) sum(sm[[role_column]] == level, na.rm = TRUE), integer(1))
  issues <- character()
  if (!counts[["sample"]]) issues <- c(issues, "No biological samples were found.")
  if (isTRUE(recipe$preprocessing$blank_filter$enabled) && !counts[["blank"]]) issues <- c(issues, "Blank filtering is enabled but no blanks were found.")
  if (isTRUE(recipe$preprocessing$qc_rsd_filter$enabled) && counts[["qc"]] < 3) issues <- c(issues, "QC RSD filtering needs at least three QCs.")

  group_column <- recipe$design$group
  group_counts <- integer()
  if (!group_column %in% names(sm)) {
    issues <- c(issues, paste("Design group column not found:", group_column))
  } else if (counts[["sample"]]) {
    sample_rows <- sm[[role_column]] == roles$levels$sample
    groups <- as.character(sm[[group_column]][sample_rows])
    if (anyNA(groups) || any(!nzchar(trimws(groups)))) {
      issues <- c(issues, "Biological sample group levels must be non-empty and non-missing.")
    } else {
      group_counts <- table(groups)
      if (length(group_counts) < 2L) issues <- c(issues, "At least two biological sample group levels are required.")
      if (any(group_counts < 3L)) issues <- c(issues, "At least one group has fewer than three samples.")
      if (isTRUE(recipe$analyses$plsda$enabled)) {
        if (any(group_counts < 4L)) issues <- c(issues, "PLS-DA needs at least four samples per group.")
        if (any(group_counts < (recipe$analyses$plsda$folds %||% 3L)))
          issues <- c(issues, "PLS-DA folds exceed the sample count in at least one group.")
      }
      for (contrast in recipe$design$contrasts %||% list()) {
        if (!all(c(contrast$numerator, contrast$denominator) %in% names(group_counts)))
          issues <- c(issues, paste("Contrast levels are absent from biological samples:", contrast$numerator, "vs", contrast$denominator))
      }
      block_column <- recipe$design$block
      if (!is.null(block_column)) {
        if (!block_column %in% names(sm)) {
          issues <- c(issues, paste("Design block column not found:", block_column))
        } else {
          blocks <- as.character(sm[[block_column]][sample_rows])
          if (anyNA(blocks) || any(!nzchar(trimws(blocks)))) {
            issues <- c(issues, "Biological sample block levels must be non-empty and non-missing.")
          } else {
            cross <- table(blocks, groups)
            if (nrow(cross) < 2L || any(rowSums(cross > 0L) < 2L))
              issues <- c(issues, "Verified blocks require at least two blocks, each spanning at least two groups.")
            full <- stats::model.matrix(~ 0 + group + block, data.frame(group = factor(groups), block = factor(blocks)))
            if (qr(full)$rank < ncol(full) || nrow(full) <= ncol(full))
              issues <- c(issues, "The group-plus-block design must have full rank and positive residual degrees of freedom.")
          }
        }
      }
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
  if (!report$valid) stop("Dataset validation failed: ", paste(issues, collapse = " "), call. = FALSE)
  report
}
