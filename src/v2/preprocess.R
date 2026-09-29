role_rows <- function(metadata, roles, names) {
  metadata[[roles$column]] %in% unlist(roles$levels[names], use.names = FALSE)
}

blank_filter_matrix <- function(x, metadata, roles, config) {
  sample_rows <- role_rows(metadata, roles, "sample")
  blank_rows <- role_rows(metadata, roles, "blank")
  if (!any(sample_rows) || !any(blank_rows)) stop("Blank filtering requires samples and blanks.", call. = FALSE)
  sample_median <- apply(x[sample_rows, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  blank_median <- apply(x[blank_rows, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  blank_prevalence <- colMeans(x[blank_rows, , drop = FALSE] > 0, na.rm = TRUE)
  ratio <- ifelse(blank_median > 0, sample_median / blank_median, Inf)
  if (any(!is.finite(sample_median) | !is.finite(blank_median) | !is.finite(blank_prevalence))) stop("Blank filtering needs finite feature summaries.", call. = FALSE)
  minimum_ratio <- as.numeric(config$minimum_sample_to_blank_ratio %||% 5)
  minimum_prevalence <- as.numeric(config$minimum_blank_prevalence %||% 0.5)
  contaminated <- blank_prevalence >= minimum_prevalence & ratio < minimum_ratio
  diagnostics <- data.frame(
    feature_id = colnames(x), sample_median = sample_median, blank_median = blank_median,
    blank_prevalence = blank_prevalence, sample_to_blank_ratio = ratio,
    removed_as_blank = contaminated, stringsAsFactors = FALSE
  )
  list(data = x[, !contaminated, drop = FALSE], diagnostics = diagnostics)
}

normalize_matrix <- function(x, method) {
  method <- tolower(method %||% "none")
  if (method == "none") return(x)
  if (method %in% c("total_sum", "median")) {
    size <- if (method == "total_sum") rowSums(x, na.rm = TRUE) else apply(x, 1, stats::median, na.rm = TRUE)
    target <- stats::median(size[size > 0], na.rm = TRUE)
    if (any(!is.finite(size) | size <= 0) || !is.finite(target) || target <= 0) stop("Normalization encountered a sample with zero or invalid size.", call. = FALSE)
    return(sweep(x, 1, size / target, "/"))
  }
  if (method == "pqn") {
    reference <- apply(x, 2, stats::median, na.rm = TRUE)
    quotients <- sweep(x, 2, reference, "/")
    quotients[!is.finite(quotients)] <- NA_real_
    factors <- apply(quotients, 1, stats::median, na.rm = TRUE)
    if (any(!is.finite(factors) | factors <= 0)) stop("PQN encountered an invalid dilution factor.", call. = FALSE)
    return(sweep(x, 1, factors, "/"))
  }
  stop("Unknown normalization method: ", method, call. = FALSE)
}

impute_matrix <- function(x, method) {
  method <- tolower(method %||% "half_minimum")
  if (any(is.infinite(x))) stop("Intensities must be finite or missing.", call. = FALSE)
  missing <- is.na(x) | x <= 0
  if (method == "none") return(x)
  if (method == "half_minimum") {
    positive <- x[is.finite(x) & x > 0]
    if (!length(positive)) stop("Cannot impute a matrix without positive values.", call. = FALSE)
    x[missing] <- min(positive) / 2
    return(x)
  }
  stop("Unknown imputation method: ", method, call. = FALSE)
}

transform_matrix <- function(x, method) {
  method <- tolower(method %||% "none")
  switch(method,
    none = x,
    log2 = log2(x),
    log10 = log10(x),
    log1p = log1p(x),
    stop("Unknown transformation method: ", method, call. = FALSE)
  )
}

scale_matrix <- function(x, method) {
  method <- tolower(method %||% "none")
  if (method == "none") return(x)
  means <- colMeans(x, na.rm = TRUE)
  sd_values <- apply(x, 2, stats::sd, na.rm = TRUE)
  denominator <- if (method == "pareto") sqrt(sd_values) else if (method == "autoscale") sd_values else stop("Unknown scaling method: ", method, call. = FALSE)
  denominator[!is.finite(denominator) | denominator == 0] <- 1
  sweep(sweep(x, 2, means, "-"), 2, denominator, "/")
}

preprocess_dataset <- function(dataset, recipe) {
  config <- recipe$preprocessing %||% list()
  roles <- recipe$roles %||% list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank"))
  x <- dataset$raw
  if (!is.matrix(x) || !is.numeric(x) || !ncol(x) || any(is.infinite(x))) stop("Raw intensities must contain finite or missing features.", call. = FALSE)
  blank_diagnostics <- data.frame()
  if (isTRUE(config$blank_filter$enabled)) {
    filtered <- blank_filter_matrix(x, dataset$sample_metadata, roles, config$blank_filter)
    x <- filtered$data
    blank_diagnostics <- filtered$diagnostics
  }
  if (!ncol(x)) stop("No features survived blank filtering.", call. = FALSE)
  after_blank <- x

  analysis_rows <- role_rows(dataset$sample_metadata, roles, "sample")
  analysis_metadata <- dataset$sample_metadata[analysis_rows, , drop = FALSE]
  analysis_x <- x[analysis_rows, , drop = FALSE]
  if (!nrow(analysis_x)) stop("No biological samples were found.", call. = FALSE)
  imputed <- impute_matrix(analysis_x, config$missing_values$method %||% "half_minimum")
  normalized <- normalize_matrix(imputed, config$normalization$method %||% "none")

  qc_diagnostics <- data.frame()
  if (isTRUE(config$qc_rsd_filter$enabled)) {
    qc_rows <- role_rows(dataset$sample_metadata, roles, "qc")
    qc_x <- x[qc_rows, , drop = FALSE]
    if (!nrow(qc_x)) stop("QC filtering requires QC samples.", call. = FALSE)
    qc_x <- impute_matrix(qc_x, config$missing_values$method %||% "half_minimum")
    qc_x <- normalize_matrix(qc_x, config$normalization$method %||% "none")
    qc_mean <- colMeans(qc_x, na.rm = TRUE)
    qc_rsd <- 100 * apply(qc_x, 2, stats::sd, na.rm = TRUE) / qc_mean
    threshold <- as.numeric(config$qc_rsd_filter$threshold_percent %||% 30)
    keep <- is.finite(qc_rsd) & qc_rsd <= threshold
    qc_diagnostics <- data.frame(feature_id = colnames(qc_x), qc_rsd_percent = qc_rsd, removed_by_qc_rsd = !keep)
    if (!any(keep)) stop("No features survived QC RSD filtering.", call. = FALSE)
    imputed <- imputed[, keep, drop = FALSE]
    normalized <- normalized[, keep, drop = FALSE]
  }

  if (!ncol(normalized)) stop("No features survived preprocessing.", call. = FALSE)
  if (any(!is.finite(normalized))) stop("Normalization produced missing or non-finite intensities.", call. = FALSE)
  transformed <- transform_matrix(normalized, config$transformation$method %||% "log2")
  if (any(!is.finite(transformed))) stop("Transformation produced missing or non-finite intensities; check imputation and transformation methods.", call. = FALSE)
  scaled <- scale_matrix(transformed, config$scaling$method %||% "pareto")
  if (any(!is.finite(scaled))) stop("Scaling produced missing or non-finite intensities.", call. = FALSE)
  feature_ids <- colnames(scaled)
  list(
    matrices = list(raw = dataset$raw, blank_filtered = after_blank, imputed = imputed, normalized = normalized, transformed = transformed, scaled = scaled),
    sample_metadata = analysis_metadata,
    variable_metadata = dataset$variable_metadata[feature_ids, , drop = FALSE],
    diagnostics = list(blank = blank_diagnostics, qc = qc_diagnostics),
    methods = config
  )
}

