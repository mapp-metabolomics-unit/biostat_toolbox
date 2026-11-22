role_rows <- function(metadata, roles, names) {
  metadata[[roles$column]] %in% unlist(roles$levels[names], use.names = FALSE)
}

blank_filter_matrix <- function(x, metadata, roles, config) {
  sample_rows <- role_rows(metadata, roles, "sample")
  blank_rows <- role_rows(metadata, roles, "blank")
  if (!any(sample_rows) || !any(blank_rows)) stop("Blank filtering requires samples and blanks.", call. = FALSE)
  sample_median <- apply(x[sample_rows, , drop = FALSE], 2, function(values)
    if (all(is.na(values))) 0 else stats::median(values, na.rm = TRUE))
  blank_median <- apply(x[blank_rows, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  blank_observed <- colSums(is.finite(x[blank_rows, , drop = FALSE]))
  blank_detected <- colSums(x[blank_rows, , drop = FALSE] > 0, na.rm = TRUE)
  blank_maximum <- apply(x[blank_rows, , drop = FALSE], 2, function(values)
    if (all(is.na(values))) NA_real_ else max(values, na.rm = TRUE))
  blank_prevalence <- colSums(x[blank_rows, , drop = FALSE] > 0, na.rm = TRUE) / pmax(blank_observed, 1L)
  blank_prevalence[blank_observed == 0L] <- NA_real_
  ratio <- ifelse(is.na(blank_median), NA_real_,
                  ifelse(blank_median > 0, sample_median / blank_median, Inf))
  if (any(!is.finite(sample_median))) stop("Blank filtering needs finite sample summaries.", call. = FALSE)
  minimum_ratio <- as.numeric(config$minimum_sample_to_blank_ratio %||% 5)
  minimum_prevalence <- as.numeric(config$minimum_blank_prevalence %||% 0.5)
  contaminated <- !is.na(blank_prevalence) & !is.na(ratio) &
    blank_prevalence >= minimum_prevalence & ratio < minimum_ratio
  diagnostics <- data.frame(
    feature_id = colnames(x), sample_median = sample_median, blank_median = blank_median,
    blank_observed = blank_observed, blank_detected = blank_detected,
    blank_maximum = blank_maximum, blank_prevalence = blank_prevalence,
    sample_to_blank_ratio = ratio, removed_as_blank = contaminated,
    blank_data_insufficient = blank_observed == 0L, stringsAsFactors = FALSE
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

impute_matrix <- function(x, method, reference = NULL) {
  method <- tolower(method %||% "half_minimum")
  if (any(is.infinite(x))) stop("Intensities must be finite or missing.", call. = FALSE)
  if (method == "none") return(x)
  if (method == "half_minimum") {
    if (!is.null(reference) && (length(reference) != ncol(x) ||
                                any(!is.finite(reference) | reference <= 0)))
      stop("Imputation requires a positive observed biological minimum for each feature.", call. = FALSE)
    for (column in seq_len(ncol(x))) {
      observed <- x[, column]
      positive <- observed[is.finite(observed) & observed > 0]
      if (!length(positive) && is.null(reference))
        stop("Cannot impute feature without an observed positive intensity: ", colnames(x)[column], call. = FALSE)
      missing <- is.na(observed) | observed <= 0
      if (any(missing)) {
        replacement <- (if (is.null(reference)) min(positive) else reference[column]) / 2
        if (!is.finite(replacement) || replacement <= 0)
          stop("Feature half-minimum underflows or is invalid: ", colnames(x)[column], call. = FALSE)
        x[missing, column] <- replacement
      }
    }
    return(x)
  }
  stop("Unknown imputation method: ", method, call. = FALSE)
}

qc_drift_correct <- function(x, metadata, roles, settings) {
  enabled <- isTRUE(settings$enabled)
  if (enabled && !identical(settings$method, "loess")) stop("QC drift method must be loess.", call. = FALSE)
  order_column <- settings$order_column %||% (if (!enabled) "injection_order" else NULL)
  if (enabled && (!is.character(order_column) || length(order_column) != 1L ||
                  is.na(order_column) || !nzchar(order_column) || !order_column %in% names(metadata)))
    stop("QC drift requires an existing order_column in sample metadata.", call. = FALSE)
  injection_order <- if (is.character(order_column) && length(order_column) == 1L &&
                         !is.na(order_column) && order_column %in% names(metadata)) metadata[[order_column]] else NULL
  valid_order <- is.numeric(injection_order) && length(injection_order) == nrow(x) &&
    all(is.finite(injection_order) & injection_order > 0) && !anyDuplicated(injection_order)
  if (enabled && !valid_order)
    stop("QC drift injection orders must be unique, finite positive numbers for all injections.", call. = FALSE)
  qc_rows <- role_rows(metadata, roles, "qc")
  valid_qc_order <- is.numeric(injection_order) && length(injection_order) == nrow(x) &&
    all(is.finite(injection_order[qc_rows]) & injection_order[qc_rows] > 0) &&
    !anyDuplicated(injection_order[qc_rows])
  minimum_qc <- settings$minimum_qc %||% 6L
  minimum_qc_fraction <- settings$minimum_qc_fraction %||% 0.8
  span <- settings$span %||% 0.75
  maximum_gap <- settings$maximum_gap_fraction %||% 0.25
  maximum_fold <- settings$maximum_fold_change %||% 5
  if (enabled) {
    if (!is.numeric(minimum_qc) || length(minimum_qc) != 1L || !is.finite(minimum_qc) ||
        minimum_qc < 6 || minimum_qc != floor(minimum_qc) ||
        !is.numeric(minimum_qc_fraction) || length(minimum_qc_fraction) != 1L ||
        !is.finite(minimum_qc_fraction) || minimum_qc_fraction < 0.5 || minimum_qc_fraction > 1 ||
        !is.numeric(span) || length(span) != 1L || !is.finite(span) || span < 0.5 || span > 1 ||
        !is.numeric(maximum_gap) || length(maximum_gap) != 1L || !is.finite(maximum_gap) ||
        maximum_gap < 0.05 || maximum_gap > 1 ||
        !is.numeric(maximum_fold) || length(maximum_fold) != 1L ||
        !is.finite(maximum_fold) || maximum_fold < 1 || maximum_fold > 100)
      stop("QC drift requires minimum_qc >= 6, minimum_qc_fraction in [0.5, 1], span in [0.5, 1], maximum_gap_fraction in [0.05, 1] and maximum_fold_change in [1, 100].", call. = FALSE)
    if (sum(qc_rows) < minimum_qc)
      stop("QC drift requires at least minimum_qc pooled QC injections.", call. = FALSE)
  }
  diagnostics <- data.frame(
    feature_id = colnames(x), qc_observed = integer(ncol(x)), qc_missing = integer(ncol(x)),
    qc_zero = integer(ncol(x)), qc_negative = integer(ncol(x)), qc_na = integer(ncol(x)),
    qc_positive_minimum = rep(NA_real_, ncol(x)),
    qc_rsd_before_percent = rep(NA_real_, ncol(x)), qc_rsd_after_percent = rep(NA_real_, ncol(x)),
    qc_order_spearman_before = rep(NA_real_, ncol(x)),
    qc_order_spearman_after = rep(NA_real_, ncol(x)),
    trend_status = rep("insufficient_qc", ncol(x)),
    status = rep(if (enabled) "insufficient_qc" else "disabled", ncol(x)),
    corrected_positive = integer(ncol(x)), outside_qc_range = integer(ncol(x)),
    stringsAsFactors = FALSE
  )
  qc_order <- if (valid_qc_order) injection_order[qc_rows] else NULL
  run_range <- if (enabled) diff(range(injection_order)) else NA_real_
  qc_summary <- function(values, order) {
    if (length(values) < 3L) return(c(rsd = NA_real_, trend = NA_real_))
    mean_value <- mean(values)
    spread <- stats::sd(values)
    rsd <- if (is.finite(mean_value) && mean_value > 0 && is.finite(spread))
      100 * (spread / mean_value) else NA_real_
    trend <- if (is.numeric(order) && length(order) == length(values) &&
                 is.finite(spread) && spread > 1e-12 * mean_value)
      stats::cor(order, values, method = "spearman") else NA_real_
    c(rsd = rsd, trend = trend)
  }
  for (column in seq_len(ncol(x))) {
    qc_values <- x[qc_rows, column]
    observed <- is.finite(qc_values) & qc_values > 0
    observed_order <- if (valid_qc_order) qc_order[observed] else NULL
    diagnostics$qc_observed[column] <- sum(observed)
    diagnostics$qc_missing[column] <- length(qc_values) - sum(observed)
    diagnostics$qc_zero[column] <- sum(!is.na(qc_values) & qc_values == 0)
    diagnostics$qc_negative[column] <- sum(!is.na(qc_values) & qc_values < 0)
    diagnostics$qc_na[column] <- sum(is.na(qc_values))
    if (!any(observed)) next
    diagnostics$qc_positive_minimum[column] <- min(qc_values[observed])
    before <- qc_summary(qc_values[observed], observed_order)
    diagnostics$qc_rsd_before_percent[column] <- diagnostics$qc_rsd_after_percent[column] <- before[["rsd"]]
    diagnostics$qc_order_spearman_before[column] <- diagnostics$qc_order_spearman_after[column] <- before[["trend"]]
    diagnostics$trend_status[column] <- if (length(qc_values[observed]) < 3L) "insufficient_qc" else
      if (!is.finite(before[["rsd"]])) "invalid_qc" else if (!valid_qc_order) "order_unavailable" else
      if (is.na(before[["trend"]])) "constant_qc" else "available"
    if (!enabled) next
    if (sum(observed) < minimum_qc || mean(observed) < minimum_qc_fraction) next
    if (!is.finite(before[["rsd"]])) {
      diagnostics$status[column] <- "invalid_qc"
      next
    }
    if (max(diff(sort(observed_order))) / run_range > maximum_gap) {
      diagnostics$status[column] <- "qc_gap"
      next
    }
    if (is.finite(before[["rsd"]]) && is.na(before[["trend"]]) && before[["rsd"]] < 1e-10) {
      diagnostics$status[column] <- "stable_qc"
      next
    }
    fit <- tryCatch(stats::loess(log_intensity ~ order,
      data = data.frame(log_intensity = log(qc_values[observed]), order = observed_order),
      span = span, degree = 1L,
      control = stats::loess.control(surface = "direct")),
      warning = function(condition) NULL, error = function(condition) NULL)
    if (is.null(fit)) {
      diagnostics$status[column] <- "fit_failed"
      next
    }
    values <- x[, column]
    eligible <- is.finite(values) & values > 0
    inside <- eligible & injection_order >= min(observed_order) & injection_order <= max(observed_order)
    diagnostics$outside_qc_range[column] <- sum(eligible & !inside)
    predicted <- tryCatch(stats::predict(fit, newdata = data.frame(order = injection_order[inside])),
      warning = function(condition) NULL, error = function(condition) NULL)
    reference_fit <- tryCatch(stats::predict(fit, newdata = data.frame(order = observed_order)),
      warning = function(condition) NULL, error = function(condition) NULL)
    if (is.null(predicted) || is.null(reference_fit) || any(!is.finite(c(predicted, reference_fit)))) {
      diagnostics$status[column] <- "fit_failed"
      next
    }
    reference <- stats::median(reference_fit)
    if (any(abs(reference - predicted) > log(maximum_fold))) {
      diagnostics$status[column] <- "excessive_correction"
      next
    }
    corrected <- exp(log(values[inside]) + reference - predicted)
    if (any(!is.finite(corrected) | corrected <= 0)) {
      diagnostics$status[column] <- "fit_failed"
      next
    }
    x[inside, column] <- corrected
    diagnostics$corrected_positive[column] <- sum(inside)
    diagnostics$status[column] <- if (diagnostics$outside_qc_range[column]) "partial_no_extrapolation" else "corrected"
    after <- qc_summary(x[qc_rows, column][observed], observed_order)
    diagnostics$qc_rsd_after_percent[column] <- after[["rsd"]]
    diagnostics$qc_order_spearman_after[column] <- after[["trend"]]
  }
  list(data = x, diagnostics = diagnostics)
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

preprocess_ordination <- function(x, config, feature_ids, biological_minimum) {
  x <- x[, feature_ids, drop = FALSE]
  imputed <- impute_matrix(x, config$missing_values$method %||% "half_minimum", biological_minimum[feature_ids])
  normalized <- normalize_matrix(imputed, config$normalization$method %||% "none")
  normalized <- normalized[, feature_ids, drop = FALSE]
  if (any(!is.finite(normalized))) stop("Ordination normalization produced missing or non-finite intensities.", call. = FALSE)
  transformed <- transform_matrix(normalized, config$transformation$method %||% "log2")
  if (any(!is.finite(transformed))) stop("Ordination transformation produced missing or non-finite intensities.", call. = FALSE)
  scaled <- scale_matrix(transformed, config$scaling$method %||% "pareto")
  if (any(!is.finite(scaled))) stop("Ordination scaling produced missing or non-finite intensities.", call. = FALSE)
  list(imputed = imputed, normalized = normalized, transformed = transformed, scaled = scaled)
}


qc_rsd_diagnostics <- function(qc_x, settings, missing_method = "half_minimum", normalization = "none") {
  if (nrow(qc_x) < 3L) stop("QC filtering requires at least three QC injections.", call. = FALSE)
  minimum_fraction <- settings$minimum_qc_fraction %||% 0.8
  threshold <- settings$threshold_percent %||% 30
  if (!is.numeric(minimum_fraction) || length(minimum_fraction) != 1L ||
      !is.finite(minimum_fraction) || minimum_fraction < 0.5 || minimum_fraction > 1 ||
      !is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold < 0)
    stop("QC RSD requires minimum_qc_fraction in [0.5, 1] and nonnegative threshold_percent.", call. = FALSE)
  observed <- is.finite(qc_x) & qc_x > 0
  qc_count <- colSums(observed)
  qc_fraction <- qc_count / nrow(qc_x)
  sufficient <- qc_count >= 3L & qc_fraction >= minimum_fraction
  qc_rsd <- rep(NA_real_, ncol(qc_x))
  if (any(sufficient)) {
    supported <- qc_x[, sufficient, drop = FALSE]
    qc_complete <- impute_matrix(supported, missing_method)
    qc_complete <- normalize_matrix(qc_complete, normalization)
    qc_complete[!observed[, sufficient, drop = FALSE]] <- NA_real_
    qc_rsd[sufficient] <- 100 * apply(qc_complete, 2, stats::sd, na.rm = TRUE) /
      colMeans(qc_complete, na.rm = TRUE)
  }
  keep <- sufficient & is.finite(qc_rsd) & qc_rsd <= threshold
  data.frame(
    feature_id = colnames(qc_x), qc_observed = qc_count, qc_support_fraction = qc_fraction,
    insufficient_qc_support = !sufficient, qc_rsd_percent = qc_rsd,
    removed_by_qc_rsd = !keep, stringsAsFactors = FALSE
  )
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
  if (!any(analysis_rows)) stop("No biological samples were found.", call. = FALSE)
  observed <- colSums(is.finite(dataset$raw[analysis_rows, , drop = FALSE]) &
                        dataset$raw[analysis_rows, , drop = FALSE] > 0)
  detection <- data.frame(
    feature_id = colnames(dataset$raw), observed_biological = observed,
    missing_biological = sum(analysis_rows) - observed,
    excluded_no_biological_detection = observed == 0L,
    removed_by_blank = colnames(dataset$raw) %in% blank_diagnostics$feature_id[blank_diagnostics$removed_as_blank],
    stringsAsFactors = FALSE
  )
  x <- x[, observed[colnames(x)] > 0L, drop = FALSE]
  if (!ncol(x)) stop("No biologically detected features survived blank filtering.", call. = FALSE)
  corrected <- qc_drift_correct(x, dataset$sample_metadata, roles, config$qc_drift %||% list())
  x <- corrected$data
  analysis_x <- x[analysis_rows, , drop = FALSE]
  biological_minimum <- apply(analysis_x, 2, function(values) min(values[is.finite(values) & values > 0]))
  imputed <- impute_matrix(analysis_x, config$missing_values$method %||% "half_minimum")
  normalized <- normalize_matrix(imputed, config$normalization$method %||% "none")

  qc_diagnostics <- data.frame()
  if (isTRUE(config$qc_rsd_filter$enabled)) {
    qc_rows <- role_rows(dataset$sample_metadata, roles, "qc")
    qc_x <- x[qc_rows, , drop = FALSE]
    qc_diagnostics <- qc_rsd_diagnostics(
      qc_x, config$qc_rsd_filter, config$missing_values$method %||% "half_minimum",
      config$normalization$method %||% "none")
    keep <- !qc_diagnostics$removed_by_qc_rsd
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
  ordination <- list()
  for (analysis in c("pca", "pcoa")) {
    selection <- recipe$analyses[[analysis]]
    if (!isTRUE(selection$enabled) ||
        !(isTRUE(selection$include_qc) || isTRUE(selection$include_blanks))) next
    if (identical(analysis, "pcoa") && !is.null(ordination$pca) &&
        identical(c(isTRUE(selection$include_qc), isTRUE(selection$include_blanks)),
                  c(isTRUE(recipe$analyses$pca$include_qc), isTRUE(recipe$analyses$pca$include_blanks)))) {
      ordination$pcoa <- ordination$pca
      next
    }
    rows <- analysis_rows
    if (isTRUE(selection$include_qc)) rows <- rows | role_rows(dataset$sample_metadata, roles, "qc")
    if (isTRUE(selection$include_blanks)) rows <- rows | role_rows(dataset$sample_metadata, roles, "blank")
    ordination[[analysis]] <- list(
      matrices = preprocess_ordination(x[rows, , drop = FALSE], config, feature_ids, biological_minimum),
      sample_metadata = dataset$sample_metadata[rows, , drop = FALSE]
    )
  }
  list(
    matrices = list(raw = dataset$raw, blank_filtered = after_blank, imputed = imputed, normalized = normalized, transformed = transformed, scaled = scaled),
    ordination = ordination,
    sample_metadata = analysis_metadata,
    variable_metadata = dataset$variable_metadata[feature_ids, , drop = FALSE],
    diagnostics = list(blank = blank_diagnostics, qc = qc_diagnostics, detection = detection,
                       drift = corrected$diagnostics),
    methods = config
  )
}

