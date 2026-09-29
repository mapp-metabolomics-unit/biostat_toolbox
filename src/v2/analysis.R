run_pca <- function(x, metadata, group_column, components = 5) {
  if (nrow(x) < 2L || !ncol(x) || any(!is.finite(x))) stop("PCA requires at least two samples and finite, nonempty features.", call. = FALSE)
  variable <- apply(x, 2, stats::sd, na.rm = TRUE) > 0
  if (!any(variable)) stop("PCA has no varying features.", call. = FALSE)
  fit <- stats::prcomp(x[, variable, drop = FALSE], center = TRUE, scale. = FALSE, rank. = min(components, nrow(x) - 1, sum(variable)))
  scores <- data.frame(sample = rownames(fit$x), fit$x, check.names = FALSE)
  scores <- cbind(scores, metadata[match(scores$sample, rownames(metadata)), , drop = FALSE])
  loadings <- data.frame(feature_id = rownames(fit$rotation), fit$rotation, check.names = FALSE)
  variance <- data.frame(component = colnames(fit$x), variance_percent = 100 * fit$sdev[seq_len(ncol(fit$x))]^2 / sum(fit$sdev^2))
  list(scores = scores, loadings = loadings, variance = variance)
}

bray_curtis <- function(x) {
  n <- nrow(x)
  result <- matrix(0, n, n, dimnames = list(rownames(x), rownames(x)))
  if (n < 2L) stop("Bray-Curtis requires at least two samples.", call. = FALSE)
  for (i in seq_len(n - 1L)) for (j in seq.int(i + 1L, n)) {
    denominator <- sum(x[i, ] + x[j, ], na.rm = TRUE)
    value <- if (denominator == 0) 0 else sum(abs(x[i, ] - x[j, ]), na.rm = TRUE) / denominator
    result[i, j] <- result[j, i] <- value
  }
  stats::as.dist(result)
}

run_pcoa <- function(x, metadata, method = "bray", components = 5) {
  if (nrow(x) < 2L || !ncol(x) || any(!is.finite(x))) stop("PCoA requires at least two samples and finite, nonempty features.", call. = FALSE)
  if (method == "bray" && any(x < 0, na.rm = TRUE)) stop("Bray-Curtis requires a nonnegative matrix.", call. = FALSE)
  distance <- if (method == "bray") bray_curtis(x) else stats::dist(x, method = method)
  fit <- stats::cmdscale(distance, k = min(components, nrow(x) - 1), eig = TRUE, add = TRUE)
  points <- if (is.list(fit)) fit$points else fit
  points <- as.matrix(points)
  rownames(points) <- rownames(x)
  colnames(points) <- paste0("PCoA", seq_len(ncol(points)))
  scores <- data.frame(sample = rownames(points), points, check.names = FALSE)
  scores <- cbind(scores, metadata[match(scores$sample, rownames(metadata)), , drop = FALSE])
  positive <- fit$eig[fit$eig > 0]
  if (!length(positive)) stop("PCoA distances have no positive eigenvalues.", call. = FALSE)
  variance <- data.frame(component = paste0("PCoA", seq_along(positive)), variance_percent = 100 * positive / sum(positive))
  list(scores = scores, variance = variance)
}

design_matrix <- function(metadata, group_column, block_column = NULL) {
  if (is.null(group_column) || !group_column %in% names(metadata)) stop("Design group column is missing.", call. = FALSE)
  group <- factor(metadata[[group_column]])
  if (anyNA(group)) stop("The design group contains missing values.", call. = FALSE)
  frame <- data.frame(group = group)
  formula <- ~ 0 + group
  if (!is.null(block_column) && nzchar(block_column)) {
    if (!block_column %in% names(metadata)) stop("Block column not found: ", block_column, call. = FALSE)
    frame$block <- factor(metadata[[block_column]])
    if (anyNA(frame$block)) stop("The block column contains missing values.", call. = FALSE)
    formula <- ~ 0 + group + block
  }
  matrix <- stats::model.matrix(formula, frame)
  colnames(matrix) <- sub("^group", "", colnames(matrix))
  if (anyDuplicated(colnames(matrix)) || qr(matrix)$rank < ncol(matrix)) stop("The group/block design is rank-deficient; planned effects cannot be estimated.", call. = FALSE)
  if (nrow(matrix) <= ncol(matrix)) stop("The design has no residual degrees of freedom.", call. = FALSE)
  list(matrix = matrix, group = group)
}

fit_contrast <- function(x, design, numerator, denominator, contrast_name, effect_label) {
  if (!all(c(numerator, denominator) %in% colnames(design))) stop("Contrast levels are absent from the design: ", numerator, ", ", denominator, call. = FALSE)
  if (identical(numerator, denominator)) stop("A contrast must compare distinct groups.", call. = FALSE)
  if (any(!is.finite(x))) stop("Differential analysis requires finite transformed intensities.", call. = FALSE)
  contrast <- numeric(ncol(design)); names(contrast) <- colnames(design)
  contrast[numerator] <- 1; contrast[denominator] <- -1
  fit <- stats::lm.fit(design, x)
  if (fit$rank < ncol(design) || nrow(design) <= fit$rank) stop("Contrast design is rank-deficient or lacks residual degrees of freedom.", call. = FALSE)
  coefficients <- fit$coefficients
  effects <- as.numeric(crossprod(contrast, coefficients))
  residual_df <- nrow(design) - fit$rank
  sigma2 <- colSums(fit$residuals^2, na.rm = TRUE) / residual_df
  covariance <- chol2inv(qr.R(qr(design)))
  standard_error <- sqrt(as.numeric(crossprod(contrast, covariance %*% contrast)) * sigma2)
  statistic <- effects / standard_error
  p_value <- 2 * stats::pt(abs(statistic), df = residual_df, lower.tail = FALSE)
  data.frame(
    feature_id = colnames(x), contrast = contrast_name, numerator = numerator, denominator = denominator,
    effect = effects, effect_scale = effect_label, standard_error = standard_error,
    statistic = statistic, degrees_of_freedom = residual_df, p_value = p_value,
    q_value = stats::p.adjust(p_value, method = "BH"), stringsAsFactors = FALSE
  )
}

run_differential <- function(x, metadata, design_config, transformation_method) {
  group_column <- design_config$group
  block_column <- design_config$block %||% NULL
  design <- design_matrix(metadata, group_column, block_column)
  contrasts <- design_config$contrasts %||% list()
  if (!length(contrasts)) return(data.frame())
  effect_label <- if (tolower(transformation_method) == "log2") "log2_fold_change" else "difference_on_transformed_scale"
  rows <- lapply(contrasts, function(contrast) {
    fit_contrast(x, design$matrix, contrast$numerator, contrast$denominator, contrast$name %||% paste(contrast$numerator, "vs", contrast$denominator, sep = "_"), effect_label)
  })
  do.call(rbind, rows)
}

run_omnibus <- function(x, metadata, group_column, block_column = NULL) {
  design_full <- design_matrix(metadata, group_column, block_column)$matrix
  if (any(!is.finite(x))) stop("Omnibus analysis requires finite transformed intensities.", call. = FALSE)
  if (is.null(block_column) || !nzchar(block_column)) {
    design_reduced <- matrix(1, nrow(x), 1)
  } else {
    design_reduced <- stats::model.matrix(~ factor(metadata[[block_column]]))
  }
  full <- stats::lm.fit(design_full, x)
  reduced <- stats::lm.fit(design_reduced, x)
  df1 <- full$rank - reduced$rank
  df2 <- nrow(x) - full$rank
  if (full$rank < ncol(design_full) || reduced$rank < ncol(design_reduced) || df1 < 1L || df2 < 1L) stop("Omnibus design is rank-deficient or lacks test degrees of freedom.", call. = FALSE)
  rss_full <- colSums(full$residuals^2)
  rss_reduced <- colSums(reduced$residuals^2)
  statistic <- (pmax(0, rss_reduced - rss_full) / df1) / (rss_full / df2)
  p_value <- stats::pf(statistic, df1, df2, lower.tail = FALSE)
  data.frame(feature_id = colnames(x), statistic = statistic, numerator_df = df1, denominator_df = df2, p_value = p_value, q_value = stats::p.adjust(p_value, "BH"))
}

run_analyses <- function(processed, recipe) {
  analyses <- recipe$analyses %||% list()
  group_column <- recipe$design$group
  result <- list()
  if (isTRUE(analyses$pca$enabled)) result$pca <- run_pca(processed$matrices$scaled, processed$sample_metadata, group_column, analyses$pca$components %||% 5)
  if (isTRUE(analyses$pcoa$enabled)) {
    stage <- analyses$pcoa$stage %||% "normalized"
    if (!stage %in% c("imputed", "normalized", "transformed", "scaled")) stop("PCoA stage must contain biological samples: imputed, normalized, transformed, or scaled.", call. = FALSE)
    result$pcoa <- run_pcoa(processed$matrices[[stage]], processed$sample_metadata, analyses$pcoa$distance %||% "bray", analyses$pcoa$components %||% 5)
  }
  if (isTRUE(analyses$omnibus$enabled)) result$omnibus <- run_omnibus(processed$matrices$transformed, processed$sample_metadata, group_column, recipe$design$block %||% NULL)
  if (isTRUE(analyses$differential$enabled)) result$differential <- run_differential(processed$matrices$transformed, processed$sample_metadata, recipe$design, processed$methods$transformation$method %||% "log2")
  result
}

