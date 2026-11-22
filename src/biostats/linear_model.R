# A deliberately restricted model specification: factor aliases and their
# interactions only. Neither metadata nor recipe strings are evaluated as R code.
run_configured_models <- function(x, metadata, config) {
  scalar_label <- function(value, label) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
      stop(label, " must be a nonempty string.", call. = FALSE)
    value
  }
  terms_from <- function(value, label, aliases) {
    if (is.null(value)) value <- character()
    if (is.list(value)) value <- if (length(value)) unlist(value, use.names = FALSE) else character()
    if (!is.character(value) || anyNA(value) || anyDuplicated(value))
      stop(label, " must be a list of distinct factor aliases or interactions.", call. = FALSE)
    for (term in value) {
      parts <- strsplit(term, ":", fixed = TRUE)[[1L]]
      if (!length(parts) || startsWith(term, ":") || endsWith(term, ":") ||
          any(!nzchar(parts)) || anyDuplicated(parts) || !all(parts %in% aliases) ||
          grepl("::", term, fixed = TRUE))
        stop("Invalid ", label, " term: ", term, call. = FALSE)
    }
    unname(value)
  }
  named_specs <- function(value, label) {
    if (is.null(value)) return(list())
    if (!is.list(value) || (length(value) && !all(vapply(value, is.list, logical(1L)))))
      stop(label, " must be a list of named specifications.", call. = FALSE)
    labels <- vapply(value, function(item) scalar_label(item$name, paste0(label, " name")), character(1L))
    if (anyDuplicated(labels)) stop(label, " names must be unique.", call. = FALSE)
    value
  }
  design <- function(terms, frame) {
    formula <- if (length(terms)) stats::reformulate(terms, intercept = TRUE) else ~ 1
    stats::model.matrix(formula, data = frame)
  }
  fit <- function(matrix, response, label) {
    fitted <- stats::lm.fit(matrix, response)
    if (fitted$rank != ncol(matrix)) stop(label, " design is rank-deficient.", call. = FALSE)
    df <- nrow(matrix) - fitted$rank
    if (df <= 0L) stop(label, " design has no residual degrees of freedom.", call. = FALSE)
    residual_squares <- fitted$residuals^2
    rss <- if (is.null(dim(residual_squares))) sum(residual_squares) else colSums(residual_squares)
    list(fitted = fitted, rss = rss, df = df)
  }

  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x)) ||
      is.null(rownames(x)) || is.null(colnames(x)) || anyNA(rownames(x)) ||
      anyNA(colnames(x)) || any(!nzchar(rownames(x))) || any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
    stop("x must be a finite numeric sample-by-feature matrix with unique sample and feature IDs.", call. = FALSE)
  if (!is.data.frame(metadata) || is.null(rownames(metadata)) || anyDuplicated(rownames(metadata)))
    stop("metadata must have unique sample row names.", call. = FALSE)
  index <- match(rownames(x), rownames(metadata))
  if (anyNA(index)) stop("Every x sample must have a matching metadata row.", call. = FALSE)
  metadata <- metadata[index, , drop = FALSE]
  if (!is.list(config)) stop("Model config must be a list.", call. = FALSE)
  factors <- config$factors
  if (!is.list(factors) || is.null(names(factors)) || anyNA(names(factors)) ||
      anyDuplicated(names(factors)) || any(!grepl("^[A-Za-z][A-Za-z0-9_.]*$", names(factors))) ||
      any(names(factors) %in% c("if", "else", "repeat", "while", "function", "for", "in", "next", "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA")))
    stop("factors must be a named list with distinct simple aliases.", call. = FALSE)
  aliases <- names(factors)
  for (alias in aliases) {
    item <- factors[[alias]]
    if (!is.list(item)) stop("Factor ", alias, " needs column and levels.", call. = FALSE)
    column <- scalar_label(item$column, paste0("Factor ", alias, " column"))
    levels <- if (is.list(item$levels)) unlist(item$levels, use.names = FALSE) else item$levels
    if (!column %in% names(metadata)) stop("Factor column not found: ", column, call. = FALSE)
    if (!is.character(levels) || length(levels) < 2L || anyNA(levels) ||
        any(!nzchar(levels)) || anyDuplicated(levels))
      stop("Factor ", alias, " needs at least two distinct ordered levels.", call. = FALSE)
    factors[[alias]]$levels <- levels
  }
  terms <- terms_from(config$terms, "Model", aliases)
  contrasts <- named_specs(config$contrasts, "Contrast")
  tests <- named_specs(config$tests, "Test")
  for (test in tests) {
    if (is.null(test$reduced_terms) || is.null(test$full_terms))
      stop("Test ", test$name, " needs reduced_terms and full_terms.", call. = FALSE)
    terms_from(test$reduced_terms, paste0("Test ", test$name, " reduced"), aliases)
    terms_from(test$full_terms, paste0("Test ", test$name, " full"), aliases)
  }
  for (contrast in contrasts) {
    weights <- contrast$coefficients
    if (is.list(weights)) weights <- unlist(weights, use.names = TRUE)
    if (!is.numeric(weights) || !length(weights) || is.null(names(weights)) ||
        anyNA(names(weights)) || any(!nzchar(names(weights))) ||
        anyDuplicated(names(weights)) || any(!is.finite(weights)) || !any(weights != 0))
      stop("Contrast ", contrast$name, " needs finite, named, nonzero coefficients.", call. = FALSE)
  }

  if (!is.null(config$paired_change)) {
    pair <- config$paired_change
    if (!is.list(pair)) stop("paired_change must specify unit, time, from and to.", call. = FALSE)
    unit <- scalar_label(pair$unit, "paired_change unit")
    time <- scalar_label(pair$time, "paired_change time")
    from <- scalar_label(pair$from, "paired_change from")
    to <- scalar_label(pair$to, "paired_change to")
    if (!unit %in% names(metadata) || !time %in% names(metadata) || identical(from, to))
      stop("paired_change needs existing unit/time columns and distinct time labels.", call. = FALSE)
    relevant <- as.character(metadata[[time]]) %in% c(from, to)
    units <- as.character(metadata[[unit]])
    if (anyNA(units[relevant]) || any(!nzchar(units[relevant])))
      stop("paired_change unit IDs must be present.", call. = FALSE)
    before <- which(as.character(metadata[[time]]) == from)
    after <- which(as.character(metadata[[time]]) == to)
    if (anyDuplicated(units[before]) || anyDuplicated(units[after]))
      stop("paired_change requires unique unit IDs at each time; duplicate pair found.", call. = FALSE)
    match_after <- match(units[before], units[after])
    before <- before[!is.na(match_after)]
    after <- after[match_after[!is.na(match_after)]]
    if (!length(before)) stop("paired_change found no complete matched pairs.", call. = FALSE)
    for (alias in aliases) {
      column <- factors[[alias]]$column
      if (anyNA(metadata[[column]][before]) || anyNA(metadata[[column]][after]) ||
          any(as.character(metadata[[column]][before]) != as.character(metadata[[column]][after])))
        stop("Factor ", alias, " changed within a pair; group/batch must be invariant.", call. = FALSE)
    }
    x <- x[after, , drop = FALSE] - x[before, , drop = FALSE]
    metadata <- metadata[before, , drop = FALSE]
  }

  frame <- data.frame(row.names = seq_len(nrow(metadata)))
  for (alias in aliases) {
    column <- factors[[alias]]$column
    observed <- as.character(metadata[[column]])
    if (anyNA(observed) || !all(observed %in% factors[[alias]]$levels))
      stop("Factor ", alias, " contains missing or unlisted levels.", call. = FALSE)
    frame[[alias]] <- factor(observed, levels = factors[[alias]]$levels)
  }
  empty_contrasts <- data.frame(feature_id = character(), contrast = character(), effect = numeric(),
                                standard_error = numeric(), statistic = numeric(), residual_df = integer(),
                                p_value = numeric(), q_value = numeric(), stringsAsFactors = FALSE)
  empty_omnibus <- data.frame(feature_id = character(), test = character(), numerator_df = integer(),
                              residual_df = integer(), partial_r2 = numeric(), statistic = numeric(),
                              p_value = numeric(), q_value = numeric(), stringsAsFactors = FALSE)
  if (length(contrasts)) {
    if (!length(terms)) stop("Model terms are required for contrasts.", call. = FALSE)
    model <- design(terms, frame)
    f <- fit(model, x, "Contrast")
    pivot <- order(f$fitted$qr$pivot)
    covariance <- chol2inv(qr.R(f$fitted$qr))[pivot, pivot, drop = FALSE]
    coefficients <- f$fitted$coefficients
    contrast_rows <- lapply(contrasts, function(item) {
      weights <- item$coefficients
      if (is.list(weights)) weights <- unlist(weights, use.names = TRUE)
      if (!all(names(weights) %in% colnames(model)))
        stop("Contrast ", item$name, " references missing model coefficients: ",
             paste(setdiff(names(weights), colnames(model)), collapse = ", "), call. = FALSE)
      vector <- numeric(ncol(model))
      vector[match(names(weights), colnames(model))] <- weights
      effect <- drop(crossprod(vector, coefficients))
      standard_error <- sqrt(pmax(0, drop(crossprod(vector, covariance %*% vector))) * f$rss / f$df)
      statistic <- effect / standard_error
      p_value <- 2 * stats::pt(abs(statistic), df = f$df, lower.tail = FALSE)
      masked <- mask_unstable_p(p_value, standard_error)
      statistic[is.na(masked$p_value)] <- NA_real_
      data.frame(feature_id = colnames(x), contrast = item$name, effect = effect,
                 standard_error = standard_error, statistic = statistic, residual_df = f$df,
                 p_value = masked$p_value, q_value = masked$q_value, stringsAsFactors = FALSE)
    })
    contrast_table <- do.call(rbind, contrast_rows)
    rownames(contrast_table) <- NULL
  } else contrast_table <- empty_contrasts

  if (length(tests)) {
    omnibus_rows <- lapply(tests, function(item) {
      reduced_terms <- terms_from(item$reduced_terms, paste0("Test ", item$name, " reduced"), aliases)
      full_terms <- terms_from(item$full_terms, paste0("Test ", item$name, " full"), aliases)
      reduced <- design(reduced_terms, frame)
      full <- design(full_terms, frame)
      f0 <- fit(reduced, x, paste0("Test ", item$name, " reduced"))
      f1 <- fit(full, x, paste0("Test ", item$name, " full"))
      numerator_df <- ncol(full) - ncol(reduced)
      if (numerator_df <= 0L ||
          any(abs(qr.resid(f1$fitted$qr, reduced)) >
              1e-8 * max(1, abs(reduced))))
        stop("Test ", item$name, " reduced model is not strictly nested in its full model.", call. = FALSE)
      gain <- pmax(0, f0$rss - f1$rss)
      statistic <- (gain / numerator_df) / (f1$rss / f1$df)
      p_value <- stats::pf(statistic, numerator_df, f1$df, lower.tail = FALSE)
      masked <- mask_unstable_p(p_value, sqrt(f1$rss / f1$df))
      statistic[is.na(masked$p_value)] <- NA_real_
      partial_r2 <- gain / f0$rss
      partial_r2[is.na(masked$p_value)] <- NA_real_
      data.frame(feature_id = colnames(x), test = item$name,
                 numerator_df = numerator_df, residual_df = f1$df,
                 partial_r2 = partial_r2, statistic = statistic,
                 p_value = masked$p_value, q_value = masked$q_value, stringsAsFactors = FALSE)
    })
    omnibus_table <- do.call(rbind, omnibus_rows)
    rownames(omnibus_table) <- NULL
  } else omnibus_table <- empty_omnibus
  list(contrasts = contrast_table, omnibus = omnibus_table)
}
