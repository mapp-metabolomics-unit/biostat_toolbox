# PLS-DA fits one-hot class responses with pls::plsr. All feature selection and
# preprocessing parameters are estimated inside each training fold.
stratified_folds <- function(labels, folds) {
  assignment <- integer(length(labels))
  for (level in levels(labels)) {
    indices <- which(labels == level)
    assignment[indices[sample.int(length(indices))]] <- rep(seq_len(folds), length.out = length(indices))
  }
  assignment
}

supervised_normalize_test <- function(x, training, method) {
  method <- tolower(method %||% "none")
  if (method == "none") return(x)
  if (method %in% c("total_sum", "median")) {
    measure <- function(z) if (method == "total_sum") rowSums(z) else apply(z, 1, stats::median)
    size <- measure(x)
    training_size <- measure(training)
    target <- stats::median(training_size)
    if (any(!is.finite(size) | size <= 0) || !is.finite(target) || target <= 0) stop("PLS-DA normalization has an invalid sample size.", call. = FALSE)
    return(sweep(x, 1, size / target, "/"))
  }
  if (method == "pqn") {
    reference <- apply(training, 2, stats::median)
    quotients <- sweep(x, 2, reference, "/")
    quotients[!is.finite(quotients)] <- NA_real_
    factors <- apply(quotients, 1, stats::median, na.rm = TRUE)
    if (any(!is.finite(factors) | factors <= 0)) stop("PLS-DA PQN has an invalid dilution factor.", call. = FALSE)
    return(sweep(x, 1, factors, "/"))
  }
  stop("Unknown normalization method: ", method, call. = FALSE)
}

supervised_preprocess <- function(dataset, recipe, sample_rows, train, holdout = integer()) {
  settings <- recipe$preprocessing %||% list()
  roles <- recipe$roles %||% list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank"))
  metadata <- dataset$sample_metadata
  raw <- dataset$raw
  keep <- rep(TRUE, ncol(raw))
  if (isTRUE(settings$blank_filter$enabled)) {
    blanks <- which(role_rows(metadata, roles, "blank"))
    filter_rows <- c(sample_rows[train], blanks)
    filtered <- blank_filter_matrix(raw[filter_rows, , drop = FALSE], metadata[filter_rows, , drop = FALSE], roles, settings$blank_filter)
    keep <- !filtered$diagnostics$removed_as_blank
  }
  if (!any(keep)) stop("PLS-DA training fold has no features after blank filtering.", call. = FALSE)
  train_x <- raw[sample_rows[train], keep, drop = FALSE]
  test_x <- raw[sample_rows[holdout], keep, drop = FALSE]
  imputation <- tolower(settings$missing_values$method %||% "half_minimum")
  if (imputation == "half_minimum") {
    positive <- train_x[is.finite(train_x) & train_x > 0]
    if (!length(positive)) stop("PLS-DA training fold has no positive intensities.", call. = FALSE)
    replacement <- min(positive) / 2
    train_x[is.na(train_x) | train_x <= 0] <- replacement
    test_x[is.na(test_x) | test_x <= 0] <- replacement
  } else if (imputation != "none") stop("Unknown imputation method: ", imputation, call. = FALSE)
  if (any(!is.finite(train_x)) || any(!is.finite(test_x))) stop("PLS-DA requires finite intensities after training-fold imputation.", call. = FALSE)

  normalization <- settings$normalization$method %||% "none"
  normalized_test <- supervised_normalize_test(test_x, train_x, normalization)
  normalized_train <- normalize_matrix(train_x, normalization)
  if (isTRUE(settings$qc_rsd_filter$enabled)) {
    qc_rows <- which(role_rows(metadata, roles, "qc"))
    if (length(qc_rows) < 3L) stop("QC RSD filtering requires at least three QCs.", call. = FALSE)
    qc_x <- raw[qc_rows, keep, drop = FALSE]
    qc_x <- normalize_matrix(impute_matrix(qc_x, imputation), normalization)
    qc_rsd <- 100 * apply(qc_x, 2, stats::sd, na.rm = TRUE) / colMeans(qc_x, na.rm = TRUE)
    stable <- is.finite(qc_rsd) & qc_rsd <= as.numeric(settings$qc_rsd_filter$threshold_percent %||% 30)
    if (!any(stable)) stop("PLS-DA training fold has no features after QC filtering.", call. = FALSE)
    normalized_train <- normalized_train[, stable, drop = FALSE]
    normalized_test <- normalized_test[, stable, drop = FALSE]
  }
  transform <- settings$transformation$method %||% "log2"
  training <- transform_matrix(normalized_train, transform)
  testing <- transform_matrix(normalized_test, transform)
  if (any(!is.finite(training)) || any(!is.finite(testing))) stop("PLS-DA transformation produced non-finite intensities.", call. = FALSE)
  scaling <- tolower(settings$scaling$method %||% "pareto")
  if (scaling != "none") {
    if (!scaling %in% c("pareto", "autoscale")) stop("Unknown scaling method: ", scaling, call. = FALSE)
    center <- colMeans(training)
    deviations <- apply(training, 2, stats::sd)
    divisor <- if (scaling == "pareto") sqrt(deviations) else deviations
    divisor[!is.finite(divisor) | divisor == 0] <- 1
    testing <- sweep(sweep(testing, 2, center, "-"), 2, divisor, "/")
    training <- sweep(sweep(training, 2, center, "-"), 2, divisor, "/")
  }
  # Constant features cannot contribute latent components and are identified only
  # from the training set, never from held-out samples.
  variable <- apply(training, 2, stats::sd) > 0
  if (!any(variable)) stop("PLS-DA training fold has no varying features.", call. = FALSE)
  list(train = training[, variable, drop = FALSE], test = testing[, variable, drop = FALSE])
}

supervised_components <- function(x, limit) {
  min(limit, ncol(x), nrow(x) - 1L, qr(sweep(x, 2, colMeans(x), "-"))$rank)
}

supervised_fit <- function(x, labels, components) {
  response <- matrix(0, nrow(x), nlevels(labels), dimnames = list(NULL, levels(labels)))
  response[cbind(seq_along(labels), as.integer(labels))] <- 1
  pls::plsr(response ~ x, data = list(response = response, x = x),
            ncomp = components, method = "oscorespls", validation = "none", scale = FALSE)
}

supervised_predict <- function(fit, x, components, classes) {
  predicted <- predict(fit, newdata = list(x = x), ncomp = seq_len(components))
  if (any(!is.finite(predicted))) stop("PLS-DA produced non-finite held-out decision values.", call. = FALSE)
  result <- matrix(NA_character_, nrow(x), components)
  for (component in seq_len(components)) {
    values <- matrix(predicted[, , component], nrow = nrow(x), ncol = length(classes))
    result[, component] <- classes[max.col(values, ties.method = "first")]
  }
  result
}

supervised_balanced_accuracy <- function(truth, predicted, classes) {
  mean(vapply(classes, function(level) mean(predicted[truth == level] == level), numeric(1)))
}

select_supervised_components <- function(dataset, recipe, sample_rows, training, labels, limit, folds) {
  assignment <- stratified_folds(labels[training], folds)
  predictions <- vector("list", folds)
  possible <- integer(folds)
  for (fold in seq_len(folds)) {
    inner_train <- training[assignment != fold]
    inner_test <- training[assignment == fold]
    processed <- supervised_preprocess(dataset, recipe, sample_rows, inner_train, inner_test)
    possible[fold] <- supervised_components(processed$train, limit)
    if (possible[fold] < 1L) stop("PLS-DA inner training fold cannot fit a latent component.", call. = FALSE)
    model <- supervised_fit(processed$train, labels[inner_train], possible[fold])
    predictions[[fold]] <- supervised_predict(model, processed$test, possible[fold], levels(labels))
  }
  candidate_count <- min(possible)
  predicted <- matrix(NA_character_, length(training), candidate_count)
  for (fold in seq_len(folds)) predicted[assignment == fold, ] <- predictions[[fold]][, seq_len(candidate_count), drop = FALSE]
  quality <- vapply(seq_len(candidate_count), function(component) {
    supervised_balanced_accuracy(as.character(labels[training]), predicted[, component], levels(labels))
  }, numeric(1))
  which.max(quality) # The smallest component count wins ties.
}

validate_supervised <- function(dataset, recipe, sample_rows, labels, options) {
  fold_predictions <- vector("list", options$repeats * options$folds)
  result_index <- 0L
  for (repeat_index in seq_len(options$repeats)) {
    assignment <- stratified_folds(labels, options$folds)
    for (fold in seq_len(options$folds)) {
      test <- which(assignment == fold)
      train <- which(assignment != fold)
      inner_folds <- min(options$folds, min(table(labels[train])))
      processed <- supervised_preprocess(dataset, recipe, sample_rows, train, test)
      limit <- supervised_components(processed$train, options$components)
      if (limit < 1L) stop("PLS-DA outer training fold cannot fit a latent component.", call. = FALSE)
      selected <- select_supervised_components(dataset, recipe, sample_rows, train, labels, limit, inner_folds)
      model <- supervised_fit(processed$train, labels[train], selected)
      predicted <- supervised_predict(model, processed$test, selected, levels(labels))[, selected]
      result_index <- result_index + 1L
      fold_predictions[[result_index]] <- data.frame(
        repeat_index = repeat_index, fold = fold, sample = rownames(dataset$raw)[sample_rows[test]],
        truth = as.character(labels[test]), predicted = predicted, components = selected,
        stringsAsFactors = FALSE
      )
    }
  }
  predictions <- do.call(rbind, fold_predictions)
  rownames(predictions) <- NULL
  repeat_quality <- vapply(seq_len(options$repeats), function(index) {
    subset <- predictions[predictions$repeat_index == index, , drop = FALSE]
    supervised_balanced_accuracy(subset$truth, subset$predicted, levels(labels))
  }, numeric(1))
  list(predictions = predictions, balanced_accuracy = mean(repeat_quality), accuracy = mean(predictions$truth == predictions$predicted))
}

run_supervised <- function(dataset, recipe) {
  options <- recipe$analyses$plsda %||% list()
  if (!isTRUE(options$enabled)) return(NULL)
  if (!requireNamespace("pls", quietly = TRUE)) stop("PLS-DA requires the CRAN 'pls' package.", call. = FALSE)
  if (!is.null(recipe$design$block)) stop("PLS-DA does not account for design.block; remove the block or disable PLS-DA.", call. = FALSE)
  group_column <- recipe$design$group
  if (is.null(group_column) || !group_column %in% names(dataset$sample_metadata)) stop("PLS-DA requires a valid design.group column.", call. = FALSE)
  roles <- recipe$roles %||% list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank"))
  sample_rows <- which(role_rows(dataset$sample_metadata, roles, "sample"))
  labels <- factor(dataset$sample_metadata[[group_column]][sample_rows])
  if (anyNA(labels)) stop("PLS-DA groups contain missing values.", call. = FALSE)
  integer_option <- function(name, default, minimum) {
    value <- options[[name]] %||% default
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value < minimum ||
        value > .Machine$integer.max || value != floor(value)) {
      stop("PLS-DA ", name, " must be an integer >= ", minimum, ".", call. = FALSE)
    }
    as.integer(value)
  }
  settings <- list(components = integer_option("components", 3L, 1L), folds = integer_option("folds", 3L, 2L),
                   repeats = integer_option("repeats", 5L, 2L), permutations = integer_option("permutations", 99L, 19L))
  empty <- data.frame()
  withheld <- function(reason) list(summary = empty, fold_predictions = empty, permutations = empty,
                                    scores = empty, status = "withheld", reason = reason)
  if (nlevels(labels) < 2L) return(withheld("PLS-DA needs at least two biological groups."))
  counts <- table(labels)
  if (any(counts < 4L)) return(withheld("PLS-DA needs at least four biological samples per group for nested validation."))
  if (settings$folds > min(counts)) return(withheld("PLS-DA folds exceed the smallest group size."))
  if (ncol(dataset$raw) < 1L) return(withheld("PLS-DA needs at least one measured feature."))
  set.seed(as.integer(recipe$seed %||% 1L))
  observed <- validate_supervised(dataset, recipe, sample_rows, labels, settings)
  null <- lapply(seq_len(settings$permutations), function(index) {
    shuffled <- factor(sample(as.character(labels)), levels = levels(labels))
    validation <- validate_supervised(dataset, recipe, sample_rows, shuffled, settings)
    data.frame(permutation = index, balanced_accuracy = validation$balanced_accuracy, accuracy = validation$accuracy)
  })
  null <- do.call(rbind, null)
  p_value <- (1 + sum(null$balanced_accuracy >= observed$balanced_accuracy)) / (settings$permutations + 1)
  all_samples <- seq_along(sample_rows)
  full_data <- supervised_preprocess(dataset, recipe, sample_rows, all_samples)
  full_limit <- supervised_components(full_data$train, settings$components)
  if (full_limit < 1L) stop("PLS-DA full training data cannot fit a latent component.", call. = FALSE)
  full_components <- select_supervised_components(dataset, recipe, sample_rows, all_samples, labels,
                                                  full_limit, settings$folds)
  full_fit <- supervised_fit(full_data$train, labels, full_components)
  if (any(!is.finite(full_fit$scores[, seq_len(full_components), drop = FALSE]))) stop("PLS-DA produced non-finite full-fit scores.", call. = FALSE)
  scores <- data.frame(sample = rownames(dataset$raw)[sample_rows],
                       full_fit$scores[, seq_len(full_components), drop = FALSE], check.names = FALSE)
  names(scores)[-1L] <- paste0("LV", seq_len(full_components))
  scores[[group_column]] <- as.character(labels)
  summary <- data.frame(model = "one_hot_pls_da", validation = "nested_repeated_stratified_cv",
                        selection = "inner_stratified_cv_balanced_accuracy",
                        preprocessing = "fold_local_blank_qc_imputation_normalization_transformation_scaling",
                        samples = length(labels), groups = nlevels(labels), folds = settings$folds,
                        repeats = settings$repeats, permutations = settings$permutations,
                        balanced_accuracy = observed$balanced_accuracy, accuracy = observed$accuracy,
                        null_mean_balanced_accuracy = mean(null$balanced_accuracy),
                        null_sd_balanced_accuracy = stats::sd(null$balanced_accuracy),
                        permutation_p_value = p_value, full_fit_components = full_components,
                        stringsAsFactors = FALSE)
  list(summary = summary, fold_predictions = observed$predictions, permutations = null,
       scores = scores, status = "validated", reason = NA_character_)
}
