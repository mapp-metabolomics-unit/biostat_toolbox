# One-way Hellinger/Euclidean PERMANOVA and bias-adjusted centroid dispersion.
# x contains normalized, nonnegative intensities for biological samples only.
run_permanova_dispersion <- function(x, metadata, config) {
  scalar_label <- function(value, name) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(trimws(value)))
      stop(name, " must be a nonempty string.", call. = FALSE)
    value
  }
  integer_setting <- function(value, name, minimum) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < minimum || value > .Machine$integer.max || value != floor(value))
      stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
    as.integer(value)
  }
  if (!is.list(config)) stop("Variance config must be a list.", call. = FALSE)
  group_column <- scalar_label(config$group_column, "group_column")
  strata_column <- scalar_label(config$strata_column, "strata_column")
  permutations <- integer_setting(config$permutations, "permutations", 99L)
  seed <- integer_setting(config$seed, "seed", 0L)
  scopes <- config$scopes
  if (!is.list(scopes) || !length(scopes))
    stop("Variance config requires one or more scopes.", call. = FALSE)
  scope_names <- vapply(scopes, function(scope) {
    if (!is.list(scope)) stop("Each scope must be an object.", call. = FALSE)
    scalar_label(scope$name, "scope name")
  }, character(1))
  if (anyDuplicated(scope_names)) stop("Variance scope names must be unique.", call. = FALSE)

  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x))
    stop("Variance x must be a nonempty numeric sample-by-feature matrix.", call. = FALSE)
  if (any(!is.finite(x))) stop("Variance matrix must contain only finite values.", call. = FALSE)
  if (any(x < 0)) stop("Variance matrix must be nonnegative.", call. = FALSE)
  totals <- rowSums(x)
  if (any(!is.finite(totals) | totals <= 0))
    stop("Variance matrix has zero or invalid sample totals.", call. = FALSE)
  if (!is.data.frame(metadata) || nrow(metadata) != nrow(x) ||
      !all(c(group_column, strata_column) %in% names(metadata)))
    stop("Variance metadata must have one row per sample and the configured group and strata columns.", call. = FALSE)
  if (!is.null(rownames(x))) {
    if (anyDuplicated(rownames(x)) || anyDuplicated(rownames(metadata)) ||
        !setequal(rownames(x), rownames(metadata)))
      stop("Variance matrix and metadata sample row names must match uniquely.", call. = FALSE)
    metadata <- metadata[rownames(x), , drop = FALSE]
  }
  labels <- as.character(metadata[[group_column]])
  strata <- as.character(metadata[[strata_column]])
  if (anyNA(labels) || any(!nzchar(trimws(labels))) ||
      anyNA(strata) || any(!nzchar(trimws(strata))))
    stop("Variance group and strata labels must not be missing or empty.", call. = FALSE)
  scope_exclusions <- lapply(scopes, function(scope) {
    exclusions <- scope$exclude_groups
    if (is.list(exclusions) && is.null(names(exclusions)))
      exclusions <- unlist(exclusions, use.names = FALSE)
    if (!is.null(exclusions) && (!is.character(exclusions) || anyNA(exclusions) ||
                                  any(!nzchar(trimws(exclusions))) ||
                                  anyDuplicated(exclusions)))
      stop("scope exclude_groups must be distinct nonempty labels.", call. = FALSE)
    if (length(setdiff(exclusions, labels)))
      stop("scope exclude_groups contains unknown group labels.", call. = FALSE)
    exclusions
  })

  # Transform once using every retained feature and the full-study sample totals;
  # subsetting groups or strata must never change the feature space.
  hellinger <- sqrt(x / totals)
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) previous_rng <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_rng) assign(".Random.seed", previous_rng, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)

  permutation_p <- function(observed, null) (1 + sum(null >= observed)) / (1 + length(null))
  variance_fit <- function(distance_sq, group) {
    n <- length(group)
    total <- sum(distance_sq) / (2 * n)
    within <- sum(vapply(unique(group), function(level) {
      indices <- which(group == level)
      sum(distance_sq[indices, indices, drop = FALSE]) / (2 * length(indices))
    }, numeric(1)))
    between <- max(0, total - within)
    groups <- length(unique(group))
    c(r2 = between / total,
      f = if (within == 0) Inf else (between / (groups - 1L)) / (within / (n - groups)))
  }
  centroid_distances <- function(a, group) {
    distances <- numeric(nrow(a))
    for (level in unique(group)) {
      indices <- which(group == level)
      n <- length(indices)
      z <- a[indices, , drop = FALSE]
      centered <- sweep(z, 2, colMeans(z), "-")
      distances[indices] <- sqrt(rowSums(centered^2)) * sqrt(n / (n - 1L))
    }
    distances
  }
  oneway_f <- function(distances, group) {
    levels <- unique(group)
    center <- mean(distances)
    means <- vapply(levels, function(level) mean(distances[group == level]), numeric(1))
    between <- sum(vapply(seq_along(levels), function(i) {
      sum(group == levels[[i]]) * (means[[i]] - center)^2
    }, numeric(1)))
    within <- sum(vapply(seq_along(levels), function(i) {
      sum((distances[group == levels[[i]]] - means[[i]])^2)
    }, numeric(1)))
    if (within == 0) return(if (between == 0) 0 else Inf)
    (between / (length(levels) - 1L)) / (within / (length(group) - length(levels)))
  }

  permanova_rows <- list()
  dispersion_rows <- list()
  group_rows <- list()
  across_rows <- list()
  for (scope_index in seq_along(scopes)) {
    scope <- scopes[[scope_index]]
    scope_name <- scope_names[[scope_index]]
    selected <- which(!labels %in% scope_exclusions[[scope_index]])
    strata_levels <- unique(strata[selected])
    if (length(strata_levels) < 2L)
      stop("Variance scope ", scope_name, " needs at least two strata for across-stratum spread.", call. = FALSE)
    distance_all <- numeric(length(selected))
    for (stratum in strata_levels) {
      positions <- which(strata[selected] == stratum)
      indices <- selected[positions]
      group <- labels[indices]
      counts <- table(group)
      n <- length(indices)
      g <- length(counts)
      if (g < 2L || any(counts < 2L) || n <= g)
        stop("Variance scope ", scope_name, "/", stratum,
             " requires at least two groups, each with >= 2 replicates.", call. = FALSE)
      a <- hellinger[indices, , drop = FALSE]
      distance_sq <- as.matrix(stats::dist(a))^2
      fit <- variance_fit(distance_sq, group)
      if (!is.finite(fit[["f"]]))
        stop("Undefined PERMANOVA: zero residual variance for ", scope_name, "/", stratum,
             ".", call. = FALSE)
      null_f <- replicate(permutations,
                          variance_fit(distance_sq, group[sample.int(n)])[["f"]])
      permanova_rows[[length(permanova_rows) + 1L]] <- data.frame(
        scope = scope_name, stratum = stratum, samples = n, groups = g,
        features = ncol(hellinger), treatment_r2 = unname(fit[["r2"]]),
        treatment_adjusted_r2 = 1 - (1 - fit[["r2"]]) * (n - 1) / (n - g),
        permanova_f = unname(fit[["f"]]),
        permanova_p = permutation_p(fit[["f"]], null_f), permutations = permutations)

      d <- centroid_distances(a, group)
      distance_all[positions] <- d
      observed <- oneway_f(d, group)
      residual <- d - ave(d, group, FUN = mean)
      null_spread <- replicate(permutations, oneway_f(residual[sample.int(n)], group))
      dispersion_rows[[length(dispersion_rows) + 1L]] <- data.frame(
        scope = scope_name, stratum = stratum, dispersion_f = observed,
        dispersion_p = permutation_p(observed, null_spread), permutations = permutations)
      for (level in names(counts)) group_rows[[length(group_rows) + 1L]] <- data.frame(
        scope = scope_name, stratum = stratum, group = level, n = unname(counts[[level]]),
        mean_centroid_distance = mean(d[group == level]),
        sd_centroid_distance = stats::sd(d[group == level]))
    }
    group <- factor(labels[selected])
    stratum <- factor(strata[selected], levels = strata_levels)
    reduced <- stats::model.matrix(~ group)
    full <- stats::model.matrix(~ group + stratum)
    if (qr(full)$rank != ncol(full) || length(selected) <= ncol(full))
      stop("Variance scope ", scope_name,
           " has a non-estimable group-adjusted strata model.", call. = FALSE)
    reduced_rss <- sum(stats::lm.fit(reduced, distance_all)$residuals^2)
    full_f <- function(stratum_labels) {
      design <- stats::model.matrix(~ group + stratum_labels)
      full_rss <- sum(stats::lm.fit(design, distance_all)$residuals^2)
      between <- max(0, reduced_rss - full_rss)
      if (full_rss == 0) return(if (between == 0) 0 else Inf)
      (between / (ncol(design) - ncol(reduced))) /
        (full_rss / (length(distance_all) - ncol(design)))
    }
    observed <- full_f(stratum)
    null_across <- replicate(permutations, {
      permuted <- stratum
      for (level in levels(group)) {
        positions <- which(group == level)
        permuted[positions] <- stratum[positions[sample.int(length(positions))]]
      }
      full_f(permuted)
    })
    across_p <- permutation_p(observed, null_across)
    for (level in strata_levels) across_rows[[length(across_rows) + 1L]] <- data.frame(
      scope = scope_name, stratum = level, samples = length(selected),
      strata = length(strata_levels), strata_f_adjusted_for_group = observed,
      strata_p_within_group_permutation = across_p, permutations = permutations,
      mean_centroid_distance = mean(distance_all[stratum == level]))
  }
  permanova <- do.call(rbind, permanova_rows)
  permanova$adjusted_r2_rank <- ave(-permanova$treatment_adjusted_r2,
                                     permanova$scope, FUN = function(z) rank(z, ties.method = "min"))
  rownames(permanova) <- NULL
  dispersion <- do.call(rbind, dispersion_rows)
  rownames(dispersion) <- NULL
  groups <- do.call(rbind, group_rows)
  rownames(groups) <- NULL
  across_strata <- do.call(rbind, across_rows)
  rownames(across_strata) <- NULL
  # Across-stratum inference permutes batch/strata labels within groups while
  # retaining estimated group-by-stratum centroids; its p value is approximate.
  list(permanova = permanova, dispersion = dispersion, groups = groups,
       across_strata = across_strata)
}
