repo_root <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- if (file.exists("src/biostats/variance.R")) "." else ".."
source(file.path(normalizePath(repo_root, mustWork = TRUE), "src", "biostats", "variance.R"))

expect_variance_error <- function(expression, pattern) {
  failure <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(failure, "error"), grepl(pattern, conditionMessage(failure), ignore.case = TRUE))
}

local({
  # Three genuinely different group centroids; spread is almost zero in the
  # first stratum and markedly greater in the second for every group.
  offsets <- c(-3, -2, -1, 1, 2, 3)
  rows <- list()
  strata <- character()
  groups <- character()
  for (stratum in c("site_one", "site_two")) for (group in c("alpha", "beta", "gamma")) {
    center <- rep(25, 3)
    dominant <- match(group, c("alpha", "beta", "gamma"))
    center[dominant] <- 100
    companion <- dominant %% 3 + 1L
    for (offset in offsets) {
      row <- center
      row[companion] <- row[companion] + offset *
        c(alpha = 0.2, beta = 1, gamma = 0.6)[[group]] *
        if (stratum == "site_one") 0.1 else 6
      rows[[length(rows) + 1L]] <- c(row, fixed_feature = 7)
      strata <- c(strata, stratum)
      groups <- c(groups, group)
    }
  }
  x <- do.call(rbind, rows)
  colnames(x) <- c("a", "b", "c", "fixed_feature")
  rownames(x) <- paste0("sample", seq_len(nrow(x)))
  metadata <- data.frame(treatment = groups, site = strata, row.names = rownames(x))
  config <- list(group_column = "treatment", strata_column = "site", permutations = 199L,
                 seed = 29L, scopes = list(list(name = "core", exclude_groups = "gamma"),
                                           list(name = "complete")))

  # The engine may change its implementation, but the public output must
  # preserve per-scope ranking and bias-adjusted within-cell distances.
  set.seed(83)
  before <- .Random.seed
  result <- run_permanova_dispersion(x, metadata, config)
  stopifnot(identical(before, .Random.seed),
            identical(names(result), c("permanova", "dispersion", "groups", "across_strata")),
            all(vapply(result, is.data.frame, logical(1))),
            nrow(result$permanova) == 4L, nrow(result$dispersion) == 4L,
            nrow(result$groups) == 10L, nrow(result$across_strata) == 4L,
            setequal(result$permanova$scope, c("core", "complete")),
            all(result$permanova$features == ncol(x)),
            all(vapply(split(result$permanova$adjusted_r2_rank, result$permanova$scope),
                       function(ranks) setequal(ranks, c(1L, 2L)), logical(1))),
            all(result$permanova$treatment_r2 > 0.5),
            all(result$permanova$permanova_p <= 0.05),
            all(result$dispersion$dispersion_p < 0.05),
            all(result$groups$n == length(offsets)),
            !any(result$groups$group == "gamma" & result$groups$scope == "core"))
  # Known two-point geometry: correction sqrt(n/(n-1)) gives each member
  # Euclidean distance between the pair divided by sqrt(2).
  selected <- c(1L, 2L, 7L, 8L, 19L, 20L, 25L, 26L)
  pair <- x[c(1L, 2L), , drop = FALSE]
  pair_hellinger <- sqrt(pair / rowSums(pair))
  expected_pair <- sqrt(sum((pair_hellinger[1L, ] - pair_hellinger[2L, ])^2) / 2)
  paired_config <- config
  paired_config$scopes <- list(list(name = "paired"))
  paired <- run_permanova_dispersion(x[selected, , drop = FALSE],
                                    metadata[selected, , drop = FALSE], paired_config)
  observed_pair <- paired$groups$mean_centroid_distance[
    paired$groups$stratum == "site_one" & paired$groups$group == "alpha"]
  stopifnot(length(observed_pair) == 1L,
            isTRUE(all.equal(observed_pair, expected_pair, tolerance = 1e-12)),
            all(result$groups$mean_centroid_distance >= 0),
            all(result$groups$sd_centroid_distance >= 0),
            all(result$groups$mean_centroid_distance[result$groups$stratum == "site_two"] >
                max(result$groups$mean_centroid_distance[result$groups$stratum == "site_one"])),
            all(result$across_strata$strata_p_within_group_permutation < 0.05),
            all(result$across_strata$mean_centroid_distance[
              result$across_strata$stratum == "site_two"] >
              max(result$across_strata$mean_centroid_distance[
                result$across_strata$stratum == "site_one"])))
  again <- run_permanova_dispersion(x, metadata, config)
  stopifnot(identical(result, again), identical(before, .Random.seed))
  altered <- config
  altered$scopes <- list(list(name = "only_core", exclude_groups = "gamma"))
  stopifnot(nrow(run_permanova_dispersion(x, metadata, altered)$permanova) == 2L)
  yaml_shaped <- config
  yaml_shaped$scopes[[1]]$exclude_groups <- list("gamma")
  stopifnot(identical(run_permanova_dispersion(x, metadata, yaml_shaped), result))

  invalid <- x
  invalid[1, 1] <- -1
  expect_variance_error(run_permanova_dispersion(invalid, metadata, config), "nonnegative")
  invalid <- x
  invalid[1, 1] <- NA_real_
  expect_variance_error(run_permanova_dispersion(invalid, metadata, config), "finite")
  invalid <- x
  invalid[1, ] <- 0
  expect_variance_error(run_permanova_dispersion(invalid, metadata, config), "zero|positive")
  no_residual <- x
  for (site in unique(metadata$site)) for (treatment in unique(metadata$treatment)) {
    members <- which(metadata$site == site & metadata$treatment == treatment)
    no_residual[members, ] <- matrix(no_residual[members[1L], ], length(members),
                                      ncol(x), byrow = TRUE)
  }
  expect_variance_error(run_permanova_dispersion(no_residual, metadata, config),
                        "zero residual")
  invalid_metadata <- metadata
  invalid_metadata$treatment[1] <- "singleton"
  expect_variance_error(run_permanova_dispersion(x, invalid_metadata, config), "replicat|singleton")
  invalid_metadata <- metadata
  invalid_metadata$treatment[invalid_metadata$site == "site_one"] <- "alpha"
  expect_variance_error(run_permanova_dispersion(x, invalid_metadata, config), "group|estimab")
  disconnected <- metadata
  at_second <- disconnected$site == "site_two"
  disconnected$treatment[at_second] <-
    paste0("other_", disconnected$treatment[at_second])
  all_groups <- config
  all_groups$scopes <- list(list(name = "all_groups"))
  expect_variance_error(run_permanova_dispersion(x, disconnected, all_groups), "estimab")
  stopifnot(identical(before, .Random.seed))
  invalid_config <- config
  invalid_config$permutations <- 98L
  expect_variance_error(run_permanova_dispersion(x, metadata, invalid_config), "permutations")
  invalid_config <- config
  invalid_config$seed <- NULL
  expect_variance_error(run_permanova_dispersion(x, metadata, invalid_config), "seed")
  invalid_config <- config
  invalid_config$scopes <- list(list(name = "first", exclude_groups = "gamma"),
                                list(name = "first"))
  expect_variance_error(run_permanova_dispersion(x, metadata, invalid_config), "scope")
})
