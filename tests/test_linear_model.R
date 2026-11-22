repo_root <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- if (file.exists("src/biostats/linear_model.R")) "." else ".."
repo_root <- normalizePath(repo_root, mustWork = TRUE)
source(file.path(repo_root, "src", "biostats", "analysis.R"))
source(file.path(repo_root, "src", "biostats", "linear_model.R"))

expect_model_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), ignore.case = TRUE))
}

local({
  # Four equally represented cells, with independently known residual sum of
  # squares 8. Treatment coding yields groupB=2, batchsecond=3, interaction=4.
  metadata <- expand.grid(rep = 1:3, batch_label = c("first", "second"),
                          treatment = c("A", "B"), KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)
  rownames(metadata) <- paste0("sample", seq_len(nrow(metadata)))
  g <- as.integer(metadata$treatment == "B")
  b <- as.integer(metadata$batch_label == "second")
  noise <- c(-1, 0, 1)[metadata$rep]
  y <- 5 + 2 * g + 3 * b + 4 * g * b + noise
  x <- cbind(signal = y, flat = rep(7, length(y)))
  rownames(x) <- rownames(metadata)
  factors <- list(group = list(column = "treatment", levels = c("A", "B")),
                  batch = list(column = "batch_label", levels = c("first", "second")))
  additive <- list(factors = factors, terms = c("group", "batch"),
                   contrasts = list(list(name = "treated_adjusted", coefficients = c(groupB = 1))))
  # Use a genuinely additive feature to test adjustment separately from the
  # deliberately interacting signal above.
  x[, "signal"] <- 5 + 2 * g + 3 * b + noise
  adjusted <- run_configured_models(x, metadata[rev(rownames(metadata)), , drop = FALSE], additive)
  stopifnot(identical(names(adjusted), c("contrasts", "omnibus")),
            identical(adjusted$contrasts$feature_id, colnames(x)),
            nrow(adjusted$omnibus) == 0L,
            identical(adjusted$contrasts$contrast, rep("treated_adjusted", 2L)),
            isTRUE(all.equal(adjusted$contrasts$effect[1], 2, tolerance = 1e-12)),
            isTRUE(all.equal(adjusted$contrasts$standard_error[1], sqrt(8 / 27), tolerance = 1e-12)),
            adjusted$contrasts$residual_df[1] == 9L,
            isTRUE(all.equal(adjusted$contrasts$p_value[1],
                             2 * pt(2 / sqrt(8 / 27), 9, lower.tail = FALSE), tolerance = 1e-12)),
            identical(adjusted$contrasts$q_value[1], adjusted$contrasts$p_value[1]),
            is.na(adjusted$contrasts$p_value[2]), is.na(adjusted$contrasts$q_value[2]))
  recipe_shape <- additive
  recipe_shape$factors$group$levels <- list("A", "B")
  recipe_shape$factors$batch$levels <- list("first", "second")
  recipe_shape$terms <- list("group", "batch")
  recipe_shape$contrasts[[1]]$coefficients <- list(groupB = 1)
  stopifnot(isTRUE(all.equal(run_configured_models(x, metadata, recipe_shape)$contrasts,
                           run_configured_models(x, metadata, additive)$contrasts)))

  x[, "signal"] <- y
  factorial <- list(factors = factors, terms = c("group", "batch", "group:batch"),
                    contrasts = list(list(name = "difference_in_differences",
                                          coefficients = c("groupB:batchsecond" = 1))),
                    tests = list(list(name = "added_interaction", reduced_terms = c("group", "batch"),
                                      full_terms = c("group", "batch", "group:batch"))))
  result <- run_configured_models(x, metadata, factorial)
  c1 <- result$contrasts[1, ]
  o1 <- result$omnibus[1, ]
  stopifnot(identical(c1$contrast, "difference_in_differences"),
            isTRUE(all.equal(c1$effect, 4, tolerance = 1e-12)),
            isTRUE(all.equal(c1$standard_error, sqrt(4 / 3), tolerance = 1e-12)),
            c1$residual_df == 8L,
            identical(o1$test, "added_interaction"), o1$numerator_df == 1L,
            o1$residual_df == 8L,
            isTRUE(all.equal(o1$statistic, 12, tolerance = 1e-12)),
            isTRUE(all.equal(o1$partial_r2, 0.6, tolerance = 1e-12)),
            isTRUE(all.equal(o1$p_value, pf(12, 1, 8, lower.tail = FALSE), tolerance = 1e-12)),
            is.na(result$omnibus$p_value[2]), is.na(result$omnibus$q_value[2]),
            is.na(result$omnibus$partial_r2[2]))
  only_test <- run_configured_models(x, metadata, list(factors = factors,
    tests = list(list(name = "batch_given_group", reduced_terms = "group",
                      full_terms = c("group", "batch")))))
  stopifnot(nrow(only_test$contrasts) == 0L, nrow(only_test$omnibus) == 2L,
            identical(only_test$omnibus$test, rep("batch_given_group", 2L)))

  # Four cells of three plants, with one additional unpaired baseline plant;
  # differences must be calculated before fitting, not modeled as 25 independent rows.
  paired_meta <- metadata
  paired_meta$unit <- paste0("plant", seq_len(nrow(metadata)))
  paired_meta$time <- "before"
  after <- paired_meta
  after$time <- "after"
  rownames(after) <- paste0("after", seq_len(nrow(after)))
  rownames(paired_meta) <- paste0("before", seq_len(nrow(paired_meta)))
  extra <- paired_meta[1, , drop = FALSE]
  extra$unit <- "unpaired"
  rownames(extra) <- "unpaired_before"
  pairs <- rbind(paired_meta, after, extra)
  before <- 10 + seq_len(nrow(metadata))
  paired_x <- matrix(c(before, before + 1 + 2 * g + 3 * b + noise, 99), ncol = 1,
                     dimnames = list(rownames(pairs), "signal"))
  paired_cfg <- list(factors = factors, terms = c("group", "batch"),
                     paired_change = list(unit = "unit", time = "time", from = "before", to = "after"),
                     contrasts = list(list(name = "infected_minus_healthy_delta", coefficients = c(groupB = 1))))
  paired_result <- run_configured_models(paired_x, pairs, paired_cfg)$contrasts
  stopifnot(isTRUE(all.equal(paired_result$effect, 2, tolerance = 1e-12)),
            isTRUE(all.equal(paired_result$standard_error, sqrt(8 / 27), tolerance = 1e-12)),
            paired_result$residual_df == 9L)
  duplicate <- rbind(pairs, pairs["before1", , drop = FALSE])
  rownames(duplicate)[nrow(duplicate)] <- "duplicate_before"
  duplicate_x <- rbind(paired_x, paired_x["before1", , drop = FALSE])
  rownames(duplicate_x)[nrow(duplicate_x)] <- "duplicate_before"
  expect_model_error(run_configured_models(duplicate_x, duplicate, paired_cfg), "duplicate|unique")
  changed <- pairs
  changed["after1", "batch_label"] <- "second"
  expect_model_error(run_configured_models(paired_x, changed, paired_cfg), "invariant|changed")
  expect_model_error(run_configured_models(x, metadata, list(factors = factors,
    terms = c("group", "batch"), contrasts = list(list(name = "invalid", coefficients = c(batchthird = 1))))),
    "coefficient")
  expect_model_error(run_configured_models(x, metadata, list(factors = factors,
    terms = c("group", "batch"), contrasts = list(list(name = "zero", coefficients = c(groupB = 0))))),
    "zero|nonzero")
  expect_model_error(run_configured_models(x, metadata, list(factors = factors,
    terms = c("group + batch"), contrasts = list())), "term")
  confounded <- metadata
  confounded$batch_label <- ifelse(confounded$treatment == "A", "first", "second")
  expect_model_error(run_configured_models(x, confounded, additive), "rank|level")
  expect_model_error(run_configured_models(x[1:3, , drop = FALSE], metadata[1:3, , drop = FALSE],
    list(factors = list(group = factors$group), terms = "group",
         contrasts = list(list(name = "B", coefficients = c(groupB = 1))))), "rank|level|residual")
  expect_model_error(run_configured_models(x[c(1, 7), , drop = FALSE],
    metadata[c(1, 7), , drop = FALSE],
    list(factors = list(group = factors$group), terms = "group",
         contrasts = list(list(name = "B", coefficients = c(groupB = 1))))), "residual")
  expect_model_error(run_configured_models(x, metadata,
    list(factors = factors, tests = list(list(name = "not_nested",
      reduced_terms = "group", full_terms = "batch")))), "nested")
  multiple <- x
  multiple[, "signal"] <- 5 + 2 * g + 3 * b + noise
  multiple[, "flat"] <- 5 + 2 * g + 3 * b + 2 * noise
  two <- run_configured_models(multiple, metadata, list(
    factors = factors, terms = c("group", "batch"),
    contrasts = list(list(name = "group_effect", coefficients = c(groupB = 1)),
                     list(name = "batch_effect", coefficients = c(batchsecond = 1)))))$contrasts
  stopifnot(identical(two$contrast, rep(c("group_effect", "batch_effect"), each = 2L)),
            isTRUE(all.equal(two$q_value[1:2], p.adjust(two$p_value[1:2], "BH"))),
            isTRUE(all.equal(two$q_value[3:4], p.adjust(two$p_value[3:4], "BH"))))
  empty <- run_configured_models(x, metadata, list(factors = factors, terms = character(0),
                                                   contrasts = list(), tests = list()))
  stopifnot(is.data.frame(empty$contrasts), is.data.frame(empty$omnibus),
            nrow(empty$contrasts) == 0L, nrow(empty$omnibus) == 0L)
})
