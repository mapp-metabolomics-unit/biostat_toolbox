repo_root <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- if (file.exists("src/biostats/preprocess.R")) "." else ".."
repo_root <- normalizePath(repo_root, mustWork = TRUE)
for (module in c("utils.R", "config.R", "preprocess.R", "annotations.R", "export.R"))
  source(file.path(repo_root, "src", "biostats", module))

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

fixture <- function() {
  order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  role <- c("QC", "sample", "QC", "blank", "QC", "sample",
            "QC", "blank", "QC", "sample", "QC", "sample")
  sample <- paste0("injection", order)
  metadata <- data.frame(sample_type = role, injection_order = order,
                         group = ifelse(role == "sample", "A", NA_character_),
                         row.names = sample)
  raw <- cbind(
    drifting = ifelse(role == "blank", 0, 100 * exp(order / 10)),
    stable = ifelse(role == "blank", 0, 40),
    sparse_qc = ifelse(role == "sample", 20, 0),
    undetected = ifelse(role == "sample", NA_real_, 50),
    weak_qc = ifelse(role == "sample", 5, 0)
  )
  raw[role == "QC", "sparse_qc"] <- c(0, 0, 0, 0, 0, 20)
  raw[role == "QC", "weak_qc"] <- c(5, 0, 5, 0, 5, 0)
  rownames(raw) <- sample
  list(raw = raw, sample_metadata = metadata,
       variable_metadata = data.frame(feature_id = colnames(raw), row.names = colnames(raw)),
       annotations = list())
}

base_recipe <- function() {
  list(kind = "study", display_group = "group",
       preprocessing = list(missing_values = list(method = "half_minimum"),
                            normalization = list(method = "none"),
                            transformation = list(method = "log2"),
                            scaling = list(method = "none")),
       analyses = list())
}

local({
  dataset <- fixture()
  recipe <- base_recipe()
  dataset$raw[c("injection2", "injection6"), "stable"] <- c(0, 40)
  dataset$raw[c("injection2", "injection6"), "sparse_qc"] <- c(NA_real_, 20)
  original <- dataset$raw
  processed <- preprocess_dataset(dataset, recipe)
  stopifnot(identical(dataset$raw, original), identical(processed$matrices$raw, original),
            !"undetected" %in% colnames(processed$matrices$imputed),
            processed$matrices$imputed["injection2", "stable"] == 20,
            processed$matrices$imputed["injection2", "sparse_qc"] == 10,
            processed$matrices$imputed["injection2", "drifting"] > 100,
            processed$diagnostics$detection$observed_biological[
              processed$diagnostics$detection$feature_id == "undetected"] == 0L,
            processed$diagnostics$detection$excluded_no_biological_detection[
              processed$diagnostics$detection$feature_id == "undetected"],
            all(processed$diagnostics$drift$status == "disabled"),
            processed$diagnostics$drift$qc_observed[
              processed$diagnostics$drift$feature_id == "sparse_qc"] == 1L,
            processed$diagnostics$drift$qc_order_spearman_before[
              processed$diagnostics$drift$feature_id == "drifting"] == 1)
  recipe$analyses$pca <- list(enabled = TRUE, include_qc = TRUE, include_blanks = TRUE)
  with_controls <- preprocess_dataset(dataset, recipe)
  stopifnot(all(is.finite(with_controls$ordination$pca$matrices$scaled)),
            with_controls$ordination$pca$matrices$imputed["injection2", "stable"] == 20,
            with_controls$ordination$pca$matrices$imputed["injection4", "sparse_qc"] == 10)
  dataset$raw[, "stable"] <- 0
  dataset$raw[dataset$sample_metadata$sample_type == "sample", "drifting"] <- NA_real_
  dataset$raw[dataset$sample_metadata$sample_type == "sample", "sparse_qc"] <- 0
  dataset$raw[dataset$sample_metadata$sample_type == "sample", "weak_qc"] <- NA_real_
  expect_error(preprocess_dataset(dataset, base_recipe()), "No biologically detected features")
})

local({
  dataset <- fixture()
  blank_rows <- dataset$sample_metadata$sample_type == "blank"
  dataset$raw[blank_rows, "stable"] <- NA_real_
  dataset$raw[blank_rows, "undetected"] <- 50
  recipe <- base_recipe()
  recipe$preprocessing$blank_filter <- list(enabled = TRUE,
    minimum_sample_to_blank_ratio = 5, minimum_blank_prevalence = 0.5)
  recipe$analyses$pca <- list(enabled = TRUE, include_blanks = TRUE)
  result <- preprocess_dataset(dataset, recipe)
  blank <- result$diagnostics$blank
  no_blanks <- blank[blank$feature_id == "stable", , drop = FALSE]
  stopifnot(no_blanks$blank_data_insufficient, is.na(no_blanks$blank_prevalence),
            !no_blanks$removed_as_blank,
            no_blanks$blank_detected == 0L, is.na(no_blanks$blank_maximum),
            blank$blank_detected[blank$feature_id == "undetected"] == 2L,
            all(is.finite(result$ordination$pca$matrices$imputed)),
            result$diagnostics$detection$excluded_no_biological_detection[
              result$diagnostics$detection$feature_id == "undetected"])
  dataset$sample_metadata$injection_order[dataset$sample_metadata$sample_type == "sample"] <- NA_real_
  incomplete_biological_order <- preprocess_dataset(dataset, recipe)
  stopifnot(incomplete_biological_order$diagnostics$drift$qc_order_spearman_before[
              incomplete_biological_order$diagnostics$drift$feature_id == "drifting"] == 1)
  dataset$sample_metadata$injection_order <- NULL
  without_order <- preprocess_dataset(dataset, recipe)
  stopifnot(all(without_order$diagnostics$drift$status == "disabled"),
            all(without_order$diagnostics$drift$trend_status == "order_unavailable" |
                  without_order$diagnostics$drift$trend_status == "insufficient_qc"))
})

local({
  dataset <- fixture()
  recipe <- base_recipe()
  recipe$preprocessing$qc_drift <- list(enabled = TRUE, method = "loess",
                                       order_column = "injection_order", span = 0.75,
                                       minimum_qc = 6L, maximum_gap_fraction = 0.25)
  validate_recipe_config(recipe)
  original <- dataset$raw
  corrected <- preprocess_dataset(dataset, recipe)
  d <- corrected$diagnostics$drift
  drift <- d[d$feature_id == "drifting", , drop = FALSE]
  stable <- d[d$feature_id == "stable", , drop = FALSE]
  sparse <- d[d$feature_id == "sparse_qc", , drop = FALSE]
  stopifnot(identical(dataset$raw, original), identical(corrected$matrices$raw, original),
            identical(corrected$matrices$blank_filtered, original),
            drift$status == "partial_no_extrapolation", drift$outside_qc_range == 1L,
            drift$qc_observed == 6L, drift$qc_missing == 0L,
            drift$qc_rsd_after_percent < drift$qc_rsd_before_percent / 10,
            abs(corrected$matrices$imputed["injection2", "drifting"] -
                corrected$matrices$imputed["injection6", "drifting"]) < 1e-5,
            corrected$matrices$imputed["injection12", "drifting"] == original["injection12", "drifting"],
            stable$status == "stable_qc", stable$corrected_positive == 0L,
            abs(corrected$matrices$imputed["injection6", "stable"] - 40) < 1e-8,
            sparse$status == "insufficient_qc", sparse$qc_observed == 1L,
            sparse$qc_missing == 5L, sparse$qc_zero == 5L, sparse$corrected_positive == 0L,
            corrected$matrices$imputed["injection2", "sparse_qc"] == 20)

  recipe$preprocessing$qc_drift$maximum_gap_fraction <- 0.1
  no_fit <- preprocess_dataset(dataset, recipe)
  stopifnot(no_fit$diagnostics$drift$status[
              no_fit$diagnostics$drift$feature_id == "drifting"] == "qc_gap",
            no_fit$matrices$imputed["injection2", "drifting"] == original["injection2", "drifting"])
  recipe$preprocessing$qc_drift$maximum_gap_fraction <- 0.25
  recipe$preprocessing$qc_drift$maximum_fold_change <- 1
  capped <- preprocess_dataset(dataset, recipe)
  stopifnot(capped$diagnostics$drift$status[
              capped$diagnostics$drift$feature_id == "drifting"] == "excessive_correction",
            capped$matrices$imputed["injection2", "drifting"] == original["injection2", "drifting"])
  recipe$preprocessing$qc_drift$maximum_fold_change <- NULL
  missing_order <- dataset
  missing_order$sample_metadata$injection_order <- NULL
  expect_error(preprocess_dataset(missing_order, recipe), "order_column")
  duplicate_order <- dataset
  duplicate_order$sample_metadata$injection_order[2] <- 1
  expect_error(preprocess_dataset(duplicate_order, recipe), "unique, finite positive")
  too_few_qcs <- dataset
  too_few_qcs$raw <- too_few_qcs$raw[rownames(too_few_qcs$raw) != "injection11", , drop = FALSE]
  too_few_qcs$sample_metadata <- too_few_qcs$sample_metadata[
    rownames(too_few_qcs$sample_metadata) != "injection11", , drop = FALSE]
  expect_error(preprocess_dataset(too_few_qcs, recipe), "at least minimum_qc")
  for (bad in list(list(span = 0.1), list(minimum_qc = 5L),
                   list(minimum_qc_fraction = 0.4), list(method = "linear"),
                   list(maximum_fold_change = 0.5))) {
    invalid <- recipe
    invalid$preprocessing$qc_drift[names(bad)] <- bad
    expect_error(validate_recipe_config(invalid), "qc_drift|QC drift")
  }
})

local({
  dataset <- fixture()
  recipe <- base_recipe()
  recipe$preprocessing$qc_rsd_filter <- list(enabled = TRUE, threshold_percent = 100)
  processed <- preprocess_dataset(dataset, recipe)
  qc <- processed$diagnostics$qc
  weak <- qc[qc$feature_id == "weak_qc", , drop = FALSE]
  sparse <- qc[qc$feature_id == "sparse_qc", , drop = FALSE]
  stopifnot(weak$qc_observed == 3L, weak$qc_support_fraction == 0.5,
            weak$insufficient_qc_support, weak$removed_by_qc_rsd,
            is.na(weak$qc_rsd_percent),
            sparse$qc_observed == 1L, sparse$removed_by_qc_rsd,
            !"weak_qc" %in% colnames(processed$matrices$scaled))
  recipe$preprocessing$qc_rsd_filter$minimum_qc_fraction <- 0.5
  less_strict <- preprocess_dataset(dataset, recipe)
  stopifnot("weak_qc" %in% colnames(less_strict$matrices$scaled))
  recipe$preprocessing$qc_rsd_filter$minimum_qc_fraction <- 0.4
  expect_error(validate_recipe_config(recipe), "QC RSD minimum_qc_fraction")

  directory <- tempfile("qc-preprocess-export-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  export_results(dataset, processed, list(), list(effective_recipe = recipe), directory)
  detection <- utils::read.delim(file.path(directory, "tables", "feature_detection.tsv"))
  drift <- utils::read.delim(file.path(directory, "tables", "qc_drift.tsv"))
  stopifnot(any(detection$feature_id == "undetected" & detection$excluded_no_biological_detection),
            identical(unique(drift$status), "disabled"),
            all(c("qc_observed", "qc_order_spearman_before", "qc_rsd_after_percent") %in% names(drift)))
})
