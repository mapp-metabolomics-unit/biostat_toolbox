repo_root <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- if (file.exists("src/biostats/runner.R")) "." else ".."
repo_root <- normalizePath(repo_root, mustWork = TRUE)
source(file.path(repo_root, "src", "biostats", "utils.R"))
source(file.path(repo_root, "src", "biostats", "runner.R"))
source_biostats_modules(repo_root)

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

local({
  design <- list(group = "group")
  expected <- c("B_vs_A", "C_vs_A", "D_vs_A", "E_vs_A", "F_vs_A",
                "C_vs_B", "D_vs_B", "E_vs_B", "F_vs_B",
                "D_vs_C", "E_vs_C", "F_vs_C", "E_vs_D", "F_vs_D", "F_vs_E")
  for (selection in list(NULL, "all")) {
    design$contrasts <- selection
    expanded <- expand_biostats_contrasts(design, c("F", "C", "A", "E", "B", "D"))
    stopifnot(identical(vapply(expanded, `[[`, character(1), "name"), expected),
              identical(vapply(expanded, `[[`, character(1), "numerator"),
                        sub("_vs_.*", "", expected)),
              identical(vapply(expanded, `[[`, character(1), "denominator"),
                        sub(".*_vs_", "", expected)),
              length(unique(vapply(expanded, function(item)
                paste(sort(c(item$numerator, item$denominator)), collapse = "/"), character(1)))) == 15L)
  }
  stopifnot(identical(expand_biostats_contrasts(list(contrasts = "all"), c("WT_O", "WT_N"))[[1]]$name,
                      "WT_O_vs_WT_N"))
  design$contrasts <- list(list(name = "custom_label", numerator = "F", denominator = "A"))
  stopifnot(identical(expand_biostats_contrasts(design, LETTERS[1:6])[[1]]$name, "custom_label"))
  expect_error(expand_biostats_contrasts(design, LETTERS[1:5]), "Contrast levels")
  for (invalid in list("ALL", "none", TRUE, list(all = TRUE), list("all"))) {
    design$contrasts <- invalid
    expect_error(expand_biostats_contrasts(design, LETTERS[1:6]), "design.contrasts|contrast")
  }
})


local({
  work <- tempfile("mapp-study-")
  batch <- file.path(work, "research-group", "docs", "mapp_project_00085", "mapp_batch_00270")
  dir.create(batch, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  samples <- paste0("sample", seq_len(18))
  groups <- rep(c("WT", "mutant", "other"), each = 6)
  parts <- rep(c("leaf", "root"), times = 9)
  metadata <- data.frame(
    filename = c(samples, "qc1", "qc2", "qc3", "blank1", "blank2"),
    sample_id = c(samples, "qc1", "qc2", "qc3", "blank1", "blank2"),
    sample_type = c(rep("sample", 18), rep("QC", 3), rep("blank", 2)),
    group = c(groups, rep(NA_character_, 5)),
    part = c(parts, rep(NA_character_, 5))
  )
  quant <- data.frame(`row ID` = paste0("f", seq_len(6)), check.names = FALSE)
  for (index in seq_along(samples)) {
    group <- match(groups[[index]], c("WT", "mutant", "other"))
    replicate <- (index - 1L) %% 6L
    quant[[paste0(samples[[index]], " Peak height")]] <-
      c(10 + 14 * (group == 2L) + replicate, 8 + 9 * (group == 3L) + replicate,
        14 + group + replicate %% 3L, 5 + 2 * group + replicate,
        17 + group * replicate, 8 + 2 * group + (replicate + group) %% 4L)
  }
  for (index in seq_len(3))
    quant[[paste0("qc", index, " Peak height")]] <- c(20, 13, 16, 10, 19, 12) + index
  for (index in seq_len(2))
    quant[[paste0("blank", index, " Peak height")]] <- rep(0, nrow(quant))
  utils::write.table(metadata, file.path(batch, "metadata.tsv"), sep = "\t",
                     row.names = FALSE, quote = TRUE, na = "")
  utils::write.csv(quant, file.path(batch, "quant.csv"), row.names = FALSE, na = "")
  dataset_path <- file.path(batch, "dataset.yaml")
  study_path <- file.path(batch, "study.yaml")
  comparison_path <- file.path(batch, "comparison.yaml")
  write_yaml_file(list(id = "mapp_project_00085/mapp_batch_00270", batch_dir = ".",
                       metadata = "metadata.tsv", quantification = "quant.csv",
                       feature_id_column = "row ID", intensity_measure = "height"),
                  dataset_path)
  study_recipe <- list(kind = "study", seed = 1L, display_group = "sample_type",
    preprocessing = list(blank_filter = list(enabled = TRUE,
                           minimum_sample_to_blank_ratio = 5, minimum_blank_prevalence = 0.5),
                         qc_rsd_filter = list(enabled = FALSE),
                         missing_values = list(method = "half_minimum"),
                         normalization = list(method = "none"),
                         transformation = list(method = "log2"),
                         scaling = list(method = "pareto")),
    analyses = list(pca = list(enabled = TRUE, components = 3L, include_qc = FALSE,
                              include_blanks = FALSE),
                    pcoa = list(enabled = TRUE, components = 3L, stage = "normalized",
                                distance = "bray", include_qc = FALSE, include_blanks = FALSE)))
  write_yaml_file(study_recipe, study_path)
  study_config <- resolve_biostats_config(dataset_path, study_path, repo_root)
  stopifnot(validate_biostats_run(study_config)$validation$valid)
  study_dir <- execute_biostats_run(study_config)
  study_hash <- basename(study_dir)
  study_manifest <- read_yaml_file(file.path(study_dir, "manifest.yaml"))
  study_scores <- sha256_file(file.path(study_dir, "tables", "pca_scores.tsv"))
  stopifnot(identical(study_manifest$kind, "study"),
            identical(study_manifest$validation$role_counts$sample, 18L),
            completed_run(study_dir, study_hash),
            identical(execute_biostats_run(study_config), study_dir),
            !file.exists(file.path(study_dir, "tables", "differential.tsv")),
            nrow(utils::read.delim(file.path(study_dir, "tables", "pca_scores.tsv"))) == 18L)

  comparison <- list(kind = "comparison", study_hash = study_hash, seed = 1L,
    sample_plan = list(group_column = "group", groups = list(
      list(source = "WT", target = "WT"),
      list(source = "mutant", target = "mutant"))),
    design = list(group = "analysis_group", contrasts = list(
      list(name = "mutant_vs_WT", numerator = "mutant", denominator = "WT"))),
    analyses = list(omnibus = list(enabled = TRUE), differential = list(enabled = TRUE),
                    plsda = list(enabled = FALSE)))
  write_yaml_file(comparison, comparison_path)
  config <- resolve_biostats_config(dataset_path, comparison_path, repo_root)
  checked <- validate_biostats_run(config)
  stopifnot(checked$validation$valid,
            identical(checked$identity$preprocess_hash, study_manifest$hashes$preprocess_hash),
            identical(as.integer(unlist(checked$validation$group_counts)), c(6L, 6L)))
  child_dir <- execute_biostats_run(config)
  first_child <- basename(child_dir)
  first_manifest <- read_yaml_file(file.path(child_dir, "manifest.yaml"))
  selected <- utils::read.delim(file.path(child_dir, "tables", "comparison_sample_metadata.tsv"))
  effects <- utils::read.delim(file.path(child_dir, "tables", "differential.tsv"))
  stopifnot(identical(first_manifest$kind, "comparison"),
            identical(first_manifest$study_hash, study_hash),
            identical(first_manifest$effective_recipe$preprocessing,
                      study_manifest$effective_recipe$preprocessing),
            identical(nrow(selected), 12L), setequal(unique(selected$analysis_group), c("WT", "mutant")),
            identical(unique(effects$contrast), "mutant_vs_WT"),
            all(is.finite(effects$effect)),
            effects$effect[effects$feature_id == "f1"] > 0,
            !file.exists(file.path(child_dir, "tables", "pca_scores.tsv")),
            identical(sha256_file(file.path(study_dir, "tables", "pca_scores.tsv")), study_scores),
            completed_run(child_dir, first_child),
            identical(execute_biostats_run(config), child_dir))
  label_only <- config
  label_only$recipe$presentation <- list(description = "Same comparison with a new display label")
  stopifnot(identical(build_run_identity(label_only)$run_hash, first_child))
  comparison$sample_plan <- list(
    group_column = "group", second_column = "part", groups = list(
      list(source = '["WT","leaf"]', target = "WT"),
      list(source = '["mutant","leaf"]', target = "mutant")))
  write_yaml_file(comparison, comparison_path)
  leaf_config <- resolve_biostats_config(dataset_path, comparison_path, repo_root)
  leaf <- validate_biostats_run(leaf_config)
  stopifnot(leaf$validation$valid,
            identical(leaf$identity$preprocess_hash, checked$identity$preprocess_hash),
            !identical(leaf$identity$run_hash, first_child),
            identical(as.integer(unlist(leaf$validation$group_counts)), c(3L, 3L)))
  leaf_dir <- execute_biostats_run(leaf_config)
  leaf_metadata <- utils::read.delim(file.path(leaf_dir, "tables", "comparison_sample_metadata.tsv"))
  stopifnot(nrow(leaf_metadata) == 6L, all(leaf_metadata$part == "leaf"),
            identical(sha256_file(file.path(study_dir, "tables", "pca_scores.tsv")), study_scores),
            completed_run(study_dir, study_hash))
  invalid <- comparison
  invalid$sample_plan$groups[[2L]]$source <- '["WT","leaf"]'
  expect_error(validate_sample_plan(invalid$sample_plan), "unique")
  bad_dataset <- read_yaml_file(dataset_path)
  bad_dataset$intensity_measure <- "area"
  bad_dataset_path <- file.path(batch, "bad_dataset.yaml")
  write_yaml_file(bad_dataset, bad_dataset_path)
  expect_error(validate_biostats_run(resolve_biostats_config(bad_dataset_path, comparison_path, repo_root)),
               "differs from its sealed study")
  published <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(repo_root, "scripts", "publish_stats_dashboard.R")),
      "--batch", shQuote(batch)), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(published, "status"))) stop(paste(published, collapse = "\n"))
  catalogue <- jsonlite::read_json(file.path(batch, "results", "stats", "dashboard", "runs.json"),
                                   simplifyVector = FALSE)
  entries <- catalogue$runs
  stopifnot(length(entries) == 3L,
            sum(vapply(entries, function(item) identical(item$kind, "study"), logical(1))) == 1L,
            sum(vapply(entries, function(item) identical(item$kind, "comparison"), logical(1))) == 2L)
  child <- Filter(function(item) identical(item$hash, first_child), entries)[[1L]]
  parent <- Filter(function(item) identical(item$hash, study_hash), entries)[[1L]]
  stopifnot(identical(parent$hashes$dataset_hash, study_manifest$hashes$dataset_hash),
            identical(parent$hashes$preprocess_hash, study_manifest$hashes$preprocess_hash),
            identical(child$hashes$dataset_hash, parent$hashes$dataset_hash),
            identical(child$hashes$preprocess_hash, parent$hashes$preprocess_hash))
  stopifnot(identical(child$study_hash, study_hash),
            identical(child$outputs[["tables/pca_scores.tsv"]], study_manifest$outputs[["tables/pca_scores.tsv"]]),
            identical(child$sources[["tables/differential.tsv"]],
                      paste("comparisons", study_hash, first_child, sep = "/")),
            is.null(child$sources[["tables/pca_scores.tsv"]]))
  writeLines("tampered", file.path(child_dir, "tables", "differential.tsv"))
  expect_error(completed_run(child_dir, first_child), "corrupt")
  stopifnot(completed_run(study_dir, study_hash))
})
local({
  set.seed(42)
  values <- matrix(rexp(40) + 1, 8, 5,
                   dimnames = list(paste0("s", seq_len(8)), paste0("f", seq_len(5))))
  values[5:8, 1:2] <- values[5:8, 1:2] + 4
  samples <- data.frame(sample_type = rep("sample", 8), group = rep(c("A", "B"), each = 4),
                        row.names = rownames(values))
  ordination <- run_pca(values, samples, "group", 3L)
  stopifnot(identical(ordination$scores$sample, rownames(values)),
            identical(ordination$loadings$feature_id, colnames(values)),
            identical(names(ordination$scores)[2:4], paste0("PC", 1:3)),
            identical(names(ordination$loadings)[2:4], paste0("PC", 1:3)),
            all(ordination$variance$variance_percent > 0),
            sum(ordination$variance$variance_percent) <= 100.0001)
  supervised_recipe <- list(seed = 42L,
    roles = list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank")),
    preprocessing = list(blank_filter = list(enabled = FALSE), qc_rsd_filter = list(enabled = FALSE),
      missing_values = list(method = "half_minimum"), normalization = list(method = "none"),
      transformation = list(method = "log2"), scaling = list(method = "pareto")),
    design = list(group = "group"),
    analyses = list(plsda = list(enabled = TRUE, components = 1L, folds = 2L,
                                repeats = 2L, permutations = 19L)))
  progress_steps <- numeric()
  progress_totals <- numeric()
  progress_phases <- character()
  model <- run_supervised(list(raw = values, sample_metadata = samples), supervised_recipe,
    progress = function(phase, completed, total) {
      progress_steps <<- c(progress_steps, completed)
      progress_totals <<- c(progress_totals, total)
      progress_phases <<- c(progress_phases, phase)
    })
  stopifnot(identical(model$status, "validated"), model$summary$groups == 2L,
            nrow(model$scores) == 8L, identical(model$scores$sample, rownames(values)),
            nrow(model$fold_predictions) == 16L,
            all(table(model$fold_predictions$sample) == 2L),
            all(model$fold_predictions$predicted %in% c("A", "B")),
            model$summary$balanced_accuracy > 0.5,
            nrow(model$permutations) == 19L,
            model$summary$permutation_p_value >= 1/20,
            model$summary$permutation_p_value <= 1,
            identical(model$loadings$feature_id, colnames(values)),
            identical(names(model$loadings)[-1L], paste0("LV", seq_len(model$summary$full_fit_components))),
            all(vapply(model$loadings[-1L], function(axis) is.numeric(axis) && all(is.finite(axis)), logical(1))),
            identical(model$vip$feature_id, model$loadings$feature_id),
            identical(names(model$vip), c("feature_id", "vip")),
            is.numeric(model$vip$vip), all(is.finite(model$vip$vip) & model$vip$vip >= 0))
  expected_fits <- model$summary$folds * model$summary$repeats *
    (model$summary$permutations + 1L) + 1L
  permutation_progress <- startsWith(progress_phases, "Permutation ")
  stopifnot(all(progress_totals == expected_fits),
            identical(progress_steps[1L], 0),
            identical(progress_steps[length(progress_steps)], as.double(expected_fits)),
            all(diff(progress_steps) >= 0),
            sum(permutation_progress) == model$summary$folds * model$summary$repeats *
              model$summary$permutations,
            length(unique(progress_phases[permutation_progress])) == model$summary$permutations,
            identical(progress_phases[length(progress_phases)], "PLS-DA complete"))
})

cat("MAPP pipeline contract checks passed.\n")
