repo_root <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- if (file.exists("src/v2/runner.R")) "." else ".."
repo_root <- normalizePath(repo_root, mustWork = TRUE)
source(file.path(repo_root, "src", "v2", "utils.R"))
source(file.path(repo_root, "src", "v2", "runner.R"))
source_v2_modules(repo_root)

expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

local({
  work <- tempfile("mapp-contract-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input_a <- file.path(work, "batch_a")
  input_b <- file.path(work, "batch_b")
  dir.create(input_a)
  dir.create(input_b)
  samples <- paste0("sample", seq_len(8))
  metadata <- data.frame(
    filename = samples, sample_id = samples,
    sample_type = c(rep("sample", 6), "QC", "blank"),
    group = c(rep("A", 3), rep("B", 3), NA, NA)
  )
  quant <- data.frame(`row ID` = c("f1", "f2", "f3"), extra = c("x", "y", "z"), check.names = FALSE)
  height <- rbind(c(2, 3, 2, 20, 22, 24, 10, 0),
                  c(12, 11, 13, 11, 12, 10, 11, 10),
                  c(6, 8, 9, 7, 8, 10, 8, 0))
  for (index in seq_along(samples)) {
    quant[[paste0(samples[index], " Peak height")]] <- height[, index]
    quant[[paste0(samples[index], " Peak area")]] <- height[, index] * 10
  }
  for (path in c(input_a, input_b)) {
    utils::write.table(metadata, file.path(path, "metadata.tsv"), sep = "\t", row.names = FALSE, quote = TRUE, na = "")
    utils::write.csv(quant, file.path(path, "quant.csv"), row.names = FALSE, na = "")
  }
  input <- function(path, measure = NULL) {
    list(id = "example/batch", batch_dir = path,
         metadata = file.path(path, "metadata.tsv"), quantification = file.path(path, "quant.csv"),
         feature_id_column = "row ID", intensity_measure = measure)
  }
  expect_error(read_mzmine_dataset(input(input_a)), "intensity_measure")
  dataset <- read_mzmine_dataset(input(input_a, "height"))
  area <- read_mzmine_dataset(input(input_a, "area"))
  stopifnot(identical(rownames(dataset$raw), samples),
            identical(colnames(dataset$raw), c("f1", "f2", "f3")),
            identical(unname(area$raw), unname(dataset$raw * 10)))

  recipe <- list(seed = 2026L,
    roles = list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank")),
    preprocessing = list(blank_filter = list(enabled = FALSE), qc_rsd_filter = list(enabled = FALSE),
      missing_values = list(method = "half_minimum"), normalization = list(method = "none"),
      transformation = list(method = "log2"), scaling = list(method = "pareto")),
    design = list(group = "group", contrasts = list(list(name = "B_vs_A", numerator = "B", denominator = "A"))),
    analyses = list(pca = list(enabled = TRUE, components = 2L),
      pcoa = list(enabled = TRUE, stage = "normalized", distance = "bray", components = 2L),
      omnibus = list(enabled = TRUE), differential = list(enabled = TRUE), plsda = list(enabled = FALSE)))
  config <- function(dataset, output_root) list(dataset = dataset, recipe = recipe,
    paths = list(output_root = output_root), repo_root = repo_root)
  first <- build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))
  relocated <- build_run_identity(config(input(input_b, "height"), file.path(work, "elsewhere")))
  alternative <- build_run_identity(config(input(input_a, "area"), file.path(work, "runs")))
  stopifnot(identical(first$run_hash, relocated$run_hash),
            !identical(first$dataset_hash, alternative$dataset_hash))
  recipe$design$contrasts[[1]]$denominator <- "absent"
  expect_error(validate_dataset(dataset, recipe), "Contrast levels")
  recipe$design$contrasts[[1]]$denominator <- "A"
  processed <- preprocess_dataset(dataset, recipe)
  differential <- run_analyses(processed, recipe)$differential
  effect <- differential$effect[differential$feature_id == "f1"]
  stopifnot(is.finite(effect), effect > 2,
            all(differential$q_value >= differential$p_value, na.rm = TRUE))

  run_config <- config(input(input_a, "height"), file.path(work, "runs"))
  output <- execute_v2_run(run_config)
  stopifnot(completed_run(output, first$run_hash),
            identical(execute_v2_run(run_config), output))
  writeLines("tampered", file.path(output, "tables", "differential.tsv"))
  expect_error(completed_run(output, first$run_hash), "corrupt")
})
cat("MAPP pipeline contract checks passed.\n")
