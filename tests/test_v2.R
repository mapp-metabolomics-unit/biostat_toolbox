repo_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo_root, "src", "v2", "utils.R"))
source(file.path(repo_root, "src", "v2", "preprocess.R"))
source(file.path(repo_root, "src", "v2", "analysis.R"))

samples <- paste0("s", 1:8)
metadata <- data.frame(
  filename = samples,
  sample_id = samples,
  sample_type = c(rep("sample", 4), rep("QC", 2), rep("blank", 2)),
  group = c("A", "A", "B", "B", NA, NA, NA, NA),
  row.names = samples
)
x <- rbind(
  s1 = c(100, 10, 20), s2 = c(110, 12, 21), s3 = c(400, 11, 19), s4 = c(420, 9, 22),
  s5 = c(250, 10, 20), s6 = c(260, 10, 21), s7 = c(1, 10, 0), s8 = c(1, 11, 0)
)
colnames(x) <- c("signal", "contaminant", "stable")
recipe <- list(
  roles = list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank")),
  preprocessing = list(
    blank_filter = list(enabled = TRUE, minimum_sample_to_blank_ratio = 5, minimum_blank_prevalence = 0.5),
    qc_rsd_filter = list(enabled = FALSE), missing_values = list(method = "half_minimum"),
    normalization = list(method = "none"), transformation = list(method = "log2"), scaling = list(method = "pareto")
  ),
  design = list(group = "group", contrasts = list(list(name = "B_vs_A", numerator = "B", denominator = "A"))),
  analyses = list(pca = list(enabled = TRUE), pcoa = list(enabled = TRUE, stage = "normalized", distance = "bray"), omnibus = list(enabled = TRUE), differential = list(enabled = TRUE))
)
dataset <- list(raw = x, sample_metadata = metadata, variable_metadata = data.frame(feature_id = colnames(x), row.names = colnames(x)))
processed <- preprocess_dataset(dataset, recipe)
stopifnot(!"contaminant" %in% colnames(processed$matrices$scaled))
stopifnot(identical(rownames(processed$matrices$scaled), samples[1:4]))
analyses <- run_analyses(processed, recipe)
effect <- analyses$differential$effect[analyses$differential$feature_id == "signal"]
stopifnot(is.finite(effect), effect > 1)
stopifnot(nrow(analyses$pca$scores) == 4, nrow(analyses$pcoa$scores) == 4)
cat("V2 unit smoke tests passed.\n")
