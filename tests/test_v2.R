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
  design <- list(group = "group")
  expected <- c("B_vs_A", "C_vs_A", "D_vs_A", "E_vs_A", "F_vs_A",
                "C_vs_B", "D_vs_B", "E_vs_B", "F_vs_B",
                "D_vs_C", "E_vs_C", "F_vs_C", "E_vs_D", "F_vs_D", "F_vs_E")
  for (selection in list(NULL, "all")) {
    design$contrasts <- selection
    expanded <- expand_v2_contrasts(design, c("F", "C", "A", "E", "B", "D"))
    stopifnot(identical(vapply(expanded, `[[`, character(1), "name"), expected),
              identical(vapply(expanded, `[[`, character(1), "numerator"),
                        sub("_vs_.*", "", expected)),
              identical(vapply(expanded, `[[`, character(1), "denominator"),
                        sub(".*_vs_", "", expected)),
              length(unique(vapply(expanded, function(item)
                paste(sort(c(item$numerator, item$denominator)), collapse = "/"), character(1)))) == 15L)
  }
  stopifnot(identical(expand_v2_contrasts(list(contrasts = "all"), c("WT_O", "WT_N"))[[1]]$name,
                      "WT_O_vs_WT_N"))
  design$contrasts <- list(list(name = "custom_label", numerator = "F", denominator = "A"))
  stopifnot(identical(expand_v2_contrasts(design, LETTERS[1:6])[[1]]$name, "custom_label"))
  expect_error(expand_v2_contrasts(design, LETTERS[1:5]), "Contrast levels")
  for (invalid in list("ALL", "none", TRUE, list(all = TRUE), list("all"))) {
    design$contrasts <- invalid
    expect_error(expand_v2_contrasts(design, LETTERS[1:6]), "design.contrasts|contrast")
  }
})


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
  horizontal <- data.frame(feature_id = c("f1", "f2", "f3"),
                           sources_number_IK2D = c(0, 2, NA),
                           sources_IK2D = c("", "gnps & sirius", ""),
                           gnps_IK2D = c("", "AAAAAAAAAAAAAA", ""),
                           sirius_IK2D = c("", "AAAAAAAAAAAAAA", ""),
                           isdb_IK2D = c("", "BBBBBBBBBBBBBB", ""),
                           isdb_libname = c("", "ISDB", ""),
                           isdb_score_taxo = c(NA, 7, NA),
                           isdb_organism_name = c("", "Ensifer meliloti", ""),
                           gnps_SMILES = c("", "CCO", ""),
                           sirius_smiles = c("", "CCN", ""), check.names = FALSE)
  for (path in c(input_a, input_b)) {
    utils::write.table(metadata, file.path(path, "metadata.tsv"), sep = "\t", row.names = FALSE, quote = TRUE, na = "")
    utils::write.csv(quant, file.path(path, "quant.csv"), row.names = FALSE, na = "")
    utils::write.table(horizontal, file.path(path, "horizontal.tsv"), sep = "\t",
                       row.names = FALSE, quote = TRUE, na = "")
  }
  input <- function(path, measure = NULL) {
    list(id = "example/batch", batch_dir = path,
         metadata = file.path(path, "metadata.tsv"), quantification = file.path(path, "quant.csv"),
         feature_id_column = "row ID", intensity_measure = measure,
         annotations = if (file.exists(file.path(path, "horizontal.tsv")))
           list(horizontal = file.path(path, "horizontal.tsv")) else NULL)
  }
  expect_error(read_mzmine_dataset(input(input_a)), "intensity_measure")
  dataset <- read_mzmine_dataset(input(input_a, "height"))
  area <- read_mzmine_dataset(input(input_a, "area"))
  stopifnot(identical(rownames(dataset$raw), samples),
            identical(colnames(dataset$raw), c("f1", "f2", "f3")),
            identical(unname(area$raw), unname(dataset$raw * 10)))
  summary <- horizontal_annotation_summary(dataset$annotations$horizontal)
  stopifnot(identical(as.character(summary$feature_id), c("f1", "f2", "f3")),
            identical(as.numeric(summary$sources_number_IK2D), c(0, 2, NA_real_)),
            identical(summary$gnps_IK2D[2], summary$sirius_IK2D[2]),
            !identical(summary$gnps_IK2D[2], summary$isdb_IK2D[2]),
            identical(as.numeric(summary$isdb_score_taxo[2]), 7),
            identical(summary$isdb_libname[2], "ISDB"),
            identical(summary$isdb_organism_name[2], "Ensifer meliloti"),
            identical(summary$gnps_SMILES[2], "CCO"),
            identical(summary$sirius_smiles[2], "CCN"))
  scored_dataset <- dataset
  scored_dataset$annotations$met_annot_enhancer <- data.frame(
    feature_id = "f2", structure_smiles = "CCO", structure_inchikey = "AAAAAAAAAAAAAA-UHFFFAOYSA-N",
    score_taxo = 7, organism_name = "Ensifer meliloti")
  scored <- normalize_annotations(scored_dataset)
  stopifnot(identical(scored$taxonomic_score[which(scored$source == "met_annot_enhancer")], 7),
            all(is.na(scored$taxonomic_score[which(!is.na(scored$source) &
                                                   scored$source != "met_annot_enhancer")])))
  duplicated_horizontal <- dataset$annotations$horizontal
  duplicated_horizontal$feature_id[2] <- "f1"
  stopifnot(identical(as.character(horizontal_annotation_summary(duplicated_horizontal)$feature_id),
                      c("f1", "f1", "f3")))
  duplicated_horizontal$feature_id[2] <- ""
  expect_error(horizontal_annotation_summary(duplicated_horizontal), "non-empty")
  duplicated_horizontal$feature_id[2] <- "f2"
  duplicated_horizontal$sources_number_IK2D[2] <- 1.5
  expect_error(horizontal_annotation_summary(duplicated_horizontal), "whole counts")
  mapped_metadata <- metadata
  names(mapped_metadata)[names(mapped_metadata) == "sample_id"] <- "mapp_sample_id"
  names(mapped_metadata)[names(mapped_metadata) == "sample_type"] <- "ATTRIBUTE_sample_type"
  utils::write.table(mapped_metadata, file.path(input_b, "mapped.tsv"), sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  mapped_input <- input(input_b, "height")
  mapped_input$metadata <- file.path(input_b, "mapped.tsv")
  mapped_input$metadata_columns <- list(sample_id = "mapp_sample_id", sample_type = "ATTRIBUTE_sample_type")
  mapped <- read_mzmine_dataset(mapped_input)
  stopifnot(identical(mapped$sample_metadata$sample_type, dataset$sample_metadata$sample_type),
            identical(mapped$sample_metadata$sample_id, dataset$sample_metadata$sample_id))
  expect_error(read_mzmine_dataset(within(mapped_input, metadata_columns <- list(sample_type = "missing"))), "Metadata lacks")

  recipe <- list(seed = 2026L,
    roles = list(column = "sample_type", levels = list(sample = "sample", qc = "QC", blank = "blank")),
    preprocessing = list(blank_filter = list(enabled = FALSE), qc_rsd_filter = list(enabled = FALSE),
      missing_values = list(method = "half_minimum"), normalization = list(method = "none"),
      transformation = list(method = "log2"), scaling = list(method = "pareto")),
    design = list(group = "group", contrasts = list(list(name = "custom_label", numerator = "B", denominator = "A"))),
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
  pinned <- input(input_a, "height")
  pinned$input_sha256 <- list(metadata = sha256_file(pinned$metadata),
                              quantification = sha256_file(pinned$quantification))
  stopifnot(identical(build_run_identity(config(pinned, file.path(work, "runs")))$run_hash,
                      first$run_hash))
  dataset_yaml <- file.path(work, "dataset.yaml")
  recipe_yaml <- file.path(work, "recipe.yaml")
  authored_recipe <- recipe
  authored_recipe$roles <- NULL
  write_yaml_file(input(input_a, "height"), dataset_yaml)
  write_yaml_file(authored_recipe, recipe_yaml)
  static <- resolve_v2_config(dataset_yaml, recipe_yaml, repo_root)
  stopifnot(identical(static$recipe$roles, default_v2_roles()),
            validate_dataset(read_mzmine_dataset(static$dataset), static$recipe)$valid)
  write_yaml_file(recipe, recipe_yaml)
  expect_error(resolve_v2_config(dataset_yaml, recipe_yaml, repo_root), "recipe.roles is obsolete")
  write_yaml_file(authored_recipe, recipe_yaml)
  obsolete_dataset <- input(input_a, "height")
  obsolete_dataset$sample_plan <- list(group_column = "group", groups = list(list(source = "A", target = "WT")))
  write_yaml_file(obsolete_dataset, dataset_yaml)
  expect_error(resolve_v2_config(dataset_yaml, recipe_yaml, repo_root), "dataset.sample_plan is obsolete")

  role_metadata <- mapped_metadata
  role_metadata$ATTRIBUTE_sample_type <- c(rep("bio", 6), "pool", "solvent")
  roles_path <- file.path(input_b, "roles.tsv")
  utils::write.table(role_metadata, roles_path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  shared <- mapped_input
  shared$metadata <- roles_path
  shared$roles <- list(column = "ATTRIBUTE_sample_type",
                       levels = list(sample = "bio", qc = "pool", blank = "solvent"))
  authored_recipe$design <- list(group = "analysis_group",
    contrasts = list(list(name = "mutant_vs_WT", numerator = "mutant", denominator = "WT")))
  authored_recipe$sample_plan <- list(group_column = "group",
    groups = list(list(source = "A", target = "WT"), list(source = "B", target = "mutant")))
  write_yaml_file(shared, dataset_yaml)
  write_yaml_file(authored_recipe, recipe_yaml)
  merged <- resolve_v2_config(dataset_yaml, recipe_yaml, repo_root)
  mapped_groups <- read_mzmine_dataset(merged$dataset, merged$recipe$sample_plan)
  stopifnot(identical(merged$recipe$roles$column, "sample_type"),
            identical(merged$recipe$roles$levels, shared$roles$levels),
            identical(mapped_groups$sample_metadata$analysis_group, c(rep("WT", 3), rep("mutant", 3), NA_character_, NA_character_)),
            identical(rownames(mapped_groups$raw), samples),
            identical(as.integer(unlist(validate_dataset(mapped_groups, merged$recipe)$group_counts)), c(3L, 3L)))
  first_mapping_hash <- build_run_identity(merged)
  retargeted_recipe <- merged$recipe
  retargeted_recipe$sample_plan$groups[[1]]$target <- "otherWT"
  retargeted_hash <- build_run_identity(within(merged, recipe <- retargeted_recipe))
  stopifnot(identical(first_mapping_hash$dataset_hash, retargeted_hash$dataset_hash),
            !identical(first_mapping_hash$preprocess_hash, retargeted_hash$preprocess_hash),
            !identical(first_mapping_hash$run_hash, retargeted_hash$run_hash))
  other_recipe <- authored_recipe
  other_recipe$design$contrasts <- list(list(name = "treated_vs_control",
                                             numerator = "treated", denominator = "control"))
  other_recipe$sample_plan$groups[[1]]$target <- "control"
  other_recipe$sample_plan$groups[[2]]$target <- "treated"
  write_yaml_file(other_recipe, recipe_yaml)
  other_mapping <- resolve_v2_config(dataset_yaml, recipe_yaml, repo_root)
  other_groups <- read_mzmine_dataset(other_mapping$dataset, other_mapping$recipe$sample_plan)
  other_hash <- build_run_identity(other_mapping)
  stopifnot(identical(other_groups$sample_metadata$analysis_group,
                      c(rep("control", 3), rep("treated", 3), NA_character_, NA_character_)),
            identical(first_mapping_hash$dataset_hash, other_hash$dataset_hash),
            !identical(first_mapping_hash$run_hash, other_hash$run_hash))
  unmapped_roles <- shared
  unmapped_roles$roles$levels$sample <- "absent"
  expect_error(read_mzmine_dataset(unmapped_roles, authored_recipe$sample_plan), "Unknown source sample roles")
  invalid_plan <- authored_recipe
  invalid_plan$sample_plan$role_column <- "sample_type"
  write_yaml_file(invalid_plan, recipe_yaml)
  expect_error(resolve_v2_config(dataset_yaml, recipe_yaml, repo_root), "Unsupported recipe.sample_plan key")
  guard <- shared
  guard$input_sha256 <- list(metadata = sha256_file(guard$metadata),
                             quantification = sha256_file(guard$quantification))
  stopifnot(identical(first_mapping_hash$run_hash,
                      build_run_identity(within(merged, dataset <- guard))$run_hash))
  changed <- guard
  changed$metadata <- file.path(input_b, "changed_roles.tsv")
  changed_metadata <- role_metadata
  changed_metadata$sample_id[1] <- "renamed"
  utils::write.table(changed_metadata, changed$metadata, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
  expect_error(read_mzmine_dataset(changed, authored_recipe$sample_plan), "Canonical metadata SHA256 mismatch")
  expect_error(build_run_identity(within(merged, dataset <- changed)), "Canonical metadata SHA256 mismatch")
  changed$input_sha256 <- NULL
  stopifnot(!identical(first_mapping_hash$run_hash,
                       build_run_identity(within(merged, dataset <- changed))$run_hash))
  recipe$design$contrasts[[1]]$denominator <- "absent"
  expect_error(validate_dataset(dataset, recipe), "Contrast levels")
  recipe$design$contrasts[[1]]$denominator <- "A"
  explicit_hash <- build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))$run_hash
  recipe$design$contrasts <- "all"
  validate_recipe_config(recipe)
  stopifnot(validate_dataset(dataset, recipe)$valid,
            !identical(build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))$run_hash,
                       explicit_hash))
  all_hash <- build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))$run_hash
  recipe$design$contrasts <- NULL
  validate_recipe_config(recipe)
  stopifnot(validate_dataset(dataset, recipe)$valid,
            !identical(build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))$run_hash,
                       explicit_hash))
  stopifnot(identical(build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))$run_hash,
                      all_hash))
  recipe$design$contrasts <- list()
  expect_error(validate_recipe_config(recipe), "nonempty design.contrasts")
  recipe$design$contrasts <- list(list(name = "custom_label", numerator = "B", denominator = "A"))
  processed <- preprocess_dataset(dataset, recipe)
  differential <- run_analyses(processed, recipe)$differential
  effect <- differential$effect[differential$feature_id == "f1"]
  stopifnot(is.finite(effect), effect > 2,
            identical(unique(differential$contrast), "custom_label"),
            all(differential$q_value >= differential$p_value, na.rm = TRUE))
  recipe$presentation <- list(description = "Biological groups A and B")
  stopifnot(identical(build_run_identity(config(input(input_a, "height"), file.path(work, "runs")))$run_hash,
                      first$run_hash))

  run_config <- config(input(input_a, "height"), file.path(work, "runs"))
  output <- execute_v2_run(run_config)
  stopifnot(completed_run(output, first$run_hash),
            identical(execute_v2_run(run_config), output))
  sealed <- read_yaml_file(file.path(output, "manifest.yaml"))
  stopifnot(identical(sealed$description, recipe$presentation$description),
            is.null(sealed$effective_recipe$presentation),
            identical(as.integer(sealed$validation$dimensions$samples), 8L))
  exported_metadata <- utils::read.delim(file.path(output, "tables", "sample_metadata.tsv"),
                                          check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(exported_metadata) == 8L,
            identical(as.integer(table(exported_metadata$sample_type)[c("sample", "QC", "blank")]), c(6L, 1L, 1L)))
  exported_horizontal <- utils::read.delim(file.path(output, "tables", "annotation_horizontal.tsv"),
                                            check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(as.numeric(exported_horizontal$sources_number_IK2D), c(0, 2, NA_real_)),
            identical(exported_horizontal$gnps_IK2D[2], exported_horizontal$sirius_IK2D[2]),
            identical(as.numeric(exported_horizontal$isdb_score_taxo[2]), 7),
            identical(exported_horizontal$isdb_libname[2], "ISDB"),
            identical(exported_horizontal$gnps_SMILES[2], "CCO"),
            identical(exported_horizontal$sirius_smiles[2], "CCN"))
  explicit_table <- utils::read.delim(file.path(output, "tables", "differential.tsv"),
                                      check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(unique(explicit_table$contrast), "custom_label"),
            nrow(explicit_table) == nrow(quant))
  writeLines("tampered", file.path(output, "tables", "differential.tsv"))
  filtered_config <- config(input(input_a, "height"), file.path(work, "filtered_runs"))
  filtered_config$recipe$preprocessing$blank_filter <- list(
    enabled = TRUE, minimum_sample_to_blank_ratio = 5, minimum_blank_prevalence = 0.5)
  filtered_run <- execute_v2_run(filtered_config)
  input_features <- utils::read.delim(file.path(filtered_run, "tables", "variable_metadata_input.tsv"),
                                     check.names = FALSE, stringsAsFactors = FALSE)
  retained_features <- utils::read.delim(file.path(filtered_run, "tables", "variable_metadata.tsv"),
                                        check.names = FALSE, stringsAsFactors = FALSE)
  raw_peaks <- utils::read.delim(file.path(filtered_run, "tables", "matrix_raw.tsv"),
                                check.names = FALSE, stringsAsFactors = FALSE)
  blank_peaks <- utils::read.delim(file.path(filtered_run, "tables", "matrix_blank_filtered.tsv"),
                                  check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot("f2" %in% input_features$feature_id,
            identical(input_features$extra[input_features$feature_id == "f2"], "y"),
            !"f2" %in% retained_features$feature_id,
            "f2" %in% names(raw_peaks), any(raw_peaks$f2 > 0),
            !"f2" %in% names(blank_peaks))
  expect_error(completed_run(output, first$run_hash), "corrupt")

  six_dir <- file.path(work, "six_groups")
  dir.create(six_dir)
  six_groups <- rep(c("F", "C", "A", "E", "B", "D"), each = 3L)
  six_samples <- paste0("six", seq_along(six_groups))
  six_metadata <- data.frame(filename = c(six_samples, "sixQC", "sixBlank"),
                             sample_id = c(six_samples, "sixQC", "sixBlank"),
                             sample_type = c(rep("sample", 18), "QC", "blank"),
                             group = c(six_groups, NA, NA))
  utils::write.table(six_metadata, file.path(six_dir, "metadata.tsv"), sep = "\t", row.names = FALSE,
                     quote = TRUE, na = "")
  six_quant <- data.frame(`row ID` = c("f1", "f2", "f3"), check.names = FALSE)
  for (index in seq_along(six_samples)) {
    level <- match(six_groups[index], LETTERS)
    replicate <- (index - 1L) %% 3L
    six_quant[[paste0(six_samples[index], " Peak height")]] <-
      c(10 + 3 * level + replicate, 20 + level + replicate^2, 5 + 2 * level + (replicate + level) %% 3L)
  }
  six_quant[["sixQC Peak height"]] <- c(15, 23, 9)
  six_quant[["sixBlank Peak height"]] <- c(0, 0, 0)
  utils::write.csv(six_quant, file.path(six_dir, "quant.csv"), row.names = FALSE)
  recipe$design$contrasts <- "all"
  six_input <- input(six_dir, "height")
  six_dataset <- read_mzmine_dataset(six_input)
  stopifnot(validate_dataset(six_dataset, recipe)$valid)
  all_output <- execute_v2_run(config(six_input, file.path(work, "six_runs")))
  all_table <- utils::read.delim(file.path(all_output, "tables", "differential.tsv"),
                                 check.names = FALSE, stringsAsFactors = FALSE)
  pairs <- unique(all_table[c("contrast", "numerator", "denominator")])
  stopifnot(nrow(pairs) == 15L, nrow(all_table) == 15L * nrow(six_quant),
            identical(pairs$contrast,
                      vapply(expand_v2_contrasts(recipe$design, six_groups), `[[`, character(1), "name")),
            !anyDuplicated(vapply(seq_len(nrow(pairs)), function(index)
              paste(sort(c(pairs$numerator[index], pairs$denominator[index])), collapse = "/"), character(1))),
            all(all_table$q_value >= all_table$p_value, na.rm = TRUE))
  recipe$design$contrasts <- NULL
  omitted_table <- run_analyses(preprocess_dataset(six_dataset, recipe), recipe)$differential
  stopifnot(identical(omitted_table$contrast, all_table$contrast),
            isTRUE(all.equal(omitted_table$effect, all_table$effect)))
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
