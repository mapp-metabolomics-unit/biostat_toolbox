repo_root <- Sys.getenv("MAPP_STATS_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) repo_root <- if (file.exists("scripts/materialize_stats_recipe.R")) "." else ".."
repo_root <- normalizePath(repo_root, mustWork = TRUE)
source(file.path(repo_root, "src", "v2", "utils.R"))
source(file.path(repo_root, "src", "v2", "runner.R"))
source_v2_modules(repo_root)

local({
  work <- tempfile("materialize-contract-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  batch <- file.path(work, "project_00001", "batch_00001")
  dir.create(batch, recursive = TRUE)
  input_dir <- file.path(batch, "metadata")
  dir.create(input_dir)
  mzmine_dir <- file.path(batch, "results", "mzmine")
  dir.create(mzmine_dir, recursive = TRUE)
  groups <- rep(c("WT_N", "WT_O", "xdhA1", "xdhA2", "xdhB1", "xdhB2"), each = 5L)
  filenames <- sprintf("sample%02d", seq_len(40L))
  metadata <- data.frame(filename = filenames, sample_id = filenames,
                         sample_type = c(rep("sample", 30L), rep("QC", 7L), rep("blank", 3L)),
                         ATTRIBUTE_mutant = c(groups, rep(NA_character_, 10L)),
                         source_taxon = c(rep("taxon_A", 15L), rep("taxon_B", 15L), rep(NA_character_, 10L)))
  metadata_path <- file.path(input_dir, "original.tsv")
  write_tsv(metadata, metadata_path)
  quant <- data.frame(`row ID` = c("feature1", "feature2", "feature3"),
                      `row m/z` = c(100.1, 200.2, 300.3), check.names = FALSE)
  for (i in seq_along(filenames)) {
    quant[[paste0(filenames[[i]], " Peak height")]] <- c(i, 100 + i, 200 + i)
    quant[[paste0(filenames[[i]], " Peak area")]] <- c(i, 100 + i, 200 + i) * 10
  }
  quant_path <- file.path(mzmine_dir, "original.csv")
  utils::write.csv(quant, quant_path, row.names = FALSE)
  annotation_path <- file.path(batch, "metadata", "annotation.tsv")
  write_tsv(data.frame(feature_id = quant[["row ID"]], category = c("a", "b", "c")), annotation_path)
  id <- "project_00001/batch_00001"
  hash <- paste(rep("a", 64L), collapse = "")
  run_dir <- file.path(batch, "results", "stats_v2", hash)
  dir.create(file.path(run_dir, "tables"), recursive = TRUE)
  write_tsv(metadata, file.path(run_dir, "tables", "sample_metadata.tsv"))
  dataset <- list(id = id, batch_dir = batch, metadata = metadata_path, quantification = quant_path,
                  feature_id_column = "row ID", intensity_measure = "height",
                  annotations = list(reference = annotation_path))
  preprocessing <- list(blank_filter = list(enabled = TRUE, minimum_sample_to_blank_ratio = 5,
                                            minimum_blank_prevalence = 0.5),
                        qc_rsd_filter = list(enabled = TRUE, threshold_percent = 30),
                        missing_values = list(method = "half_minimum"),
                        normalization = list(method = "pqn"), transformation = list(method = "log2"),
                        scaling = list(method = "pareto"))
  analyses <- list(pca = list(enabled = TRUE, components = 3),
                   pcoa = list(enabled = TRUE, components = 3, stage = "normalized", distance = "bray"),
                   omnibus = list(enabled = TRUE), differential = list(enabled = TRUE),
                   plsda = list(enabled = FALSE, components = 2, folds = 3, repeats = 2, permutations = 19))
  manifest <- list(schema_version = "1.0.0", status = "complete", dataset_id = id,
                   hashes = list(run_hash = hash),
                   inputs = list(metadata = sha256_file(metadata_path), quantification = sha256_file(quant_path),
                                 annotations = list(reference = sha256_file(annotation_path))),
                   effective_dataset = dataset,
                   effective_recipe = list(seed = 42, roles = list(column = "sample_type",
                       levels = list(sample = "sample", qc = "QC", blank = "blank")),
                       design = list(group = "ATTRIBUTE_mutant"), preprocessing = preprocessing,
                       analyses = analyses),
                   outputs = list(`tables/sample_metadata.tsv` = sha256_file(file.path(run_dir, "tables", "sample_metadata.tsv"))))
  write_yaml_file(manifest, file.path(run_dir, "manifest.yaml"))
  write_yaml_file(list(run_hash = hash, content_sha256 = sha256_object(run_contents(run_dir))),
                  file.path(run_dir, "COMPLETE"))
  script <- file.path(repo_root, "scripts", "materialize_stats_recipe.R")
  invoke <- function(plan, success) {
    plan_file <- file.path(work, "plan.json")
    write_json(plan, plan_file)
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                c(shQuote(script), "--batch", shQuote(batch), "--plan", shQuote(plan_file)),
                stdout = TRUE, stderr = TRUE))
    if (success && !is.null(attr(output, "status")))
      stop("Materialization of ", plan$slug, " failed: ", paste(output, collapse = "\n"), call. = FALSE)
    if (!success) stopifnot(!is.null(attr(output, "status")))
    output
  }
  mapping <- list(list(source = "WT_N", target = "WT"), list(source = "WT_O", target = "WT"),
                  list(source = "xdhA1", target = "xdhA"), list(source = "xdhA2", target = "xdhA"),
                  list(source = "xdhB1", target = "xdhB"), list(source = "xdhB2", target = "xdhB"))
  plan <- list(schema_version = 1L, source_run_hash = hash, slug = "all_merged",
               group_column = "ATTRIBUTE_mutant", groups = mapping, seed = 42L,
               preprocessing = preprocessing, analyses = analyses, contrasts = "all")
  original_hashes <- c(metadata = sha256_file(metadata_path), quantification = sha256_file(quant_path))
  for (exclude in c(FALSE, TRUE)) {
    next_plan <- plan
    if (exclude) {
      next_plan$slug <- "exclude_xdhb2"
      next_plan$groups <- mapping[-6L]
    } else {
      next_plan$description <- "WT, xdhA, and xdhB merged contrasts"
    }
    invoke(next_plan, TRUE)
    out <- file.path(batch, "recipes", next_plan$slug)
    generated <- resolve_v2_config(file.path(out, "dataset.yaml"), file.path(out, "recipe.yaml"), repo_root)
    parsed <- read_mzmine_dataset(generated$dataset, generated$recipe$sample_plan)
    report <- validate_dataset(parsed, generated$recipe)
    kept <- if (exclude) filenames[-seq.int(26L, 30L)] else filenames
    generated_dataset <- read_yaml_file(file.path(out, "dataset.yaml"))
    generated_recipe <- read_yaml_file(file.path(out, "recipe.yaml"))
    stopifnot(identical(sort(list.files(out)), c("dataset.yaml", "recipe.yaml")),
              !dir.exists(file.path(batch, "results", "stats_v2_inputs", next_plan$slug)),
              identical(generated_dataset$metadata, "metadata/original.tsv"),
              identical(generated_dataset$quantification, "results/mzmine/original.csv"),
              is.null(generated_dataset$sample_plan), is.null(generated_dataset$roles),
              identical(generated_recipe$sample_plan$group_column, "ATTRIBUTE_mutant"),
              identical(generated_recipe$sample_plan$groups, next_plan$groups),
              is.null(generated_recipe$roles), is.null(generated_recipe$output_root),
              identical(generated_dataset$input_sha256$metadata, original_hashes[["metadata"]]),
              identical(generated_dataset$input_sha256$quantification, original_hashes[["quantification"]]),
              identical(c(metadata = sha256_file(metadata_path), quantification = sha256_file(quant_path)),
                        original_hashes),
              report$valid, nrow(parsed$sample_metadata) == length(kept),
              identical(sort(parsed$sample_metadata$filename), sort(kept)),
              identical(sort(rownames(parsed$raw)), sort(kept)),
              identical(as.numeric(parsed$raw[kept, "feature1"]), as.numeric(match(kept, filenames))),
              identical(parsed$variable_metadata[["row m/z"]], quant[["row m/z"]]),
              identical(parsed$annotations$reference$category, c("a", "b", "c")),
              identical(as.integer(report$role_counts$qc), 7L),
              identical(as.integer(report$role_counts$blank), 3L),
              identical(as.integer(report$role_counts$sample), if (exclude) 25L else 30L))
    if (!exclude) {
      stopifnot(identical(as.integer(unlist(report$group_counts)), c(10L, 10L, 10L)),
                identical(generated$recipe$presentation$description, next_plan$description),
                identical(parsed$sample_metadata$analysis_group[1:30],
                          rep(c("WT", "xdhA", "xdhB"), each = 10L)),
                all(is.na(parsed$sample_metadata$analysis_group[31:40])))
      repeated <- resolve_v2_config(file.path(out, "dataset.yaml"),
                                    file.path(out, "recipe.yaml"), repo_root)
      stopifnot(identical(build_run_identity(generated)$run_hash,
                          build_run_identity(repeated)$run_hash))
      alternate_description <- generated
      alternate_description$recipe$presentation$description <- "Different display label"
      stopifnot(identical(build_run_identity(generated)$run_hash,
                          build_run_identity(alternate_description)$run_hash))
      alternate_mapping <- generated
      alternate_mapping$recipe$sample_plan$groups[[1L]]$target <- "new_WT"
      stopifnot(!identical(build_run_identity(generated)$run_hash,
                           build_run_identity(alternate_mapping)$run_hash))
    }
  }
  alternate <- plan
  alternate$slug <- "other_metadata_group"
  alternate$group_column <- "source_taxon"
  alternate$groups <- list(list(source = "taxon_A", target = "A"),
                           list(source = "taxon_B", target = "B"))
  invoke(alternate, TRUE)
  alternate_dir <- file.path(batch, "recipes", alternate$slug)
  alternate_config <- resolve_v2_config(file.path(alternate_dir, "dataset.yaml"),
                                        file.path(alternate_dir, "recipe.yaml"), repo_root)
  alternate_report <- validate_dataset(read_mzmine_dataset(alternate_config$dataset,
                                                             alternate_config$recipe$sample_plan),
                                        alternate_config$recipe)
  stopifnot(identical(as.integer(unlist(alternate_report$group_counts)), c(15L, 15L)),
            identical(alternate_config$recipe$design$group, "analysis_group"))
  invalid_description <- plan
  invalid_description$slug <- "bad_description"
  for (description in c("  padded", "line\nbreak", "tab\tcharacter", strrep("x", 161L))) {
    invalid_description$description <- description
    stopifnot(any(grepl("description", invoke(invalid_description, FALSE))))
  }
  stopifnot(!dir.exists(file.path(batch, "recipes", invalid_description$slug)))
  existing <- plan
  existing$slug <- "all_merged"
  stopifnot(any(grepl("Refusing to overwrite", invoke(existing, FALSE))))
  linked <- plan
  linked$slug <- "linked_plan"
  link <- file.path(batch, "recipes", linked$slug)
  stopifnot(file.symlink(run_dir, link))
  stopifnot(any(grepl("Refusing to overwrite", invoke(linked, FALSE))))
  unlink(link)
  duplicate <- plan
  duplicate$slug <- "duplicate_source"
  duplicate$groups[[2L]]$source <- "WT_N"
  stopifnot(any(grepl("unique", invoke(duplicate, FALSE))),
            !dir.exists(file.path(batch, "recipes", duplicate$slug)))
  tampered <- plan
  tampered$slug <- "tampered_input"
  write_tsv(metadata[-1L, , drop = FALSE], metadata_path)
  stopifnot(any(grepl("SHA256 mismatch", invoke(tampered, FALSE))),
            !dir.exists(file.path(batch, "recipes", tampered$slug)))
  derived_config <- resolve_v2_config(file.path(batch, "recipes", "exclude_xdhb2", "dataset.yaml"),
                                      file.path(batch, "recipes", "exclude_xdhb2", "recipe.yaml"), repo_root)
  mismatch <- tryCatch({ read_mzmine_dataset(derived_config$dataset,
                                             derived_config$recipe$sample_plan); NULL }, error = identity)
  stopifnot(inherits(mismatch, "error"), grepl("Canonical metadata SHA256 mismatch", conditionMessage(mismatch)))
  write_tsv(metadata, metadata_path)
  changed_quant <- quant
  changed_quant[["sample01 Peak height"]][1L] <- 999
  utils::write.csv(changed_quant, quant_path, row.names = FALSE)
  mismatch <- tryCatch({ validate_v2_run(derived_config); NULL }, error = identity)
  stopifnot(inherits(mismatch, "error"), grepl("Canonical quantification SHA256 mismatch", conditionMessage(mismatch)))
  utils::write.csv(quant, quant_path, row.names = FALSE)
  stopifnot(identical(c(metadata = sha256_file(metadata_path), quantification = sha256_file(quant_path)),
                      original_hashes))
  manifest$inputs$quantification <- NULL
  write_yaml_file(manifest, file.path(run_dir, "manifest.yaml"))
  write_yaml_file(list(run_hash = hash, content_sha256 = sha256_object(run_contents(run_dir))),
                  file.path(run_dir, "COMPLETE"))
  missing_hash <- plan
  missing_hash$slug <- "missing_hash"
  stopifnot(any(grepl("Invalid quantification SHA256", invoke(missing_hash, FALSE))),
            !dir.exists(file.path(batch, "recipes", missing_hash$slug)))
  manifest$inputs$quantification <- original_hashes[["quantification"]]
  for (field in c("effective_recipe", "effective_dataset")) {
    altered <- manifest
    altered[[field]]$sample_plan <- list(group_column = "ATTRIBUTE_mutant", groups = mapping)
    write_yaml_file(altered, file.path(run_dir, "manifest.yaml"))
    invalid_source <- plan
    invalid_source$slug <- paste0("planned_", field)
    stopifnot(any(grepl("corrupt or has an incompatible COMPLETE marker",
                       invoke(invalid_source, FALSE))))
    write_yaml_file(list(run_hash = hash, content_sha256 = sha256_object(run_contents(run_dir))),
                    file.path(run_dir, "COMPLETE"))
    stopifnot(any(grepl("canonical full-batch run without a sample plan", invoke(invalid_source, FALSE))),
              !dir.exists(file.path(batch, "recipes", invalid_source$slug)))
  }
  aliased <- metadata
  names(aliased)[names(aliased) == "sample_type"] <- "raw_role"
  aliased$raw_role[aliased$raw_role == "sample"] <- "specimen"
  aliased$raw_role[aliased$raw_role == "QC"] <- "pool"
  write_tsv(aliased, metadata_path)
  manifest$inputs$metadata <- sha256_file(metadata_path)
  manifest$effective_dataset$metadata_columns <- list(sample_type = "raw_role")
  manifest$effective_recipe$roles <- list(column = "raw_role",
                                           levels = list(sample = "specimen", qc = "pool", blank = "blank"))
  write_yaml_file(manifest, file.path(run_dir, "manifest.yaml"))
  write_yaml_file(list(run_hash = hash, content_sha256 = sha256_object(run_contents(run_dir))),
                  file.path(run_dir, "COMPLETE"))
  aliased_plan <- plan
  aliased_plan$slug <- "aliased_roles"
  invoke(aliased_plan, TRUE)
  alias_dir <- file.path(batch, "recipes", aliased_plan$slug)
  alias_dataset <- read_yaml_file(file.path(alias_dir, "dataset.yaml"))
  alias_recipe <- read_yaml_file(file.path(alias_dir, "recipe.yaml"))
  alias_config <- resolve_v2_config(file.path(alias_dir, "dataset.yaml"),
                                    file.path(alias_dir, "recipe.yaml"), repo_root)
  alias_data <- read_mzmine_dataset(alias_config$dataset, alias_config$recipe$sample_plan)
  alias_report <- validate_dataset(alias_data, alias_config$recipe)
  stopifnot(identical(alias_dataset$roles,
                      list(column = "sample_type",
                           levels = list(sample = "specimen", qc = "pool", blank = "blank"))),
            identical(alias_dataset$input_sha256$metadata, manifest$inputs$metadata),
            is.null(alias_recipe$roles), is.null(alias_recipe$output_root),
            identical(alias_data$sample_metadata$analysis_group[1:30],
                      rep(c("WT", "xdhA", "xdhB"), each = 10L)),
            alias_report$valid, identical(as.integer(alias_report$role_counts$qc), 7L))
  manifest$effective_recipe$roles <- NULL
  manifest$effective_dataset$roles <- list(column = "raw_role",
                                            levels = list(sample = "specimen", qc = "pool", blank = "blank"))
  write_yaml_file(manifest, file.path(run_dir, "manifest.yaml"))
  write_yaml_file(list(run_hash = hash, content_sha256 = sha256_object(run_contents(run_dir))),
                  file.path(run_dir, "COMPLETE"))
  fallback <- plan
  fallback$slug <- "dataset_roles"
  invoke(fallback, TRUE)
  fallback_dir <- file.path(batch, "recipes", fallback$slug)
  stopifnot(identical(read_yaml_file(file.path(fallback_dir, "dataset.yaml"))$roles,
                      alias_dataset$roles))
})
