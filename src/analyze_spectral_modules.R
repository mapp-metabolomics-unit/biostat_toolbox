#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(rlang)
  library(stringr)
  library(tibble)
  library(digest)
  library(igraph)
  library(MAPPstructToolbox)
})

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("--file=", args_full)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}

source(file.path(script_dir, "helpers.r"), local = TRUE)

option_list <- list(
  make_option(c("-p", "--params"), default = "params/params.yaml", help = "Path to params.yaml [default %default]"),
  make_option(c("-u", "--params-user"), default = "params/params_user.yaml", help = "Path to params_user.yaml [default %default]"),
  make_option(c("-H", "--hash"), default = NULL, help = "Hash of the result directory under the configured stats output path"),
  make_option(c("-r", "--results-dir"), default = NULL, help = "Override directory that contains DE.rds and foldchange_pvalues.csv"),
  make_option(c("-n", "--network-file"), default = NULL, help = "Override path to filtered_pairs.tsv"),
  make_option(c("-o", "--output-dir"), default = NULL, help = "Directory used to save spectral-module outputs [default results_dir/spectral_modules]"),
  make_option(c("--pvalue-column"), default = NULL, help = "Name of the p-value column to use (default: first *_p_value column)"),
  make_option(c("--module-mode"), default = "neutral", help = "Module discovery mode: neutral or discriminant [default %default]"),
  make_option(c("--module-method"), default = "auto", help = "Module definition: auto, louvain, or component [default %default]"),
  make_option(c("--discriminant-stat"), default = "t", help = "Feature-level discriminant statistic for discriminant mode: t or signed_log10p [default %default]"),
  make_option(c("--same-direction-only"), default = "FALSE", help = "For discriminant mode, keep only edges between features with the same effect direction: TRUE or FALSE [default %default]"),
  make_option(c("--min-module-size"), type = "integer", default = 3, help = "Minimum number of features per retained module [default %default]"),
  make_option(c("--min-cosine"), type = "double", default = 0, help = "Optional minimum cosine threshold for network edges [default %default]"),
  make_option(c("--n-permutations"), type = "integer", default = 999, help = "Number of label permutations for module p-values [default %default]"),
  make_option(c("--top-n-plots"), type = "integer", default = 12, help = "Number of top modules to include in the summary plot [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

normalize_option_name <- function(opt, underscore_name, hyphen_name) {
  if (is.null(opt[[underscore_name]]) && !is.null(opt[[hyphen_name]])) {
    opt[[underscore_name]] <- opt[[hyphen_name]]
  }
  opt
}

opt <- normalize_option_name(opt, "network_file", "network-file")
opt <- normalize_option_name(opt, "output_dir", "output-dir")
opt <- normalize_option_name(opt, "pvalue_column", "pvalue-column")
opt <- normalize_option_name(opt, "module_mode", "module-mode")
opt <- normalize_option_name(opt, "module_method", "module-method")
opt <- normalize_option_name(opt, "discriminant_stat", "discriminant-stat")
opt <- normalize_option_name(opt, "same_direction_only", "same-direction-only")
opt <- normalize_option_name(opt, "min_module_size", "min-module-size")
opt <- normalize_option_name(opt, "min_cosine", "min-cosine")
opt <- normalize_option_name(opt, "n_permutations", "n-permutations")
opt <- normalize_option_name(opt, "top_n_plots", "top-n-plots")

raw_cli_args <- commandArgs(trailingOnly = TRUE)
extract_flag_value <- function(args, flag_long, flag_short = NULL) {
  check_flag <- function(flg) {
    if (is.null(flg)) {
      return(NULL)
    }
    eq_pattern <- paste0("^", flg, "=")
    idx_eq <- grep(eq_pattern, args)
    if (length(idx_eq)) {
      return(sub(eq_pattern, "", args[idx_eq[1]]))
    }
    idx_space <- which(args == flg)
    if (length(idx_space)) {
      pos <- idx_space[1]
      if (pos < length(args)) {
        return(args[pos + 1])
      }
    }
    NULL
  }
  val <- check_flag(flag_long)
  if (is.null(val)) {
    val <- check_flag(flag_short)
  }
  val
}

explicit_results_dir <- extract_flag_value(raw_cli_args, "--results-dir", "-r")
explicit_params <- extract_flag_value(raw_cli_args, "--params", "-p")
explicit_params_user <- extract_flag_value(raw_cli_args, "--params-user", "-u")
explicit_hash <- extract_flag_value(raw_cli_args, "--hash", "-H")
explicit_network_file <- extract_flag_value(raw_cli_args, "--network-file", "-n")
explicit_output_dir <- extract_flag_value(raw_cli_args, "--output-dir", "-o")
explicit_pvalue_column <- extract_flag_value(raw_cli_args, "--pvalue-column")
explicit_module_mode <- extract_flag_value(raw_cli_args, "--module-mode")
explicit_module_method <- extract_flag_value(raw_cli_args, "--module-method")
explicit_discriminant_stat <- extract_flag_value(raw_cli_args, "--discriminant-stat")
explicit_same_direction_only <- extract_flag_value(raw_cli_args, "--same-direction-only")
explicit_min_module_size <- extract_flag_value(raw_cli_args, "--min-module-size")
explicit_min_cosine <- extract_flag_value(raw_cli_args, "--min-cosine")
explicit_n_permutations <- extract_flag_value(raw_cli_args, "--n-permutations")
explicit_top_n_plots <- extract_flag_value(raw_cli_args, "--top-n-plots")

if (!is.null(explicit_results_dir)) {
  opt$results_dir <- explicit_results_dir
}
if (!is.null(explicit_params)) {
  opt$params <- explicit_params
}
if (!is.null(explicit_params_user)) {
  opt$params_user <- explicit_params_user
}
if (!is.null(explicit_hash)) {
  opt$hash <- explicit_hash
}
if (!is.null(explicit_network_file)) {
  opt$network_file <- explicit_network_file
}
if (!is.null(explicit_output_dir)) {
  opt$output_dir <- explicit_output_dir
}
if (!is.null(explicit_pvalue_column)) {
  opt$pvalue_column <- explicit_pvalue_column
}
if (!is.null(explicit_module_mode)) {
  opt$module_mode <- explicit_module_mode
}
if (!is.null(explicit_module_method)) {
  opt$module_method <- explicit_module_method
}
if (!is.null(explicit_discriminant_stat)) {
  opt$discriminant_stat <- explicit_discriminant_stat
}
if (!is.null(explicit_same_direction_only)) {
  opt$same_direction_only <- explicit_same_direction_only
}
if (!is.null(explicit_min_module_size)) {
  opt$min_module_size <- as.integer(explicit_min_module_size)
}
if (!is.null(explicit_min_cosine)) {
  opt$min_cosine <- as.numeric(explicit_min_cosine)
}
if (!is.null(explicit_n_permutations)) {
  opt$n_permutations <- as.integer(explicit_n_permutations)
}
if (!is.null(explicit_top_n_plots)) {
  opt$top_n_plots <- as.integer(explicit_top_n_plots)
}

trim_option_value <- function(value) {
  if (is.character(value) && length(value)) {
    trimmed <- trimws(value)
    return(trimmed[nzchar(trimmed)])
  }
  value
}

has_value <- function(value) {
  if (is.null(value) || !length(value)) {
    return(FALSE)
  }
  val <- as.character(value[1])
  if (is.na(val)) {
    return(FALSE)
  }
  nzchar(trimws(val))
}

opt$params <- trim_option_value(opt$params)
opt$params_user <- trim_option_value(opt$params_user)
opt$hash <- trim_option_value(opt$hash)
opt$results_dir <- trim_option_value(opt$results_dir)
opt$network_file <- trim_option_value(opt$network_file)
opt$output_dir <- trim_option_value(opt$output_dir)
opt$pvalue_column <- trim_option_value(opt$pvalue_column)
opt$module_mode <- trim_option_value(opt$module_mode)
opt$module_method <- trim_option_value(opt$module_method)
opt$discriminant_stat <- trim_option_value(opt$discriminant_stat)
opt$same_direction_only <- trim_option_value(opt$same_direction_only)

if (is.null(opt$params) || !nzchar(opt$params)) {
  opt$params <- "params/params.yaml"
}
if (is.null(opt$params_user) || !nzchar(opt$params_user)) {
  opt$params_user <- "params/params_user.yaml"
}

resolve_relative_path <- function(path_value, fallback_dir) {
  if (is.null(path_value) || !nzchar(path_value)) {
    return(path_value)
  }
  if (grepl("^/", path_value)) {
    return(path_value)
  }
  if (file.exists(path_value)) {
    return(normalizePath(path_value))
  }
  candidate <- file.path(fallback_dir, path_value)
  if (file.exists(candidate)) {
    return(normalizePath(candidate))
  }
  normalizePath(path_value, mustWork = FALSE)
}

coerce_path_scalar <- function(value, label) {
  if (is.null(value) || !length(value)) {
    stop(sprintf("No value provided for %s.", label))
  }
  value <- as.character(value)
  if (!nzchar(value[1])) {
    stop(sprintf("Empty path received for %s.", label))
  }
  value[1]
}

opt$params <- coerce_path_scalar(resolve_relative_path(opt$params, script_dir), "--params")
opt$params_user <- coerce_path_scalar(resolve_relative_path(opt$params_user, script_dir), "--params-user")
if (has_value(opt$network_file)) {
  opt$network_file <- coerce_path_scalar(resolve_relative_path(opt$network_file, script_dir), "--network-file")
}

ensure_yaml_exists <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("YAML file %s not found.", path))
  }
}

ensure_yaml_exists(opt$params)
ensure_yaml_exists(opt$params_user)

params <- yaml.load_file(opt$params)
params_user <- yaml.load_file(opt$params_user)

params$paths$docs <- params_user$paths$docs
params$paths$output <- params_user$paths$output
params$operating_system$system <- params_user$operating_system$system
params$operating_system$pandoc <- params_user$operating_system$pandoc
params$target$sample_metadata_header <- tolower(params$target$sample_metadata_header)

valid_module_modes <- c("neutral", "discriminant")
module_mode <- tolower(opt$module_mode)
if (!module_mode %in% valid_module_modes) {
  stop(sprintf("Invalid module mode '%s'. Choose one of: %s", opt$module_mode, paste(valid_module_modes, collapse = ", ")))
}

valid_module_methods <- c("auto", "louvain", "component")
module_method <- tolower(opt$module_method)
if (!module_method %in% valid_module_methods) {
  stop(sprintf("Invalid module method '%s'. Choose one of: %s", opt$module_method, paste(valid_module_methods, collapse = ", ")))
}

valid_discriminant_stats <- c("t", "signed_log10p")
discriminant_stat <- tolower(opt$discriminant_stat)
if (!discriminant_stat %in% valid_discriminant_stats) {
  stop(sprintf("Invalid discriminant stat '%s'. Choose one of: %s", opt$discriminant_stat, paste(valid_discriminant_stats, collapse = ", ")))
}

parse_bool_string <- function(x, label) {
  if (is.logical(x) && length(x) == 1) {
    return(isTRUE(x))
  }
  x <- tolower(trimws(as.character(x[1])))
  if (x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("%s must be TRUE or FALSE.", label))
}

same_direction_only <- parse_bool_string(opt$same_direction_only, "--same-direction-only")

if (is.na(opt$min_module_size) || opt$min_module_size < 2) {
  stop("--min-module-size must be an integer >= 2.")
}
if (is.na(opt$min_cosine) || opt$min_cosine < 0 || opt$min_cosine > 1) {
  stop("--min-cosine must be between 0 and 1.")
}
if (is.na(opt$n_permutations) || opt$n_permutations < 99) {
  stop("--n-permutations must be >= 99.")
}
if (is.na(opt$top_n_plots) || opt$top_n_plots < 1) {
  stop("--top-n-plots must be >= 1.")
}

config_hash <- convert_yaml_to_single_row_df_with_hash(params)$hash

resolve_results_base_dir <- function(params) {
  if (!is.null(params$paths$output) && nzchar(params$paths$output)) {
    return(params$paths$output)
  }
  file.path(params$paths$docs, params$mapp_project, params$mapp_batch, "results", "stats")
}

resolve_results_dir <- function(params, inferred_hash, hash_override, override) {
  if (has_value(override)) {
    return(trimws(as.character(override[1])))
  }
  if (has_value(hash_override)) {
    selected_hash <- trimws(as.character(hash_override[1]))
  } else {
    selected_hash <- inferred_hash
  }
  file.path(resolve_results_base_dir(params), selected_hash)
}

results_base_dir <- resolve_results_base_dir(params)
results_dir <- resolve_results_dir(params, config_hash, opt$hash, opt$results_dir)
if (has_value(opt$results_dir)) {
  message(sprintf("Using user-specified results directory: %s", results_dir))
} else if (has_value(opt$hash)) {
  message(sprintf("Using results directory from user-specified hash %s: %s", opt$hash, results_dir))
} else {
  message(sprintf("Using inferred results directory: %s", results_dir))
}
if (!dir.exists(results_dir)) {
  available_hashes <- character()
  if (dir.exists(results_base_dir)) {
    available_hashes <- basename(list.dirs(results_base_dir, full.names = TRUE, recursive = FALSE))
  }
  if (!has_value(opt$results_dir) && length(available_hashes)) {
    stop(sprintf(
      paste(
        "Results directory %s does not exist.",
        "Available hashes under %s: %s",
        "Use --hash <hash> to select one of them or --results-dir to pass a full path.",
        sep = " "
      ),
      results_dir,
      results_base_dir,
      paste(available_hashes, collapse = ", ")
    ))
  }
  stop(sprintf("Results directory %s does not exist.", results_dir))
}

output_dir <- opt$output_dir
if (is.null(output_dir) || !nzchar(output_dir)) {
  output_dir <- file.path(results_dir, "spectral_modules")
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

working_directory <- file.path(params$paths$docs, params$mapp_project, params$mapp_batch)
network_file <- opt$network_file
if (is.null(network_file) || !nzchar(network_file)) {
  network_file <- file.path(
    working_directory,
    "results",
    "met_annot_enhancer",
    params$gnps_job_id,
    "nf_output",
    "networking",
    "filtered_pairs.tsv"
  )
}

de_path <- file.path(results_dir, "DE.rds")
foldchange_path <- file.path(results_dir, "foldchange_pvalues.csv")
if (!file.exists(de_path)) {
  stop(sprintf("Missing DE.rds at %s", de_path))
}
if (!file.exists(foldchange_path)) {
  stop(sprintf("Missing foldchange_pvalues.csv at %s", foldchange_path))
}
if (!file.exists(network_file)) {
  stop(sprintf("Missing filtered_pairs.tsv at %s", network_file))
}

message("Loading DE object...")
DE <- readRDS(de_path)
data_matrix <- as.data.frame(DE$data)
sample_meta <- DE$sample_meta
if (!nrow(data_matrix) || !ncol(data_matrix)) {
  stop("DE$data is empty.")
}
if (!nrow(sample_meta)) {
  stop("DE$sample_meta is empty.")
}

sample_ids <- rownames(data_matrix)
if (is.null(sample_ids) || !length(sample_ids)) {
  stop("DE$data does not contain sample row names.")
}
if (!all(sample_ids %in% rownames(sample_meta))) {
  missing_samples <- setdiff(sample_ids, rownames(sample_meta))
  stop(sprintf("The following samples are missing from DE$sample_meta: %s", paste(missing_samples, collapse = ", ")))
}
sample_meta <- sample_meta[sample_ids, , drop = FALSE]

target_factor <- params$target$sample_metadata_header
if (!target_factor %in% colnames(sample_meta)) {
  stop(sprintf("Sample metadata header %s is not present in DE$sample_meta.", target_factor))
}

group_factor <- factor(sample_meta[[target_factor]])
if (nlevels(group_factor) < 2) {
  stop(sprintf("Target factor %s must contain at least two groups.", target_factor))
}

message("Loading fold-change/p-value table...")
foldchange_tbl <- readr::read_csv(foldchange_path, show_col_types = FALSE) %>%
  mutate(feature_id = as.character(feature_id))

pvalue_cols <- grep("_p_value$", names(foldchange_tbl), value = TRUE)
if (!length(pvalue_cols)) {
  stop("No columns ending with _p_value found in foldchange_pvalues.csv")
}
if (!is.null(opt$pvalue_column)) {
  if (!opt$pvalue_column %in% pvalue_cols) {
    stop(sprintf("Requested p-value column %s was not found. Available options: %s", opt$pvalue_column, paste(pvalue_cols, collapse = ", ")))
  }
  p_value_column <- opt$pvalue_column
} else if (length(pvalue_cols) == 1) {
  p_value_column <- pvalue_cols
} else {
  message(sprintf("Multiple p-value columns detected (%s); defaulting to %s. Use --pvalue-column to override.", paste(pvalue_cols, collapse = ", "), pvalue_cols[1]))
  p_value_column <- pvalue_cols[1]
}

comparison_prefix <- sub("_p_value$", "", p_value_column)
log2_fc_column <- paste0(comparison_prefix, "_fold_change_log2")
fold_change_column <- paste0(comparison_prefix, "_fold_change")
comparison_groups <- strsplit(comparison_prefix, "_vs_", fixed = TRUE)[[1]]
if (length(comparison_groups) == 2 && all(comparison_groups %in% unique(as.character(sample_meta[[target_factor]])))) {
  group_factor <- factor(sample_meta[[target_factor]], levels = comparison_groups)
} else {
  group_factor <- factor(sample_meta[[target_factor]])
}
feature_p_threshold <- params$posthoc$p_value
if (is.null(feature_p_threshold) || is.na(feature_p_threshold)) {
  feature_p_threshold <- 0.05
}

compute_feature_discriminant_scores <- function(data_matrix, group_factor, stat_method, foldchange_tbl, p_value_column, log2_fc_column) {
  feature_ids <- colnames(data_matrix)
  score_tbl <- tibble(
    feature_id = feature_ids,
    discriminant_score = 0,
    abs_discriminant_score = 0,
    discriminant_direction = 0
  )

  if (stat_method == "signed_log10p") {
    if (!all(c(p_value_column, log2_fc_column) %in% colnames(foldchange_tbl))) {
      stop(sprintf("Columns %s and %s are required for --discriminant-stat signed_log10p.", p_value_column, log2_fc_column))
    }
    score_tbl <- foldchange_tbl %>%
      filter(feature_id %in% feature_ids) %>%
      transmute(
        feature_id = feature_id,
        discriminant_score = sign(.data[[log2_fc_column]]) * (-log10(pmax(.data[[p_value_column]], 1e-300))),
        abs_discriminant_score = abs(discriminant_score),
        discriminant_direction = sign(discriminant_score)
      )
  } else {
    if (nlevels(group_factor) != 2) {
      stop("Discriminant mode with --discriminant-stat t currently requires exactly two target groups.")
    }
    group_one_idx <- group_factor == levels(group_factor)[1]
    group_two_idx <- group_factor == levels(group_factor)[2]
    group_one <- as.matrix(data_matrix[group_one_idx, , drop = FALSE])
    group_two <- as.matrix(data_matrix[group_two_idx, , drop = FALSE])
    mean_one <- colMeans(group_one, na.rm = TRUE)
    mean_two <- colMeans(group_two, na.rm = TRUE)
    var_one <- apply(group_one, 2, stats::var, na.rm = TRUE)
    var_two <- apply(group_two, 2, stats::var, na.rm = TRUE)
    denom <- sqrt((var_one / nrow(group_one)) + (var_two / nrow(group_two)))
    t_stat <- (mean_one - mean_two) / denom
    mean_diff <- mean_one - mean_two
    t_stat[!is.finite(t_stat)] <- mean_diff[!is.finite(t_stat)]
    t_stat[!is.finite(t_stat)] <- 0
    score_tbl <- tibble(
      feature_id = feature_ids,
      discriminant_score = as.numeric(t_stat),
      abs_discriminant_score = abs(as.numeric(t_stat)),
      discriminant_direction = sign(as.numeric(t_stat))
    )
  }

  abs_vals <- score_tbl$abs_discriminant_score
  keep_vals <- is.finite(abs_vals)
  if (any(keep_vals)) {
    vals <- abs_vals[keep_vals]
    if (length(unique(vals)) == 1) {
      score_tbl$discriminant_relevance <- ifelse(keep_vals, 1, 0.25)
    } else {
      score_tbl$discriminant_relevance <- 0.25
      score_tbl$discriminant_relevance[keep_vals] <- 0.25 + 0.75 * ((vals - min(vals)) / (max(vals) - min(vals)))
    }
  } else {
    score_tbl$discriminant_relevance <- 0.25
  }

  score_tbl
}

feature_score_tbl <- compute_feature_discriminant_scores(
  data_matrix = data_matrix,
  group_factor = group_factor,
  stat_method = discriminant_stat,
  foldchange_tbl = foldchange_tbl,
  p_value_column = p_value_column,
  log2_fc_column = log2_fc_column
)

message("Loading GNPS spectral network...")
edge_tbl <- readr::read_tsv(network_file, show_col_types = FALSE)
required_edge_cols <- c("CLUSTERID1", "CLUSTERID2", "ComponentIndex", "Cosine")
missing_edge_cols <- setdiff(required_edge_cols, names(edge_tbl))
if (length(missing_edge_cols)) {
  stop(sprintf("Network file is missing required columns: %s", paste(missing_edge_cols, collapse = ", ")))
}

available_features <- colnames(data_matrix)
edge_tbl <- edge_tbl %>%
  transmute(
    from = as.character(CLUSTERID1),
    to = as.character(CLUSTERID2),
    component_index = as.character(ComponentIndex),
    cosine = as.numeric(Cosine),
    delta_mz = if ("DeltaMZ" %in% names(edge_tbl)) as.numeric(DeltaMZ) else NA_real_
  ) %>%
  filter(from != to, !is.na(cosine), cosine >= opt$min_cosine) %>%
  filter(from %in% available_features, to %in% available_features)

if (!nrow(edge_tbl)) {
  stop("No spectral-network edges remain after matching to DE$data and applying the cosine threshold.")
}

edge_tbl <- edge_tbl %>%
  left_join(
    feature_score_tbl %>%
      select(feature_id, discriminant_score, abs_discriminant_score, discriminant_direction, discriminant_relevance) %>%
      rename(
        from = feature_id,
        from_discriminant_score = discriminant_score,
        from_abs_discriminant_score = abs_discriminant_score,
        from_discriminant_direction = discriminant_direction,
        from_discriminant_relevance = discriminant_relevance
      ),
    by = "from"
  ) %>%
  left_join(
    feature_score_tbl %>%
      select(feature_id, discriminant_score, abs_discriminant_score, discriminant_direction, discriminant_relevance) %>%
      rename(
        to = feature_id,
        to_discriminant_score = discriminant_score,
        to_abs_discriminant_score = abs_discriminant_score,
        to_discriminant_direction = discriminant_direction,
        to_discriminant_relevance = discriminant_relevance
      ),
    by = "to"
  )

if (module_mode == "discriminant") {
  if (same_direction_only) {
    edge_tbl <- edge_tbl %>%
      filter(
        from_discriminant_direction != 0,
        to_discriminant_direction != 0,
        from_discriminant_direction == to_discriminant_direction
      )
  }
  edge_tbl <- edge_tbl %>%
    mutate(
      cluster_weight = cosine * sqrt(from_discriminant_relevance * to_discriminant_relevance)
    )
} else {
  edge_tbl <- edge_tbl %>%
    mutate(cluster_weight = cosine)
}

edge_tbl <- edge_tbl %>%
  filter(is.finite(cluster_weight), cluster_weight > 0)

if (!nrow(edge_tbl)) {
  stop("No spectral-network edges remain after applying the selected module mode and directionality constraints.")
}

vertex_component_tbl <- bind_rows(
  edge_tbl %>% transmute(feature_id = from, component_index = component_index),
  edge_tbl %>% transmute(feature_id = to, component_index = component_index)
) %>%
  distinct(feature_id, component_index)

network_graph <- graph_from_data_frame(
  edge_tbl %>% select(from, to, cosine, cluster_weight, component_index, delta_mz),
  directed = FALSE,
  vertices = vertex_component_tbl %>% transmute(name = feature_id)
)
network_graph <- simplify(
  network_graph,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = list(
    cosine = "max",
    cluster_weight = "max",
    component_index = "first",
    delta_mz = "first"
  )
)

node_strength_tbl <- tibble(
  feature_id = names(strength(network_graph, vids = V(network_graph), weights = E(network_graph)$cluster_weight)),
  weighted_degree = as.numeric(strength(network_graph, vids = V(network_graph), weights = E(network_graph)$cluster_weight)),
  spectral_weighted_degree = as.numeric(strength(network_graph, vids = V(network_graph), weights = E(network_graph)$cosine)),
  degree = as.numeric(degree(network_graph, mode = "all"))
)

first_non_missing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x)) & x != "NA"]
  if (!length(x)) {
    return(NA_character_)
  }
  counts <- sort(table(x), decreasing = TRUE)
  names(counts)[1]
}

get_first_available_value <- function(tbl, cols, default = NA_character_) {
  existing_cols <- cols[cols %in% colnames(tbl)]
  if (!length(existing_cols) || !nrow(tbl)) {
    return(default)
  }
  for (col_name in existing_cols) {
    val <- tbl[[col_name]][1]
    if (!is.null(val) && length(val) && !is.na(val) && nzchar(trimws(as.character(val)))) {
      return(as.character(val))
    }
  }
  default
}

safe_stat <- function(x, fn, default = NA_real_) {
  vals <- x[is.finite(x)]
  if (!length(vals)) {
    return(default)
  }
  fn(vals)
}

get_module_assignments <- function(graph, vertex_component_tbl, method, min_size) {
  if (method == "component") {
    module_tbl <- vertex_component_tbl %>%
      mutate(module_index = component_index)
  } else {
    membership <- cluster_louvain(graph, weights = E(graph)$cluster_weight)$membership
    module_tbl <- tibble(
      feature_id = V(graph)$name,
      module_index = as.character(unname(membership))
    ) %>%
      left_join(vertex_component_tbl, by = "feature_id")
  }

  module_tbl %>%
    mutate(module_id = paste0("component_", component_index, "_module_", module_index)) %>%
    add_count(module_id, name = "module_size") %>%
    filter(module_size >= min_size) %>%
    arrange(component_index, module_id, feature_id)
}

module_tbl <- NULL
effective_method <- module_method
if (module_method == "auto") {
  louvain_tbl <- get_module_assignments(network_graph, vertex_component_tbl, "louvain", opt$min_module_size)
  if (nrow(louvain_tbl)) {
    module_tbl <- louvain_tbl
    effective_method <- "louvain"
  } else {
    module_tbl <- get_module_assignments(network_graph, vertex_component_tbl, "component", opt$min_module_size)
    effective_method <- "component"
  }
} else {
  module_tbl <- get_module_assignments(network_graph, vertex_component_tbl, module_method, opt$min_module_size)
}

if (!nrow(module_tbl)) {
  stop(sprintf(
    "No modules met the minimum size threshold (%s features) using method '%s'.",
    opt$min_module_size,
    module_method
  ))
}

message(sprintf(
  "Retained %s spectral modules using method '%s' (minimum size %s).",
  dplyr::n_distinct(module_tbl$module_id),
  effective_method,
  opt$min_module_size
))

module_edge_tbl <- edge_tbl %>%
  inner_join(module_tbl %>% select(feature_id, module_id) %>% rename(from = feature_id), by = "from") %>%
  inner_join(module_tbl %>% select(feature_id, module_id) %>% rename(to = feature_id, module_id_to = module_id), by = "to") %>%
  filter(module_id == module_id_to) %>%
  select(module_id, from, to, component_index, cosine, cluster_weight, delta_mz, starts_with("from_discriminant"), starts_with("to_discriminant"))

module_feature_tbl <- module_tbl %>%
  left_join(node_strength_tbl, by = "feature_id") %>%
  left_join(feature_score_tbl, by = "feature_id") %>%
  left_join(foldchange_tbl, by = "feature_id")

safe_name <- function(x) {
  x <- gsub("[^[:alnum:]]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

rescale_numeric <- function(x, to = c(0, 1), default = mean(to)) {
  x <- as.numeric(x)
  out <- rep(default, length(x))
  keep <- is.finite(x)
  if (!any(keep)) {
    return(out)
  }
  vals <- x[keep]
  if (length(unique(vals)) == 1) {
    out[keep] <- mean(to)
    return(out)
  }
  rng <- range(vals, na.rm = TRUE)
  out[keep] <- to[1] + (vals - rng[1]) * (to[2] - to[1]) / (rng[2] - rng[1])
  out
}

make_diverging_colors <- function(values, low = "#3b4cc0", mid = "#f7f7f7", high = "#b40426", na_color = "#bdbdbd") {
  values <- as.numeric(values)
  out <- rep(na_color, length(values))
  keep <- is.finite(values)
  if (!any(keep)) {
    return(out)
  }
  vals <- values[keep]
  max_abs <- max(abs(vals), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) {
    out[keep] <- mid
    return(out)
  }
  palette_fn <- grDevices::colorRampPalette(c(low, mid, high))
  palette_vals <- palette_fn(201)
  scaled <- round(((vals + max_abs) / (2 * max_abs)) * 200) + 1
  scaled <- pmax(1, pmin(201, scaled))
  out[keep] <- palette_vals[scaled]
  out
}

choose_annotation_label <- function(tbl) {
  feature_id <- get_first_available_value(tbl, c("feature_id"), default = "NA")
  label <- get_first_available_value(
    tbl,
    c("sirius_chebiasciiname", "sirius_name", "canopus_npc_class", "feature_id"),
    default = NA_character_
  )
  if (is.na(label) || !nzchar(trimws(label)) || identical(label, feature_id)) {
    return(feature_id)
  }
  paste0(
    feature_id,
    "\n",
    stringr::str_trunc(label, width = 50)
  )
}

compute_module_score <- function(sample_feature_matrix) {
  sample_feature_matrix <- as.matrix(sample_feature_matrix)
  if (!ncol(sample_feature_matrix)) {
    stop("Cannot compute a module score with zero features.")
  }
  keep_cols <- apply(sample_feature_matrix, 2, function(x) {
    vals <- x[is.finite(x)]
    length(vals) > 1 && stats::sd(vals) > 0
  })
  sample_feature_matrix <- sample_feature_matrix[, keep_cols, drop = FALSE]
  if (!ncol(sample_feature_matrix)) {
    zero_score <- rep(0, nrow(sample_feature_matrix))
    names(zero_score) <- rownames(sample_feature_matrix)
    return(list(
      score = zero_score,
      loading_tbl = tibble(feature_id = character(), loading = numeric()),
      variance_explained = 0,
      score_method = "zero_variance_module"
    ))
  }

  if (ncol(sample_feature_matrix) == 1) {
    score_vec <- as.numeric(scale(sample_feature_matrix[, 1]))
    if (all(!is.finite(score_vec))) {
      score_vec <- rep(0, nrow(sample_feature_matrix))
    }
    names(score_vec) <- rownames(sample_feature_matrix)
    return(list(
      score = score_vec,
      loading_tbl = tibble(feature_id = colnames(sample_feature_matrix), loading = 1),
      variance_explained = 1,
      score_method = "single_feature_zscore"
    ))
  }

  scaled_matrix <- scale(sample_feature_matrix)
  scaled_matrix[!is.finite(scaled_matrix)] <- 0
  pca_fit <- stats::prcomp(scaled_matrix, center = FALSE, scale. = FALSE)
  module_score <- pca_fit$x[, 1]
  loadings <- pca_fit$rotation[, 1]
  mean_profile <- rowMeans(scaled_matrix, na.rm = TRUE)
  cor_sign <- suppressWarnings(stats::cor(module_score, mean_profile, use = "complete.obs"))
  if (!is.na(cor_sign) && cor_sign < 0) {
    module_score <- -module_score
    loadings <- -loadings
  }

  list(
    score = stats::setNames(as.numeric(module_score), rownames(sample_feature_matrix)),
    loading_tbl = tibble(feature_id = names(loadings), loading = as.numeric(loadings)),
    variance_explained = (pca_fit$sdev[1]^2) / sum(pca_fit$sdev^2),
    score_method = "pc1_eigenmetabolite"
  )
}

compute_group_statistic <- function(score, group_factor) {
  keep <- is.finite(score) & !is.na(group_factor)
  score <- score[keep]
  group_factor <- droplevels(group_factor[keep])
  if (nlevels(group_factor) < 2) {
    stop("At least two groups with finite module scores are required.")
  }

  if (nlevels(group_factor) == 2) {
    level_names <- levels(group_factor)
    group_one <- score[group_factor == level_names[1]]
    group_two <- score[group_factor == level_names[2]]
    mean_diff <- mean(group_one) - mean(group_two)
    stat_denom <- sqrt(stats::var(group_one) / length(group_one) + stats::var(group_two) / length(group_two))
    if (!is.finite(stat_denom) || stat_denom == 0) {
      statistic <- abs(mean_diff)
    } else {
      statistic <- abs(mean_diff / stat_denom)
    }
    pooled_sd <- sqrt(((length(group_one) - 1) * stats::var(group_one) + (length(group_two) - 1) * stats::var(group_two)) /
      max(length(group_one) + length(group_two) - 2, 1))
    effect_size <- if (is.finite(pooled_sd) && pooled_sd > 0) mean_diff / pooled_sd else mean_diff
    effect_name <- paste0("cohens_d_", safe_name(level_names[1]), "_minus_", safe_name(level_names[2]))
    group_means <- stats::setNames(
      c(mean(group_one), mean(group_two)),
      paste0("mean_", safe_name(level_names))
    )
    effect_values <- stats::setNames(mean_diff, paste0("delta_", safe_name(level_names[1]), "_minus_", safe_name(level_names[2])))
    list(
      statistic = statistic,
      effect_name = effect_name,
      effect_size = effect_size,
      group_means = group_means,
      group_effects = effect_values
    )
  } else {
    fit <- stats::lm(score ~ group_factor)
    anova_tbl <- stats::anova(fit)
    statistic <- anova_tbl[["F value"]][1]
    if (!is.finite(statistic)) {
      statistic <- 0
    }
    group_means <- tapply(score, group_factor, mean)
    named_means <- stats::setNames(as.numeric(group_means), paste0("mean_", safe_name(names(group_means))))
    list(
      statistic = as.numeric(statistic),
      effect_name = "anova_r_squared",
      effect_size = summary(fit)$r.squared,
      group_means = named_means,
      group_effects = numeric()
    )
  }
}

permutation_p_value <- function(score, group_factor, n_permutations) {
  observed <- compute_group_statistic(score, group_factor)
  perm_stats <- numeric(n_permutations)
  for (i in seq_len(n_permutations)) {
    perm_stats[i] <- compute_group_statistic(score, sample(group_factor))$statistic
  }
  p_value <- (sum(perm_stats >= observed$statistic) + 1) / (n_permutations + 1)
  list(
    observed = observed,
    p_value = p_value
  )
}

module_ids <- unique(module_tbl$module_id)
module_score_list <- vector("list", length(module_ids))
module_summary_list <- vector("list", length(module_ids))
module_member_list <- vector("list", length(module_ids))

for (i in seq_along(module_ids)) {
  module_id <- module_ids[i]
  this_module_features <- module_tbl %>%
    filter(module_id == !!module_id) %>%
    pull(feature_id) %>%
    unique()

  module_matrix <- data_matrix[, this_module_features, drop = FALSE]
  score_fit <- compute_module_score(module_matrix)
  module_score_vec <- score_fit$score
  perm_fit <- permutation_p_value(module_score_vec, group_factor, opt$n_permutations)

  score_tbl <- tibble(
    sample_id = names(module_score_vec),
    module_id = module_id,
    module_score = as.numeric(module_score_vec),
    !!target_factor := sample_meta[names(module_score_vec), target_factor]
  )

  module_scores_this <- score_tbl
  module_scores_this$score_method <- score_fit$score_method
  module_scores_this$variance_explained <- score_fit$variance_explained
  module_score_list[[i]] <- module_scores_this

  module_members_this <- module_feature_tbl %>%
    filter(module_id == !!module_id) %>%
    left_join(score_fit$loading_tbl, by = "feature_id") %>%
    mutate(
      loading = dplyr::coalesce(loading, 0),
      abs_loading = abs(loading)
    ) %>%
    arrange(desc(weighted_degree), desc(abs_loading), feature_id)
  module_member_list[[i]] <- module_members_this

  edge_stats_this <- module_edge_tbl %>%
    filter(module_id == !!module_id)
  hub_feature <- module_members_this %>%
    arrange(desc(weighted_degree), desc(abs_loading), feature_id) %>%
    slice(1)
  feature_p_values <- module_members_this[[p_value_column]]
  log2_fc_values <- if (log2_fc_column %in% colnames(module_members_this)) module_members_this[[log2_fc_column]] else numeric()

  group_mean_values <- perm_fit$observed$group_means
  group_effect_values <- perm_fit$observed$group_effects

  summary_tbl <- tibble(
    module_id = module_id,
    component_index = first(module_members_this$component_index),
    module_mode = module_mode,
    discriminant_stat = if (module_mode == "discriminant") discriminant_stat else NA_character_,
    same_direction_only = if (module_mode == "discriminant") same_direction_only else NA,
    module_size = nrow(module_members_this),
    n_edges = nrow(edge_stats_this),
    mean_cosine = if (nrow(edge_stats_this)) mean(edge_stats_this$cosine, na.rm = TRUE) else NA_real_,
    max_cosine = if (nrow(edge_stats_this)) max(edge_stats_this$cosine, na.rm = TRUE) else NA_real_,
    mean_cluster_weight = if (nrow(edge_stats_this)) mean(edge_stats_this$cluster_weight, na.rm = TRUE) else NA_real_,
    max_cluster_weight = if (nrow(edge_stats_this)) max(edge_stats_this$cluster_weight, na.rm = TRUE) else NA_real_,
    module_statistic = perm_fit$observed$statistic,
    module_p_value = perm_fit$p_value,
    score_method = score_fit$score_method,
    variance_explained = score_fit$variance_explained,
    effect_metric = perm_fit$observed$effect_name,
    effect_size = perm_fit$observed$effect_size,
    mean_abs_discriminant_score = safe_stat(abs(module_members_this$discriminant_score), mean),
    mean_signed_discriminant_score = safe_stat(module_members_this$discriminant_score, mean),
    n_feature_level_hits = sum(feature_p_values < feature_p_threshold, na.rm = TRUE),
    min_feature_p_value = safe_stat(feature_p_values, min),
    mean_feature_p_value = safe_stat(feature_p_values, mean),
    mean_abs_feature_log2_fc = safe_stat(abs(log2_fc_values), mean),
    mean_signed_feature_log2_fc = safe_stat(log2_fc_values, mean),
    hub_feature_id = hub_feature$feature_id,
    hub_feature_name = get_first_available_value(hub_feature, c("sirius_chebiasciiname", "sirius_name", "feature_id")),
    top_sirius_name = if ("sirius_chebiasciiname" %in% colnames(module_members_this)) first_non_missing(module_members_this$sirius_chebiasciiname) else NA_character_,
    top_canopus_pathway = if ("canopus_npc_pathway" %in% colnames(module_members_this)) first_non_missing(module_members_this$canopus_npc_pathway) else NA_character_,
    top_canopus_superclass = if ("canopus_npc_superclass" %in% colnames(module_members_this)) first_non_missing(module_members_this$canopus_npc_superclass) else NA_character_,
    top_canopus_class = if ("canopus_npc_class" %in% colnames(module_members_this)) first_non_missing(module_members_this$canopus_npc_class) else NA_character_
  )

  if (length(group_mean_values)) {
    for (nm in names(group_mean_values)) {
      summary_tbl[[nm]] <- unname(group_mean_values[[nm]])
    }
  }
  if (length(group_effect_values)) {
    for (nm in names(group_effect_values)) {
      summary_tbl[[nm]] <- unname(group_effect_values[[nm]])
    }
  }

  module_summary_list[[i]] <- summary_tbl
}

module_scores_tbl <- bind_rows(module_score_list)
module_members_tbl <- bind_rows(module_member_list) %>%
  arrange(module_id, desc(weighted_degree), desc(abs_loading), feature_id)
module_summary_tbl <- bind_rows(module_summary_list) %>%
  mutate(
    module_q_value = p.adjust(module_p_value, method = "BH")
  ) %>%
  arrange(module_q_value, module_p_value, desc(abs(effect_size)))

summary_path <- file.path(output_dir, "spectral_module_summary.csv")
members_path <- file.path(output_dir, "spectral_module_members.csv")
scores_path <- file.path(output_dir, "spectral_module_scores.csv")

readr::write_csv(module_summary_tbl, summary_path)
readr::write_csv(module_members_tbl, members_path)
readr::write_csv(module_scores_tbl, scores_path)

plot_module_ids <- module_summary_tbl %>%
  slice_head(n = opt$top_n_plots) %>%
  pull(module_id)

plot_path <- file.path(output_dir, "spectral_module_scores_top.png")
if (length(plot_module_ids)) {
  plot_df <- module_scores_tbl %>%
    filter(module_id %in% plot_module_ids) %>%
    mutate(module_id = factor(module_id, levels = plot_module_ids))

  plot_obj <- ggplot(plot_df, aes(x = .data[[target_factor]], y = module_score, fill = .data[[target_factor]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_point(position = position_jitter(width = 0.15), alpha = 0.8, size = 2) +
    facet_wrap(~ module_id, scales = "free_y") +
    labs(
      x = target_factor,
      y = "Spectral module score",
      title = "Top spectrally coherent modules associated with the target factor",
      subtitle = sprintf("Method: %s | p-values from %s label permutations", effective_method, opt$n_permutations)
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 9)
    )

  ggsave(plot_path, plot = plot_obj, width = 12, height = 8)
}

graphml_dir <- file.path(output_dir, "graphml")
network_plot_dir <- file.path(output_dir, "network_plots")
dir.create(graphml_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_plot_dir, recursive = TRUE, showWarnings = FALSE)

if (length(plot_module_ids)) {
  for (module_id in plot_module_ids) {
    module_summary_row <- module_summary_tbl %>%
      filter(module_id == !!module_id) %>%
      slice(1)
    module_vertices <- module_members_tbl %>%
      filter(module_id == !!module_id) %>%
      mutate(
        feature_p_value = if (p_value_column %in% colnames(.)) .data[[p_value_column]] else NA_real_,
        feature_fold_change = if (fold_change_column %in% colnames(.)) .data[[fold_change_column]] else NA_real_,
        feature_log2_fc = if (log2_fc_column %in% colnames(.)) .data[[log2_fc_column]] else NA_real_,
        node_label = vapply(seq_len(n()), function(i) choose_annotation_label(.[i, , drop = FALSE]), character(1)),
        is_hub = feature_id == module_summary_row$hub_feature_id[1]
      )
    module_edges <- module_edge_tbl %>%
      filter(module_id == !!module_id) %>%
      mutate(abs_delta_mz = abs(delta_mz))

    vertex_export_tbl <- module_vertices %>%
      transmute(
        name = feature_id,
        feature_id = feature_id,
        module_id = module_id,
        component_index = component_index,
        weighted_degree = weighted_degree,
        spectral_weighted_degree = spectral_weighted_degree,
        degree = degree,
        discriminant_score = discriminant_score,
        abs_discriminant_score = abs_discriminant_score,
        discriminant_direction = discriminant_direction,
        discriminant_relevance = discriminant_relevance,
        loading = loading,
        abs_loading = abs_loading,
        feature_p_value = feature_p_value,
        feature_fold_change = feature_fold_change,
        feature_log2_fc = feature_log2_fc,
        annotation_label = node_label,
        is_hub = is_hub,
        feature_id_full = if ("feature_id_full" %in% colnames(module_vertices)) feature_id_full else NA_character_,
        sirius_chebiasciiname = if ("sirius_chebiasciiname" %in% colnames(module_vertices)) sirius_chebiasciiname else NA_character_,
        sirius_name = if ("sirius_name" %in% colnames(module_vertices)) sirius_name else NA_character_,
        canopus_npc_pathway = if ("canopus_npc_pathway" %in% colnames(module_vertices)) canopus_npc_pathway else NA_character_,
        canopus_npc_superclass = if ("canopus_npc_superclass" %in% colnames(module_vertices)) canopus_npc_superclass else NA_character_,
        canopus_npc_class = if ("canopus_npc_class" %in% colnames(module_vertices)) canopus_npc_class else NA_character_
      )

    edge_export_tbl <- module_edges %>%
      transmute(
        from = from,
        to = to,
        module_id = module_id,
        component_index = component_index,
        cosine = cosine,
        cluster_weight = cluster_weight,
        from_discriminant_score = from_discriminant_score,
        to_discriminant_score = to_discriminant_score,
        from_discriminant_direction = from_discriminant_direction,
        to_discriminant_direction = to_discriminant_direction,
        delta_mz = delta_mz,
        abs_delta_mz = abs_delta_mz
      )

    module_graph <- graph_from_data_frame(
      d = edge_export_tbl,
      directed = FALSE,
      vertices = vertex_export_tbl
    )

    layout_coords <- if (vcount(module_graph) > 1) {
      layout_with_fr(module_graph, weights = E(module_graph)$cosine)
    } else {
      matrix(c(0, 0), ncol = 2)
    }
    V(module_graph)$x <- layout_coords[, 1]
    V(module_graph)$y <- layout_coords[, 2]

    module_graph <- set_graph_attr(module_graph, "module_id", module_id)
    module_graph <- set_graph_attr(module_graph, "component_index", module_summary_row$component_index[1])
    module_graph <- set_graph_attr(module_graph, "module_mode", module_mode)
    module_graph <- set_graph_attr(module_graph, "discriminant_stat", if (module_mode == "discriminant") discriminant_stat else "")
    module_graph <- set_graph_attr(module_graph, "same_direction_only", same_direction_only)
    module_graph <- set_graph_attr(module_graph, "module_p_value", module_summary_row$module_p_value[1])
    module_graph <- set_graph_attr(module_graph, "module_q_value", module_summary_row$module_q_value[1])
    module_graph <- set_graph_attr(module_graph, "effect_metric", module_summary_row$effect_metric[1])
    module_graph <- set_graph_attr(module_graph, "effect_size", module_summary_row$effect_size[1])
    module_graph <- set_graph_attr(module_graph, "target_factor", target_factor)
    module_graph <- set_graph_attr(module_graph, "comparison_prefix", comparison_prefix)

    graphml_path <- file.path(graphml_dir, paste0(module_id, ".graphml"))
    write_graph(module_graph, file = graphml_path, format = "graphml")

    node_colors <- make_diverging_colors(V(module_graph)$feature_log2_fc)
    node_sizes <- rescale_numeric(V(module_graph)$weighted_degree, to = c(14, 28), default = 18)
    edge_widths <- rescale_numeric(E(module_graph)$cosine, to = c(1, 6), default = 2)
    edge_alpha <- rescale_numeric(E(module_graph)$cosine, to = c(0.25, 0.9), default = 0.5)
    edge_colors <- vapply(edge_alpha, function(alpha_val) {
      grDevices::adjustcolor("#4d4d4d", alpha.f = alpha_val)
    }, character(1))
    label_cex <- if (vcount(module_graph) <= 8) 0.8 else if (vcount(module_graph) <= 15) 0.65 else 0.5

    preview_path <- file.path(network_plot_dir, paste0(module_id, ".png"))
    grDevices::png(filename = preview_path, width = 1600, height = 1200, res = 160)
    old_par <- par(no.readonly = TRUE)
    par(mar = c(2, 2, 5, 2))
    plot(
      module_graph,
      layout = layout_coords,
      vertex.size = node_sizes,
      vertex.color = node_colors,
      vertex.label = V(module_graph)$annotation_label,
      vertex.label.cex = label_cex,
      vertex.label.family = "sans",
      vertex.label.color = "#1f1f1f",
      vertex.frame.color = ifelse(V(module_graph)$is_hub, "#111111", "#666666"),
      vertex.frame.width = ifelse(V(module_graph)$is_hub, 2.5, 1),
      edge.width = edge_widths,
      edge.color = edge_colors,
      main = sprintf("%s | %s | q=%.3g | effect=%.3g", module_id, module_mode, module_summary_row$module_q_value[1], module_summary_row$effect_size[1])
    )
    mtext(
      sprintf(
        "Node color: %s | range %.2f to %.2f | Edge width: cosine %.2f to %.2f | clustering weight %.2f to %.2f",
        log2_fc_column,
        safe_stat(V(module_graph)$feature_log2_fc, min, default = 0),
        safe_stat(V(module_graph)$feature_log2_fc, max, default = 0),
        safe_stat(E(module_graph)$cosine, min, default = 0),
        safe_stat(E(module_graph)$cosine, max, default = 0),
        safe_stat(E(module_graph)$cluster_weight, min, default = 0),
        safe_stat(E(module_graph)$cluster_weight, max, default = 0)
      ),
      side = 3,
      line = 0.5,
      cex = 0.75
    )
    mtext("Hub nodes have thicker borders; labels come from Sirius/Canopus when available.", side = 1, line = 0.25, cex = 0.7)
    par(old_par)
    grDevices::dev.off()
  }
}

message(sprintf("Saved module summary to %s", summary_path))
message(sprintf("Saved module members to %s", members_path))
message(sprintf("Saved module scores to %s", scores_path))
if (length(plot_module_ids)) {
  message(sprintf("Saved top-module plot to %s", plot_path))
  message(sprintf("Saved GraphML module networks to %s", graphml_dir))
  message(sprintf("Saved network preview plots to %s", network_plot_dir))
}
