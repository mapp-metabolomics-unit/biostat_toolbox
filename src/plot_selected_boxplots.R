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
  library(viridisLite)
  library(wesanderson)
  library(digest)
  library(structToolbox)
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
  make_option(c("-r", "--results-dir"), default = NULL, help = "Override directory that contains DE.rds and foldchange_pvalues.csv"),
  make_option(c("-f", "--features"), default = NULL, help = "Comma separated list of feature IDs to plot"),
  make_option(c("-F", "--features-file"), default = NULL, help = "Optional file that lists feature IDs (one per line)"),
  make_option(c("-o", "--output-dir"), default = NULL, help = "Directory used to save plots and data (default results_dir/selected_boxplots_cli)"),
  make_option(c("--pvalue-column"), default = NULL, help = "Name of the p-value column to use (default: first *_p_value column)"),
  make_option(c("--data-output"), default = NULL, help = "Optional path to save the long-format table (default output_dir/selected_boxplot_data.csv)"),
  make_option(c("-t", "--plot-type"), default = "box", help = "Plot geometry: box, violin, or violin_box [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

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
explicit_plot_type <- extract_flag_value(raw_cli_args, "--plot-type", "-t")

if (!is.null(explicit_results_dir)) {
  opt$results_dir <- explicit_results_dir
}
if (!is.null(explicit_params)) {
  opt$params <- explicit_params
}
if (!is.null(explicit_params_user)) {
  opt$params_user <- explicit_params_user
}
if (!is.null(explicit_plot_type)) {
  opt$plot_type <- explicit_plot_type
}

trim_option_value <- function(value) {
  if (is.character(value) && length(value)) {
    trimmed <- trimws(value)
    trimmed[nzchar(trimmed)]
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
opt$results_dir <- trim_option_value(opt$results_dir)
opt$features <- trim_option_value(opt$features)
opt$features_file <- trim_option_value(opt$features_file)
opt$output_dir <- trim_option_value(opt$output_dir)
opt$pvalue_column <- trim_option_value(opt$pvalue_column)
opt$data_output <- trim_option_value(opt$data_output)
opt$plot_type <- trim_option_value(opt$plot_type)

if (is.null(opt$params) || !nzchar(opt$params)) {
  opt$params <- "params/params.yaml"
}
if (is.null(opt$params_user) || !nzchar(opt$params_user)) {
  opt$params_user <- "params/params_user.yaml"
}
if (is.null(opt$plot_type) || !nzchar(opt$plot_type)) {
  opt$plot_type <- "box"
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

collect_feature_ids <- function(opt) {
  ids <- character()
  if (!is.null(opt$features) && nzchar(opt$features)) {
    ids <- c(ids, trimws(unlist(strsplit(opt$features, ","))))
  }
  if (!is.null(opt$features_file) && file.exists(opt$features_file)) {
    ids <- c(ids, trimws(readr::read_lines(opt$features_file)))
  } else if (!is.null(opt$features_file) && !file.exists(opt$features_file)) {
    stop(sprintf("Feature file %s does not exist.", opt$features_file))
  }
  ids <- unique(ids[nzchar(ids)])
  if (!length(ids)) {
    stop("No feature ids provided via --features or --features-file.")
  }
  ids
}

feature_ids <- collect_feature_ids(opt)

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

valid_plot_types <- c("box", "violin", "violin_box")
plot_type <- tolower(opt$plot_type)
if (!plot_type %in% valid_plot_types) {
  stop(sprintf("Invalid plot type '%s'. Choose one of: %s", opt$plot_type, paste(valid_plot_types, collapse = ", ")))
}

config_hash <- convert_yaml_to_single_row_df_with_hash(params)$hash

resolve_results_dir <- function(params, hash, override) {
  if (has_value(override)) {
    return(trimws(as.character(override[1])))
  }
  if (!is.null(params$paths$output) && nzchar(params$paths$output)) {
    return(file.path(params$paths$output, hash))
  }
  file.path(params$paths$docs, params$mapp_project, params$mapp_batch, "results", "stats", hash)
}

results_dir <- resolve_results_dir(params, config_hash, opt$results_dir)
if (has_value(opt$results_dir)) {
  message(sprintf("Using user-specified results directory: %s", results_dir))
} else {
  message(sprintf("Using inferred results directory: %s", results_dir))
}
if (!dir.exists(results_dir)) {
  stop(sprintf("Results directory %s does not exist.", results_dir))
}

output_dir <- opt$output_dir
if (is.null(output_dir) || !nzchar(output_dir)) {
  output_dir <- file.path(results_dir, "selected_boxplots_cli")
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

data_output_path <- opt$data_output
if (is.null(data_output_path) || !nzchar(data_output_path)) {
  data_output_path <- file.path(output_dir, "selected_boxplot_data.csv")
}

de_path <- file.path(results_dir, "DE.rds")
foldchange_path <- file.path(results_dir, "foldchange_pvalues.csv")

if (!file.exists(de_path)) {
  stop(sprintf("Missing DE.rds at %s", de_path))
}
if (!file.exists(foldchange_path)) {
  stop(sprintf("Missing foldchange_pvalues.csv at %s", foldchange_path))
}

message("Loading DE object...")
DE <- readRDS(de_path)
data_matrix <- as.data.frame(DE$data)
data_feature_names <- colnames(data_matrix)
if (!length(data_feature_names)) {
  stop("DE$data does not contain any features.")
}

message("Loading fold-change/p-value table...")
foldchange_tbl <- readr::read_csv(foldchange_path, show_col_types = FALSE)

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

if (!"feature_id" %in% names(foldchange_tbl)) {
  stop("Column feature_id is missing from foldchange_pvalues.csv")
}

foldchange_tbl <- foldchange_tbl %>%
  mutate(feature_id = as.character(feature_id))

missing_in_stats <- setdiff(feature_ids, foldchange_tbl$feature_id)
if (length(missing_in_stats)) {
  warning(sprintf(
    "The following requested feature IDs were not found in foldchange_pvalues.csv and will be skipped: %s",
    paste(missing_in_stats, collapse = ", ")
  ))
}

selected_features <- unique(feature_ids[feature_ids %in% foldchange_tbl$feature_id])
if (!length(selected_features)) {
  stop("None of the requested feature IDs were found in foldchange_pvalues.csv.")
}

missing_in_data <- setdiff(selected_features, data_feature_names)
if (length(missing_in_data)) {
  warning(sprintf(
    "The following features were found in foldchange_pvalues.csv but not in the DE object and will be skipped: %s",
    paste(missing_in_data, collapse = ", ")
  ))
}

selected_features <- setdiff(selected_features, missing_in_data)
if (!length(selected_features)) {
  stop("None of the requested feature IDs are available in the DE object.")
}

foldchange_selected <- foldchange_tbl %>%
  filter(feature_id %in% selected_features)

target_factor <- params$target$sample_metadata_header
if (!target_factor %in% colnames(DE$sample_meta)) {
  stop(sprintf("Sample metadata header %s is not present in DE$sample_meta.", target_factor))
}

if (!all(selected_features %in% colnames(data_matrix))) {
  missing_in_data <- setdiff(selected_features, colnames(data_matrix))
  stop(sprintf("The following feature IDs are missing from DE$data: %s", paste(missing_in_data, collapse = ", ")))
}

data_subset <- data_matrix[, selected_features, drop = FALSE]
colnames(data_subset) <- paste0("X", colnames(data_subset))

merged_df <- merge(DE$sample_meta, data_subset, by = "row.names", sort = FALSE)

merged_df[[target_factor]] <- factor(merged_df[[target_factor]])

data_subset_for_boxplots <- merged_df %>%
  select(Row.names, all_of(target_factor), starts_with("X")) %>%
  rename_with(~ sub("^X", "", .x), starts_with("X"))

plot_data <- data_subset_for_boxplots %>%
  pivot_longer(
    cols = -all_of(c("Row.names", target_factor)),
    names_to = "variable",
    values_to = "value"
  ) %>%
  select(-Row.names) %>%
  left_join(foldchange_selected, by = c("variable" = "feature_id"))

if (!nrow(plot_data)) {
  stop("No rows available after joining measurement data with metadata.")
}

build_custom_colors <- function(sample_meta, color_params) {
  factor_levels <- unique(sample_meta)
  if (isTRUE(color_params$continuous)) {
    sorted_levels <- as.character(sort(as.numeric(factor_levels)))
    return(stats::setNames(viridisLite::viridis(length(sorted_levels)), sorted_levels))
  }

  color_keys <- color_params$all$key
  color_values <- color_params$all$value

  if (length(color_keys)) {
    if (!all(color_keys %in% factor_levels)) {
      missing_cols <- color_keys[!color_keys %in% factor_levels]
      factor_str <- paste(unique(factor_levels), collapse = ", ")
      stop(sprintf(
        "Values %s defined under params$colors$all$key are missing from the metadata. Available groups: %s",
        paste(missing_cols, collapse = ", "),
        factor_str
      ))
    }
    return(stats::setNames(color_values, color_keys))
  }

  palette_values <- unique(unlist(wes_palettes[names(wes_palettes)]))
  if (length(palette_values) < length(factor_levels)) {
    palette_values <- grDevices::colorRampPalette(palette_values)(length(factor_levels))
  } else {
    palette_values <- sample(palette_values, length(factor_levels))
  }
  names(palette_values) <- factor_levels
  palette_values
}

custom_colors <- build_custom_colors(merged_df[[target_factor]], params$colors)

feature_value_table <- as.data.frame(t(data_matrix[, selected_features, drop = FALSE]))
feature_value_table <- tibble::rownames_to_column(feature_value_table, var = "feature_id")
readr::write_csv(feature_value_table, data_output_path)
message(sprintf("Saved feature intensity table to %s", data_output_path))

safe_filename <- function(x) {
  slug <- stringr::str_replace_all(x, "[^[:alnum:]]+", "_")
  slug <- stringr::str_replace_all(slug, "_+", "_")
  slug <- stringr::str_replace(slug, "^_+", "")
  slug <- stringr::str_replace(slug, "_+$", "")
  if (!nzchar(slug)) {
    slug <- "feature"
  }
  slug
}

target_sym <- sym(target_factor)

for (feature_id in selected_features) {
  data_for_plot <- plot_data %>% filter(variable == feature_id)
  if (!nrow(data_for_plot)) {
    warning(sprintf("No data available to plot feature %s", feature_id))
    next
  }

  rounded_p_value <- round(first(data_for_plot[[p_value_column]]), 5)
  feature_name_col <- "sirius_chebiasciiname"
  feature_details_col <- "feature_id_full"
  feature_name <- first(data_for_plot[[feature_name_col]])
  if (is.null(feature_name) || is.na(feature_name) || length(feature_name) == 0) {
    feature_name <- "NA"
  } else if (!is.character(feature_name)) {
    feature_name <- as.character(feature_name)
  }
  feature_details <- first(data_for_plot[[feature_details_col]])
  if (is.null(feature_details) || is.na(feature_details) || length(feature_details) == 0) {
    feature_details <- "NA"
  } else if (!is.character(feature_details)) {
    feature_details <- as.character(feature_details)
  }
  feature_title_label <- feature_id

  group_stats <- data_for_plot %>%
    group_by(!!target_sym) %>%
    summarise(
      count = n(),
      mean_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  count_lookup <- setNames(group_stats$count, group_stats[[target_factor]])
  fc_text <- "Fold-change unavailable (requires ≥2 groups)"
  if (nrow(group_stats) >= 2) {
    g1 <- levels(data_for_plot[[target_factor]])[1]
    g2 <- levels(data_for_plot[[target_factor]])[2]
    if (!is.na(g1) && !is.na(g2) && g1 %in% group_stats[[target_factor]] && g2 %in% group_stats[[target_factor]]) {
      mean_g1 <- group_stats$mean_value[group_stats[[target_factor]] == g1]
      mean_g2 <- group_stats$mean_value[group_stats[[target_factor]] == g2]
      if (length(mean_g1) && length(mean_g2) && mean_g2 != 0) {
        fc_value <- mean_g1 / mean_g2
        log2_fc <- log2(fc_value)
        fc_text <- paste0(
          "Fold-change (", g1, "/", g2, ") = ",
          signif(fc_value, 3),
          " (log2=", signif(log2_fc, 3), ")"
        )
      }
    }
  }
  caption_text <- paste(
    paste("p-value (", p_value_column, ") =", rounded_p_value),
    fc_text,
    paste("Group counts:", paste(
      paste(names(count_lookup), paste0("n=", count_lookup)),
      collapse = "; "
    )),
    sep = " | "
  )
  label_fn <- function(x) {
    counts <- count_lookup[x]
    counts[is.na(counts)] <- 0
    paste0(x, " (n=", counts, ")")
  }

  palette_for_plot <- custom_colors[names(custom_colors) %in% levels(merged_df[[target_factor]])]
  missing_levels <- setdiff(levels(merged_df[[target_factor]]), names(palette_for_plot))
  if (length(missing_levels)) {
    stop(sprintf("Missing colors for levels: %s", paste(missing_levels, collapse = ", ")))
  }

  plot_obj <- ggplot(data_for_plot, aes(x = !!target_sym, y = value, fill = !!target_sym))

  if (plot_type == "violin") {
    plot_obj <- plot_obj +
      geom_violin(trim = FALSE, alpha = 0.7)
  } else if (plot_type == "violin_box") {
    plot_obj <- plot_obj +
      geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
      geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8)
  } else {
    plot_obj <- plot_obj +
      geom_boxplot()
  }

  plot_obj <- plot_obj +
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
    labs(
      x = target_factor,
      y = "Normalized Intensity",
      title = paste("Compared intensities for feature:", feature_title_label),
      subtitle = paste(
        "\n",
        "Compound name: ", feature_name, "\n",
        "Feature details: ", feature_details
      ),
      caption = caption_text
    ) +
    theme(
      plot.caption = element_text(hjust = 0, face = "italic"),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.position = "bottom",
      legend.box = "vertical"
    ) +
    scale_fill_manual(
      name = "Groups",
      values = palette_for_plot,
      labels = label_fn
    ) +
    guides(
      fill = guide_legend(
        nrow = if (length(levels(merged_df[[target_factor]])) <= 3) 1 else ceiling(length(levels(merged_df[[target_factor]])) / 3),
        byrow = TRUE
      )
    )

  file_prefix <- switch(
    plot_type,
    "violin" = "violinplot_",
    "violin_box" = "violinboxplot_",
    "boxplot_"
  )
  file_name <- file.path(output_dir, paste0(file_prefix, safe_filename(feature_id), ".png"))
  ggsave(plot = plot_obj, filename = file_name, width = 8, height = 8)
  message(sprintf("Saved %s plot for %s to %s", plot_type, feature_id, file_name))
}
