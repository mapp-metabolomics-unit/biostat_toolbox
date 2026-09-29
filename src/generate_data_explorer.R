#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(htmltools)
  library(janitor)
  library(jsonlite)
  library(optparse)
  library(plotly)
  library(readr)
  library(yaml)
})

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
if (length(script_path)) {
  script_path <- normalizePath(script_path[1])
} else {
  script_path <- normalizePath(file.path(getwd(), "src", "generate_data_explorer.R"), mustWork = FALSE)
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

option_list <- list(
  make_option(c("-p", "--params"), default = file.path(repo_root, "params", "params.yaml"), help = "Path to params.yaml"),
  make_option(c("-u", "--params-user"), default = file.path(repo_root, "params", "params_user.yaml"), help = "Path to params_user.yaml"),
  make_option(c("-o", "--output-dir"), default = NULL, help = "Directory where data_explorer.html is written"),
  make_option(c("--max-features"), default = 5000, type = "integer", help = "Maximum number of features exported to the browser payload"),
  make_option(c("--chunk-size"), default = 250, type = "integer", help = "Number of features per lazy-loaded intensity chunk")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

normalize_option_name <- function(opt, underscore_name, hyphen_name) {
  if (is.null(opt[[underscore_name]]) && !is.null(opt[[hyphen_name]])) {
    opt[[underscore_name]] <- opt[[hyphen_name]]
  }
  opt
}

opt <- normalize_option_name(opt, "params_user", "params-user")
opt <- normalize_option_name(opt, "output_dir", "output-dir")
opt <- normalize_option_name(opt, "max_features", "max-features")
opt <- normalize_option_name(opt, "chunk_size", "chunk-size")

has_value <- function(value) {
  !is.null(value) && length(value) && !is.na(value[1]) && nzchar(trimws(as.character(value[1])))
}

resolve_path <- function(path_value, fallback_dir = getwd(), must_work = FALSE) {
  if (!has_value(path_value)) {
    return(NULL)
  }
  path_value <- trimws(as.character(path_value[1]))
  if (grepl("^/", path_value)) {
    return(normalizePath(path_value, mustWork = must_work))
  }
  if (file.exists(path_value)) {
    return(normalizePath(path_value, mustWork = must_work))
  }
  normalizePath(file.path(fallback_dir, path_value), mustWork = must_work)
}

json_for_script <- function(value) {
  json <- jsonlite::toJSON(value, dataframe = "rows", auto_unbox = TRUE, na = "null", digits = 10)
  gsub("</", "<\\/", as.character(json), fixed = TRUE)
}

first_existing_column <- function(data, candidates) {
  found <- intersect(candidates, colnames(data))
  if (!length(found)) {
    return(NULL)
  }
  found[1]
}

read_params <- function(params_path, params_user_path) {
  params <- yaml.load_file(params_path)
  params_user <- yaml.load_file(params_user_path)
  params$paths$docs <- params_user$paths$docs
  params$paths$output <- params_user$paths$output
  params
}

load_quant_matrix <- function(working_directory, params) {
  quant_name <- if (identical(as.character(params$actions$run_with_gap_filled), "TRUE")) {
    paste0(params$mapp_batch, "_gf_quant.csv")
  } else {
    paste0(params$mapp_batch, "_quant.csv")
  }
  quant_path <- file.path(working_directory, "results", "mzmine", quant_name)
  if (!file.exists(quant_path)) {
    stop(sprintf("Quant table not found: %s", quant_path))
  }
  feature_table <- readr::read_delim(quant_path, delim = ",", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  feature_table <- feature_table %>%
    rename(
      feature_id = `row ID`,
      feature_mz = `row m/z`,
      feature_rt = `row retention time`
    )
  feature_table$feature_id <- as.character(feature_table$feature_id)
  feature_table$feature_id_full <- paste(feature_table$feature_id, round(feature_table$feature_mz, 2), round(feature_table$feature_rt, 1), sep = "_")

  intensity_columns <- grep(" Peak area$| Peak height$", colnames(feature_table), value = TRUE)
  if (!length(intensity_columns)) {
    stop("No intensity columns ending in ' Peak area' or ' Peak height' were found in the quant table.")
  }
  if (any(grepl(" Peak area$", colnames(feature_table))) && any(grepl(" Peak height$", colnames(feature_table)))) {
    intensity_columns <- grep(" Peak area$", colnames(feature_table), value = TRUE)
  }
  intensity_names <- gsub(" Peak area$| Peak height$", "", intensity_columns)
  intensity_matrix <- as.data.frame(t(as.matrix(feature_table[, intensity_columns, drop = FALSE])))
  rownames(intensity_matrix) <- intensity_names
  colnames(intensity_matrix) <- feature_table$feature_id
  intensity_matrix[] <- lapply(intensity_matrix, function(column) suppressWarnings(as.numeric(column)))

  feature_metadata <- feature_table %>%
    select(feature_id, feature_id_full, feature_mz, feature_rt)

  list(data = intensity_matrix, feature_metadata = feature_metadata)
}

load_sample_metadata <- function(working_directory, params) {
  metadata_path <- file.path(working_directory, "metadata", "treated", paste(params$mapp_batch, "metadata.tsv", sep = "_"))
  if (!file.exists(metadata_path)) {
    stop(sprintf("Sample metadata not found: %s", metadata_path))
  }
  sample_metadata <- readr::read_delim(metadata_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  sample_metadata <- as.data.frame(sample_metadata) %>% janitor::clean_names(case = "snake")
  if (!"filename" %in% colnames(sample_metadata)) {
    stop("Sample metadata must contain a filename column.")
  }
  rownames(sample_metadata) <- sample_metadata$filename
  sample_metadata
}

load_canopus_metadata <- function(working_directory) {
  canopus_path <- if (file.exists(file.path(working_directory, "results", "sirius", "canopus_structure_summary.tsv"))) {
    file.path(working_directory, "results", "sirius", "canopus_structure_summary.tsv")
  } else {
    file.path(working_directory, "results", "sirius", "canopus_compound_summary.tsv")
  }
  if (!file.exists(canopus_path)) {
    warning(sprintf("CANOPUS table not found: %s", canopus_path))
    return(data.frame(feature_id = character()))
  }
  canopus <- readr::read_delim(canopus_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  canopus <- as.data.frame(canopus) %>% janitor::clean_names(case = "snake")
  feature_id_column <- first_existing_column(canopus, c("mapping_feature_id", "feature_id", "id"))
  if (is.null(feature_id_column)) {
    warning("CANOPUS table has no recognized feature id column.")
    return(data.frame(feature_id = character()))
  }
  canopus$feature_id <- as.character(canopus[[feature_id_column]])
  npc_pathway_column <- first_existing_column(canopus, c("npc_pathway", "npc_number_pathway"))
  npc_superclass_column <- first_existing_column(canopus, c("npc_superclass", "npc_number_superclass"))
  npc_class_column <- first_existing_column(canopus, c("npc_class", "npc_number_class"))
  npc_pathway_probability_column <- first_existing_column(canopus, c("npc_pathway_probability", "npc_number_pathway_probability"))
  npc_superclass_probability_column <- first_existing_column(canopus, c("npc_superclass_probability", "npc_number_superclass_probability"))
  npc_class_probability_column <- first_existing_column(canopus, c("npc_class_probability", "npc_number_class_probability"))
  canopus$npc_pathway <- if (!is.null(npc_pathway_column)) as.character(canopus[[npc_pathway_column]]) else NA_character_
  canopus$npc_superclass <- if (!is.null(npc_superclass_column)) as.character(canopus[[npc_superclass_column]]) else NA_character_
  canopus$npc_class <- if (!is.null(npc_class_column)) as.character(canopus[[npc_class_column]]) else NA_character_
  canopus$npc_pathway_probability <- if (!is.null(npc_pathway_probability_column)) suppressWarnings(as.numeric(canopus[[npc_pathway_probability_column]])) else NA_real_
  canopus$npc_superclass_probability <- if (!is.null(npc_superclass_probability_column)) suppressWarnings(as.numeric(canopus[[npc_superclass_probability_column]])) else NA_real_
  canopus$npc_class_probability <- if (!is.null(npc_class_probability_column)) suppressWarnings(as.numeric(canopus[[npc_class_probability_column]])) else NA_real_
  keep_columns <- intersect(
    c(
      "feature_id",
      "npc_pathway",
      "npc_superclass",
      "npc_class",
      "npc_pathway_probability",
      "npc_superclass_probability",
      "npc_class_probability",
      "molecular_formula",
      "adduct"
    ),
    colnames(canopus)
  )
  canopus[, keep_columns, drop = FALSE]
}

load_archived_volcano_runs <- function(output_dir) {
  archive_root <- file.path(output_dir, "reprocess_by_extraction_phase")
  if (!dir.exists(archive_root)) {
    return(list())
  }
  run_dirs <- list.dirs(archive_root, recursive = TRUE, full.names = TRUE)
  run_dirs <- run_dirs[
    file.exists(file.path(run_dirs, "session_info.txt")) &
      file.exists(file.path(run_dirs, "params.yaml")) &
      file.exists(file.path(run_dirs, "foldchange_pvalues.csv"))
  ]
  archived_runs <- list()
  for (run_dir in run_dirs) {
    run_params <- tryCatch(yaml::yaml.load_file(file.path(run_dir, "params.yaml")), error = function(error) NULL)
    if (is.null(run_params)) next
    stats_table <- tryCatch(
      readr::read_csv(file.path(run_dir, "foldchange_pvalues.csv"), show_col_types = FALSE),
      error = function(error) NULL
    )
    if (is.null(stats_table) || !nrow(stats_table)) next
    feature_column <- first_existing_column(stats_table, c("feature_id", "row_id"))
    if (is.null(feature_column)) next
    p_columns <- grep("_p_value$", colnames(stats_table), value = TRUE)
    if (!length(p_columns)) next
    sample_file <- file.path(run_dir, "formatted_sample_metadata.tsv")
    sample_ids <- character(0)
    if (file.exists(sample_file)) {
      run_samples <- tryCatch(readr::read_tsv(sample_file, show_col_types = FALSE), error = function(error) NULL)
      if (!is.null(run_samples)) {
        sample_column <- first_existing_column(run_samples, c("filename", "sample_id"))
        if (!is.null(sample_column)) sample_ids <- as.character(run_samples[[sample_column]])
      }
    }
    group_column <- janitor::make_clean_names(as.character(run_params$target$sample_metadata_header), case = "snake")
    phase <- basename(dirname(run_dir))
    run_hash <- basename(run_dir)
    for (p_column in p_columns) {
      contrast <- sub("_p_value$", "", p_column)
      fold_column <- paste0(contrast, "_fold_change_log2")
      if (!fold_column %in% colnames(stats_table)) next
      groups <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
      if (length(groups) != 2) next
      archived_runs[[length(archived_runs) + 1]] <- list(
        id = paste(phase, run_hash, contrast, sep = "/"),
        phase = phase,
        run_hash = run_hash,
        contrast = contrast,
        group_column = group_column,
        groups = as.list(groups),
        sample_ids = as.list(sort(unique(sample_ids))),
        configured_p_value = if (!is.null(run_params$posthoc$p_value)) as.numeric(run_params$posthoc$p_value) else NA_real_,
        feature_ids = as.list(as.character(stats_table[[feature_column]])),
        p_values = as.list(as.numeric(stats_table[[p_column]])),
        log2_fold_changes = as.list(as.numeric(stats_table[[fold_column]]))
      )
    }
  }
  archived_runs
}

load_component_metadata <- function(working_directory, params) {
  candidates <- c(
    file.path(
      working_directory,
      "results",
      "tmp",
      paste0(params$mapp_batch, "_met_annot_unified_horizontal.tsv")
    ),
    file.path(
      working_directory,
      "results",
      "met_annot_enhancer",
      params$mapp_batch,
      paste0(params$mapp_batch, "_spectral_match_results_repond_flat.tsv")
    ),
    file.path(
      working_directory,
      "results",
      "tmp",
      paste0(params$mapp_batch, "_met_annot_unified_vertical.tsv")
    )
  )
  component_path <- candidates[file.exists(candidates)][1]
  if (is.na(component_path)) {
    warning("No component metadata table found. Component index will default to -1.")
    return(data.frame(feature_id = character(), component_id = character()))
  }
  component_table <- readr::read_delim(component_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
  component_table <- as.data.frame(component_table) %>% janitor::clean_names(case = "snake")
  feature_id_column <- first_existing_column(component_table, c("feature_id", "mapping_feature_id", "id"))
  component_id_column <- first_existing_column(component_table, c("gnps_mn_component", "component_id", "gnps_component_id", "isdb_component_id"))
  if (is.null(feature_id_column) || is.null(component_id_column)) {
    warning(sprintf("Component table has no recognized feature/component columns: %s", component_path))
    return(data.frame(feature_id = character(), component_id = character()))
  }
  score_columns <- setdiff(
    grep("(probability|score|confidence|p_value|q_value|pvalue|qvalue|fdr)", colnames(component_table), ignore.case = TRUE, value = TRUE),
    c(feature_id_column, component_id_column)
  )
  score_columns <- score_columns[vapply(score_columns, function(column) {
    any(is.finite(suppressWarnings(as.numeric(component_table[[column]]))))
  }, logical(1))]
  structure_columns <- setdiff(
    grep("(smiles|inchi|inchikey)", colnames(component_table), ignore.case = TRUE, value = TRUE),
    c(feature_id_column, component_id_column, score_columns)
  )
  structure_columns <- structure_columns[vapply(structure_columns, function(column) {
    values <- as.character(component_table[[column]])
    any(!is.na(values) & nzchar(values) & toupper(values) != "NA")
  }, logical(1))]
  component_table <- component_table %>%
    transmute(
      feature_id = as.character(.data[[feature_id_column]]),
      component_id = as.character(.data[[component_id_column]]),
      across(all_of(score_columns), ~ suppressWarnings(as.numeric(.x))),
      across(all_of(structure_columns), as.character)
    ) %>%
    filter(!is.na(feature_id), nzchar(feature_id)) %>%
    mutate(component_id = ifelse(is.na(component_id) | !nzchar(component_id), "-1", component_id)) %>%
    distinct(feature_id, .keep_all = TRUE)
  component_table
}

metadata_columns_for_browser <- function(sample_metadata) {
  colnames(sample_metadata)[vapply(sample_metadata, function(column) {
    values <- unique(as.character(column))
    values <- values[!is.na(values) & nzchar(values)]
    length(values) > 0 && length(values) <= 250
  }, logical(1))]
}

get_fixed_npc_palettes <- function() {
  list(
    reference_pathway_colors = c(
      "Alkaloids" = "#514300",
      "Alkaloids x Amino acids and Peptides" = "#715e00",
      "Alkaloids x Terpenoids" = "#756101",
      "Amino acids and Peptides" = "#ca5a04",
      "Amino acids and Peptides x Polyketides" = "#d37f3e",
      "Amino acids and Peptides x Shikimates and Phenylpropanoids" = "#ca9f04",
      "Carbohydrates" = "#485f2f",
      "Fatty acids" = "#612ece",
      "Polyketides" = "#865993",
      "Polyketides x Terpenoids" = "#6a5c8a",
      "Shikimates and Phenylpropanoids" = "#6ba148",
      "Terpenoids" = "#63acf5",
      "Other" = "#B7B7B7",
      "Unclassified" = "#B7B7B7"
    ),
    reference_superclass_pathway = c(
      "Diterpenoids" = "Terpenoids",
      "Steroids" = "Terpenoids",
      "Triterpenoids" = "Terpenoids",
      "Fatty acyls" = "Fatty acids",
      "Fatty amides" = "Fatty acids",
      "Fatty esters" = "Fatty acids",
      "Glycerolipids" = "Fatty acids",
      "Glycerophospholipids" = "Fatty acids",
      "Sphingolipids" = "Fatty acids",
      "Linear polyketides" = "Polyketides",
      "Macrolides" = "Polyketides",
      "Pseudoalkaloids" = "Alkaloids x Terpenoids",
      "Pseudoalkaloids x" = "Alkaloids x Terpenoids"
    ),
    pathway_family = c(
      "Terpenoids" = "micro_cvd_purple",
      "Fatty acids" = "micro_cvd_blue",
      "Polyketides" = "micro_cvd_orange",
      "Alkaloids" = "micro_cvd_green",
      "Shikimates and Phenylpropanoids" = "micro_cvd_turquoise",
      "Amino acids and Peptides" = "micro_orange",
      "Carbohydrates" = "micro_purple",
      "Other" = "micro_cvd_gray",
      "Unclassified" = "micro_cvd_gray"
    ),
    shade_hex = list(
      micro_cvd_gray = c("#616161", "#8B8B8B", "#B7B7B7", "#D6D6D6", "#F5F5F5"),
      micro_cvd_purple = c("#7D3560", "#A1527F", "#CC79A7", "#E794C1", "#EFB6D6"),
      micro_cvd_blue = c("#098BD9", "#56B4E9", "#7DCCFF", "#BCE1FF", "#E7F4FF"),
      micro_cvd_orange = c("#9D654C", "#C17754", "#F09163", "#FCB076", "#FFD5AF"),
      micro_cvd_green = c("#4E7705", "#6D9F06", "#97CE2F", "#BDEC6F", "#DDFFA0"),
      micro_cvd_turquoise = c("#148F77", "#009E73", "#43BA8F", "#48C9B0", "#A3E4D7"),
      micro_orange = c("#ff7f00", "#fe9929", "#fdae6b", "#fec44f", "#feeda0"),
      micro_purple = c("#6a51a3", "#807dba", "#9e9ac8", "#bcbddc", "#dadaeb")
    )
  )
}

build_fixed_npc_color_maps <- function(feature_metadata) {
  palettes <- get_fixed_npc_palettes()
  ranked <- feature_metadata %>%
    mutate(
      npc_pathway = dplyr::coalesce(na_if(npc_pathway, ""), "Unclassified"),
      npc_superclass = dplyr::coalesce(na_if(npc_superclass, ""), "Other"),
      npc_class = dplyr::coalesce(na_if(npc_class, ""), "Other")
    )

  superclass_rank <- ranked %>%
    count(npc_pathway, npc_superclass, name = "abundance") %>%
    group_by(npc_pathway) %>%
    arrange(desc(abundance), npc_superclass, .by_group = TRUE) %>%
    mutate(shade_index = pmin(row_number(), 5L)) %>%
    mutate(shade_index = ifelse(row_number() > 4L, 5L, shade_index)) %>%
    ungroup() %>%
    mutate(
      palette_family = dplyr::coalesce(unname(palettes$pathway_family[npc_pathway]), "micro_cvd_gray"),
      hex = purrr::map2_chr(palette_family, shade_index, ~ palettes$shade_hex[[.x]][.y])
    )

  microshades_pathway_colors <- ranked %>%
    distinct(npc_pathway) %>%
    mutate(
      palette_family = dplyr::coalesce(unname(palettes$pathway_family[npc_pathway]), "micro_cvd_gray"),
      hex = purrr::map_chr(palette_family, ~ palettes$shade_hex[[.x]][1])
    ) %>%
    select(npc_pathway, hex) %>%
    tibble::deframe()

  microshades_superclass_colors <- superclass_rank %>%
    distinct(npc_superclass, hex) %>%
    tibble::deframe()

  microshades_class_colors <- ranked %>%
    distinct(npc_pathway, npc_superclass, npc_class) %>%
    left_join(
      superclass_rank %>% select(npc_pathway, npc_superclass, hex),
      by = c("npc_pathway", "npc_superclass")
    ) %>%
    mutate(hex = dplyr::coalesce(hex, microshades_pathway_colors[npc_pathway], "#F5F5F5")) %>%
    distinct(npc_class, hex) %>%
    tibble::deframe()

  reference_pathway_colors <- palettes$reference_pathway_colors
  reference_superclass_colors <- ranked %>%
    distinct(npc_pathway, npc_superclass) %>%
    mutate(
      reference_pathway = dplyr::coalesce(
        unname(palettes$reference_superclass_pathway[npc_superclass]),
        ifelse(npc_pathway %in% names(reference_pathway_colors), npc_pathway, NA_character_),
        "Other"
      ),
      hex = unname(reference_pathway_colors[reference_pathway])
    ) %>%
    distinct(npc_superclass, hex) %>%
    tibble::deframe()

  reference_class_colors <- ranked %>%
    distinct(npc_pathway, npc_superclass, npc_class) %>%
    mutate(
      reference_pathway = dplyr::coalesce(
        unname(palettes$reference_superclass_pathway[npc_superclass]),
        ifelse(npc_pathway %in% names(reference_pathway_colors), npc_pathway, NA_character_),
        "Other"
      ),
      hex = unname(reference_pathway_colors[reference_pathway])
    ) %>%
    distinct(npc_class, hex) %>%
    tibble::deframe()

  list(
    pathway = as.list(reference_pathway_colors),
    superclass = as.list(reference_superclass_colors),
    class = as.list(reference_class_colors),
    microshades_pathway = as.list(microshades_pathway_colors),
    microshades_superclass = as.list(microshades_superclass_colors),
    microshades_class = as.list(microshades_class_colors),
    reference_superclass_pathway = as.list(palettes$reference_superclass_pathway)
  )
}

params_path <- resolve_path(opt$params, repo_root, must_work = TRUE)
params_user_path <- resolve_path(opt$params_user, repo_root, must_work = TRUE)
params <- read_params(params_path, params_user_path)
working_directory <- file.path(params$paths$docs, params$mapp_project, params$mapp_batch)

output_dir <- if (has_value(opt$output_dir)) {
  resolve_path(opt$output_dir, getwd(), must_work = FALSE)
} else if (has_value(params$paths$output)) {
  normalizePath(params$paths$output, mustWork = FALSE)
} else {
  file.path(working_directory, "results", "stats")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

archive_stats_dir <- if (has_value(params$paths$output) && dir.exists(params$paths$output)) {
  normalizePath(params$paths$output)
} else {
  output_dir
}
archived_volcano_runs <- load_archived_volcano_runs(archive_stats_dir)
message(sprintf("Loaded %d archived Volcano contrast(s).", length(archived_volcano_runs)))
quant <- load_quant_matrix(working_directory, params)
sample_metadata <- load_sample_metadata(working_directory, params)
canopus_metadata <- load_canopus_metadata(working_directory)
component_metadata <- load_component_metadata(working_directory, params)

common_samples <- intersect(rownames(quant$data), rownames(sample_metadata))
if (!length(common_samples)) {
  stop("No quant table samples match metadata filenames.")
}
intensity_matrix <- quant$data[common_samples, , drop = FALSE]
sample_metadata <- sample_metadata[common_samples, , drop = FALSE]

feature_metadata <- quant$feature_metadata %>%
  left_join(canopus_metadata, by = "feature_id") %>%
  left_join(component_metadata, by = "feature_id")
feature_ids <- intersect(colnames(intensity_matrix), feature_metadata$feature_id)
if (!length(feature_ids)) {
  stop("No feature metadata rows match quant table feature ids.")
}
if (length(feature_ids) > opt$max_features) {
  feature_totals <- colSums(intensity_matrix[, feature_ids, drop = FALSE], na.rm = TRUE)
  feature_ids <- names(sort(feature_totals, decreasing = TRUE))[seq_len(opt$max_features)]
  warning(sprintf("Exporting top %d features by total intensity. Increase --max-features to export more.", opt$max_features))
}
intensity_matrix <- intensity_matrix[, feature_ids, drop = FALSE]
feature_metadata <- feature_metadata[match(feature_ids, feature_metadata$feature_id), , drop = FALSE]
feature_metadata$component_id[is.na(feature_metadata$component_id) | !nzchar(as.character(feature_metadata$component_id))] <- "-1"
feature_metadata$feature_total_intensity <- as.numeric(colSums(intensity_matrix[, feature_ids, drop = FALSE], na.rm = TRUE))
feature_metadata$feature_mean_intensity <- as.numeric(colMeans(intensity_matrix[, feature_ids, drop = FALSE], na.rm = TRUE))
feature_metadata$feature_max_intensity <- as.numeric(apply(intensity_matrix[, feature_ids, drop = FALSE], 2, max, na.rm = TRUE))

feature_label <- as.character(feature_metadata$feature_id)
has_class <- !is.na(feature_metadata$npc_class) & nzchar(as.character(feature_metadata$npc_class))
feature_label[has_class] <- paste(feature_metadata$feature_id[has_class], feature_metadata$npc_class[has_class], sep = " | ")
feature_metadata$feature_label <- feature_label
feature_metadata <- as.data.frame(feature_metadata)
character_columns <- colnames(feature_metadata)[vapply(feature_metadata, is.character, logical(1))]
for (column in character_columns) {
  feature_metadata[[column]][is.na(feature_metadata[[column]])] <- ""
}
feature_metadata_columns <- lapply(colnames(feature_metadata), function(column) {
  values <- feature_metadata[[column]]
  if (is.numeric(values)) {
    values[!is.finite(values)] <- NA_real_
  }
  as.list(values)
})
names(feature_metadata_columns) <- colnames(feature_metadata)

intensity_payload <- lapply(feature_ids, function(feature_id) {
  as.numeric(intensity_matrix[, feature_id])
})
names(intensity_payload) <- feature_ids

assets_dir <- file.path(output_dir, "data_explorer_assets")
chunks_dir <- file.path(assets_dir, "intensity_chunks")
if (dir.exists(chunks_dir)) {
  unlink(Sys.glob(file.path(chunks_dir, "*.js")))
}
dir.create(chunks_dir, recursive = TRUE, showWarnings = FALSE)

chunk_size <- max(1, as.integer(opt$chunk_size))
asset_version <- paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "-", Sys.getpid())
feature_chunks <- split(feature_ids, ceiling(seq_along(feature_ids) / chunk_size))
chunk_manifest <- lapply(seq_along(feature_chunks), function(index) {
  chunk_id <- sprintf("chunk_%04d", index)
  chunk_features <- feature_chunks[[index]]
  chunk_payload <- intensity_payload[chunk_features]
  chunk_file <- file.path(chunks_dir, paste0(chunk_id, ".js"))
  writeLines(
    paste0(
      "window.MAPP_DATA_EXPLORER_INTENSITY_CHUNKS = window.MAPP_DATA_EXPLORER_INTENSITY_CHUNKS || {};",
      "window.MAPP_DATA_EXPLORER_INTENSITY_CHUNKS[", jsonlite::toJSON(chunk_id, auto_unbox = TRUE), "] = ",
      json_for_script(chunk_payload),
      ";"
    ),
    con = chunk_file,
    useBytes = TRUE
  )
  list(
    id = chunk_id,
    path = paste0(file.path("data_explorer_assets", "intensity_chunks", paste0(chunk_id, ".js")), "?v=", asset_version),
    features = as.list(chunk_features)
  )
})
feature_chunk_map <- stats::setNames(
  rep(vapply(chunk_manifest, `[[`, character(1), "id"), lengths(feature_chunks)),
  unlist(feature_chunks, use.names = FALSE)
)

payload <- list(
  title = paste("Data explorer for", params$mapp_batch),
  mapp_project = params$mapp_project,
  mapp_batch = params$mapp_batch,
  asset_version = asset_version,
  default_group = janitor::make_clean_names(as.character(params$target$sample_metadata_header), case = "snake"),
  archived_volcano_runs = archived_volcano_runs,
  sample_ids = common_samples,
  sample_metadata = sample_metadata,
  metadata_columns = metadata_columns_for_browser(sample_metadata),
  feature_metadata_columns = feature_metadata_columns,
  npc_color_maps = build_fixed_npc_color_maps(feature_metadata),
  sample_totals = as.numeric(rowSums(intensity_matrix[, feature_ids, drop = FALSE], na.rm = TRUE)),
  intensity_chunks = chunk_manifest,
  feature_chunk_map = as.list(feature_chunk_map)
)

payload_file <- file.path(assets_dir, "data_explorer_payload.js")
css_file <- file.path(assets_dir, "data_explorer.css")
js_file <- file.path(assets_dir, "data_explorer.js")
html_file <- file.path(output_dir, "data_explorer.html")
writeLines(
  paste0("window.MAPP_DATA_EXPLORER = ", json_for_script(payload), ";"),
  con = payload_file,
  useBytes = TRUE
)

dummy_plotly <- plotly::plot_ly(x = 1, y = 1, type = "scatter", mode = "markers") %>%
  plotly::layout(width = 1, height = 1, margin = list(l = 0, r = 0, t = 0, b = 0)) %>%
  plotly::config(displaylogo = FALSE)

data_explorer_css <- r"(
:root {
  --rail-bg: #20242a;
  --rail-bg-soft: #2b3037;
  --chrome: #eef0f2;
  --chrome-strong: #e1e4e7;
  --panel: #f8f9fa;
  --panel-strong: #ffffff;
  --canvas: #dfe3e7;
  --viewport: #f5f6f7;
  --field: #ffffff;
  --line: #cfd4d9;
  --line-strong: #aeb6bf;
  --ink: #16191d;
  --ink-soft: #454b53;
  --ink-muted: #747d87;
  --accent: #1f6f78;
  --accent-2: #b7791f;
  --accent-soft: #dcebed;
  --shadow: 0 10px 28px rgba(20, 24, 28, 0.12);
  --shadow-soft: 0 1px 2px rgba(20, 24, 28, 0.08);
  --radius-sm: 3px;
  --radius-md: 5px;
  --rail-width: 56px;
  --sidebar-width: 340px;
  color-scheme: light;
}

* { box-sizing: border-box; }
html, body { width: 100%; height: 100%; margin: 0; overflow: hidden; }
body {
  background: var(--canvas);
  color: var(--ink);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
}

select, input, textarea {
  width: 100%;
  min-width: 0;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: var(--field);
  color: var(--ink);
  font: inherit;
  font-size: 0.78rem;
}
select, input {
  height: 32px;
  padding: 0 0.55rem;
}
textarea {
  min-height: 78px;
  resize: vertical;
  padding: 0.45rem 0.55rem;
  line-height: 1.35;
}
select[multiple] {
  height: auto;
  min-height: 78px;
  padding: 0.35rem 0.45rem;
}
select:focus, input:focus, textarea:focus {
  outline: 2px solid rgba(31, 111, 120, 0.24);
  border-color: var(--accent);
}
label { min-width: 0; }

.app-container {
  display: flex;
  width: 100vw;
  height: 100vh;
  overflow: hidden;
}

.app-rail {
  width: var(--rail-width);
  flex: 0 0 var(--rail-width);
  background: var(--rail-bg);
  color: #f8fafc;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 0.7rem;
  padding: 0.7rem 0.45rem;
  border-right: 1px solid rgba(255, 255, 255, 0.08);
}
.rail-mark {
  width: 34px;
  height: 34px;
  display: grid;
  place-items: center;
  border-radius: 7px;
  background: linear-gradient(135deg, var(--accent), #2d9aa4);
  font-weight: 800;
  font-size: 0.92rem;
  box-shadow: inset 0 0 0 1px rgba(255,255,255,0.18);
}
.rail-chip {
  width: 38px;
  height: 38px;
  display: grid;
  place-items: center;
  border-radius: 6px;
  background: transparent;
  color: #b9c1ca;
  font-size: 0.68rem;
  font-weight: 700;
  writing-mode: vertical-rl;
  transform: rotate(180deg);
  letter-spacing: 0.04em;
}
.rail-chip.is-active {
  background: var(--rail-bg-soft);
  color: #ffffff;
  box-shadow: inset 3px 0 0 var(--accent-2);
}
.rail-spacer { flex: 1; }
.rail-dot {
  width: 9px;
  height: 9px;
  border-radius: 50%;
  background: #34d399;
  box-shadow: 0 0 0 3px rgba(52, 211, 153, 0.15);
}

.app-sidebar {
  width: var(--sidebar-width);
  flex: 0 0 var(--sidebar-width);
  background: var(--chrome);
  border-right: 1px solid var(--line);
  display: flex;
  flex-direction: column;
  min-width: 0;
  transition: width 150ms ease, flex-basis 150ms ease;
}
.app-sidebar.is-collapsed {
  width: 44px;
  flex-basis: 44px;
}

.sidebar-header, .app-header-bar {
  height: 50px;
  flex-shrink: 0;
  border-bottom: 1px solid var(--line);
  display: flex;
  align-items: center;
  padding: 0 0.9rem;
}

.sidebar-header { justify-content: space-between; align-items: center; }
.sidebar-toggle {
  width: 28px;
  height: 26px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: #ffffff;
  color: var(--ink-soft);
  cursor: pointer;
  font: inherit;
  font-size: 0.76rem;
  font-weight: 820;
  flex: 0 0 28px;
}
.app-sidebar.is-collapsed .sidebar-header {
  padding: 0 0.45rem;
  justify-content: center;
}
.app-sidebar.is-collapsed .sidebar-header > div,
.app-sidebar.is-collapsed .sidebar-content {
  display: none;
}
.brand-title { font-size: 0.86rem; font-weight: 800; letter-spacing: 0.01em; }
.brand-subtitle { margin-top: 0.12rem; font-size: 0.68rem; color: var(--ink-muted); }

.sidebar-content {
  flex: 1;
  overflow-y: auto;
  padding: 0.55rem;
  display: flex;
  flex-direction: column;
  gap: 0.55rem;
}

.sidebar-card {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow-soft);
  padding: 0.65rem;
}
.card-header {
  display: flex;
  align-items: center;
  gap: 0.45rem;
  min-height: 22px;
  margin: -0.65rem -0.65rem 0.6rem;
  padding: 0.38rem 0.65rem;
  border-bottom: 1px solid var(--line);
  background: var(--chrome-strong);
  color: var(--ink-soft);
  text-transform: uppercase;
  letter-spacing: 0.055em;
  font-size: 0.66rem;
  font-weight: 800;
}
.form-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; }
.form-stack { display: grid; gap: 0.52rem; }
.sample-filter-row {
  display: grid;
  grid-template-columns: minmax(0, 0.95fr) minmax(0, 1.2fr);
  gap: 0.42rem;
  padding-bottom: 0.48rem;
  border-bottom: 1px dashed var(--line);
}
.sample-filter-row:last-child {
  padding-bottom: 0;
  border-bottom: none;
}
.numeric-filter-panel {
  display: grid;
  gap: 0.58rem;
}
.numeric-filter-panel:empty { display: none; }
.numeric-filter-row {
  display: grid;
  gap: 0.24rem;
}
.numeric-filter-head {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 0.4rem;
}
.numeric-filter-name {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  color: var(--ink-muted);
  font-size: 0.64rem;
  font-weight: 760;
  letter-spacing: 0.03em;
  text-transform: uppercase;
}
.numeric-filter-value {
  flex-shrink: 0;
  color: var(--ink-soft);
  font-size: 0.68rem;
  font-weight: 780;
}
input[type="range"] {
  height: 18px;
  padding: 0;
}
.form-group { display: grid; gap: 0.26rem; min-width: 0; }
.form-group > span, .checkbox-row span:first-child {
  font-size: 0.64rem;
  font-weight: 760;
  color: var(--ink-muted);
  text-transform: uppercase;
  letter-spacing: 0.03em;
}
.checkbox-row {
  display: flex;
  align-items: center;
  gap: 0.45rem;
  color: var(--ink-soft);
  font-size: 0.76rem;
}
.checkbox-row input { width: 14px; height: 14px; }

.app-main {
  flex: 1;
  min-width: 0;
  display: flex;
  flex-direction: column;
  overflow: hidden;
}

.app-header-bar {
  background: var(--panel-strong);
  justify-content: space-between;
  gap: 0.75rem;
  box-shadow: var(--shadow-soft);
  z-index: 3;
}
.header-title { min-width: 0; }
.header-title h1 { margin: 0; font-size: 0.9rem; line-height: 1.2; font-weight: 800; }
.header-summary { margin-top: 0.08rem; font-size: 0.68rem; color: var(--ink-muted); }
.header-search { width: min(460px, 44vw); }
.header-actions {
  display: flex;
  align-items: center;
  gap: 0.45rem;
  flex-shrink: 0;
}

.workspace {
  flex: 1;
  min-height: 0;
  overflow: hidden;
  padding: 0.9rem;
  display: flex;
  gap: 0.8rem;
  align-items: stretch;
  background:
    linear-gradient(rgba(255,255,255,0.45), rgba(255,255,255,0.45)),
    radial-gradient(circle at 1px 1px, rgba(120,128,136,0.26) 1px, transparent 0);
  background-size: auto, 18px 18px;
}
.workspace-content {
  flex: 1;
  min-width: 0;
  overflow-y: auto;
  padding-right: 0.1rem;
}
.plot-grid {
  display: grid;
  grid-template-columns: repeat(var(--grid-columns, 4), minmax(0, 1fr));
  gap: 0.78rem;
  align-content: start;
}
.workspace-tabs {
  display: flex;
  align-items: center;
  gap: 0.35rem;
  margin-bottom: 0.75rem;
  padding: 0.35rem;
  width: fit-content;
  max-width: 100%;
  background: rgba(248, 249, 250, 0.96);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow-soft);
}
.workspace-tab {
  min-width: 92px;
  max-width: 280px;
  height: 30px;
  border: 1px solid transparent;
  border-radius: var(--radius-sm);
  background: transparent;
  color: var(--ink-soft);
  padding: 0 0.62rem;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  cursor: pointer;
  font: inherit;
  font-size: 0.74rem;
  font-weight: 780;
}
.workspace-tab:hover {
  background: #eef1f3;
  color: var(--ink);
}
.workspace-tab.is-active {
  background: var(--rail-bg);
  color: #ffffff;
  border-color: var(--rail-bg);
}
.workspace-tab[hidden] { display: none; }
.tab-panel[hidden] { display: none; }
.shared-legend {
  position: sticky;
  top: -0.9rem;
  z-index: 5;
  display: block;
  margin-bottom: 0.75rem;
  padding: 0.55rem 0.7rem;
  background: rgba(248, 249, 250, 0.96);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow);
  backdrop-filter: blur(8px);
}
.shared-legend[hidden] { display: none; }
.shared-legend-header {
  display: flex;
  align-items: center;
  gap: 0.55rem;
  min-width: 0;
}
.shared-legend-title {
  color: var(--ink-muted);
  font-size: 0.68rem;
  font-weight: 820;
  letter-spacing: 0.045em;
  text-transform: uppercase;
}
.shared-legend-count {
  color: var(--ink-muted);
  font-size: 0.68rem;
  font-weight: 740;
  white-space: nowrap;
}
.shared-legend-spacer { flex: 1; }
.shared-legend-group {
  display: inline-flex;
  align-items: center;
  gap: 0.32rem;
  color: var(--ink-muted);
  font-size: 0.68rem;
  font-weight: 760;
  white-space: nowrap;
}
.shared-legend-group select {
  height: 26px;
  min-width: 138px;
  max-width: 220px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: #ffffff;
  color: var(--ink-soft);
  padding: 0 0.38rem;
  font: inherit;
  font-size: 0.7rem;
  font-weight: 720;
}
.shared-legend-toggle {
  height: 26px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: #ffffff;
  color: var(--ink-soft);
  padding: 0 0.58rem;
  cursor: pointer;
  font: inherit;
  font-size: 0.7rem;
  font-weight: 780;
}
.shared-legend-toggle:hover {
  background: #eef1f3;
  color: var(--ink);
}
.shared-legend-items {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 0.45rem 0.8rem;
  margin-top: 0.55rem;
  max-height: 180px;
  overflow: auto;
}
.legend-item {
  display: inline-flex;
  align-items: center;
  gap: 0.32rem;
  color: var(--ink-soft);
  font-size: 0.74rem;
  font-weight: 680;
  padding: 0.1rem 0.18rem;
  border-radius: var(--radius-sm);
}
.legend-item:hover {
  background: rgba(20, 24, 28, 0.06);
}
.legend-item.is-hidden {
  opacity: 0.42;
  text-decoration: line-through;
}
.legend-visible-toggle {
  width: 16px;
  height: 16px;
  margin: 0;
  flex: 0 0 16px;
  cursor: pointer;
}
.legend-swatch {
  width: 24px;
  height: 24px;
  padding: 0;
  border-radius: 5px;
  border: 1px solid rgba(20, 24, 28, 0.18);
  flex: 0 0 24px;
  cursor: pointer;
  overflow: hidden;
}
.legend-color-code {
  width: 82px;
  height: 24px;
  padding: 0 0.35rem;
  font-size: 0.68rem;
  font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;
}
.pagination-bar {
  position: sticky;
  bottom: 0;
  z-index: 2;
  margin-top: 0.8rem;
  padding: 0.6rem 0.7rem;
  display: flex;
  align-items: center;
  justify-content: center;
  gap: 0.55rem;
  background: rgba(238, 240, 242, 0.92);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow-soft);
  backdrop-filter: blur(8px);
}
.pagination-bar[hidden] { display: none; }
.pagination-label {
  min-width: 150px;
  text-align: center;
  color: var(--ink-soft);
  font-size: 0.76rem;
  font-weight: 760;
}
.pagination-count {
  color: var(--ink-muted);
  font-size: 0.72rem;
}
.plot-card {
  background: var(--panel-strong);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  display: flex;
  flex-direction: column;
  min-width: 0;
  overflow: hidden;
  box-shadow: var(--shadow-soft);
  transition: border-color 120ms ease, box-shadow 120ms ease, transform 120ms ease;
}
.plot-card:hover {
  border-color: var(--line-strong);
  box-shadow: var(--shadow);
}
.plot-card.is-selected {
  border-color: var(--accent);
  box-shadow: 0 0 0 2px rgba(31, 111, 120, 0.16), var(--shadow-soft);
}
.plot-card-header {
  min-height: 36px;
  padding: 0.42rem 0.58rem;
  background: #f1f3f5;
  border-bottom: 1px solid var(--line);
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 0.5rem;
}
.plot-card-title {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  font-size: 0.75rem;
  font-weight: 760;
}
.plot-card-meta {
  flex-shrink: 0;
  color: var(--ink-muted);
  font-size: 0.64rem;
  padding: 0.12rem 0.36rem;
  border-radius: 999px;
  background: #e7eaed;
}
.plot-card-actions {
  display: inline-flex;
  align-items: center;
  gap: 0.35rem;
  flex-shrink: 0;
}
.plot-card-action {
  height: 24px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: var(--field);
  color: var(--ink-soft);
  padding: 0 0.48rem;
  cursor: pointer;
  font: inherit;
  font-size: 0.68rem;
  font-weight: 780;
}
.plot-card-drilldown {
  height: 24px;
  max-width: 122px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: var(--field);
  color: var(--ink-soft);
  padding: 0 0.28rem;
  cursor: pointer;
  font: inherit;
  font-size: 0.68rem;
  font-weight: 780;
}
.plot-card-action:hover,
.plot-card-drilldown:hover {
  border-color: var(--line-strong);
  color: var(--ink);
  background: #eef1f3;
}
.plot {
  width: 100%;
  height: min(66vh, 680px);
  min-height: 500px;
}
.plot.mini {
  height: 300px;
  min-height: 300px;
}

.drilldown-panel {
  background: rgba(248, 249, 250, 0.96);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow-soft);
  overflow: hidden;
}
.composition-panel {
  display: grid;
  gap: 0.78rem;
}
.composition-toolbar {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(130px, 1fr));
  gap: 0.58rem;
  align-items: end;
  padding: 0.7rem;
  background: rgba(248, 249, 250, 0.96);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow-soft);
}
.composition-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 0.78rem;
}
.composition-grid.is-volcano {
  grid-template-columns: minmax(0, 1fr);
}
.composition-grid.is-volcano .composition-plot {
  height: min(76vh, 820px);
  min-height: 600px;
}
.composition-card {
  background: var(--panel-strong);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  overflow: hidden;
  box-shadow: var(--shadow-soft);
  min-width: 0;
}
.composition-card-header {
  min-height: 38px;
  padding: 0.48rem 0.65rem;
  background: #f1f3f5;
  border-bottom: 1px solid var(--line);
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 0.5rem;
}
.composition-card-title {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  font-size: 0.78rem;
  font-weight: 820;
}
.composition-card-meta {
  flex-shrink: 0;
  color: var(--ink-muted);
  font-size: 0.66rem;
}
.composition-plot {
  width: 100%;
  height: min(72vh, 760px);
  min-height: 560px;
}
.composition-empty {
  padding: 1rem;
  color: var(--ink-muted);
  font-size: 0.8rem;
}
.drilldown-header {
  min-height: 42px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 0.7rem;
  padding: 0.52rem 0.7rem;
  border-bottom: 1px solid var(--line);
  background: #f1f3f5;
}
.drilldown-title {
  min-width: 0;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
  font-size: 0.78rem;
  font-weight: 820;
}
.drilldown-grid {
  display: grid;
  grid-template-columns: repeat(var(--grid-columns, 4), minmax(0, 1fr));
  gap: 0.78rem;
  padding: 0.78rem;
}

.selection-panel {
  width: 330px;
  flex: 0 0 330px;
  min-width: 0;
  max-height: 100%;
  overflow: hidden;
  background: rgba(248, 249, 250, 0.96);
  border: 1px solid var(--line);
  border-radius: var(--radius-md);
  box-shadow: var(--shadow-soft);
  display: flex;
  flex-direction: column;
  transition: flex-basis 150ms ease, width 150ms ease;
}
.selection-panel.is-collapsed {
  width: 42px;
  flex-basis: 42px;
}
.inspector-header {
  min-height: 38px;
  padding: 0.45rem 0.58rem;
  display: flex;
  align-items: center;
  gap: 0.5rem;
  border-bottom: 1px solid var(--line);
  background: #f1f3f5;
}
.inspector-heading {
  min-width: 0;
  flex: 1;
  color: var(--ink-soft);
  font-size: 0.72rem;
  font-weight: 820;
  letter-spacing: 0.045em;
  text-transform: uppercase;
  white-space: nowrap;
}
.inspector-toggle {
  width: 28px;
  height: 26px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: #ffffff;
  color: var(--ink-soft);
  cursor: pointer;
  font: inherit;
  font-size: 0.76rem;
  font-weight: 820;
}
.selection-panel.is-collapsed .inspector-heading,
.selection-panel.is-collapsed .inspector-body {
  display: none;
}
.inspector-body {
  overflow-y: auto;
  padding: 0.65rem 0.8rem;
}
.inspector-config {
  display: grid;
  gap: 0.48rem;
  margin-bottom: 0.65rem;
  padding-bottom: 0.65rem;
  border-bottom: 1px solid var(--line);
}
.inspector-config:empty { display: none; }
.inspector-config select {
  width: 100%;
  min-width: 0;
}
.selection-panel .inspector-card {
  border-bottom: none;
  padding-bottom: 0;
  margin-bottom: 0;
}
.selection-panel .inspector-card + .inspector-card {
  margin-top: 0.65rem;
}
.selection-panel .feature-table {
  width: 100%;
  table-layout: fixed;
  border-collapse: collapse;
}
.selection-panel .feature-table tr,
.selection-panel .feature-table th,
.selection-panel .feature-table td {
  border-bottom: 1px solid rgba(180, 185, 190, 0.45);
  padding: 0.32rem 0;
  vertical-align: top;
}
.selection-panel .feature-table th {
  width: 42%;
  padding-right: 0.55rem;
  font-size: 0.66rem;
}
.selection-panel .feature-table td {
  font-size: 0.76rem;
  overflow-wrap: anywhere;
}
.inspector-empty {
  min-height: 44px;
  color: var(--ink-muted);
}
.inspector-card {
  border-bottom: 1px solid var(--line);
  padding-bottom: 1rem;
  margin-bottom: 1rem;
}
.inspector-card:last-child { border-bottom: none; }
.inspector-title {
  font-size: 0.9rem;
  font-weight: 800;
  line-height: 1.25;
  overflow-wrap: anywhere;
}
.inspector-subtitle {
  margin-top: 0.3rem;
  color: var(--ink-muted);
  font-size: 0.72rem;
}
.structure-card {
  display: grid;
  gap: 0.45rem;
}
.structure-card img {
  width: 100%;
  max-height: 260px;
  object-fit: contain;
  background: #ffffff;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
}
.structure-smiles {
  color: var(--ink-muted);
  font-size: 0.68rem;
  overflow-wrap: anywhere;
}
.feature-table {
  width: 100%;
  table-layout: fixed;
  border-collapse: collapse;
  font-size: 0.76rem;
}
.feature-table th, .feature-table td {
  padding: 0.35rem 0;
  border-bottom: 1px solid var(--line);
  text-align: left;
  vertical-align: top;
  overflow-wrap: anywhere;
}
.feature-table th {
  width: 42%;
  color: var(--ink-muted);
  padding-right: 0.5rem;
}

.btn {
  height: 30px;
  border: 1px solid var(--line);
  border-radius: var(--radius-sm);
  background: var(--field);
  color: var(--ink-soft);
  padding: 0 0.7rem;
  cursor: pointer;
  font: inherit;
  font-size: 0.74rem;
  font-weight: 760;
}
.btn:hover { border-color: var(--line-strong); color: var(--ink); background: #eef1f3; }
.btn:disabled {
  cursor: not-allowed;
  opacity: 0.48;
  background: #eef1f3;
}
.muted { color: var(--ink-muted); font-size: 0.76rem; }
.hidden-dependency { display: none; }

@media (max-width: 1100px) {
  .app-rail { display: none; }
  .app-sidebar { position: fixed; inset: 0 auto 0 0; z-index: 30; box-shadow: var(--shadow); transform: translateX(-100%); }
  .app-sidebar.is-open { transform: translateX(0); }
  .header-search { width: 44vw; }
  .plot-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .composition-toolbar { grid-template-columns: repeat(2, minmax(0, 1fr)); }
  .composition-grid { grid-template-columns: 1fr; }
}
@media (max-width: 720px) {
  .plot-grid { grid-template-columns: 1fr; }
  .composition-toolbar { grid-template-columns: 1fr; }
  .header-search { width: 100%; }
  .app-header-bar { height: auto; min-height: 52px; flex-wrap: wrap; padding: 0.6rem 0.8rem; }
}
)"

data_explorer_js <- r"(
document.addEventListener('DOMContentLoaded', function () {
  const state = window.MAPP_DATA_EXPLORER;
  if (!state) {
    document.body.innerHTML = '<main class="workspace"><section class="plot-card"><div class="plot-card-header"><div class="plot-card-title">MAPP data explorer</div></div><p class="muted" style="padding:1rem">Payload missing: data_explorer_assets/data_explorer_payload.js</p></section></main>';
    return;
  }
  if (!state.feature_metadata && state.feature_metadata_columns) {
    const featureColumns = state.feature_metadata_columns;
    const featureColumnNames = Object.keys(featureColumns);
    const featureCount = featureColumnNames.length ? featureColumns[featureColumnNames[0]].length : 0;
    state.feature_metadata = Array.from({ length: featureCount }, function(_, index) {
      const row = {};
      featureColumnNames.forEach(function(column) {
        row[column] = featureColumns[column][index];
      });
      return row;
    });
    delete state.feature_metadata_columns;
  }

  const levelSelect = document.querySelector('[data-level]');
  const searchInput = document.querySelector('[data-search]');
  const globalSearchInput = document.querySelector('[data-global-search]');
	  const filterPanel = document.querySelector('[data-filter-panel]');
	  const filterToggle = document.querySelector('[data-filter-toggle]');
	  const itemSelect = document.querySelector('[data-item]');
	  const groupSelect = document.querySelector('[data-group]');
	  const facetSelect = document.querySelector('[data-facet]');
	  const colorSelect = document.querySelector('[data-color]');
	  const sampleFilterRows = Array.from(document.querySelectorAll('[data-sample-filter]'));
	  const groupOrderInput = document.querySelector('[data-group-order]');
	  const groupColorsInput = document.querySelector('[data-group-colors]');
	  const useCurrentGroupsButton = document.querySelector('[data-use-current-groups]');
	  const pathwaySelect = document.querySelector('[data-pathway]');
	  const superclassSelect = document.querySelector('[data-superclass]');
	  const classSelect = document.querySelector('[data-class]');
	  const componentSelect = document.querySelector('[data-component]');
	  const dropSingletonsToggle = document.querySelector('[data-drop-singletons]');
	  const numericFilterPanel = document.querySelector('[data-numeric-filters]');
	  const valueModeSelect = document.querySelector('[data-value-mode]');
	  const plotTypeSelect = document.querySelector('[data-plot-type]');
	  const sortModeSelect = document.querySelector('[data-sort-mode]');
	  const legendModeSelect = document.querySelector('[data-legend-mode]');
	  const gridToggle = document.querySelector('[data-grid-toggle]');
  const pointsToggle = document.querySelector('[data-points-toggle]');
  const gridColumnsInput = document.querySelector('[data-grid-columns]');
	  const gridRowsInput = document.querySelector('[data-grid-rows]');
	  const overviewTab = document.querySelector('[data-tab-overview]');
	  const drilldownTab = document.querySelector('[data-tab-drilldown]');
	  const compositionTab = document.querySelector('[data-tab-composition]');
	  const overviewPanel = document.querySelector('[data-overview-panel]');
	  const compositionPanel = document.querySelector('[data-composition-panel]');
	  const compositionGroupSelect = document.querySelector('[data-composition-group]');
	  const compositionASelect = document.querySelector('[data-composition-a]');
	  const compositionBSelect = document.querySelector('[data-composition-b]');
	  const compositionDepthSelect = document.querySelector('[data-composition-depth]');
	  const compositionValueModeSelect = document.querySelector('[data-composition-value-mode]');
	  const compositionViewSelect = document.querySelector('[data-composition-view]');
	  const compositionVolcanoSourceSelect = document.querySelector('[data-composition-volcano-source]');
	  const compositionPValueInput = document.querySelector('[data-composition-p-value]');
	  const compositionFoldChangeInput = document.querySelector('[data-composition-fold-change]');
	  const compositionTreemapControls = Array.from(document.querySelectorAll('[data-composition-treemap-control]'));
	  const compositionVolcanoControls = Array.from(document.querySelectorAll('[data-composition-volcano-control]'));
	  const compositionGrid = document.querySelector('[data-composition-grid]');
	  const compositionCardA = document.querySelector('[data-composition-card-a]');
	  const compositionCardB = document.querySelector('[data-composition-card-b]');
	  const compositionPlotA = document.querySelector('[data-composition-plot-a]');
	  const compositionPlotB = document.querySelector('[data-composition-plot-b]');
	  const compositionTitleA = document.querySelector('[data-composition-title-a]');
	  const compositionTitleB = document.querySelector('[data-composition-title-b]');
	  const compositionMetaA = document.querySelector('[data-composition-meta-a]');
	  const compositionMetaB = document.querySelector('[data-composition-meta-b]');
	  const compositionResetButton = document.querySelector('[data-composition-reset]');
	  const sharedLegend = document.querySelector('[data-shared-legend]');
	  const plotGrid = document.querySelector('[data-plot-grid]');
	  const paginationBar = document.querySelector('[data-pagination]');
	  const paginationLabel = document.querySelector('[data-page-label]');
	  const paginationCount = document.querySelector('[data-page-count]');
	  const previousPageButton = document.querySelector('[data-page-previous]');
	  const nextPageButton = document.querySelector('[data-page-next]');
	  const drilldownPanel = document.querySelector('[data-drilldown]');
	  const drilldownTitle = document.querySelector('[data-drilldown-title]');
	  const drilldownGrid = document.querySelector('[data-drilldown-grid]');
	  const drilldownCloseButton = document.querySelector('[data-drilldown-close]');
	  const inspectorPanel = document.querySelector('[data-inspector-panel]');
	  const inspectorToggle = document.querySelector('[data-inspector-toggle]');
	  const smilesColumnSelect = document.querySelector('[data-smiles-column]');
	  const info = document.querySelector('[data-info]');
  const countLabel = document.querySelector('[data-count]');
  const titleLabel = document.querySelector('[data-title]');
  const summaryLabel = document.querySelector('[data-summary]');

	  titleLabel.textContent = state.title;
	  summaryLabel.textContent = state.sample_ids.length + ' samples - ' + state.feature_metadata.length + ' features';
	  let currentPage = 1;
	  let compositionFocusId = 'root';
	  let compositionRenderToken = 0;

	  const palette = ['#2563EB','#DC2626','#059669','#7C3AED','#D97706','#0891B2','#BE123C','#4B5563','#EA580C','#0F766E','#9333EA','#65A30D'];
	  const featuresById = {};
	  state.feature_metadata.forEach(function(feature) { featuresById[feature.feature_id] = feature; });
	  let sharedLegendCollapsed = true;
	  let sharedLegendInitialized = false;
	  const hiddenLegendValues = {};
  function valueText(value) { return value === null || value === undefined || value === '' ? 'NA' : String(value); }
  function escapeHtml(value) {
    return valueText(value).replace(/[&<>"']/g, function(char) {
      return ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' })[char];
    });
  }
  function prettyColumn(column) { return column ? column.replace(/^attribute_/, '').replace(/_/g, ' ') : 'None'; }
  function uniqueSorted(values) { return Array.from(new Set(values.map(valueText))).sort(function(a,b){ return a.localeCompare(b); }); }
  function addOption(select, value, label) {
    const option = document.createElement('option');
    option.value = value;
    option.textContent = label || value;
    select.appendChild(option);
  }
	  function selectedOptions(select) { return Array.from(select.selectedOptions).map(function(option){ return option.value; }); }
	  function enableClickToggleMultiSelect(select) {
	    select.addEventListener('mousedown', function(event) {
	      if (!event.target || event.target.tagName !== 'OPTION') return;
	      event.preventDefault();
	      event.target.selected = !event.target.selected;
	      select.focus();
	      select.dispatchEvent(new Event('change', { bubbles: true }));
	    });
	  }
		  function metadataValue(sample, column) { return column ? valueText(sample[column]) : ''; }
	  function featureValue(feature, column) { return valueText(feature[column]); }
	  function smilesColumns() {
	    if (!state.feature_metadata.length) return [];
	    const columns = Object.keys(state.feature_metadata[0]).filter(function(column) {
	      if (!/smiles/i.test(column)) return false;
	      return state.feature_metadata.some(function(feature) {
	        return feature[column] && valueText(feature[column]) !== 'NA';
	      });
	    });
	    return columns.sort(function(a, b) {
	      if (a === 'sirius_smiles') return -1;
	      if (b === 'sirius_smiles') return 1;
	      return prettyColumn(a).localeCompare(prettyColumn(b));
	    });
	  }
	  function populateSmilesColumnSelect() {
	    const columns = smilesColumns();
	    smilesColumnSelect.innerHTML = '';
	    addOption(smilesColumnSelect, '', 'Auto detect');
	    columns.forEach(function(column) {
	      addOption(smilesColumnSelect, column, prettyColumn(column));
	    });
	    if (columns.indexOf('sirius_smiles') !== -1) {
	      smilesColumnSelect.value = 'sirius_smiles';
	    } else if (columns.length) {
	      smilesColumnSelect.value = columns[0];
	    }
	  }
	  function featureSmiles(feature) {
	    if (!feature) return '';
	    if (smilesColumnSelect.value && feature[smilesColumnSelect.value] && valueText(feature[smilesColumnSelect.value]) !== 'NA') {
	      return valueText(feature[smilesColumnSelect.value]);
	    }
	    const priority = ['sirius_smiles', 'gnps_smiles', 'isdb_smiles', 'smiles', 'canonical_smiles', 'structure_smiles', 'met_annot_smiles'];
	    for (let i = 0; i < priority.length; i++) {
	      if (feature[priority[i]] && valueText(feature[priority[i]]) !== 'NA') return valueText(feature[priority[i]]);
	    }
	    const key = Object.keys(feature).find(function(column) {
	      return /smiles/i.test(column) && feature[column] && valueText(feature[column]) !== 'NA';
	    });
	    return key ? valueText(feature[key]) : '';
	  }
	  function numericFeatureValue(feature, column) {
	    const value = Number(feature[column]);
	    return Number.isFinite(value) ? value : null;
	  }
	  function formatSliderValue(value) {
	    const numeric = Number(value);
	    if (!Number.isFinite(numeric)) return 'NA';
	    if (Math.abs(numeric) < 1 && numeric !== 0) return numeric.toFixed(3);
	    if (Math.abs(numeric) >= 1000) return numeric.toExponential(2);
	    return numeric.toFixed(2).replace(/\.?0+$/, '');
	  }
	  function parseListInput(text) {
	    return String(text || '').split(/[\n,]+/).map(function(value) { return value.trim(); }).filter(Boolean);
	  }
	  function parseColorMap(text) {
	    const map = {};
	    String(text || '').split(/\n+/).forEach(function(line) {
	      const clean = line.trim();
	      if (!clean) return;
	      const parts = clean.split(/\s*[:=]\s*/);
	      if (parts.length < 2) return;
	      const key = parts.shift().trim();
	      const color = parts.join(':').trim();
	      if (key && color) map[key] = color;
	    });
	    return map;
	  }
	  function orderedValues(values, order) {
	    if (!order.length) return values;
	    const present = new Set(values);
	    const used = new Set();
	    const ordered = order.filter(function(value) {
	      const keep = present.has(value) && !used.has(value);
	      if (keep) used.add(value);
	      return keep;
	    });
	    return ordered.concat(values.filter(function(value) { return !used.has(value); }));
	  }
	  function groupOrder() {
	    return parseListInput(groupOrderInput.value);
	  }
	  function colorFor(value, index) {
	    const customColors = parseColorMap(groupColorsInput.value);
	    return customColors[value] || palette[index % palette.length];
	  }
	  function hexColorOrFallback(color, fallback) {
	    const clean = String(color || '').trim();
	    return /^#[0-9a-fA-F]{6}$/.test(clean) ? clean : fallback;
	  }
	  function setCustomColor(group, color) {
	    const cleanColor = String(color || '').trim();
	    if (!/^#[0-9a-fA-F]{6}$/.test(cleanColor)) return;
	    const colorMap = parseColorMap(groupColorsInput.value);
	    colorMap[group] = cleanColor;
	    const order = currentColorValues();
	    const orderedKeys = order.concat(Object.keys(colorMap).filter(function(key) { return order.indexOf(key) === -1; }));
	    groupColorsInput.value = orderedKeys.map(function(key) { return key + ' = ' + colorMap[key]; }).join('\n');
	  }
	  function metadataColumnPriority(column) {
	    const preferred = [
	      'taxon_16s', 'taxon_full_genome', 'attribute_source_taxon', 'source_taxon',
	      'attribute_treatment', 'treatment', 'sample_type', 'attribute_bioactivity',
	      'bioactivity', 'isolation_method'
	    ];
	    const lower = String(column || '').toLowerCase();
	    const preferredIndex = preferred.indexOf(column);
	    if (preferredIndex !== -1) return preferredIndex;
	    if (/taxon|taxonomy|organism|species|genus/.test(lower)) return 20;
	    if (/treatment|condition|group|class|type|source/.test(lower)) return 40;
	    if (/file|filename|sample[_-]?(id|name)?|mzml|usi/.test(lower)) return 1000;
	    return 100 + state.metadata_columns.indexOf(column);
	  }
	  function orderedMetadataColumns() {
	    return state.metadata_columns.slice().sort(function(a, b) {
	      const diff = metadataColumnPriority(a) - metadataColumnPriority(b);
	      return diff || prettyColumn(a).localeCompare(prettyColumn(b));
	    });
	  }
	  function metadataDistinctValues(column) {
	    if (!column) return [];
	    return uniqueSorted(state.sample_metadata.map(function(sample) {
	      return metadataValue(sample, column);
	    }).filter(function(value) { return value !== 'NA'; }));
	  }
	  function isIdentifierColumn(column) {
	    return /file|filename|sample[_-]?(id|name)?|mzml|usi/i.test(column || '');
	  }
	  function defaultGroupingColumn() {
	    if (state.default_group && state.metadata_columns.indexOf(state.default_group) !== -1 && metadataDistinctValues(state.default_group).length > 1) {
	      return state.default_group;
	    }
	    const preferred = [
	      'taxon_16s', 'taxon_full_genome', 'attribute_source_taxon', 'source_taxon',
	      'attribute_treatment', 'treatment', 'sample_type', 'attribute_bioactivity',
	      'bioactivity', 'isolation_method'
	    ];
	    for (let i = 0; i < preferred.length; i++) {
	      if (state.metadata_columns.indexOf(preferred[i]) === -1) continue;
	      if (metadataDistinctValues(preferred[i]).length > 1) return preferred[i];
	    }
	    const fallback = orderedMetadataColumns().find(function(column) {
	      const count = metadataDistinctValues(column).length;
	      return count > 1 && count <= 80 && !isIdentifierColumn(column);
	    });
	    if (fallback) return fallback;
	    const anyNonId = orderedMetadataColumns().find(function(column) {
	      return metadataDistinctValues(column).length > 1 && !isIdentifierColumn(column);
	    });
	    if (anyNonId) return anyNonId;
	    return state.metadata_columns[0] || '';
	  }

	  function populateMetadataSelect(select, includeNone) {
    select.innerHTML = '';
    if (includeNone) addOption(select, '', 'None');
    orderedMetadataColumns().forEach(function(column) { addOption(select, column, prettyColumn(column)); });
  }
	  function populateTermSelect(select, column) {
	    select.innerHTML = '';
	    uniqueSorted(state.feature_metadata.map(function(feature) { return featureValue(feature, column); }).filter(function(value) { return value !== 'NA'; }))
	      .forEach(function(value) { addOption(select, value, value); });
	  }
	  function npcColor(level, value, fallback) {
	    const maps = state.npc_color_maps || {};
	    const levelMap = maps[level] || {};
	    return levelMap[value] || fallback || '#D6D6D6';
	  }
	  function referenceNpcColor(pathway, superclass, fallback) {
	    const maps = state.npc_color_maps || {};
	    const pathwayMap = maps.pathway || {};
	    const superclassToPathway = maps.reference_superclass_pathway || {};
	    const referencePathway = superclassToPathway[superclass] || (pathwayMap[pathway] ? pathway : null);
	    return pathwayMap[referencePathway] || fallback || pathwayMap.Other || '#B7B7B7';
	  }
	  function numericFilterColumns() {
	    if (!state.feature_metadata.length) return [];
	    const first = state.feature_metadata[0];
	    return Object.keys(first).filter(function(column) {
	      if (!/(probability|score|confidence|p_value|q_value|pvalue|qvalue|fdr)/i.test(column)) return false;
	      const values = state.feature_metadata.map(function(feature) { return numericFeatureValue(feature, column); }).filter(function(value) { return value !== null; });
	      return values.length > 0;
	    });
	  }
	  function populateNumericFilters() {
	    numericFilterPanel.innerHTML = '';
	    numericFilterColumns().forEach(function(column) {
	      const values = state.feature_metadata.map(function(feature) { return numericFeatureValue(feature, column); }).filter(function(value) { return value !== null; });
	      const min = Math.min.apply(null, values);
	      const max = Math.max.apply(null, values);
	      if (!Number.isFinite(min) || !Number.isFinite(max) || min === max) return;
	      const row = document.createElement('label');
	      row.className = 'numeric-filter-row';
	      row.dataset.numericFilter = column;
	      const head = document.createElement('div');
	      head.className = 'numeric-filter-head';
	      const name = document.createElement('span');
	      name.className = 'numeric-filter-name';
	      name.textContent = prettyColumn(column);
	      const valueLabel = document.createElement('span');
	      valueLabel.className = 'numeric-filter-value';
	      valueLabel.dataset.numericValue = '';
	      const slider = document.createElement('input');
	      slider.type = 'range';
	      slider.min = String(min);
	      slider.max = String(max);
	      slider.step = max <= 1 ? '0.01' : String(Math.max((max - min) / 200, 0.001));
	      slider.value = String(min);
	      slider.dataset.numericColumn = column;
	      slider.dataset.numericMin = String(min);
	      valueLabel.textContent = '>= ' + formatSliderValue(slider.value);
	      slider.addEventListener('input', function() {
	        valueLabel.textContent = '>= ' + formatSliderValue(slider.value);
	        resetPage();
	        render();
	      });
	      head.appendChild(name);
	      head.appendChild(valueLabel);
	      row.appendChild(head);
	      row.appendChild(slider);
	      numericFilterPanel.appendChild(row);
	    });
	  }
	  function passesNumericFilters(feature) {
	    return Array.from(numericFilterPanel.querySelectorAll('[data-numeric-column]')).every(function(slider) {
	      const min = Number(slider.dataset.numericMin);
	      const threshold = Number(slider.value);
	      if (!Number.isFinite(threshold) || threshold <= min) return true;
	      const value = numericFeatureValue(feature, slider.dataset.numericColumn);
	      return value !== null && value >= threshold;
	    });
	  }
	  function populateFilterValues(row) {
	    const columnSelect = row.querySelector('[data-filter-column]');
	    const valuesSelect = row.querySelector('[data-filter-values]');
	    const previous = selectedOptions(valuesSelect);
	    valuesSelect.innerHTML = '';
	    const column = columnSelect.value;
	    if (!column) return;
	    uniqueSorted(state.sample_metadata.map(function(sample) { return metadataValue(sample, column); }))
	      .forEach(function(value) { addOption(valuesSelect, value, value); });
	    Array.from(valuesSelect.options).forEach(function(option) {
	      option.selected = previous.indexOf(option.value) !== -1;
	    });
	  }
	  function passesTermFilter(feature) {
	    const pathways = selectedOptions(pathwaySelect);
	    const superclasses = selectedOptions(superclassSelect);
	    const classes = selectedOptions(classSelect);
	    const components = selectedOptions(componentSelect);
	    if (pathways.length && pathways.indexOf(featureValue(feature, 'npc_pathway')) === -1) return false;
	    if (superclasses.length && superclasses.indexOf(featureValue(feature, 'npc_superclass')) === -1) return false;
	    if (classes.length && classes.indexOf(featureValue(feature, 'npc_class')) === -1) return false;
	    if (components.length && components.indexOf(featureValue(feature, 'component_id')) === -1) return false;
	    if (dropSingletonsToggle.checked && featureValue(feature, 'component_id') === '-1') return false;
	    if (!passesNumericFilters(feature)) return false;
	    return true;
	  }
  function matchingFeatures() {
    return state.feature_metadata.filter(passesTermFilter);
  }
  function combinedQuery() {
    return [searchInput.value, globalSearchInput.value].join(' ').trim().toLowerCase();
  }
  function buildEntities() {
    const level = levelSelect.value;
    const query = combinedQuery();
    const features = matchingFeatures();
    if (level === 'feature') {
      return features.map(function(feature) {
        return { id: feature.feature_id, label: feature.feature_label, level: 'feature', featureIds: [feature.feature_id], meta: feature };
      }).filter(function(entity) {
        return !query || JSON.stringify(entity.meta).toLowerCase().indexOf(query) !== -1;
      });
    }
	    const column = level === 'pathway' ? 'npc_pathway' : (level === 'superclass' ? 'npc_superclass' : (level === 'component' ? 'component_id' : 'npc_class'));
	    const byTerm = {};
	    features.forEach(function(feature) {
	      const term = featureValue(feature, column);
	      if (term === 'NA') return;
	      if (!byTerm[term]) byTerm[term] = [];
	      byTerm[term].push(feature.feature_id);
	    });
	    return Object.keys(byTerm).sort(function(a,b){ return a.localeCompare(b); }).map(function(term) {
	      const label = level === 'component' ? 'Component ' + term + ' (' + byTerm[term].length + ' features)' : term + ' (' + byTerm[term].length + ' features)';
	      return { id: term, label: label, level: level, featureIds: byTerm[term], meta: { level: level, term: term, n_features: byTerm[term].length } };
	    }).filter(function(entity) {
	      return !query || entity.label.toLowerCase().indexOf(query) !== -1;
	    });
  }
		  function populateItems() {
	    const current = itemSelect.value;
	    const entities = sortEntities(buildEntities());
    itemSelect.innerHTML = '';
    entities.slice(0, 1000).forEach(function(entity) { addOption(itemSelect, entity.id, entity.label); });
    if (current && Array.from(itemSelect.options).some(function(option) { return option.value === current; })) itemSelect.value = current;
	    countLabel.textContent = entities.length + ' matching ' + levelSelect.value + (entities.length === 1 ? '' : 's');
	  }
	  function resetPage() {
	    currentPage = 1;
	  }
	  function selectedSamplesBase() {
	    const filters = sampleFilterRows.map(function(row) {
	      const column = row.querySelector('[data-filter-column]').value;
	      const values = selectedOptions(row.querySelector('[data-filter-values]'));
	      return { column: column, values: values };
	    }).filter(function(filter) { return filter.column && filter.values.length; });
	    return state.sample_metadata.map(function(sample, index) { return { sample: sample, index: index }; })
	      .filter(function(item) {
	        return filters.every(function(filter) {
	          return filter.values.indexOf(metadataValue(item.sample, filter.column)) !== -1;
	        });
	      });
	  }
	  function selectedSamples() {
	    const hidden = currentHiddenLegendSet();
	    const legendColumn = currentLegendColumn();
	    return selectedSamplesBase().filter(function(item) {
	      return !legendColumn || !hidden.has(metadataValue(item.sample, legendColumn));
	    });
	  }
	  function populateCompositionValues() {
	    const previousA = compositionASelect.value;
	    const previousB = compositionBSelect.value;
	    compositionASelect.innerHTML = '';
	    compositionBSelect.innerHTML = '';
	    const column = compositionGroupSelect.value;
	    const values = orderedValues(uniqueSorted(selectedSamples().map(function(item) {
	      return metadataValue(item.sample, column);
	    })), column === groupSelect.value ? groupOrder() : []);
	    values.forEach(function(value) {
	      addOption(compositionASelect, value, value);
	      addOption(compositionBSelect, value, value);
	    });
	    if (previousA && values.indexOf(previousA) !== -1) {
	      compositionASelect.value = previousA;
	    } else if (values.length) {
	      compositionASelect.value = values[0];
	    }
	    if (previousB && values.indexOf(previousB) !== -1) {
	      compositionBSelect.value = previousB;
	    } else if (values.length > 1) {
	      compositionBSelect.value = values[1];
	    } else if (values.length) {
	      compositionBSelect.value = values[0];
	    }
	  }
	  function currentGroupValues() {
	    return orderedValues(uniqueSorted(selectedSamples().map(function(item) {
	      return metadataValue(item.sample, groupSelect.value);
	    })), groupOrder());
	  }
	  function currentColorValues() {
	    const colorColumn = colorSelect.value === '__group__' ? groupSelect.value : colorSelect.value;
	    const values = uniqueSorted(selectedSamplesBase().map(function(item) {
	      return metadataValue(item.sample, colorColumn);
	    }));
	    return orderedValues(values, colorColumn === groupSelect.value ? groupOrder() : []);
	  }
	  function currentColorIndex() {
	    const index = {};
	    currentColorValues().forEach(function(value, position) {
	      index[value] = position;
	    });
	    return index;
	  }
	  function colorLegendTitle() {
	    return colorSelect.value === '__group__' ? prettyColumn(groupSelect.value) : prettyColumn(colorSelect.value);
	  }
	  function currentLegendColumn() {
	    return colorSelect.value === '__group__' ? groupSelect.value : colorSelect.value;
	  }
	  function hiddenLegendSet(column) {
	    if (!hiddenLegendValues[column]) hiddenLegendValues[column] = new Set();
	    return hiddenLegendValues[column];
	  }
	  function currentHiddenLegendSet() {
	    return hiddenLegendSet(currentLegendColumn());
	  }
	  function toggleLegendValue(value) {
	    const hidden = currentHiddenLegendSet();
	    if (hidden.has(value)) {
	      hidden.delete(value);
	    } else {
	      hidden.add(value);
	    }
	    resetPage();
	    render();
	  }
	  function clearLegendHiddenValues() {
	    currentHiddenLegendSet().clear();
	    resetPage();
	    render();
	  }
	  function syncGroupingColumn(column) {
	    const previousGroup = groupSelect.value;
	    const previousCompositionGroup = compositionGroupSelect.value;
	    groupSelect.value = column;
	    if (!previousCompositionGroup || previousCompositionGroup === previousGroup) {
	      compositionGroupSelect.value = column;
	      populateCompositionValues();
	    }
	    resetPage();
	    render();
	  }
	  function renderSharedLegend() {
	    const useSharedLegend = legendModeSelect.value === 'shared';
	    const compositionActive = compositionPanel && !compositionPanel.hidden;
	    sharedLegend.hidden = !useSharedLegend || compositionActive;
	    if (!useSharedLegend || compositionActive) {
	      sharedLegend.innerHTML = '';
	      return;
	    }
	    const values = currentColorValues();
	    const hidden = currentHiddenLegendSet();
	    const hiddenCount = values.filter(function(value) { return hidden.has(value); }).length;
	    if (!sharedLegendInitialized) {
	      sharedLegendCollapsed = values.length > 8;
	      sharedLegendInitialized = true;
	    }
	    sharedLegend.classList.toggle('is-collapsed', sharedLegendCollapsed);
	    const header = document.createElement('div');
	    header.className = 'shared-legend-header';
	    const title = document.createElement('span');
	    title.className = 'shared-legend-title';
	    title.textContent = colorLegendTitle();
	    const count = document.createElement('span');
	    count.className = 'shared-legend-count';
	    count.textContent = values.length + (values.length === 1 ? ' group' : ' groups') + (hiddenCount ? ' - ' + hiddenCount + ' hidden' : '');
	    const spacer = document.createElement('span');
	    spacer.className = 'shared-legend-spacer';
	    const groupingControl = document.createElement('label');
	    groupingControl.className = 'shared-legend-group';
	    const groupingLabel = document.createElement('span');
	    groupingLabel.textContent = 'Group by';
	    const groupingSelect = document.createElement('select');
	    orderedMetadataColumns().forEach(function(column) {
	      addOption(groupingSelect, column, prettyColumn(column));
	    });
	    groupingSelect.value = groupSelect.value;
	    groupingSelect.addEventListener('change', function() {
	      syncGroupingColumn(groupingSelect.value);
	    });
	    groupingControl.appendChild(groupingLabel);
	    groupingControl.appendChild(groupingSelect);
	    const toggle = document.createElement('button');
	    toggle.type = 'button';
	    toggle.className = 'shared-legend-toggle';
	    toggle.textContent = sharedLegendCollapsed ? 'Show legend' : 'Hide legend';
	    toggle.addEventListener('click', function() {
	      sharedLegendCollapsed = !sharedLegendCollapsed;
	      renderSharedLegend();
	    });
	    sharedLegend.innerHTML = '';
	    header.appendChild(title);
	    header.appendChild(count);
	    header.appendChild(spacer);
	    header.appendChild(groupingControl);
	    if (hiddenCount) {
	      const clearHidden = document.createElement('button');
	      clearHidden.type = 'button';
	      clearHidden.className = 'shared-legend-toggle';
	      clearHidden.textContent = 'Show all';
	      clearHidden.title = 'Restore hidden legend groups';
	      clearHidden.addEventListener('click', clearLegendHiddenValues);
	      header.appendChild(clearHidden);
	    }
	    header.appendChild(toggle);
	    sharedLegend.appendChild(header);
	    if (sharedLegendCollapsed) return;
	    const items = document.createElement('div');
	    items.className = 'shared-legend-items';
	    values.forEach(function(value, index) {
	      const item = document.createElement('span');
	      item.className = 'legend-item';
	      item.classList.toggle('is-hidden', hidden.has(value));
	      item.title = hidden.has(value) ? 'Hidden from plots' : 'Visible in plots';
	      const visibleToggle = document.createElement('input');
	      visibleToggle.className = 'legend-visible-toggle';
	      visibleToggle.type = 'checkbox';
	      visibleToggle.checked = !hidden.has(value);
	      visibleToggle.title = hidden.has(value) ? 'Show ' + value : 'Hide ' + value;
	      visibleToggle.addEventListener('change', function() {
	        toggleLegendValue(value);
	      });
	      const fallbackColor = palette[index % palette.length];
	      const currentColor = hexColorOrFallback(colorFor(value, index), fallbackColor);
	      const swatch = document.createElement('input');
	      swatch.className = 'legend-swatch';
	      swatch.type = 'color';
	      swatch.value = currentColor;
	      swatch.title = 'Pick color for ' + value;
	      swatch.addEventListener('click', function(event) { event.stopPropagation(); });
	      const label = document.createElement('span');
	      label.textContent = value;
	      const code = document.createElement('input');
	      code.className = 'legend-color-code';
	      code.type = 'text';
	      code.value = currentColor;
	      code.title = 'Color code for ' + value;
	      code.addEventListener('click', function(event) { event.stopPropagation(); });
	      swatch.addEventListener('input', function() {
	        setCustomColor(value, swatch.value);
	        code.value = swatch.value;
	        render();
	      });
	      code.addEventListener('change', function() {
	        const nextColor = hexColorOrFallback(code.value, swatch.value);
	        code.value = nextColor;
	        swatch.value = nextColor;
	        setCustomColor(value, nextColor);
	        render();
	      });
	      item.appendChild(visibleToggle);
	      item.appendChild(swatch);
	      item.appendChild(label);
	      item.appendChild(code);
	      items.appendChild(item);
	    });
	    sharedLegend.appendChild(items);
	  }
	  function seedCurrentGroups() {
	    const groups = currentGroupValues();
	    groupOrderInput.value = groups.join('\n');
	    const existingColors = parseColorMap(groupColorsInput.value);
	    groupColorsInput.value = groups.map(function(group, index) {
	      return group + ' = ' + (existingColors[group] || palette[index % palette.length]);
	    }).join('\n');
	  }
	  const intensityChunksById = {};
	  (state.intensity_chunks || []).forEach(function(chunk) { intensityChunksById[chunk.id] = chunk; });
	  const intensityCache = {};
	  const chunkPromises = {};
	  function loadScript(src) {
	    return new Promise(function(resolve, reject) {
	      const script = document.createElement('script');
	      script.src = src;
	      script.async = true;
	      script.onload = resolve;
	      script.onerror = function() { reject(new Error('Could not load ' + src)); };
	      document.head.appendChild(script);
	    });
	  }
	  function loadChunk(chunkId) {
	    if (!chunkId) return Promise.resolve();
	    if (chunkPromises[chunkId]) return chunkPromises[chunkId];
	    const chunk = intensityChunksById[chunkId];
	    if (!chunk) return Promise.reject(new Error('Unknown intensity chunk ' + chunkId));
	    chunkPromises[chunkId] = loadScript(chunk.path).then(function() {
	      const chunkData = (window.MAPP_DATA_EXPLORER_INTENSITY_CHUNKS || {})[chunkId] || {};
	      Object.keys(chunkData).forEach(function(featureId) {
	        intensityCache[featureId] = chunkData[featureId];
	      });
	      delete (window.MAPP_DATA_EXPLORER_INTENSITY_CHUNKS || {})[chunkId];
	    });
	    return chunkPromises[chunkId];
	  }
	  function loadFeatureIntensities(featureIds) {
	    const chunkIds = Array.from(new Set(featureIds.map(function(featureId) {
	      return state.feature_chunk_map ? state.feature_chunk_map[featureId] : null;
	    }).filter(Boolean)));
	    return Promise.all(chunkIds.map(loadChunk));
	  }
  async function rawEntityValues(entity) {
    await loadFeatureIntensities(entity.featureIds);
    const values = state.sample_ids.map(function() { return 0; });
    entity.featureIds.forEach(function(featureId) {
      const featureValues = intensityCache[featureId] || [];
      featureValues.forEach(function(value, index) {
        const numeric = Number(value);
        if (Number.isFinite(numeric)) values[index] += numeric;
      });
    });
    return values;
  }
  const totals = state.sample_totals || state.sample_ids.map(function() { return 0; });
	  async function transformedValues(entity) {
    const rawValues = await rawEntityValues(entity);
    return rawValues.map(function(value, index) {
      if (valueModeSelect.value === 'log10') return Math.log10(value + 1);
      if (valueModeSelect.value === 'percent_total') return totals[index] > 0 ? 100 * value / totals[index] : null;
      return value;
    });
  }
	  function yTitle() {
    if (valueModeSelect.value === 'log10') return 'log10 intensity + 1';
    if (valueModeSelect.value === 'percent_total') return 'Percent of sample total intensity';
    return 'Raw intensity';
  }
	  function compositionValueTitle() {
	    if (compositionValueModeSelect.value === 'percent_selected') return 'Percent of displayed composition';
	    if (compositionValueModeSelect.value === 'percent_sample_total') return 'Percent of total sample intensity';
	    return 'Summed raw intensity';
	  }
	  function formatValue(value) {
	    if (value === null || value === undefined || !Number.isFinite(Number(value))) return 'NA';
	    if (valueModeSelect.value === 'raw') return Number(value).toExponential(3);
	    return Number(value).toPrecision(4);
	  }
	  function entitySummaryScore(entity, mode) {
	    const featureIds = entity.featureIds || [];
	    const values = featureIds.map(function(featureId) {
	      const feature = featuresById[featureId] || {};
	      const column = mode === 'mean_desc' ? 'feature_mean_intensity' : 'feature_max_intensity';
	      const value = Number(feature[column]);
	      return Number.isFinite(value) ? value : 0;
	    });
	    if (!values.length) return -Infinity;
	    if (mode === 'mean_desc') return values.reduce(function(sum, value) { return sum + value; }, 0) / values.length;
	    return Math.max.apply(null, values);
	  }
	  function sortEntities(entities) {
	    const mode = sortModeSelect.value;
	    if (mode === 'label') {
	      return entities.slice();
	    }
	    return entities.slice().sort(function(a, b) {
	      const aScore = entitySummaryScore(a, mode);
	      const bScore = entitySummaryScore(b, mode);
	      if (bScore !== aScore) return bScore - aScore;
	      return a.label.localeCompare(b.label);
	    });
	  }
	  function wrapPlotLabel(value, width, maxLines) {
	    const text = valueText(value);
	    if (text.length <= width) return text;
	    const words = text.split(/[\s_/:-]+/).filter(Boolean);
	    const lines = [];
	    let line = '';
	    words.forEach(function(word) {
	      const next = line ? line + ' ' + word : word;
	      if (next.length > width && line) {
	        lines.push(line);
	        line = word;
	      } else {
	        line = next;
	      }
	    });
	    if (line) lines.push(line);
	    const clipped = lines.slice(0, maxLines);
	    if (lines.length > maxLines) clipped[maxLines - 1] = clipped[maxLines - 1].replace(/\s+$/, '') + '...';
	    return clipped.join('<br>');
	  }
	  function axisSuffix(index) {
	    return index === 0 ? '' : String(index + 1);
	  }
	  async function plotEntity(entity, element) {
	    element.innerHTML = '<div class="composition-empty">Loading intensities...</div>';
	    const values = await transformedValues(entity);
	    element.innerHTML = '';
	    const sampleItems = selectedSamples();
    const rows = sampleItems.map(function(item) {
      const sample = item.sample;
      const colorColumn = colorSelect.value === '__group__' ? groupSelect.value : colorSelect.value;
      return {
        value: values[item.index],
        sample: sample,
        group: metadataValue(sample, groupSelect.value),
        color: metadataValue(sample, colorColumn),
        facet: facetSelect.value ? metadataValue(sample, facetSelect.value) : ''
      };
    }).filter(function(row) { return row.value !== null && row.value !== undefined && Number.isFinite(Number(row.value)); });
	    const order = groupOrder();
	    const facets = facetSelect.value ? uniqueSorted(rows.map(function(row) { return row.facet; })) : [''];
	    const groups = orderedValues(uniqueSorted(rows.map(function(row) { return row.group; })), order);
	    const colorValues = colorSelect.value === '__group__' ? groups : orderedValues(uniqueSorted(rows.map(function(row) { return row.color; })), colorSelect.value === groupSelect.value ? order : []);
    const colorIndex = currentColorIndex();
	    const traces = [];
	    const legendShown = new Set();
	    const annotations = [];
	    const isMiniPlot = element.classList.contains('mini');
	    const facetColumns = facets.length <= 1 ? 1 : Math.min(facets.length, isMiniPlot ? 2 : 3);
	    const facetRows = Math.ceil(facets.length / facetColumns);
	    const facetGapX = facets.length > 1 ? 0.075 : 0;
	    const facetGapY = facets.length > 1 ? 0.15 : 0;
	    const facetWidth = (1 - facetGapX * (facetColumns - 1)) / facetColumns;
	    const facetHeight = (1 - facetGapY * (facetRows - 1)) / facetRows;
	    if (facets.length > 1) {
	      element.style.height = (isMiniPlot ? Math.max(320, facetRows * 250) : Math.max(520, facetRows * 340)) + 'px';
	    } else {
	      element.style.height = '';
	    }
	    facets.forEach(function(facet, facetIndex) {
	      const facetRows = rows.filter(function(row) { return row.facet === facet; });
	      const suffix = axisSuffix(facetIndex);
	      colorValues.forEach(function(colorValue) {
	        const traceGroups = colorSelect.value === '__group__' ? [colorValue] : groups;
	        traceGroups.forEach(function(group) {
	          const traceRows = facetRows.filter(function(row) { return row.color === colorValue && row.group === group; });
	          if (!traceRows.length) return;
	          const trace = {
	            x: traceRows.map(function(row) { return row.group; }),
	            y: traceRows.map(function(row) { return row.value; }),
            text: traceRows.map(function(row) {
              return 'Sample: ' + metadataValue(row.sample, 'sample_id') +
                '<br>' + prettyColumn(groupSelect.value) + ': ' + row.group +
                '<br>Value: ' + formatValue(row.value);
            }),
            hoverinfo: 'text',
            name: colorValue,
            legendgroup: colorValue,
	            showlegend: legendModeSelect.value === 'per_plot' && !legendShown.has(colorValue),
            marker: { color: colorFor(colorValue, colorIndex[colorValue] !== undefined ? colorIndex[colorValue] : 0), size: 6, opacity: 0.82 },
	            line: { color: colorFor(colorValue, colorIndex[colorValue] !== undefined ? colorIndex[colorValue] : 0) }
	          };
          legendShown.add(colorValue);
          if (plotTypeSelect.value === 'violin') {
            trace.type = 'violin';
            trace.box = { visible: true };
            trace.meanline = { visible: true };
            trace.points = pointsToggle.checked ? 'all' : false;
          } else if (plotTypeSelect.value === 'points') {
            trace.type = 'scatter';
            trace.mode = 'markers';
          } else {
            trace.type = 'box';
            trace.boxpoints = pointsToggle.checked ? 'all' : false;
          }
	          trace.xaxis = 'x' + suffix;
	          trace.yaxis = 'y' + suffix;
	          traces.push(trace);
	        });
	      });
	      if (facet) {
	        const facetRow = Math.floor(facetIndex / facetColumns);
	        const facetColumn = facetIndex % facetColumns;
	        const x0 = facetColumn * (facetWidth + facetGapX);
	        const x1 = x0 + facetWidth;
	        const y1 = 1 - facetRow * (facetHeight + facetGapY);
	        annotations.push({
	          text: wrapPlotLabel(facet, isMiniPlot ? 16 : 24, 2),
	          xref: 'paper',
	          yref: 'paper',
	          x: (x0 + x1) / 2,
	          y: Math.min(1.08, y1 + 0.055),
	          showarrow: false,
	          font: { size: isMiniPlot ? 9 : 10, color: '#30343a' },
	          align: 'center',
	          bgcolor: 'rgba(255,255,255,0.82)',
	          borderpad: 2
	        });
	      }
	    });
    const yValues = rows.map(function(row) { return Number(row.value); });
    const minY = Math.min.apply(null, yValues);
    const maxY = Math.max.apply(null, yValues);
    const pad = Number.isFinite(minY) && Number.isFinite(maxY) && maxY !== minY ? (maxY - minY) * 0.06 : 1;
    const yRange = Number.isFinite(minY) && Number.isFinite(maxY) ? [Math.max(0, minY - pad), maxY + pad] : undefined;
    const yAxis = { title: yTitle(), gridcolor: '#E5E7EB', zeroline: false, range: yRange, tickformat: valueModeSelect.value === 'raw' ? '.2e' : '.3f' };
	    const layout = {
	      title: { text: '', font: { size: 12 } },
	      font: { family: 'Arial, sans-serif', size: 10, color: '#14161a' },
		      margin: { l: 58, r: legendModeSelect.value === 'per_plot' ? 120 : 16, t: facets.length > 1 ? 58 : 24, b: 58 },
	      paper_bgcolor: 'white',
	      plot_bgcolor: 'white',
	      xaxis: { title: prettyColumn(groupSelect.value), tickangle: -20, zeroline: false, automargin: true },
	      yaxis: yAxis,
	      boxmode: 'group',
	      violinmode: 'group',
      annotations: annotations,
	      showlegend: legendModeSelect.value === 'per_plot',
	      legend: { title: { text: colorLegendTitle() }, x: 1.02, y: 1 }
	    };
	    if (facets.length > 1) {
	      facets.forEach(function(facet, index) {
	        const suffix = axisSuffix(index);
	        const facetRow = Math.floor(index / facetColumns);
	        const facetColumn = index % facetColumns;
	        const x0 = facetColumn * (facetWidth + facetGapX);
	        const x1 = x0 + facetWidth;
	        const y1 = 1 - facetRow * (facetHeight + facetGapY);
	        const y0 = y1 - facetHeight;
	        const isLeftColumn = facetColumn === 0;
	        const isBottomRow = facetRow === facetRows - 1;
	        layout['xaxis' + suffix] = {
	          domain: [x0, x1],
	          title: isBottomRow ? prettyColumn(groupSelect.value) : '',
	          tickangle: -20,
	          zeroline: false,
	          automargin: true,
	          matches: index === 0 ? undefined : 'x'
	        };
	        layout['yaxis' + suffix] = Object.assign({}, yAxis, {
	          domain: [y0, y1],
	          title: index === 0 ? yTitle() : '',
	          showticklabels: isLeftColumn,
	          ticks: isLeftColumn ? 'outside' : '',
	          matches: index === 0 ? undefined : 'y'
	        });
	      });
	    }
	    return Plotly.react(element, traces, layout, { displaylogo: false, responsive: true, scrollZoom: false });
  }
	  function renderInfo(entity) {
	    if (!entity) {
	      info.innerHTML = '<div class="inspector-empty"><div>No selection</div><p>Select an entity or use grid mode.</p></div>';
	      return;
	    }
	    Array.from(document.querySelectorAll('.plot-card')).forEach(function(card) {
	      card.classList.toggle('is-selected', card.dataset.entityId === String(entity.id));
	    });
	    const rows = Object.keys(entity.meta).map(function(key) {
	      return '<tr><th>' + escapeHtml(key) + '</th><td>' + escapeHtml(entity.meta[key]) + '</td></tr>';
	    }).join('');
	    let structureHtml = '';
	    if ((entity.level || entity.meta.level) === 'feature' && entity.featureIds.length === 1) {
	      const feature = featuresById[entity.featureIds[0]] || entity.meta;
	      const smiles = featureSmiles(feature);
	      if (smiles) {
	        const depictionUrl = 'https://api.naturalproducts.net/latest/depict/2D?smiles=' + encodeURIComponent(smiles) + '&width=480&height=360';
	        structureHtml = '<section class="inspector-card structure-card"><div class="inspector-title">Structure</div><img alt="2D chemical structure" loading="lazy" src="' + depictionUrl + '"><div class="structure-smiles">' + escapeHtml(smiles) + '</div></section>';
	      }
	    }
    info.innerHTML = '<section class="inspector-card"><div class="inspector-title">' + escapeHtml(entity.label) + '</div><div class="inspector-subtitle">' + entity.featureIds.length + ' feature(s)</div></section>' + structureHtml + '<section class="inspector-card"><table class="feature-table">' + rows + '</table></section>';
	  }
	  function activateWorkspaceTab(tabName) {
	    const showDrilldown = tabName === 'drilldown' && !drilldownTab.hidden;
	    const showComposition = tabName === 'composition';
	    overviewTab.classList.toggle('is-active', !showDrilldown && !showComposition);
	    drilldownTab.classList.toggle('is-active', showDrilldown);
	    compositionTab.classList.toggle('is-active', showComposition);
	    overviewPanel.hidden = showDrilldown || showComposition;
	    compositionPanel.hidden = !showComposition;
	    drilldownPanel.hidden = !showDrilldown;
	    renderSharedLegend();
	    window.requestAnimationFrame(function() {
	      if (showComposition) {
	        renderComposition();
	        [compositionPlotA, compositionPlotB].forEach(function(plot) { if (plot) Plotly.Plots.resize(plot); });
	      } else {
	        Array.from(document.querySelectorAll(showDrilldown ? '[data-drilldown-grid] .plot' : '[data-overview-panel] .plot')).forEach(resetPlotZoom);
	      }
	    });
	  }
	  function featureEntitiesFrom(entity) {
	    return sortEntities(entity.featureIds.map(function(featureId) {
	      const feature = featuresById[featureId];
	      if (!feature) return null;
	      return {
	        id: feature.feature_id,
	        label: feature.feature_label,
	        level: 'feature',
	        featureIds: [feature.feature_id],
	        meta: feature
	      };
	    }).filter(Boolean));
	  }
	  function drilldownTargetOptions(entity) {
	    if (!entity || entity.featureIds.length <= 1) return [];
	    const level = entity.level || entity.meta.level || 'feature';
	    if (level === 'pathway') return [
	      { value: 'superclass', label: 'NPC superclass' },
	      { value: 'class', label: 'NPC class' },
	      { value: 'feature', label: 'Feature' }
	    ];
	    if (level === 'superclass') return [
	      { value: 'class', label: 'NPC class' },
	      { value: 'feature', label: 'Feature' }
	    ];
	    if (level === 'class') return [
	      { value: 'feature', label: 'Feature' }
	    ];
	    if (level === 'component') return [
	      { value: 'superclass', label: 'NPC superclass' },
	      { value: 'class', label: 'NPC class' },
	      { value: 'feature', label: 'Feature' }
	    ];
	    return [];
	  }
	  function drilldownEntitiesFrom(entity, targetLevel) {
	    if (targetLevel === 'feature') return featureEntitiesFrom(entity);
	    const column = targetLevel === 'superclass' ? 'npc_superclass' : 'npc_class';
	    const byTerm = {};
	    entity.featureIds.forEach(function(featureId) {
	      const feature = featuresById[featureId];
	      if (!feature) return;
	      const term = featureValue(feature, column);
	      if (term === 'NA') return;
	      if (!byTerm[term]) byTerm[term] = [];
	      byTerm[term].push(featureId);
	    });
	    return sortEntities(Object.keys(byTerm).sort(function(a, b) { return a.localeCompare(b); }).map(function(term) {
	      return {
	        id: entity.id + '|' + targetLevel + '|' + term,
	        label: term + ' (' + byTerm[term].length + ' features)',
	        level: targetLevel,
	        featureIds: byTerm[term],
	        meta: { level: targetLevel, term: term, n_features: byTerm[term].length, parent: entity.label }
	      };
	    }));
	  }
	  function clearDrilldown(activateOverview) {
	    if (activateOverview === undefined) activateOverview = true;
	    drilldownTab.hidden = true;
	    drilldownTab.classList.remove('is-active');
	    drilldownTitle.textContent = '';
	    drilldownTab.textContent = 'Exploded features';
	    drilldownGrid.innerHTML = '';
	    if (activateOverview) activateWorkspaceTab('overview');
	  }
	  function samplesForComposition(value) {
	    const column = compositionGroupSelect.value;
	    return selectedSamples().filter(function(item) {
	      return metadataValue(item.sample, column) === value;
	    });
	  }
	  function aggregateCachedFeatureForSamples(featureId, sampleItems) {
	    const values = intensityCache[featureId] || [];
	    let total = 0;
	    sampleItems.forEach(function(item) {
	      const numeric = Number(values[item.index]);
	      if (Number.isFinite(numeric)) total += numeric;
	    });
	    return total;
	  }
	  function sampleTotalForSamples(sampleItems) {
	    return sampleItems.reduce(function(sum, item) {
	      const value = Number(totals[item.index]);
	      return sum + (Number.isFinite(value) ? value : 0);
	    }, 0);
	  }
	  async function buildCompositionTree(groupValue) {
	    const sampleItems = samplesForComposition(groupValue);
	    const depth = compositionDepthSelect.value;
	    const features = matchingFeatures().filter(function(feature) {
	      return featureValue(feature, 'npc_pathway') !== 'NA';
	    });
	    await loadFeatureIntensities(features.map(function(feature) { return feature.feature_id; }));
	    const nodes = {};
	    let focusLabel = groupValue || 'Selected samples';
	    function addNode(id, label, parent, value, level, filterColumn, filterValue, color) {
	      if (!nodes[id]) {
	        nodes[id] = {
	          id: id,
	          label: label,
	          parent: parent,
	          value: 0,
	          level: level,
	          filterColumn: filterColumn || '',
	          filterValue: filterValue || '',
	          color: color || '#D6D6D6',
	          features: new Set()
	        };
	      }
	      nodes[id].value += value;
	      return nodes[id];
	    }
	    const root = addNode('root', groupValue || 'Selected samples', '', 0, 'root', '', '', '#F5F5F5');
	    features.forEach(function(feature) {
	      const raw = aggregateCachedFeatureForSamples(feature.feature_id, sampleItems);
	      if (!Number.isFinite(raw) || raw <= 0) return;
	      const pathway = featureValue(feature, 'npc_pathway') === 'NA' ? 'Unclassified' : featureValue(feature, 'npc_pathway');
	      const superclass = featureValue(feature, 'npc_superclass') === 'NA' ? 'Other' : featureValue(feature, 'npc_superclass');
	      const npcClass = featureValue(feature, 'npc_class') === 'NA' ? 'Other' : featureValue(feature, 'npc_class');
	      const taxonomyColor = referenceNpcColor(pathway, superclass, '#B7B7B7');
	      const pathwayId = 'pathway|' + pathway;
	      const superclassId = pathwayId + '|superclass|' + superclass;
	      const classId = superclassId + '|class|' + npcClass;
	      root.value += raw;
	      [root].forEach(function(node) { node.features.add(feature.feature_id); });
	      addNode(pathwayId, pathway, 'root', raw, 'pathway', 'npc_pathway', pathway, taxonomyColor).features.add(feature.feature_id);
	      if (['superclass', 'class', 'feature'].indexOf(depth) !== -1) {
	        addNode(superclassId, superclass, pathwayId, raw, 'superclass', 'npc_superclass', superclass, taxonomyColor).features.add(feature.feature_id);
	      }
	      if (['class', 'feature'].indexOf(depth) !== -1) {
	        addNode(classId, npcClass, superclassId, raw, 'class', 'npc_class', npcClass, taxonomyColor).features.add(feature.feature_id);
	      }
	      if (depth === 'feature') {
	        const featureId = classId + '|feature|' + feature.feature_id;
	        addNode(featureId, feature.feature_label, classId, raw, 'feature', 'feature_id', feature.feature_id, taxonomyColor).features.add(feature.feature_id);
	      }
	    });
	    if (compositionFocusId !== 'root' && nodes[compositionFocusId]) {
	      focusLabel = nodes[compositionFocusId].label;
	      Object.keys(nodes).forEach(function(id) {
	        if (id !== compositionFocusId && id.indexOf(compositionFocusId + '|') !== 0) {
	          delete nodes[id];
	        }
	      });
	      nodes[compositionFocusId].parent = '';
	      nodes[compositionFocusId].label = focusLabel;
	    } else if (compositionFocusId !== 'root') {
	      focusLabel = compositionFocusId.split('|').pop() || focusLabel;
	      Object.keys(nodes).forEach(function(id) { delete nodes[id]; });
	      addNode(compositionFocusId, focusLabel, '', 0, 'focus', '', '', '#F5F5F5');
	    }
	    const selectedTotal = nodes[compositionFocusId] ? nodes[compositionFocusId].value : root.value;
	    const sampleGrandTotal = sampleTotalForSamples(sampleItems);
	    const divisor = compositionValueModeSelect.value === 'percent_sample_total' ? sampleGrandTotal : selectedTotal;
	    Object.keys(nodes).forEach(function(id) {
	      if (compositionValueModeSelect.value !== 'raw') {
	        nodes[id].value = divisor > 0 ? 100 * nodes[id].value / divisor : 0;
	      }
	      nodes[id].n_features = nodes[id].features.size;
	    });
	    return {
	      nodes: Object.values(nodes).filter(function(node) { return node.id === 'root' || node.value > 0; }),
	      sampleCount: sampleItems.length,
	      featureCount: features.length,
	      focusLabel: focusLabel,
	      selectedTotal: selectedTotal,
	      sampleGrandTotal: sampleGrandTotal
	    };
	  }
	  function treemapHoverText(node) {
	    return '<b>' + node.label + '</b>' +
	      '<br>Level: ' + node.level +
	      '<br>Value: ' + (compositionValueModeSelect.value === 'raw' ? Number(node.value).toExponential(3) : Number(node.value).toFixed(3) + '%') +
	      '<br>Features: ' + node.n_features +
	      '<extra></extra>';
	  }
	  function prepareCompositionPlot(element, message) {
	    if (element && (element.data || element._fullLayout)) {
	      try { Plotly.purge(element); } catch (error) { console.warn('Could not purge previous composition plot', error); }
	    }
	    element.innerHTML = '<div class="composition-empty">' + escapeHtml(message) + '</div>';
	  }
	  function showCompositionError(element, error, renderToken) {
	    if (renderToken !== compositionRenderToken) return;
	    try { Plotly.purge(element); } catch (purgeError) { console.warn('Could not purge failed composition plot', purgeError); }
	    const detail = error && error.message ? error.message : String(error || 'Unknown rendering error');
	    element.innerHTML = '<div class="composition-empty"><b>Could not render this comparison.</b><br>' + escapeHtml(detail) + '</div>';
	    console.error('Composition rendering failed', error);
	  }
	  async function renderCompositionPlot(groupValue, element, titleElement, metaElement, renderToken) {
	    prepareCompositionPlot(element, 'Loading composition...');
	    const tree = await buildCompositionTree(groupValue);
	    if (renderToken !== compositionRenderToken) return;
	    if (!tree.nodes.length || tree.nodes.every(function(node) { return !node.value; })) {
	      element.innerHTML = '<div class="composition-empty">No classified signal for this group with the current filters.</div>';
	      titleElement.textContent = groupValue || 'No group';
	      metaElement.textContent = '0 features';
	      return;
	    }
	    titleElement.textContent = groupValue || 'Selected samples';
	    metaElement.textContent = tree.sampleCount + ' sample(s), ' + tree.featureCount + ' candidate feature(s)' + (compositionFocusId === 'root' ? '' : ', focus: ' + tree.focusLabel);
	    const trace = {
	      type: 'treemap',
	      ids: tree.nodes.map(function(node) { return node.id; }),
	      labels: tree.nodes.map(function(node) { return node.label; }),
	      parents: tree.nodes.map(function(node) { return node.parent; }),
	      values: tree.nodes.map(function(node) { return node.value; }),
	      branchvalues: 'total',
	      maxdepth: compositionDepthSelect.value === 'feature' ? 4 : 3,
	      textinfo: 'label+percent parent',
	      hovertemplate: tree.nodes.map(treemapHoverText),
	      pathbar: { visible: false },
	      root: { color: '#F5F5F5' },
	      marker: {
	        colors: tree.nodes.map(function(node) { return node.color; }),
	        line: { width: 1, color: '#FFFFFF' }
	      },
	      tiling: { packing: 'squarify', pad: 2 }
	    };
	    const layout = {
	      margin: { l: 8, r: 8, t: 8, b: 8 },
	      paper_bgcolor: 'white',
	      font: { family: 'Arial, sans-serif', size: 12, color: '#14161a' },
	      uniformtext: { minsize: 9, mode: 'hide' }
	    };
	    element.innerHTML = '';
	    await Plotly.react(element, [trace], layout, { displaylogo: false, responsive: true, scrollZoom: false });
	    if (renderToken !== compositionRenderToken) return;
	    if (typeof element.removeAllListeners === 'function') element.removeAllListeners('plotly_click');
	    element.on('plotly_click', function(event) {
	      const point = event.points && event.points[0];
	      if (!point) return;
	      const nodeId = point.id || (point.data && point.data.ids ? point.data.ids[point.pointNumber] : null);
	      const node = tree.nodes.find(function(item) { return item.id === nodeId; });
	      if (!node || !node.filterColumn || !node.filterValue) return;
	      const originalEvent = event.event || {};
	      if (!originalEvent.metaKey && !originalEvent.ctrlKey && !originalEvent.altKey) {
	        compositionFocusId = node.id;
	        renderComposition();
	        return;
	      }
	      if (node.filterColumn === 'npc_pathway') {
	        levelSelect.value = 'pathway';
	        pathwaySelect.selectedIndex = -1;
	        Array.from(pathwaySelect.options).forEach(function(option) { option.selected = option.value === node.filterValue; });
	      } else if (node.filterColumn === 'npc_superclass') {
	        levelSelect.value = 'superclass';
	        superclassSelect.selectedIndex = -1;
	        Array.from(superclassSelect.options).forEach(function(option) { option.selected = option.value === node.filterValue; });
	      } else if (node.filterColumn === 'npc_class') {
	        levelSelect.value = 'class';
	        classSelect.selectedIndex = -1;
	        Array.from(classSelect.options).forEach(function(option) { option.selected = option.value === node.filterValue; });
	      } else if (node.filterColumn === 'feature_id') {
	        levelSelect.value = 'feature';
	        searchInput.value = node.filterValue;
	      }
	      resetPage();
	      render();
	      activateWorkspaceTab('overview');
	    });
	  }
	  function meanAndVariance(values) {
	    const clean = values.filter(function(value) { return Number.isFinite(value); });
	    if (!clean.length) return { n: 0, mean: 0, variance: 0 };
	    const mean = clean.reduce(function(sum, value) { return sum + value; }, 0) / clean.length;
	    const variance = clean.length > 1 ? clean.reduce(function(sum, value) {
	      const delta = value - mean;
	      return sum + delta * delta;
	    }, 0) / (clean.length - 1) : 0;
	    return { n: clean.length, mean: mean, variance: variance };
	  }
	  function logGamma(value) {
	    const coefficients = [676.5203681218851,-1259.1392167224028,771.3234287776531,-176.6150291621406,12.507343278686905,-0.13857109526572012,9.984369578019572e-6,1.5056327351493116e-7];
	    if (value < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * value)) - logGamma(1 - value);
	    let x = 0.9999999999998099;
	    const shifted = value - 1;
	    coefficients.forEach(function(coefficient, index) { x += coefficient / (shifted + index + 1); });
	    const t = shifted + coefficients.length - 0.5;
	    return 0.5 * Math.log(2 * Math.PI) + (shifted + 0.5) * Math.log(t) - t + Math.log(x);
	  }
	  function betaContinuedFraction(a, b, x) {
	    const maxIterations = 200;
	    const epsilon = 3e-10;
	    const floor = 1e-30;
	    let qab = a + b;
	    let qap = a + 1;
	    let qam = a - 1;
	    let c = 1;
	    let d = 1 - qab * x / qap;
	    if (Math.abs(d) < floor) d = floor;
	    d = 1 / d;
	    let h = d;
	    for (let m = 1; m <= maxIterations; m += 1) {
	      const m2 = 2 * m;
	      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
	      d = 1 + aa * d;
	      if (Math.abs(d) < floor) d = floor;
	      c = 1 + aa / c;
	      if (Math.abs(c) < floor) c = floor;
	      d = 1 / d;
	      h *= d * c;
	      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
	      d = 1 + aa * d;
	      if (Math.abs(d) < floor) d = floor;
	      c = 1 + aa / c;
	      if (Math.abs(c) < floor) c = floor;
	      d = 1 / d;
	      const delta = d * c;
	      h *= delta;
	      if (Math.abs(delta - 1) < epsilon) break;
	    }
	    return h;
	  }
	  function regularizedBeta(x, a, b) {
	    if (x <= 0) return 0;
	    if (x >= 1) return 1;
	    const factor = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
	    if (x < (a + 1) / (a + b + 2)) return factor * betaContinuedFraction(a, b, x) / a;
	    return 1 - factor * betaContinuedFraction(b, a, 1 - x) / b;
	  }
	  function welchPValue(leftValues, rightValues) {
	    const left = meanAndVariance(leftValues);
	    const right = meanAndVariance(rightValues);
	    if (left.n < 2 || right.n < 2) return 1;
	    const leftTerm = left.variance / left.n;
	    const rightTerm = right.variance / right.n;
	    const standardErrorSquared = leftTerm + rightTerm;
	    if (!(standardErrorSquared > 0)) return left.mean === right.mean ? 1 : 0;
	    const statistic = Math.abs(right.mean - left.mean) / Math.sqrt(standardErrorSquared);
	    const denominator = leftTerm * leftTerm / (left.n - 1) + rightTerm * rightTerm / (right.n - 1);
	    const degreesFreedom = denominator > 0 ? standardErrorSquared * standardErrorSquared / denominator : left.n + right.n - 2;
	    const x = degreesFreedom / (degreesFreedom + statistic * statistic);
	    return Math.max(0, Math.min(1, regularizedBeta(x, degreesFreedom / 2, 0.5)));
	  }
	  async function buildVolcanoRows(leftValue, rightValue) {
	    const leftSamples = samplesForComposition(leftValue);
	    const rightSamples = samplesForComposition(rightValue);
	    const features = matchingFeatures();
	    await loadFeatureIntensities(features.map(function(feature) { return feature.feature_id; }));
	    return {
	      leftCount: leftSamples.length,
	      rightCount: rightSamples.length,
	      rows: features.map(function(feature) {
	        const values = intensityCache[feature.feature_id] || [];
	        const leftRaw = leftSamples.map(function(item) { return Number(values[item.index]); }).filter(Number.isFinite);
	        const rightRaw = rightSamples.map(function(item) { return Number(values[item.index]); }).filter(Number.isFinite);
	        const leftMean = meanAndVariance(leftRaw).mean;
	        const rightMean = meanAndVariance(rightRaw).mean;
	        const leftLog = leftRaw.map(function(value) { return Math.log2(Math.max(0, value) + 1); });
	        const rightLog = rightRaw.map(function(value) { return Math.log2(Math.max(0, value) + 1); });
	        const pValue = welchPValue(leftLog, rightLog);
	        return {
	          feature: feature,
	          featureId: feature.feature_id,
	          log2FoldChange: Math.log2((rightMean + 1) / (leftMean + 1)),
	          minusLog10P: -Math.log10(Math.max(pValue, 1e-300)),
	          pValue: pValue,
	          leftMean: leftMean,
	          rightMean: rightMean
	        };
	      })
	    };
	  }
	  function sameStringSet(left, right) {
	    const a = Array.from(new Set((left || []).map(String))).sort();
	    const b = Array.from(new Set((right || []).map(String))).sort();
	    return a.length === b.length && a.every(function(value, index) { return value === b[index]; });
	  }
	  function selectedComparisonSampleIds(leftValue, rightValue) {
	    const column = compositionGroupSelect.value;
	    return selectedSamples().filter(function(item) {
	      const value = metadataValue(item.sample, column);
	      return value === leftValue || value === rightValue;
	    }).map(function(item) { return state.sample_ids[item.index]; });
	  }
	  function matchingArchivedVolcanoRun(leftValue, rightValue) {
	    const selectedIds = selectedSamples().map(function(item) { return state.sample_ids[item.index]; });
	    const selectedPairIds = selectedComparisonSampleIds(leftValue, rightValue);
	    const requestedGroups = [leftValue, rightValue].map(String).sort();
	    const candidates = (state.archived_volcano_runs || []).filter(function(run) {
	      return run.group_column === compositionGroupSelect.value &&
	        (run.groups || []).map(String).sort().join('\u0000') === requestedGroups.join('\u0000');
	    });
	    return candidates.find(function(run) { return sameStringSet(run.sample_ids, selectedIds); }) ||
	      candidates.find(function(run) { return sameStringSet(run.sample_ids, selectedPairIds); }) || null;
	  }
	  function buildArchivedVolcanoRows(run, leftValue, rightValue) {
	    const allowedFeatures = new Set(matchingFeatures().map(function(feature) { return String(feature.feature_id); }));
	    const rows = [];
	    (run.feature_ids || []).forEach(function(featureId, index) {
	      featureId = String(featureId);
	      const feature = featuresById[featureId];
	      const rawPValue = run.p_values[index];
	      const rawFoldChange = run.log2_fold_changes[index];
	      if (rawPValue === null || rawPValue === undefined || rawFoldChange === null || rawFoldChange === undefined) return;
	      const pValue = Number(rawPValue);
	      const log2FoldChange = Number(rawFoldChange);
	      if (!feature || !allowedFeatures.has(featureId) || !Number.isFinite(pValue) || !Number.isFinite(log2FoldChange)) return;
	      rows.push({
	        feature: feature,
	        featureId: featureId,
	        log2FoldChange: log2FoldChange,
	        minusLog10P: -Math.log10(Math.max(pValue, 1e-300)),
	        pValue: pValue,
	        leftMean: null,
	        rightMean: null
	      });
	    });
	    return {
	      leftCount: samplesForComposition(leftValue).length,
	      rightCount: samplesForComposition(rightValue).length,
	      rows: rows,
	      archivedRun: run
	    };
	  }
	  async function renderVolcanoPlot(leftValue, rightValue, renderToken) {
	    prepareCompositionPlot(compositionPlotA, 'Calculating volcano plot...');
	    const requestedArchived = compositionVolcanoSourceSelect.value === 'archived';
	    const archivedRun = requestedArchived ? matchingArchivedVolcanoRun(leftValue, rightValue) : null;
	    const result = archivedRun ? buildArchivedVolcanoRows(archivedRun, leftValue, rightValue) : await buildVolcanoRows(leftValue, rightValue);
	    if (renderToken !== compositionRenderToken) return;
	    const pCutoff = Math.max(1e-300, Math.min(1, Number(compositionPValueInput.value) || 0.05));
	    const foldCutoff = Math.max(0, Number(compositionFoldChangeInput.value) || 0);
	    const colors = result.rows.map(function(row) {
	      if (row.pValue <= pCutoff && row.log2FoldChange >= foldCutoff) return '#DC2626';
	      if (row.pValue <= pCutoff && row.log2FoldChange <= -foldCutoff) return '#2563EB';
	      return '#9CA3AF';
	    });
	    const significantCount = colors.filter(function(color) { return color !== '#9CA3AF'; }).length;
	    const sourceLabel = archivedRun ? 'Archived statistics: ' + archivedRun.phase + '/' + archivedRun.run_hash : (requestedArchived ? 'Dynamic Welch fallback: no exact archived sample set' : 'Dynamic Welch');
	    compositionTitleA.textContent = archivedRun ? 'Archived contrast: ' + archivedRun.contrast : (rightValue || 'Right') + ' versus ' + (leftValue || 'Left');
	    compositionMetaA.textContent = sourceLabel + ' · ' + result.leftCount + ' left, ' + result.rightCount + ' right sample(s), ' + significantCount + ' highlighted feature(s)';
	    if (!result.rows.length || !result.leftCount || !result.rightCount) {
	      compositionPlotA.innerHTML = '<div class="composition-empty">Both groups need samples and the current filters must retain features.</div>';
	      return;
	    }
	    const trace = {
	      type: 'scattergl',
	      mode: 'markers',
	      x: result.rows.map(function(row) { return row.log2FoldChange; }),
	      y: result.rows.map(function(row) { return row.minusLog10P; }),
	      customdata: result.rows.map(function(row) { return row.featureId; }),
	      text: result.rows.map(function(row) {
	        const meanText = archivedRun ? '' :
	          '<br>Left mean: ' + row.leftMean.toExponential(3) +
	          '<br>Right mean: ' + row.rightMean.toExponential(3);
	        return '<b>' + escapeHtml(row.feature.feature_label) + '</b>' +
	          '<br>NPC pathway: ' + escapeHtml(featureValue(row.feature, 'npc_pathway')) +
	          '<br>NPC superclass: ' + escapeHtml(featureValue(row.feature, 'npc_superclass')) +
	          '<br>NPC class: ' + escapeHtml(featureValue(row.feature, 'npc_class')) +
	          meanText +
	          '<br>' + (archivedRun ? 'Archived log2 fold change' : 'log2(Right / Left)') + ': ' + row.log2FoldChange.toFixed(3) +
	          '<br>' + (archivedRun ? 'Archived p-value' : 'Welch p-value') + ': ' + row.pValue.toExponential(3);
	      }),
	      hovertemplate: '%{text}<extra></extra>',
	      marker: { color: colors, size: 7, opacity: 0.72, line: { color: 'rgba(255,255,255,0.45)', width: 0.5 } }
	    };
	    const yCutoff = -Math.log10(pCutoff);
	    const layout = {
	      margin: { l: 68, r: 25, t: 20, b: 62 },
	      paper_bgcolor: 'white',
	      plot_bgcolor: '#FAFAFA',
	      font: { family: 'Arial, sans-serif', size: 12, color: '#14161a' },
	      xaxis: { title: archivedRun ? 'Archived log2 fold change (' + archivedRun.contrast + ')' : 'log2 mean intensity fold change (Right / Left)', zeroline: true, zerolinecolor: '#6B7280', gridcolor: '#E5E7EB' },
	      yaxis: { title: archivedRun ? '-log10 archived p-value' : '-log10 Welch p-value', gridcolor: '#E5E7EB' },
	      shapes: [
	        { type: 'line', x0: -foldCutoff, x1: -foldCutoff, y0: 0, y1: 1, yref: 'paper', line: { color: '#6B7280', dash: 'dot', width: 1 } },
	        { type: 'line', x0: foldCutoff, x1: foldCutoff, y0: 0, y1: 1, yref: 'paper', line: { color: '#6B7280', dash: 'dot', width: 1 } },
	        { type: 'line', x0: 0, x1: 1, xref: 'paper', y0: yCutoff, y1: yCutoff, line: { color: '#6B7280', dash: 'dot', width: 1 } }
	      ]
	    };
	    compositionPlotA.innerHTML = '';
	    await Plotly.react(compositionPlotA, [trace], layout, { displaylogo: false, responsive: true, scrollZoom: true });
	    if (renderToken !== compositionRenderToken) return;
	    if (typeof compositionPlotA.removeAllListeners === 'function') compositionPlotA.removeAllListeners('plotly_click');
	    compositionPlotA.on('plotly_click', function(event) {
	      const point = event.points && event.points[0];
	      const featureId = point && point.customdata;
	      if (!featureId) return;
	      levelSelect.value = 'feature';
	      searchInput.value = featureId;
	      resetPage();
	      render();
	      activateWorkspaceTab('overview');
	    });
	  }
	  function updateCompositionViewControls() {
	    const volcano = compositionViewSelect.value === 'volcano';
	    compositionTreemapControls.forEach(function(control) { control.hidden = volcano; });
	    compositionVolcanoControls.forEach(function(control) { control.hidden = !volcano; });
	    compositionGrid.classList.toggle('is-volcano', volcano);
	    compositionCardB.hidden = volcano;
	  }
	  function resetCompositionPlots() {
	    compositionFocusId = 'root';
	    renderComposition();
	  }
	  function renderComposition() {
	    populateCompositionValues();
	    updateCompositionViewControls();
	    const renderToken = ++compositionRenderToken;
	    if (compositionViewSelect.value === 'volcano') {
	      return renderVolcanoPlot(compositionASelect.value, compositionBSelect.value, renderToken).catch(function(error) {
	        showCompositionError(compositionPlotA, error, renderToken);
	      });
	    }
	    return Promise.all([
	      renderCompositionPlot(compositionASelect.value, compositionPlotA, compositionTitleA, compositionMetaA, renderToken).catch(function(error) {
	        showCompositionError(compositionPlotA, error, renderToken);
	      }),
	      renderCompositionPlot(compositionBSelect.value, compositionPlotB, compositionTitleB, compositionMetaB, renderToken).catch(function(error) {
	        showCompositionError(compositionPlotB, error, renderToken);
	      })
	    ]);
	  }
	  function renderDrilldown(entity, targetLevel) {
	    if (!entity || entity.featureIds.length <= 1) {
	      clearDrilldown();
	      return;
	    }
	    const drillEntities = drilldownEntitiesFrom(entity, targetLevel || 'feature');
	    if (!drillEntities.length) {
	      clearDrilldown();
	      return;
	    }
	    const gridColumns = Math.max(1, Number(gridColumnsInput.value) || 4);
	    drilldownGrid.style.setProperty('--grid-columns', String(gridColumns));
	    const targetLabel = targetLevel === 'superclass' ? 'NPC superclass' : (targetLevel === 'class' ? 'NPC class' : 'feature');
	    drilldownTitle.textContent = 'Drill down to ' + targetLabel + ': ' + entity.label;
	    drilldownTab.textContent = 'Drill down: ' + entity.label.replace(/\s*\([0-9]+\s+features\)\s*$/, '');
	    drilldownGrid.innerHTML = '';
	    drilldownTab.hidden = false;
	    activateWorkspaceTab('drilldown');
	    const pendingPlots = [];
	    drillEntities.forEach(function(drillEntity, index) {
	      const card = document.createElement('article');
	      card.className = 'plot-card';
	      card.dataset.entityId = drillEntity.id;
	      card.addEventListener('click', function() { renderInfo(drillEntity); });
	      const header = document.createElement('div');
	      header.className = 'plot-card-header';
	      const title = document.createElement('div');
	      title.className = 'plot-card-title';
	      title.textContent = drillEntity.label;
	      const meta = document.createElement('div');
	      meta.className = 'plot-card-meta';
	      meta.textContent = drillEntity.featureIds.length === 1 ? 'feature' : drillEntity.featureIds.length + ' f.';
	      const actions = document.createElement('div');
	      actions.className = 'plot-card-actions';
	      actions.appendChild(meta);
	      const drillOptions = drilldownTargetOptions(drillEntity);
	      if (drillOptions.length) {
	        const drillSelect = document.createElement('select');
	        drillSelect.className = 'plot-card-drilldown';
	        drillSelect.title = 'Open a nested drill-down dashboard tab';
	        addOption(drillSelect, '', 'Drill down');
	        drillOptions.forEach(function(option) {
	          addOption(drillSelect, option.value, option.label);
	        });
	        drillSelect.addEventListener('click', function(event) {
	          event.stopPropagation();
	        });
	        drillSelect.addEventListener('change', function(event) {
	          event.stopPropagation();
	          if (!drillSelect.value) return;
	          renderInfo(drillEntity);
	          renderDrilldown(drillEntity, drillSelect.value);
	          drillSelect.value = '';
	        });
	        actions.appendChild(drillSelect);
	      }
	      header.appendChild(title);
	      header.appendChild(actions);
	      const plot = document.createElement('div');
	      plot.className = 'plot mini';
	      plot.id = 'drilldown_plot_' + index;
	      card.appendChild(header);
	      card.appendChild(plot);
	      drilldownGrid.appendChild(card);
	      pendingPlots.push({ entity: drillEntity, plot: plot });
	    });
	    window.requestAnimationFrame(function() {
	      pendingPlots.forEach(function(item) {
	        Promise.resolve(plotEntity(item.entity, item.plot)).then(function() {
	          resetPlotZoom(item.plot);
	        });
	      });
	    });
	  }
		  function render() {
			    populateItems();
			    const entities = sortEntities(buildEntities());
		    renderSharedLegend();
		    const selected = entities.find(function(entity) { return entity.id === itemSelect.value; }) || entities[0];
	    const gridColumns = Math.max(1, Number(gridColumnsInput.value) || 4);
	    const gridRows = Math.max(1, Number(gridRowsInput.value) || 10);
	    const gridMaxPlots = gridColumns * gridRows;
	    const totalPages = gridToggle.checked ? Math.max(1, Math.ceil(entities.length / gridMaxPlots)) : 1;
	    currentPage = Math.min(Math.max(1, currentPage), totalPages);
	    const pageStart = (currentPage - 1) * gridMaxPlots;
	    const pageEnd = pageStart + gridMaxPlots;
	    plotGrid.style.setProperty('--grid-columns', gridToggle.checked ? String(gridColumns) : '1');
	    const toPlot = gridToggle.checked ? entities.slice(pageStart, pageEnd) : (selected ? [selected] : []);
	    paginationBar.hidden = !gridToggle.checked || entities.length <= gridMaxPlots;
	    paginationLabel.textContent = 'Page ' + currentPage + ' / ' + totalPages;
	    paginationCount.textContent = entities.length ? (pageStart + 1) + '-' + Math.min(pageEnd, entities.length) + ' of ' + entities.length : '0 of 0';
	    previousPageButton.disabled = currentPage <= 1;
	    nextPageButton.disabled = currentPage >= totalPages;
	    plotGrid.innerHTML = '';
	    clearDrilldown(false);
    const pendingPlots = [];
	    toPlot.forEach(function(entity, index) {
	      const card = document.createElement('article');
	      card.className = 'plot-card';
	      card.dataset.entityId = entity.id;
	      card.addEventListener('click', function() {
	        itemSelect.value = entity.id;
	        renderInfo(entity);
	      });
      const header = document.createElement('div');
      header.className = 'plot-card-header';
      const title = document.createElement('div');
      title.className = 'plot-card-title';
      title.textContent = entity.label;
	      const meta = document.createElement('div');
	      meta.className = 'plot-card-meta';
	      meta.textContent = entity.featureIds.length + ' f.';
	      const actions = document.createElement('div');
	      actions.className = 'plot-card-actions';
	      actions.appendChild(meta);
	      const drillOptions = drilldownTargetOptions(entity);
	      if (drillOptions.length) {
	        const drillSelect = document.createElement('select');
	        drillSelect.className = 'plot-card-drilldown';
	        drillSelect.title = 'Open a drill-down dashboard tab';
	        addOption(drillSelect, '', 'Drill down');
	        drillOptions.forEach(function(option) {
	          addOption(drillSelect, option.value, option.label);
	        });
	        drillSelect.addEventListener('click', function(event) {
	          event.stopPropagation();
	        });
	        drillSelect.addEventListener('change', function(event) {
	          event.stopPropagation();
	          if (!drillSelect.value) return;
	          itemSelect.value = entity.id;
	          renderInfo(entity);
	          renderDrilldown(entity, drillSelect.value);
	          drillSelect.value = '';
	        });
	        actions.appendChild(drillSelect);
	      }
	      header.appendChild(title);
	      header.appendChild(actions);
      const plot = document.createElement('div');
      plot.className = 'plot' + (toPlot.length > 1 ? ' mini' : '');
      plot.id = 'plot_' + index;
      card.appendChild(header);
      card.appendChild(plot);
      plotGrid.appendChild(card);
      pendingPlots.push({ entity: entity, plot: plot });
    });
	    window.requestAnimationFrame(function() {
	      pendingPlots.forEach(function(item) {
	        Promise.resolve(plotEntity(item.entity, item.plot)).then(function() {
	          resetPlotZoom(item.plot);
	        });
	      });
	    });
		    renderInfo(selected);
		    if (!compositionPanel.hidden) renderComposition();
		  }

  populateMetadataSelect(groupSelect, false);
	  populateMetadataSelect(facetSelect, true);
	  populateMetadataSelect(colorSelect, false);
	  populateMetadataSelect(compositionGroupSelect, false);
	  addOption(colorSelect, '__group__', 'Same as x-axis grouping');
	  colorSelect.value = '__group__';
		  sampleFilterRows.forEach(function(row) {
		    populateMetadataSelect(row.querySelector('[data-filter-column]'), true);
		    populateFilterValues(row);
		  });
		  populateTermSelect(pathwaySelect, 'npc_pathway');
		  populateTermSelect(superclassSelect, 'npc_superclass');
		  populateTermSelect(classSelect, 'npc_class');
		  populateTermSelect(componentSelect, 'component_id');
		  Array.from(document.querySelectorAll('select[multiple]')).forEach(enableClickToggleMultiSelect);
		  populateNumericFilters();
		  populateSmilesColumnSelect();
	  groupSelect.value = defaultGroupingColumn();
	  compositionGroupSelect.value = groupSelect.value;
	  compositionViewSelect.value = 'treemap';
	  compositionVolcanoSourceSelect.value = 'archived';
	  levelSelect.value = 'pathway';
	  valueModeSelect.value = 'raw';
	  sortModeSelect.value = 'max_desc';
	  gridToggle.checked = true;
	  populateItems();
	  populateCompositionValues();
	  const compositionOnlyControls = [
	    compositionGroupSelect, compositionASelect, compositionBSelect, compositionDepthSelect,
	    compositionValueModeSelect, compositionViewSelect, compositionVolcanoSourceSelect, compositionPValueInput, compositionFoldChangeInput
	  ];
		  [
	    levelSelect, searchInput, globalSearchInput, itemSelect, groupSelect, facetSelect, colorSelect,
	    pathwaySelect, superclassSelect, classSelect, valueModeSelect,
	    componentSelect, dropSingletonsToggle, plotTypeSelect, sortModeSelect, legendModeSelect, gridToggle, pointsToggle, gridColumnsInput, gridRowsInput, groupOrderInput, groupColorsInput,
	    compositionGroupSelect, compositionASelect, compositionBSelect, compositionDepthSelect, compositionValueModeSelect,
	    compositionViewSelect, compositionVolcanoSourceSelect, compositionPValueInput, compositionFoldChangeInput
		  ].forEach(function(control) {
		    control.addEventListener('change', function() {
		      if (compositionOnlyControls.indexOf(control) !== -1) {
		        if ([compositionGroupSelect, compositionASelect, compositionBSelect, compositionDepthSelect, compositionViewSelect].indexOf(control) !== -1) compositionFocusId = 'root';
		        renderComposition();
		        return;
		      }
		      if (control !== itemSelect) resetPage();
		      render();
		    });
		    control.addEventListener('input', function() {
		      if (control === searchInput || control === globalSearchInput || control === gridColumnsInput || control === gridRowsInput || control === groupOrderInput || control === groupColorsInput) {
		        resetPage();
		        render();
		      }
		    });
		  });
	  sampleFilterRows.forEach(function(row) {
	    const columnSelect = row.querySelector('[data-filter-column]');
	    const valuesSelect = row.querySelector('[data-filter-values]');
	    columnSelect.addEventListener('change', function() {
	      populateFilterValues(row);
	      resetPage();
	      render();
	    });
	    valuesSelect.addEventListener('change', function() {
	      resetPage();
	      render();
	    });
	  });
	  previousPageButton.addEventListener('click', function() {
	    currentPage = Math.max(1, currentPage - 1);
	    render();
	  });
	  nextPageButton.addEventListener('click', function() {
	    currentPage += 1;
	    render();
	  });
	  document.querySelector('[data-clear-npc-filters]').addEventListener('click', function() {
	    pathwaySelect.selectedIndex = -1;
	    superclassSelect.selectedIndex = -1;
	    classSelect.selectedIndex = -1;
	    componentSelect.selectedIndex = -1;
	    dropSingletonsToggle.checked = false;
	    resetPage();
	    render();
	  });
	  overviewTab.addEventListener('click', function() { activateWorkspaceTab('overview'); });
	  drilldownTab.addEventListener('click', function() { activateWorkspaceTab('drilldown'); });
	  compositionTab.addEventListener('click', function() { activateWorkspaceTab('composition'); });
	  drilldownCloseButton.addEventListener('click', clearDrilldown);
	  filterToggle.addEventListener('click', function() {
	    const collapsed = !filterPanel.classList.contains('is-collapsed');
	    filterPanel.classList.toggle('is-collapsed', collapsed);
	    filterToggle.textContent = collapsed ? '›' : '‹';
	    filterToggle.title = collapsed ? 'Expand filters' : 'Collapse filters';
	    window.requestAnimationFrame(function() {
	      Array.from(document.querySelectorAll('.plot')).forEach(function(plot) {
	        if (plot) Plotly.Plots.resize(plot);
	      });
	    });
	  });
	  inspectorToggle.addEventListener('click', function() {
	    const collapsed = !inspectorPanel.classList.contains('is-collapsed');
	    inspectorPanel.classList.toggle('is-collapsed', collapsed);
	    inspectorToggle.textContent = collapsed ? '‹' : '›';
	    inspectorToggle.title = collapsed ? 'Expand selection details' : 'Collapse selection details';
	    window.requestAnimationFrame(function() {
	      Array.from(document.querySelectorAll('.plot')).forEach(function(plot) {
	        if (plot) Plotly.Plots.resize(plot);
	      });
	    });
	  });
	  smilesColumnSelect.addEventListener('change', function() {
	    const entities = sortEntities(buildEntities());
	    const selected = entities.find(function(entity) { return entity.id === itemSelect.value; }) || entities[0];
	    renderInfo(selected);
	  });
	  compositionResetButton.addEventListener('click', resetCompositionPlots);
	  useCurrentGroupsButton.addEventListener('click', function() {
	    seedCurrentGroups();
	    resetPage();
	    render();
	  });
	  function resetPlotZoom(plot) {
	    if (!plot || !plot.layout) return;
	    const update = {};
	    Object.keys(plot.layout).forEach(function(key) {
	      if (/^[xy]axis[0-9]*$/.test(key)) {
	        update[key + '.autorange'] = true;
	      }
	    });
	    Plotly.relayout(plot, update).then(function() {
	      Plotly.Plots.resize(plot);
	    });
	  }
	  document.querySelector('[data-fit-plots]').addEventListener('click', function() {
	    Array.from(document.querySelectorAll('.plot')).forEach(resetPlotZoom);
	  });
	  document.querySelector('[data-clear-filters]').addEventListener('click', function() {
	    searchInput.value = '';
	    globalSearchInput.value = '';
	    pathwaySelect.selectedIndex = -1;
	    superclassSelect.selectedIndex = -1;
	    classSelect.selectedIndex = -1;
	    componentSelect.selectedIndex = -1;
	    dropSingletonsToggle.checked = false;
	    numericFilterPanel.querySelectorAll('[data-numeric-column]').forEach(function(slider) {
	      slider.value = slider.dataset.numericMin;
	      const valueLabel = slider.closest('[data-numeric-filter]').querySelector('[data-numeric-value]');
	      if (valueLabel) valueLabel.textContent = '>= ' + formatSliderValue(slider.value);
	    });
	    sampleFilterRows.forEach(function(row) {
	      row.querySelector('[data-filter-column]').value = '';
	      row.querySelector('[data-filter-values]').innerHTML = '';
	    });
	    groupOrderInput.value = '';
	    groupColorsInput.value = '';
	    Object.keys(hiddenLegendValues).forEach(function(column) { hiddenLegendValues[column].clear(); });
	    compositionFocusId = 'root';
	    resetPage();
	    render();
	  });
  render();
});
)"

writeLines(data_explorer_css, con = css_file, useBytes = TRUE)
writeLines(data_explorer_js, con = js_file, useBytes = TRUE)

dashboard <- htmltools::browsable(tags$html(
  tags$head(
    tags$title("MAPP data explorer"),
    tags$link(rel = "stylesheet", href = paste0("data_explorer_assets/data_explorer.css?v=", asset_version)),
    tags$script(src = paste0("data_explorer_assets/data_explorer_payload.js?v=", asset_version))
  ),
  tags$body(
    tags$div(
      class = "app-container",
      tags$nav(
        class = "app-rail",
        tags$div(class = "rail-mark", "M"),
        tags$div(class = "rail-chip is-active", "Explore"),
        tags$div(class = "rail-spacer"),
        tags$div(class = "rail-dot")
      ),
      tags$aside(
        class = "app-sidebar",
        `data-filter-panel` = "",
        tags$div(
          class = "sidebar-header",
          tags$div(tags$div(class = "brand-title", "MAPP Data Explorer"), tags$div(class = "brand-subtitle", "Feature and NPC-level exploration")),
          tags$button(class = "sidebar-toggle", type = "button", `data-filter-toggle` = "", title = "Collapse filters", "‹")
        ),
        tags$div(
          class = "sidebar-content",
          tags$section(
            class = "sidebar-card",
            tags$div(class = "card-header", "Entity"),
            tags$div(
              class = "form-stack",
	              tags$label(class = "form-group", tags$span("Level"), tags$select(`data-level` = "", tags$option(value = "feature", "Feature"), tags$option(value = "component", "Component"), tags$option(value = "class", "NPC class"), tags$option(value = "superclass", "NPC superclass"), tags$option(value = "pathway", selected = "selected", "NPC pathway"))),
              tags$label(class = "form-group", tags$span("Search"), tags$input(`data-search` = "", type = "search", placeholder = "Ceramides, feature id, pathway...")),
              tags$label(class = "form-group", tags$span("Entity"), tags$select(`data-item` = ""))
            )
          ),
          tags$section(
            class = "sidebar-card",
            tags$div(class = "card-header", "Display"),
            tags$div(
              class = "form-stack",
	              tags$label(class = "form-group", tags$span("Value"), tags$select(`data-value-mode` = "", tags$option(value = "raw", selected = "selected", "Raw intensity"), tags$option(value = "log10", "log10 raw intensity + 1"), tags$option(value = "percent_total", "Relative to sample total (%)"))),
		              tags$label(class = "form-group", tags$span("Plot type"), tags$select(`data-plot-type` = "", tags$option(value = "box", "Box"), tags$option(value = "violin", "Violin"), tags$option(value = "points", "Points"))),
		              tags$label(class = "form-group", tags$span("Plot order"), tags$select(`data-sort-mode` = "", tags$option(value = "label", "Name / current order"), tags$option(value = "max_desc", selected = "selected", "Decreasing max intensity"), tags$option(value = "mean_desc", "Decreasing mean intensity"))),
		              tags$label(class = "form-group", tags$span("Legend"), tags$select(`data-legend-mode` = "", tags$option(value = "shared", "Shared dashboard legend"), tags$option(value = "per_plot", "Legend in each plot"))),
		              tags$label(class = "checkbox-row", tags$input(type = "checkbox", `data-grid-toggle` = "", checked = "checked"), tags$span("Grid mode")),
	              tags$label(class = "checkbox-row", tags$input(type = "checkbox", `data-points-toggle` = "", checked = "checked"), tags$span("Show points"))
            )
          ),
          tags$section(
            class = "sidebar-card",
            tags$div(class = "card-header", "Plot mapping"),
	            tags$div(
	              class = "form-stack",
	              tags$label(class = "form-group", tags$span("Group on x-axis"), tags$select(`data-group` = "")),
	              tags$label(class = "form-group", tags$span("Facet by"), tags$select(`data-facet` = "")),
	              tags$label(class = "form-group", tags$span("Color by"), tags$select(`data-color` = "")),
	              tags$button(class = "btn", `data-use-current-groups` = "", "Use visible groups"),
	              tags$label(
	                class = "form-group",
	                tags$span("Group order"),
	                tags$textarea(`data-group-order` = "", placeholder = "Untreated\nPET 0.5\nPS 0.5\nPET 100\nPS 100")
	              ),
	              tags$label(
	                class = "form-group",
	                tags$span("Group colors"),
	                tags$textarea(`data-group-colors` = "", placeholder = "Untreated = #4B5563\nPET 0.5 = #2563EB\nPS 0.5 = #059669")
	              )
	            )
	          ),
	          tags$section(
	            class = "sidebar-card",
	            tags$div(class = "card-header", "Samples"),
	            tags$div(
	              class = "form-stack",
	              tags$div(
	                class = "sample-filter-row",
	                `data-sample-filter` = "",
	                tags$label(class = "form-group", tags$span("Metadata"), tags$select(`data-filter-column` = "")),
	                tags$label(class = "form-group", tags$span("Values"), tags$select(`data-filter-values` = "", multiple = "multiple"))
	              ),
	              tags$div(
	                class = "sample-filter-row",
	                `data-sample-filter` = "",
	                tags$label(class = "form-group", tags$span("Metadata"), tags$select(`data-filter-column` = "")),
	                tags$label(class = "form-group", tags$span("Values"), tags$select(`data-filter-values` = "", multiple = "multiple"))
	              ),
	              tags$div(
	                class = "sample-filter-row",
	                `data-sample-filter` = "",
	                tags$label(class = "form-group", tags$span("Metadata"), tags$select(`data-filter-column` = "")),
	                tags$label(class = "form-group", tags$span("Values"), tags$select(`data-filter-values` = "", multiple = "multiple"))
	              )
	            )
	          ),
          tags$section(
            class = "sidebar-card",
            tags$div(class = "card-header", "NPC filters"),
            tags$div(
              class = "form-stack",
	              tags$label(class = "form-group", tags$span("Pathway"), tags$select(`data-pathway` = "", multiple = "multiple")),
		              tags$label(class = "form-group", tags$span("Superclass"), tags$select(`data-superclass` = "", multiple = "multiple")),
		              tags$label(class = "form-group", tags$span("Class"), tags$select(`data-class` = "", multiple = "multiple")),
		              tags$label(class = "form-group", tags$span("Component index"), tags$select(`data-component` = "", multiple = "multiple")),
		              tags$label(class = "checkbox-row", tags$input(type = "checkbox", `data-drop-singletons` = ""), tags$span("Drop singletons")),
		              tags$button(class = "btn", `data-clear-npc-filters` = "", "Clear NPC filters")
		            )
		          ),
		          tags$section(
		            class = "sidebar-card",
		            tags$div(class = "card-header", "Annotation scores"),
		            tags$div(class = "numeric-filter-panel", `data-numeric-filters` = "")
		          ),
	          tags$section(
            class = "sidebar-card",
            tags$div(class = "card-header", "Grid"),
            tags$div(
              class = "form-grid",
              tags$label(class = "form-group", tags$span("Columns"), tags$input(`data-grid-columns` = "", type = "number", min = "1", max = "8", value = "4")),
              tags$label(class = "form-group", tags$span("Rows"), tags$input(`data-grid-rows` = "", type = "number", min = "1", max = "30", value = "10"))
            )
          )
        )
      ),
      tags$main(
        class = "app-main",
        tags$header(
          class = "app-header-bar",
          tags$div(class = "header-title", tags$h1(`data-title` = "", "MAPP data explorer"), tags$div(class = "header-summary", `data-summary` = ""), tags$div(class = "header-summary", `data-count` = "")),
          tags$input(class = "header-search", `data-global-search` = "", type = "search", placeholder = "Global search"),
          tags$div(class = "header-actions", tags$button(class = "btn", `data-fit-plots` = "", "Fit plots"), tags$button(class = "btn", `data-clear-filters` = "", "Clear filters"))
        ),
		        tags$section(
			          class = "workspace",
			          tags$div(
			            class = "workspace-content",
			            tags$div(class = "shared-legend", `data-shared-legend` = ""),
		          tags$div(
		            class = "workspace-tabs",
		            tags$button(class = "workspace-tab is-active", `data-tab-overview` = "", "Overview"),
		            tags$button(class = "workspace-tab", `data-tab-composition` = "", "Composition"),
		            tags$button(class = "workspace-tab", `data-tab-drilldown` = "", hidden = "hidden", "Exploded features")
		          ),
		          tags$section(
		            class = "tab-panel",
		            `data-overview-panel` = "",
		            tags$div(class = "plot-grid", `data-plot-grid` = ""),
		            tags$div(
		              class = "pagination-bar",
		              `data-pagination` = "",
		              hidden = "hidden",
		              tags$button(class = "btn", `data-page-previous` = "", "Previous"),
		              tags$span(class = "pagination-label", `data-page-label` = "", "Page 1 / 1"),
		              tags$span(class = "pagination-count", `data-page-count` = "", "0 of 0"),
		              tags$button(class = "btn", `data-page-next` = "", "Next")
		            )
		          ),
		          tags$section(
		            class = "composition-panel tab-panel",
		            `data-composition-panel` = "",
		            hidden = "hidden",
		            tags$div(
		              class = "composition-toolbar",
			              tags$label(class = "form-group", tags$span("Compare by"), tags$select(`data-composition-group` = "")),
			              tags$label(class = "form-group", tags$span("Left"), tags$select(`data-composition-a` = "")),
			              tags$label(class = "form-group", tags$span("Right"), tags$select(`data-composition-b` = "")),
			              tags$label(
			                class = "form-group",
			                tags$span("View"),
			                tags$select(
			                  `data-composition-view` = "",
			                  tags$option(value = "treemap", "Treemaps"),
			                  tags$option(value = "volcano", "Volcano plot")
			                )
			              ),
			              tags$label(
			                class = "form-group",
			                `data-composition-treemap-control` = "",
			                tags$span("Depth"),
		                tags$select(
		                  `data-composition-depth` = "",
		                  tags$option(value = "class", "Class"),
		                  tags$option(value = "superclass", "Superclass"),
		                  tags$option(value = "pathway", "Pathway"),
		                  tags$option(value = "feature", "Feature")
		                )
		              ),
			              tags$label(
			                class = "form-group",
			                `data-composition-treemap-control` = "",
			                tags$span("Value"),
		                tags$select(
		                  `data-composition-value-mode` = "",
		                  tags$option(value = "percent_selected", "Displayed %"),
		                  tags$option(value = "percent_sample_total", "Sample total %"),
			                  tags$option(value = "raw", "Raw")
			                )
			              ),
			              tags$label(
			                class = "form-group",
			                `data-composition-volcano-control` = "",
			                hidden = "hidden",
			                tags$span("Volcano source"),
			                tags$select(
			                  `data-composition-volcano-source` = "",
			                  tags$option(value = "archived", "Archived statistics (exact)"),
			                  tags$option(value = "dynamic", "Dynamic Welch")
			                )
			              ),
			              tags$label(
			                class = "form-group",
			                `data-composition-volcano-control` = "",
			                hidden = "hidden",
			                tags$span("P-value cutoff"),
			                tags$input(`data-composition-p-value` = "", type = "number", min = "0.000001", max = "1", step = "0.01", value = "0.05")
			              ),
			              tags$label(
			                class = "form-group",
			                `data-composition-volcano-control` = "",
			                hidden = "hidden",
			                tags$span("|log2 fold| cutoff"),
			                tags$input(`data-composition-fold-change` = "", type = "number", min = "0", step = "0.25", value = "1")
			              ),
			              tags$button(class = "btn", type = "button", `data-composition-reset` = "", "Reset view")
			            ),
			            tags$div(
			              class = "composition-grid",
			              `data-composition-grid` = "",
		              tags$article(
		                class = "composition-card",
		                `data-composition-card-a` = "",
		                tags$div(
		                  class = "composition-card-header",
		                  tags$div(class = "composition-card-title", `data-composition-title-a` = "", "Left group"),
		                  tags$div(
		                    class = "plot-card-actions",
		                    tags$div(class = "composition-card-meta", `data-composition-meta-a` = "")
		                  )
		                ),
		                tags$div(class = "composition-plot", `data-composition-plot-a` = "")
		              ),
			              tags$article(
			                class = "composition-card",
			                `data-composition-card-b` = "",
		                tags$div(
		                  class = "composition-card-header",
		                  tags$div(class = "composition-card-title", `data-composition-title-b` = "", "Right group"),
		                  tags$div(
		                    class = "plot-card-actions",
		                    tags$div(class = "composition-card-meta", `data-composition-meta-b` = "")
		                  )
		                ),
		                tags$div(class = "composition-plot", `data-composition-plot-b` = "")
		              )
		            )
		          ),
		          tags$section(
		            class = "drilldown-panel tab-panel",
		            `data-drilldown` = "",
		            hidden = "hidden",
	            tags$div(
	              class = "drilldown-header",
	              tags$div(class = "drilldown-title", `data-drilldown-title` = ""),
	              tags$button(class = "btn", `data-drilldown-close` = "", "Close")
		            ),
		            tags$div(class = "drilldown-grid", `data-drilldown-grid` = "")
		          )
		          ),
		          tags$aside(
		            class = "selection-panel",
		            `data-inspector-panel` = "",
		            tags$div(
		              class = "inspector-header",
		              tags$div(class = "inspector-heading", "Selection details"),
		              tags$button(class = "inspector-toggle", type = "button", `data-inspector-toggle` = "", title = "Collapse selection details", "›")
		            ),
		            tags$div(
		              class = "inspector-body",
		              tags$section(
		                class = "inspector-config",
		                tags$label(class = "form-group", tags$span("SMILES column"), tags$select(`data-smiles-column` = ""))
		              ),
		              tags$div(`data-info` = "")
		            )
		          )
		        )
      ),
      tags$div(class = "hidden-dependency", dummy_plotly),
      tags$script(src = paste0("data_explorer_assets/data_explorer.js?v=", asset_version))
    )
  )
))

htmltools::save_html(dashboard, file = html_file, libdir = "data_explorer_files")
message(sprintf("Data explorer written to %s", html_file))
message(sprintf("Payload written to %s", payload_file))
