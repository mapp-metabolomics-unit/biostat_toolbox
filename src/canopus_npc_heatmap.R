#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(iheatmapr)
  library(htmlwidgets)
})

option_list <- list(
  make_option(c("-c", "--canopus-file"), default = NULL, help = "Path to canopus_structure_summary.tsv"),
  make_option(c("-q", "--quant-file"), default = NULL, help = "Optional override for the MZmine quant CSV"),
  make_option(c("-m", "--metadata-file"), default = NULL, help = "Optional override for the treated metadata TSV"),
  make_option(c("-o", "--output-html"), default = NULL, help = "Output HTML file [default inferred next to the CANOPUS file]"),
  make_option(c("--output-table"), default = NULL, help = "Optional TSV with the plotted matrix and annotations [default inferred next to the HTML]"),
  make_option(c("--sample-type"), default = "sample", help = "Metadata sample_type value to keep; use 'all' to keep every column [default %default]"),
  make_option(c("--sample-annotation"), default = "ATTRIBUTE_condition", help = "Metadata column used as the main sample annotation [default %default]"),
  make_option(c("--top-n"), type = "integer", default = 300, help = "Top N annotated features ranked by total intensity; use 0 for all [default %default]"),
  make_option(c("--transform"), default = "log10", help = "Intensity transform: log10 or none [default %default]"),
  make_option(c("--scale"), default = "row_zscore", help = "Heatmap scaling: row_zscore or none [default %default]"),
  make_option(c("--cluster-rows"), default = "FALSE", help = "Cluster rows: TRUE or FALSE [default %default]"),
  make_option(c("--cluster-cols"), default = "FALSE", help = "Cluster columns: TRUE or FALSE [default %default]"),
  make_option(c("--width"), type = "integer", default = 2200, help = "Widget width in pixels [default %default]"),
  make_option(c("--height"), type = "integer", default = 3200, help = "Widget height in pixels [default %default]")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

normalize_option_name <- function(opt, preferred_name, alternative_names) {
  current_value <- opt[[preferred_name]]
  if (!is.null(current_value) && length(current_value)) {
    return(opt)
  }
  for (alt_name in alternative_names) {
    alt_value <- opt[[alt_name]]
    if (!is.null(alt_value) && length(alt_value)) {
      opt[[preferred_name]] <- alt_value
      return(opt)
    }
  }
  opt
}

trim_option_value <- function(value) {
  if (is.character(value) && length(value)) {
    trimmed <- trimws(value)
    return(trimmed[nzchar(trimmed)])
  }
  value
}

as_flag <- function(value, label) {
  if (is.null(value) || !length(value)) {
    return(FALSE)
  }
  if (is.logical(value)) {
    return(value[1])
  }
  value <- tolower(trimws(as.character(value[1])))
  if (value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("Invalid logical value '%s' for %s.", value, label))
}

opt <- normalize_option_name(opt, "canopus_file", c("canopus-file"))
opt <- normalize_option_name(opt, "quant_file", c("quant-file"))
opt <- normalize_option_name(opt, "metadata_file", c("metadata-file"))
opt <- normalize_option_name(opt, "output_html", c("output-html"))
opt <- normalize_option_name(opt, "output_table", c("output-table"))
opt <- normalize_option_name(opt, "sample_type", c("sample-type"))
opt <- normalize_option_name(opt, "sample_annotation", c("sample-annotation"))
opt <- normalize_option_name(opt, "top_n", c("top-n"))
opt <- normalize_option_name(opt, "cluster_rows", c("cluster-rows", "cluster.rows"))
opt <- normalize_option_name(opt, "cluster_cols", c("cluster-cols", "cluster.cols"))

coerce_scalar_path <- function(value, label) {
  if (is.null(value) || !length(value)) {
    stop(sprintf("No value provided for %s.", label))
  }
  value <- as.character(value[1])
  if (!nzchar(trimws(value))) {
    stop(sprintf("Empty path provided for %s.", label))
  }
  value
}

ensure_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path))
  }
}

infer_batch_dir <- function(canopus_file) {
  dirname(dirname(dirname(canopus_file)))
}

infer_quant_file <- function(canopus_file) {
  batch_dir <- infer_batch_dir(canopus_file)
  batch_name <- basename(batch_dir)
  file.path(batch_dir, "results", "mzmine", paste0(batch_name, "_quant.csv"))
}

infer_metadata_file <- function(canopus_file) {
  batch_dir <- infer_batch_dir(canopus_file)
  batch_name <- basename(batch_dir)
  file.path(batch_dir, "metadata", "treated", paste0(batch_name, "_metadata.tsv"))
}

infer_output_html <- function(canopus_file) {
  file.path(dirname(canopus_file), "canopus_npc_feature_heatmap.html")
}

default_output_table <- function(output_html) {
  file.path(dirname(output_html), paste0(tools::file_path_sans_ext(basename(output_html)), "_data.tsv"))
}

make_named_palette <- function(values, palette_name = "Set 3") {
  values <- unique(values[!is.na(values)])
  values <- sort(as.character(values))
  if (!length(values)) {
    return(setNames(character(), character()))
  }
  colors <- grDevices::hcl.colors(length(values), palette = palette_name)
  setNames(colors, values)
}

get_fixed_npc_palettes <- function() {
  list(
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

build_fixed_npc_color_maps <- function(feature_meta) {
  palettes <- get_fixed_npc_palettes()

  feature_meta_ranked <- feature_meta %>%
    mutate(
      npc_pathway = dplyr::coalesce(npc_pathway, "Unclassified"),
      npc_superclass = dplyr::coalesce(npc_superclass, "Other"),
      npc_class = dplyr::coalesce(npc_class, "Other")
    )

  superclass_rank_df <- feature_meta_ranked %>%
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

  pathway_colors <- feature_meta_ranked %>%
    distinct(npc_pathway) %>%
    mutate(
      palette_family = dplyr::coalesce(unname(palettes$pathway_family[npc_pathway]), "micro_cvd_gray"),
      hex = purrr::map_chr(palette_family, ~ palettes$shade_hex[[.x]][1])
    ) %>%
    select(npc_pathway, hex) %>%
    tibble::deframe()

  superclass_colors <- superclass_rank_df %>%
    distinct(npc_superclass, hex) %>%
    tibble::deframe()

  class_colors <- feature_meta_ranked %>%
    distinct(npc_pathway, npc_superclass, npc_class) %>%
    left_join(
      superclass_rank_df %>% select(npc_pathway, npc_superclass, hex),
      by = c("npc_pathway", "npc_superclass")
    ) %>%
    mutate(hex = dplyr::coalesce(hex, pathway_colors[npc_pathway], "#F5F5F5")) %>%
    distinct(npc_class, hex) %>%
    tibble::deframe()

  row_track_colors <- feature_meta_ranked %>%
    left_join(
      superclass_rank_df %>% select(npc_pathway, npc_superclass, hex),
      by = c("npc_pathway", "npc_superclass")
    ) %>%
    mutate(
      row_track = paste(npc_pathway, npc_superclass, sep = " - "),
      hex = dplyr::coalesce(hex, pathway_colors[npc_pathway], "#F5F5F5")
    ) %>%
    distinct(row_track, hex) %>%
    tibble::deframe()

  list(
    pathway = pathway_colors,
    superclass = superclass_colors,
    class = class_colors,
    row_track = row_track_colors
  )
}

build_annotation_spec <- function(df, palette_lookup) {
  keep_cols <- names(df)[vapply(df, function(col) length(unique(col[!is.na(col)])) > 1, logical(1))]
  if (!length(keep_cols)) {
    return(NULL)
  }
  list(
    data = as.data.frame(df[, keep_cols, drop = FALSE], check.names = FALSE),
    colors = palette_lookup[keep_cols]
  )
}

clean_intensity_name <- function(column_name) {
  str_remove(column_name, " Peak height$")
}

row_zscore <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- apply(x, 1, stats::sd, na.rm = TRUE)
  scaled <- sweep(x, 1, row_means, FUN = "-")
  scaled <- sweep(scaled, 1, row_sds, FUN = "/")
  scaled[is.na(scaled)] <- 0
  scaled[is.infinite(scaled)] <- 0
  scaled
}

build_hover_matrix <- function(display_mat, raw_mat, feature_meta, sample_meta) {
  hover <- matrix("", nrow = nrow(display_mat), ncol = ncol(display_mat))
  for (i in seq_len(nrow(display_mat))) {
    for (j in seq_len(ncol(display_mat))) {
      hover[i, j] <- paste0(
        "Feature ID: ", feature_meta$feature_id[i],
        "<br>Sample: ", sample_meta$sample_label[j],
        "<br>Filename: ", sample_meta$filename[j],
        "<br>Main annotation: ", sample_meta$sample_annotation[j],
        "<br>Sample type: ", sample_meta$sample_type[j],
        "<br>NPC pathway: ", feature_meta$npc_pathway[i],
        "<br>NPC superclass: ", feature_meta$npc_superclass[i],
        "<br>NPC class: ", feature_meta$npc_class[i],
        "<br>m/z: ", sprintf("%.5f", feature_meta$ion_mass[i]),
        "<br>RT (min): ", sprintf("%.3f", feature_meta$retention_time_min[i]),
        "<br>Raw intensity: ", format(raw_mat[i, j], scientific = TRUE, digits = 4),
        "<br>Displayed value: ", sprintf("%.3f", display_mat[i, j])
      )
    }
  }
  hover
}

opt$canopus_file <- trim_option_value(opt$canopus_file)
opt$quant_file <- trim_option_value(opt$quant_file)
opt$metadata_file <- trim_option_value(opt$metadata_file)
opt$output_html <- trim_option_value(opt$output_html)
opt$output_table <- trim_option_value(opt$output_table)
opt$sample_type <- trim_option_value(opt$sample_type)
opt$sample_annotation <- trim_option_value(opt$sample_annotation)
opt$transform <- tolower(trim_option_value(opt$transform))
opt$scale <- tolower(trim_option_value(opt$scale))
opt$cluster_rows <- as_flag(opt$cluster_rows, "--cluster-rows")
opt$cluster_cols <- as_flag(opt$cluster_cols, "--cluster-cols")

if (is.null(opt$canopus_file)) {
  stop("Please provide --canopus-file.")
}

canopus_file <- normalizePath(coerce_scalar_path(opt$canopus_file, "--canopus-file"), mustWork = FALSE)
quant_file <- if (!is.null(opt$quant_file)) normalizePath(coerce_scalar_path(opt$quant_file, "--quant-file"), mustWork = FALSE) else infer_quant_file(canopus_file)
metadata_file <- if (!is.null(opt$metadata_file)) normalizePath(coerce_scalar_path(opt$metadata_file, "--metadata-file"), mustWork = FALSE) else infer_metadata_file(canopus_file)
output_html <- if (!is.null(opt$output_html)) normalizePath(coerce_scalar_path(opt$output_html, "--output-html"), mustWork = FALSE) else infer_output_html(canopus_file)
output_table <- if (!is.null(opt$output_table)) normalizePath(coerce_scalar_path(opt$output_table, "--output-table"), mustWork = FALSE) else default_output_table(output_html)

ensure_exists(canopus_file, "CANOPUS file")
ensure_exists(quant_file, "Quant file")
ensure_exists(metadata_file, "Metadata file")

dir.create(dirname(output_html), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)

if (!opt$transform %in% c("log10", "none")) {
  stop("Invalid --transform. Choose 'log10' or 'none'.")
}
if (!opt$scale %in% c("row_zscore", "none")) {
  stop("Invalid --scale. Choose 'row_zscore' or 'none'.")
}
if (is.null(opt$top_n) || is.na(opt$top_n) || opt$top_n < 0) {
  stop("--top-n must be a non-negative integer.")
}

canopus_df <- readr::read_tsv(
  canopus_file,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  transmute(
    feature_id = as.integer(mappingFeatureId),
    aligned_feature_id = as.character(alignedFeatureId),
    ion_mass = as.numeric(ionMass),
    retention_time_min = as.numeric(retentionTimeInMinutes),
    npc_pathway = coalesce(`NPC#pathway`, "Unclassified"),
    npc_superclass = coalesce(`NPC#superclass`, "Unclassified"),
    npc_class = coalesce(`NPC#class`, "Unclassified")
  ) %>%
  filter(!is.na(feature_id)) %>%
  distinct(feature_id, .keep_all = TRUE)

quant_df <- readr::read_csv(
  quant_file,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "unique_quiet"
)

intensity_cols <- names(quant_df)[str_detect(names(quant_df), " Peak height$")]
if (!length(intensity_cols)) {
  stop("No intensity columns ending with ' Peak height' were found in the quant file.")
}

metadata_df <- readr::read_tsv(
  metadata_file,
  show_col_types = FALSE,
  progress = FALSE
)

if (!"filename" %in% names(metadata_df)) {
  stop("The metadata file must contain a 'filename' column.")
}
if (!"sample_type" %in% names(metadata_df)) {
  stop("The metadata file must contain a 'sample_type' column.")
}
if (!opt$sample_annotation %in% names(metadata_df)) {
  stop(sprintf("Sample annotation column '%s' was not found in the metadata file.", opt$sample_annotation))
}

sample_manifest <- tibble(intensity_col = intensity_cols) %>%
  mutate(filename = vapply(intensity_col, clean_intensity_name, FUN.VALUE = character(1))) %>%
  left_join(metadata_df, by = "filename")

if (any(is.na(sample_manifest$sample_type))) {
  missing_files <- sample_manifest$filename[is.na(sample_manifest$sample_type)]
  stop(sprintf(
    "Metadata rows are missing for %d quant columns. First missing filename: %s",
    length(missing_files),
    missing_files[1]
  ))
}

if (!tolower(opt$sample_type) %in% c("all", "*")) {
  sample_manifest <- sample_manifest %>%
    filter(tolower(sample_type) == tolower(opt$sample_type))
  if (!nrow(sample_manifest)) {
    stop(sprintf("No columns remain after filtering sample_type == '%s'.", opt$sample_type))
  }
}

quant_selected <- quant_df %>%
  transmute(
    feature_id = as.integer(`row ID`),
    !!!rlang::syms(sample_manifest$intensity_col)
  )

feature_table <- canopus_df %>%
  inner_join(quant_selected, by = "feature_id") %>%
  mutate(
    total_intensity = rowSums(across(all_of(sample_manifest$intensity_col)), na.rm = TRUE)
  ) %>%
  arrange(npc_pathway, npc_superclass, npc_class, desc(total_intensity), feature_id)

if (opt$top_n > 0) {
  feature_table <- feature_table %>%
    slice_head(n = min(opt$top_n, nrow(feature_table)))
}

if (!nrow(feature_table)) {
  stop("No annotated features remained after joining CANOPUS and the quant table.")
}

raw_mat <- feature_table %>%
  select(all_of(sample_manifest$intensity_col)) %>%
  as.matrix()
rownames(raw_mat) <- as.character(feature_table$feature_id)
colnames(raw_mat) <- sample_manifest$filename
storage.mode(raw_mat) <- "numeric"
raw_mat[is.na(raw_mat)] <- 0

display_mat <- raw_mat
if (opt$transform == "log10") {
  display_mat <- log10(display_mat + 1)
}
if (opt$scale == "row_zscore") {
  display_mat <- row_zscore(display_mat)
}

sample_manifest <- sample_manifest %>%
  mutate(
    sample_type = dplyr::coalesce(as.character(sample_type), "Missing"),
    sample_annotation = dplyr::coalesce(as.character(.data[[opt$sample_annotation]]), "Missing"),
    sample_label = dplyr::coalesce(if ("sample_id" %in% names(sample_manifest)) as.character(sample_id) else NA_character_, filename)
  ) %>%
  arrange(sample_annotation, sample_type, sample_label)

display_mat <- display_mat[, sample_manifest$filename, drop = FALSE]
raw_mat <- raw_mat[, sample_manifest$filename, drop = FALSE]

feature_meta <- feature_table %>%
  transmute(
    feature_id = as.character(feature_id),
    aligned_feature_id = aligned_feature_id,
    ion_mass = ion_mass,
    retention_time_min = retention_time_min,
    npc_pathway = npc_pathway,
    npc_superclass = npc_superclass,
    npc_class = npc_class,
    row_track = paste(npc_pathway, npc_superclass, sep = " - "),
    row_label = paste0("F", feature_id, " | ", npc_class)
  )

hover_mat <- build_hover_matrix(display_mat, raw_mat, feature_meta, sample_manifest)

sample_annotation_palette <- make_named_palette(sample_manifest$sample_annotation, "Set 2")
sample_type_palette <- make_named_palette(sample_manifest$sample_type, "Dark 3")
npc_color_maps <- build_fixed_npc_color_maps(feature_meta)

row_annotation_spec <- build_annotation_spec(
  data.frame(
    Pathway = feature_meta$npc_pathway,
    Superclass = feature_meta$npc_superclass,
    Class = feature_meta$npc_class,
    check.names = FALSE
  ),
  list(
    Pathway = npc_color_maps$pathway,
    Superclass = npc_color_maps$superclass,
    Class = npc_color_maps$class
  )
)

col_annotation_spec <- build_annotation_spec(
  data.frame(
    SampleType = sample_manifest$sample_type,
    Group = sample_manifest$sample_annotation,
    check.names = FALSE
  ),
  list(
    SampleType = sample_type_palette,
    Group = sample_annotation_palette
  )
)

legend_title <- if (opt$scale == "row_zscore") {
  "Intensity (row z-score)"
} else if (opt$transform == "log10") {
  "Intensity (log10 peak height + 1)"
} else {
  "Intensity"
}

heatmap_title <- paste(
  "CANOPUS NPC feature heatmap",
  paste("Features:", nrow(display_mat)),
  paste("Samples:", ncol(display_mat)),
  paste("Source:", basename(canopus_file)),
  sep = "<br>"
)

widget_height <- max(opt$height, 900 + nrow(display_mat) * 8)

iheatmap <- iheatmapr::main_heatmap(
  display_mat,
  name = legend_title,
  colors = "RdBu",
  show_colorbar = TRUE,
  text = hover_mat,
  layout = list(
    width = opt$width,
    height = widget_height,
    title = list(text = heatmap_title, x = 0.05),
    margin = list(t = 110, r = 240, b = 140, l = 240),
    hovermode = "closest"
  )
) %>%
  add_row_labels(
    tickvals = NULL,
    ticktext = feature_meta$row_label,
    side = "left",
    buffer = 0.01,
    size = 0.35,
    font = list(size = 8)
  )

if (!is.null(row_annotation_spec)) {
  iheatmap <- iheatmap %>%
    add_row_annotation(
      row_annotation_spec$data,
      side = "right",
      buffer = 0.03,
      colors = row_annotation_spec$colors
    )
}

if (!is.null(col_annotation_spec)) {
  iheatmap <- iheatmap %>%
    add_col_annotation(
      col_annotation_spec$data,
      side = "top",
      buffer = 0.03,
      colors = col_annotation_spec$colors
    )
}

iheatmap <- iheatmap %>%
  add_col_labels(
    tickvals = NULL,
    ticktext = sample_manifest$sample_label,
    textangle = -90,
    size = 0.45,
    font = list(size = 10)
  )

if (opt$cluster_rows) {
  iheatmap <- iheatmap %>% add_row_clustering(side = "right")
}
if (opt$cluster_cols) {
  iheatmap <- iheatmap %>% add_col_clustering()
}

save_iheatmap(iheatmap, file = output_html)

output_table_df <- feature_table %>%
  select(
    feature_id,
    aligned_feature_id,
    ion_mass,
    retention_time_min,
    npc_pathway,
    npc_superclass,
    npc_class,
    all_of(sample_manifest$intensity_col)
  )

names(output_table_df) <- c(
  "feature_id",
  "aligned_feature_id",
  "ion_mass",
  "retention_time_min",
  "npc_pathway",
  "npc_superclass",
  "npc_class",
  sample_manifest$filename
)

readr::write_tsv(output_table_df, output_table)

message("Saved heatmap HTML to: ", output_html)
message("Saved plotted data to: ", output_table)
