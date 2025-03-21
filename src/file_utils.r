#!/usr/bin/env Rscript

#' File utilities module for biostat toolbox
#' This module handles file operations, path management, and naming conventions

#' Generate output filenames based on parameters
#' @param params Configuration parameters
#' @return List of filenames for various outputs
generate_filenames <- function(params) {
  # Generate file prefix based on parameters
  file_prefix <- ""
  
  # Define all output filenames
  filenames <- list(
    box_plots = paste0(file_prefix, "Boxplots.pdf"),
    box_plots_interactive = paste0(file_prefix, "Boxplots_interactive.html"),
    DE = paste0(file_prefix, "DE.rds"),
    DE_original = paste0(file_prefix, "DE_original.rds"),
    DE_description = paste0(file_prefix, "DE_description.txt"),
    DE_original_description = paste0(file_prefix, "DE_original_description.txt"),
    foldchange_pvalues = paste0(file_prefix, "foldchange_pvalues.csv"),
    formatted_peak_table = paste0(file_prefix, "formatted_peak_table.csv"),
    formatted_sample_data_table = paste0(file_prefix, "formatted_sample_data_table.csv"),
    formatted_sample_metadata = paste0(file_prefix, "formatted_sample_metadata.tsv"),
    formatted_variable_metadata = paste0(file_prefix, "formatted_variable_metadata.csv"),
    graphml = paste0(file_prefix, "graphml.graphml"),
    heatmap_pval = paste0(file_prefix, "Heatmap_pval.html"),
    heatmap_rf = paste0(file_prefix, "Heatmap_rf.html"),
    interactive_table = paste0(file_prefix, "interactive_table.html"),
    metaboverse_table = paste0(file_prefix, "metaboverse_table.tsv"),
    PCA = paste0(file_prefix, "PCA.pdf"),
    PCA3D = paste0(file_prefix, "PCA3D.html"),
    PCoA = paste0(file_prefix, "PCoA.pdf"),
    PCoA3D = paste0(file_prefix, "PCoA3D.html"),
    PLSDA = paste0(file_prefix, "PLSDA.pdf"),
    PLSDA_loadings = paste0(file_prefix, "PLSDA_loadings.tsv"),
    PLSDA_VIP_plot = paste0(file_prefix, "PLSDA_VIP.pdf"),
    PLSDA_VIP_table = paste0(file_prefix, "PLSDA_VIP.tsv"),
    R_script = paste0(file_prefix, "R_script_backup.R"),
    random_forest = paste0(file_prefix, "RF_importance.html"),
    random_forest_model = paste0(file_prefix, "RF_model.txt"),
    session_info = paste0(file_prefix, "session_info.txt"),
    summary_stats_table_full = paste0(file_prefix, "summary_stats_table_full.csv"),
    summary_stats_table_selected = paste0(file_prefix, "summary_stats_table_selected.csv"),
    summary_stat_output_selected_cytoscape = paste0(file_prefix, "summary_stats_table_selected_cytoscape.csv"),
    treemap = paste0(file_prefix, "Treemap_interactive.html"),
    volcano = paste0(file_prefix, "Volcano.pdf"),
    volcano_interactive = paste0(file_prefix, "Volcano_interactive.html")
  )
  
  return(filenames)
}

#' Setup output directory structure
#' @param params Configuration parameters
#' @param hash Configuration hash
#' @return Path to output directory
setup_output_directory <- function(params, hash) {
  if (params$paths$output != "") {
    output_directory <- file.path(params$paths$output, hash)
  } else {
    working_directory <- file.path(params$paths$docs, params$mapp_project, params$mapp_batch)
    output_directory <- file.path(working_directory, "results", "stats", hash)
  }
  
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Directory created:", output_directory, "\n")
  } else {
    message("Directory already exists:", output_directory, "\n")
  }
  
  return(output_directory)
}

#' Save configuration files to output directory
#' @param params_path Path to main params file
#' @param params_user_path Path to user params file
#' @param output_directory Output directory path
#' @param filenames List of output filenames
save_config_files <- function(params_path, params_user_path, output_directory, filenames) {
  message("Writing params.yaml ...")
  file.copy(params_path, file.path(output_directory, filenames$params), overwrite = TRUE)
  file.copy(params_user_path, file.path(output_directory, filenames$params_user), overwrite = TRUE)
}

#' String sanitization function
#' @param string String to sanitize
#' @return Sanitized string
sanitize_string <- function(string) {
  string <- gsub("_{2,}", "_", string)
  string <- gsub("^_|_$", "", string)
  return(string)
}

#' Format filter status
#' @param filter Filter configuration
#' @return Formatted filter status string
formatted_filter_status <- function(filter) {
  return(paste(filter$mode, filter$factor_name, paste(filter$levels, collapse = "_"), sep = "_"))
} 