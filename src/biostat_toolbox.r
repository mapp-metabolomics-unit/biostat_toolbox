#!/usr/bin/env Rscript

#' Main pipeline script for biostat toolbox
#' This script orchestrates the workflow by calling functions from the modularized components

# Source all required modules
source("src/package_management.r")
source("src/config_handler.r")
source("src/file_utils.r")
source("src/data_loader.r")

# Initialize packages
initialize_packages()

# Set paths to configuration files
path_to_params <- "./params/params.yaml"
path_to_params_user <- "./params/params_user.yaml"

# Load configuration
params <- load_configuration(path_to_params, path_to_params_user)

# Setup working directory
working_directory <- file.path(params$paths$docs, params$mapp_project, params$mapp_batch)

# Convert params to dataframe and get hash
new_row_df <- convert_yaml_to_single_row_df_with_hash(params)
common_tsv_path <- file.path(params$paths$output, "params_log.tsv")
append_to_common_df_and_save(new_row_df, common_tsv_path)

# Setup output directory
output_directory <- setup_output_directory(params, new_row_df$hash)

# Generate filenames
filenames <- generate_filenames(params)

# Save configuration files
save_config_files(path_to_params, path_to_params_user, output_directory, filenames)

# Load and process feature table
feature_table <- load_feature_table(working_directory, params)
feature_intensities <- extract_feature_intensities(feature_table)
feature_metadata <- extract_feature_metadata(feature_table)

# Load annotations
sirius_annotations_filename <- "compound_identifications.tsv"
canopus_annotations_filename <- "canopus_compound_summary.tsv"

data_sirius <- load_sirius_annotations(working_directory, sirius_annotations_filename)
data_canopus <- load_canopus_annotations(working_directory, canopus_annotations_filename)
data_met_annot <- load_met_annot_data(working_directory, params)

# Check GNPS job type
gnps2_job <- check_gnps2_job(working_directory, params)

# TODO: Add statistical analysis functions
# TODO: Add visualization functions
# TODO: Add export functions

message("Biostat toolbox pipeline completed successfully!")