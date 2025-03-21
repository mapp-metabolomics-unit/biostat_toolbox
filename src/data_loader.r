#!/usr/bin/env Rscript

#' Data loader module for biostat toolbox
#' This module handles loading and initial processing of different data types

library(readr)
library(dplyr)

#' Load feature table from MZmine output
#' @param working_directory Working directory path
#' @param params Configuration parameters
#' @return Processed feature table
load_feature_table <- function(working_directory, params) {
  # Determine input file based on gap filling setting
  if (params$actions$run_with_gap_filled == "TRUE") {
    input_file <- paste0(params$mapp_batch, "_gf_quant.csv")
  } else {
    input_file <- paste0(params$mapp_batch, "_quant.csv")
  }
  
  # Read the feature table
  feature_table <- read_delim(
    file.path(working_directory, "results", "mzmine", input_file),
    delim = ",",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  # Process the feature table
  feature_table <- feature_table %>%
    rename(
      "feature_id" = "row ID",
      "feature_mz" = "row m/z",
      "feature_rt" = "row retention time"
    )
  
  # Create feature_id_full column
  feature_table$feature_id_full <- paste(
    feature_table$feature_id,
    round(feature_table$feature_mz, digits = 2),
    round(feature_table$feature_rt, digits = 1),
    sep = "_"
  )
  
  return(feature_table)
}

#' Extract feature intensities from feature table
#' @param feature_table Full feature table
#' @return Matrix of feature intensities
extract_feature_intensities <- function(feature_table) {
  # Extract intensity data
  feature_table_intensities <- feature_table %>%
    select(feature_id, contains(" Peak height")) %>%
    rename_with(~ gsub(" Peak height", "", .x)) %>%
    column_to_rownames(var = "feature_id") %>%
    as.data.frame() %>%
    t()
  
  # Order by rownames and column names
  feature_table_intensities <- feature_table_intensities[order(row.names(feature_table_intensities)), ]
  feature_table_intensities <- feature_table_intensities[, order(colnames(feature_table_intensities))]
  
  return(as.data.frame(feature_table_intensities))
}

#' Extract feature metadata from feature table
#' @param feature_table Full feature table
#' @return Feature metadata dataframe
extract_feature_metadata <- function(feature_table) {
  feature_metadata <- feature_table %>%
    select(feature_id_full, feature_id, feature_mz, feature_rt)
  return(feature_metadata)
}

#' Load Sirius annotations
#' @param working_directory Working directory path
#' @param sirius_annotations_filename Filename for Sirius annotations
#' @return Processed Sirius annotations
load_sirius_annotations <- function(working_directory, sirius_annotations_filename) {
  chebied_file <- file.path(working_directory, "results", "sirius", paste("chebied", sirius_annotations_filename, sep = "_"))
  
  if (file.exists(chebied_file)) {
    data_sirius <- read_delim(chebied_file,
      delim = "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    )
  } else {
    data_sirius <- read_delim(
      file.path(working_directory, "results", "sirius", sirius_annotations_filename),
      delim = "\t",
      escape_double = FALSE,
      trim_ws = TRUE
    )
    
    # Process Sirius data
    for_chembiid_smiles <- unique(data_sirius$smiles)
    
    print("Getting ChEBI IDs from smiles ...")
    
    if (ping_service("chebi") == FALSE) {
      print("The ChEBI service is down. We will issue an empty DF. Please try again later.")
      chebi_ids <- data.frame(query = for_chembiid_smiles, chebiid = NA, chebiasciiname = NA)
    } else {
      print("The ChEBI service is up.")
      chebi_ids <- get_chebiid(for_chembiid_smiles, from = "smiles", to = "chebiid", match = "best")
    }
    
    # Merge and process data
    data_sirius <- merge(data_sirius, chebi_ids, by.x = "smiles", by.y = "query")
    colnames(data_sirius) <- paste("sirius", colnames(data_sirius), sep = "_")
    data_sirius$feature_id <- as.numeric(data_sirius$sirius_featureId)
    
    # Save processed data
    write.table(data_sirius, chebied_file, sep = "\t", row.names = FALSE)
  }
  
  return(data_sirius)
}

#' Load CANOPUS annotations
#' @param working_directory Working directory path
#' @param canopus_annotations_filename Filename for CANOPUS annotations
#' @return Processed CANOPUS annotations
load_canopus_annotations <- function(working_directory, canopus_annotations_filename) {
  data_canopus <- read_delim(
    file.path(working_directory, "results", "sirius", canopus_annotations_filename),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  colnames(data_canopus) <- paste("canopus", colnames(data_canopus), sep = "_")
  data_canopus$feature_id <- as.numeric(data_canopus$canopus_featureId)
  
  write.table(
    data_canopus,
    file.path(working_directory, "results", "sirius", paste("featured", canopus_annotations_filename, sep = "_")),
    sep = "\t",
    row.names = FALSE
  )
  
  return(data_canopus)
}

#' Load MetAnnot data
#' @param working_directory Working directory path
#' @param params Configuration parameters
#' @return Processed MetAnnot data
load_met_annot_data <- function(working_directory, params) {
  data_met_annot <- read_delim(
    file.path(working_directory, "results", "met_annot_enhancer",
      params$met_annot_enhancer_folder,
      paste0(params$met_annot_enhancer_folder, "_spectral_match_results_repond.tsv")
    ),
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
  )
  
  colnames(data_met_annot) <- paste("met_annot", colnames(data_met_annot), sep = "_")
  data_met_annot$feature_id <- as.numeric(data_met_annot$met_annot_feature_id)
  
  return(data_met_annot)
}

#' Check if GNPS job is from GNPS2 or legacy interface
#' @param working_directory Working directory path
#' @param params Configuration parameters
#' @return Boolean indicating if job is from GNPS2
check_gnps2_job <- function(working_directory, params) {
  gnps2_job <- file.exists(
    Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id))
  )
  
  if (gnps2_job) {
    print("This is a GNPS2 job.")
  } else {
    print("This is a job from the GNPS legacy interface.")
  }
  
  return(gnps2_job)
} 