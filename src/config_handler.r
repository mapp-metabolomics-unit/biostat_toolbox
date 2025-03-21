#!/usr/bin/env Rscript

#' Configuration handler module for biostat toolbox
#' This module handles YAML configuration loading and parameter management

library(yaml)
library(digest)
library(dplyr)

#' Ensure a newline at the end of YAML files
#' @param file_path Path to the YAML file
ensure_newline <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) > 0 && nchar(tail(lines, 1)) != 0) {
    lines <- c(lines, "")
    writeLines(lines, file_path)
  }
}

#' Flatten list with alphabetical sorting and dot notation for nested keys
#' @param x List to flatten
#' @param prefix Prefix for nested keys
flatten_list <- function(x, prefix = "") {
  flattened <- vector("list")
  
  for (name in names(x)) {
    current_path <- if (prefix == "") name else paste(prefix, name, sep = ".")
    
    if (is.list(x[[name]])) {
      flattened <- c(flattened, flatten_list(x[[name]], current_path))
    } else {
      value <- x[[name]]
      if (is.atomic(value) && length(value) > 1) {
        value <- paste(sort(value), collapse = "|")
      }
      if (is.logical(value)) {
        value <- as.character(value)
      }
      flattened[[current_path]] <- value
    }
  }
  return(flattened)
}

#' Convert YAML to single row dataframe with hash
#' @param yaml_content YAML content to convert
convert_yaml_to_single_row_df_with_hash <- function(yaml_content) {
  flattened_list <- flatten_list(yaml_content)
  content_hash <- digest(flattened_list, algo = "md5")
  df <- as.data.frame(t(unlist(flattened_list)), stringsAsFactors = FALSE)
  colnames(df) <- names(flattened_list)
  df$hash <- content_hash
  df$timestamp <- Sys.time()
  return(df)
}

#' Append configuration to common dataframe and save
#' @param new_row_df New configuration row
#' @param common_tsv_path Path to save the common TSV file
append_to_common_df_and_save <- function(new_row_df, common_tsv_path) {
  if (file.exists(common_tsv_path)) {
    common_df <- read.table(common_tsv_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    common_df$timestamp <- as.character(common_df$timestamp)
  } else {
    common_df <- data.frame(hash = character(), timestamp = character(), stringsAsFactors = FALSE)
  }

  if (!all(c("hash", "timestamp") %in% names(new_row_df))) {
    stop("The new_row_df must contain 'hash' and 'timestamp' columns.")
  }

  new_row_df$timestamp <- as.character(new_row_df$timestamp)

  for (col in intersect(names(common_df), names(new_row_df))) {
    common_df[[col]] <- type.convert(as.character(common_df[[col]]), as.is = TRUE)
    new_row_df[[col]] <- type.convert(as.character(new_row_df[[col]]), as.is = TRUE)
  }

  if (nrow(common_df) == 0) {
    common_df <- new_row_df
  } else {
    if (new_row_df$hash %in% common_df$hash) {
      index <- which(common_df$hash == new_row_df$hash)
      common_df$timestamp[index] <- new_row_df$timestamp
      message("Content hash already exists. Timestamp updated.")
    } else {
      common_df <- bind_rows(common_df, new_row_df)
    }
  }

  common_df <- common_df %>%
    select(hash, timestamp, everything()) %>%
    select(hash, timestamp, sort(names(.)[-c(1, 2)]))

  write.table(common_df, common_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Load and process configuration files
#' @param params_path Path to main params file
#' @param params_user_path Path to user params file
#' @return List containing processed parameters
load_configuration <- function(params_path, params_user_path) {
  # Ensure newlines at the end of YAML files
  ensure_newline(params_path)
  ensure_newline(params_user_path)
  
  # Load params files
  params <- yaml.load_file(params_path)
  params_user <- yaml.load_file(params_user_path)
  
  # Update params with user settings
  params$paths$docs <- params_user$paths$docs
  params$paths$output <- params_user$paths$output
  params$operating_system$system <- params_user$operating_system$system
  params$operating_system$pandoc <- params_user$operating_system$pandoc
  
  return(params)
} 