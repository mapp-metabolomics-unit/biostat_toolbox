# A set of helpers functions to be used in the main script

# Revised function to ensure a newline at the end of YAML files
ensure_newline <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) > 0 && nchar(tail(lines, 1)) != 0) {
    lines <- c(lines, "")
    writeLines(lines, file_path)
  }
}


# Flatten list with alphabetical sorting of lists and dot notation for nested keys
flatten_list <- function(x, prefix = "") {
  flattened <- vector("list")
  
  for (name in names(x)) {
    current_path <- if (prefix == "") name else paste(prefix, name, sep = ".")
    
    if (is.list(x[[name]])) {
      # Recursive call for lists
      flattened <- c(flattened, flatten_list(x[[name]], current_path))
    } else {
      # Handling atomic vectors and single elements, with sorting for vectors
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

convert_yaml_to_single_row_df_with_hash <- function(yaml_content) {
  # Convert the YAML content into a structured list
  flattened_list <- flatten_list(yaml_content)
  # Compute the hash of the flattened list
  content_hash <- digest(flattened_list, algo = "md5")
  # Convert the list to a dataframe row
  df <- as.data.frame(t(unlist(flattened_list)), stringsAsFactors = FALSE)
  colnames(df) <- names(flattened_list)
  # Append the hash to the dataframe
  df$hash <- content_hash
  # Append the timestamp to the dataframe
  df$timestamp <- Sys.time()
  return(df)
}


append_to_common_df_and_save <- function(new_row_df, common_df_path, common_tsv_path) {
  if (file.exists(common_df_path)) {
    common_df <- readRDS(common_df_path)
  } else {
    common_df <- data.frame(matrix(nrow = 0, ncol = 0))
  }

  # Ensure 'hash' column is at the beginning
  new_row_df <- bind_cols(new_row_df %>% select(hash, timestamp), new_row_df %>% select(-hash, -timestamp))

  # Check for unique hash before appending
  if (!(new_row_df$hash %in% common_df$hash)) {
    # Append the new row with timestamp
    common_df <- bind_rows(common_df, new_row_df)
    
    # Sort columns alphabetically
    common_df <- common_df %>%
      select(hash, timestamp, everything()) %>%
      select(hash, timestamp, sort(names(.)[-c(1, 2)]))

    # Save the updated common dataframe
    saveRDS(common_df, common_df_path)
    write.table(common_df, common_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    message("Content hash already exists in the dataframe. No new row added.")
  }
}

# String sanitization function

sanitize_string <- function(string) {
  # We make sure that no multiple _ exists in the filter_variable_metadata_status string
  string <- gsub("_{2,}", "_", string)
  # We also make sure that the string doesn not start or finish with an underscore
  string <- gsub("^_|_$", "", string)
  return(string)
}


# The previous lines are functionalized in

formatted_filter_status <- function(filter) {
  return(paste(filter$mode, filter$factor_name, paste(filter$levels, collapse = "_"), sep = "_"))
}
