# A set of helpers functions to be used in the main script

range01 <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }

  min_value <- min(x, na.rm = TRUE)
  max_value <- max(x, na.rm = TRUE)

  if (identical(min_value, max_value)) {
    return(rep(0, length(x)))
  }

  (x - min_value) / (max_value - min_value)
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
  # Add a human-readable description without making it part of the hash.
  df$description <- describe_params_for_hash_table(yaml_content)
  return(df)
}

param_has_value <- function(value) {
  !is.null(value) && length(value) && !all(is.na(value)) && any(nzchar(trimws(as.character(unlist(value)))))
}

format_param_vector <- function(value) {
  if (!param_has_value(value)) {
    return("")
  }
  values <- trimws(as.character(unlist(value, use.names = FALSE)))
  values <- values[!is.na(values) & nzchar(values)]
  paste(values, collapse = ", ")
}

describe_filter_param <- function(filter_param, label) {
  if (is.null(filter_param) || !param_has_value(filter_param$mode)) {
    return(character(0))
  }
  mode <- as.character(filter_param$mode[1])
  if (!mode %in% c("include", "exclude", "above", "below")) {
    return(character(0))
  }
  factor_name <- as.character(filter_param$factor_name[1])
  levels <- format_param_vector(filter_param$levels)
  level <- format_param_vector(filter_param$level)

  if (mode %in% c("include", "exclude")) {
    if (!nzchar(factor_name) || !nzchar(levels)) {
      return(character(0))
    }
    return(sprintf("%s %s %s in %s", label, mode, factor_name, levels))
  }
  if (!nzchar(factor_name) || !nzchar(level)) {
    return(character(0))
  }
  sprintf("%s %s %s %s", label, factor_name, mode, level)
}

describe_params_for_hash_table <- function(params) {
  target <- if (param_has_value(params$target$sample_metadata_header)) {
    as.character(params$target$sample_metadata_header[1])
  } else {
    "selected target"
  }
  compared_groups <- format_param_vector(params$colors$all$key)
  comparison <- if (nzchar(compared_groups)) {
    sprintf("Compare %s across %s", compared_groups, target)
  } else {
    sprintf("Compare groups across %s", target)
  }

  filters <- c(
    describe_filter_param(params$filter_sample_type, "samples"),
    describe_filter_param(params$filter_sample_metadata_one, "samples"),
    describe_filter_param(params$filter_sample_metadata_two, "samples"),
    describe_filter_param(params$filter_variable_metadata_one, "features"),
    describe_filter_param(params$filter_variable_metadata_two, "features"),
    describe_filter_param(params$filter_variable_metadata_annotated, "features"),
    describe_filter_param(params$filter_variable_metadata_num, "features")
  )
  filter_text <- if (length(filters)) {
    paste("Filters:", paste(filters, collapse = "; "))
  } else {
    "No explicit sample/feature filters"
  }

  scaling <- if (param_has_value(params$actions$scale_method)) {
    paste("Scaling:", as.character(params$actions$scale_method[1]))
  } else {
    "Scaling: unspecified"
  }

  npc_terms <- c(
    if (param_has_value(params$npc_summed_intensity$pathway)) paste("NPC pathway", format_param_vector(params$npc_summed_intensity$pathway)),
    if (param_has_value(params$npc_summed_intensity$superclass)) paste("NPC superclass", format_param_vector(params$npc_summed_intensity$superclass)),
    if (param_has_value(params$npc_summed_intensity$class)) paste("NPC class", format_param_vector(params$npc_summed_intensity$class))
  )
  npc_text <- if (length(npc_terms)) {
    paste("NPC plots:", paste(npc_terms, collapse = "; "))
  } else {
    "NPC plots: none"
  }

  paste(comparison, filter_text, scaling, npc_text, sep = ". ")
}


# append_to_common_df_and_save <- function(new_row_df, common_df_path, common_tsv_path) {
#   if (file.exists(common_df_path)) {
#     common_df <- readRDS(common_df_path)
#   } else {
#     common_df <- data.frame(matrix(nrow = 0, ncol = 0))
#   }

#   # Ensure 'hash' column is at the beginning
#   new_row_df <- bind_cols(new_row_df %>% select(hash, timestamp), new_row_df %>% select(-hash, -timestamp))

#   # Check for unique hash before appending
#   if (!(new_row_df$hash %in% common_df$hash)) {
#     # Append the new row with timestamp
#     common_df <- bind_rows(common_df, new_row_df)
    
#     # Sort columns alphabetically
#     common_df <- common_df %>%
#       select(hash, timestamp, everything()) %>%
#       select(hash, timestamp, sort(names(.)[-c(1, 2)]))

#     # Save the updated common dataframe
#     saveRDS(common_df, common_df_path)
#     write.table(common_df, common_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
#   } else {
#     message("Content hash already exists in the dataframe. No new row added.")
#   }
# }

append_to_common_df_and_save <- function(new_row_df, common_tsv_path) {
  if (file.exists(common_tsv_path)) {
    # Read the existing common dataframe from TSV
    common_df <- read.table(common_tsv_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")

    # Convert the timestamp column to character for consistency
    common_df$timestamp <- as.character(common_df$timestamp)
  } else {
    # Create an empty dataframe with appropriate columns if the TSV doesn't exist
    common_df <- data.frame(hash = character(), timestamp = character(), stringsAsFactors = FALSE)
  }

  # Ensure 'new_row_df' has the required columns
  if (!("hash" %in% names(new_row_df)) || !("timestamp" %in% names(new_row_df))) {
    stop("The new_row_df must contain 'hash' and 'timestamp' columns.")
  }

  # Ensure 'timestamp' is character for compatibility
  new_row_df$timestamp <- as.character(new_row_df$timestamp)

  # Align column types between common_df and new_row_df
  for (col in intersect(names(common_df), names(new_row_df))) {
    common_df[[col]] <- type.convert(as.character(common_df[[col]]), as.is = TRUE)
    new_row_df[[col]] <- type.convert(as.character(new_row_df[[col]]), as.is = TRUE)
  }

  # Check if the hash already exists in common_df
  if (nrow(common_df) == 0) {
    common_df <- new_row_df
  } else {
    if (new_row_df$hash %in% common_df$hash) {
      # Update the existing row so descriptions and newly added columns can be refreshed.
      index <- which(common_df$hash == new_row_df$hash)
      missing_cols <- setdiff(names(new_row_df), names(common_df))
      for (col in missing_cols) {
        common_df[[col]] <- NA
      }
      missing_cols <- setdiff(names(common_df), names(new_row_df))
      for (col in missing_cols) {
        new_row_df[[col]] <- NA
      }
      common_df[index, names(new_row_df)] <- new_row_df[1, names(new_row_df)]
      message("Content hash already exists. Timestamp updated.")
    } else {
      # Append the new row if the hash is unique
      common_df <- bind_rows(common_df, new_row_df)
    }
  }

  # Sort columns alphabetically while keeping the human-facing columns first.
  common_df <- common_df %>%
    select(hash, timestamp, any_of("description"), everything()) %>%
    select(hash, timestamp, any_of("description"), sort(names(.)[!names(.) %in% c("hash", "timestamp", "description")]))

  # Save the updated common dataframe as a TSV file
  write.table(common_df, common_tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)
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
