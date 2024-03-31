library(yaml)
library(dplyr)
library(digest)
library(digest)
library(plyr)



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


common_df_path <- "common_dataframe.rds"
common_tsv_path <- "common_dataframe.tsv"



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


# Assuming the YAML content is loaded and converted to a single-row

# Load the YAML content from a file or a string
your_yaml_string <- "
mapp_project: mapp_project_00017
mapp_batch: mapp_batch_00044
dataset_experiment:
  name: \"mapp_batch_00044 LCMS metabolomics dataset\"
  description: \"Simon Blanchoud - Tunicate metabolomics\"
actions:
  scale_data: TRUE
  filter_sample_type: TRUE
  filter_sample_metadata_one: TRUE
  test_that: TRUE
  export_data: TRUE
options:
  gnps_column_for_boxplots:
    factor_name:
        - 'process'
        - 'sample_type'
        - 'sample_group'
filter_sample_type:
  mode: 'include'
  factor_name: 'sample_type'
  levels:
   - 'BK'
   - 'BLANK'
"
yaml_content <- yaml.load(your_yaml_string)

# Convert the YAML content to a dataframe row
new_row_df <- convert_yaml_to_single_row_df_with_hash(yaml_content)

# Append this row to the common dataframe and save
append_to_common_df_and_save(new_row_df, common_df_path, common_tsv_path)