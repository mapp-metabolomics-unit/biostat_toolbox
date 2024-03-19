# A set of helpers functions to be used in the main script

generate_hash_from_yaml <- function(yaml_file_path) {
  yaml_content <- readLines(yaml_file_path)
  yaml_hash <- digest(paste(yaml_content, collapse = ""), algo = "md5")
  return(yaml_hash)
}



# update_mapping_file <- function(hash, description, mapping_file_path) {
#   line <- paste(hash, description, sep = ": ")
#   if (!file.exists(mapping_file_path)) {
#     write(line, mapping_file_path)
#   } else {
#     # Append without repeating the same hash
#     existing_hashes <- readLines(mapping_file_path)
#     if (!any(grepl(hash, existing_hashes))) {
#       write(line, mapping_file_path, append = TRUE)
#     }
#   }
# }


# Updated function to write parameters and their hash to a TSV file
update_mapping_file <- function(params, yaml_hash, mapping_file_path) {
  # Extract parameters
  params_values <- c(
    params$mapp_project,
    params$mapp_batch,
    params$polarity,
    params$actions$run_with_gap_filled,
    params$actions$scale_data,
    params$actions$filter_sample_type,
    params$actions$filter_sample_metadata_one,
    params$actions$filter_sample_metadata_two,
    params$actions$filter_variable_metadata_one,
    params$actions$filter_variable_metadata_two,
    params$actions$filter_variable_metadata_annotated,
    params$actions$filter_variable_metadata_num,
    params$actions$filter_features,
    # Assuming the following is a correction to avoid repetition and correctly list intended parameters
    params$actions$short_path,
    params$actions$run_cytoscape_connector
  )
  
  # Create a line with the hash and the parameters separated by tabs
  line <- paste(yaml_hash, paste(params_values, collapse = "\t"), sep = "\t")
  
  # Check if the mapping file exists, if not, create it with headers
  if (!file.exists(mapping_file_path)) {
    headers <- c("hash", "mapp_project", "mapp_batch", "polarity", "run_with_gap_filled", 
                 "scale_data", "filter_sample_type", "filter_sample_metadata_one", 
                 "filter_sample_metadata_two", "filter_variable_metadata_one", 
                 "filter_variable_metadata_two", "filter_variable_metadata_annotated", 
                 "filter_variable_metadata_num", "filter_features", "short_path", 
                 "run_cytoscape_connector")
    write(paste(headers, collapse = "\t"), mapping_file_path)
  }
  
  # Append the new line to the mapping file
  write(line, mapping_file_path, append = TRUE)
}