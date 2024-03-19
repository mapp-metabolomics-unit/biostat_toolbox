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


# # Updated function to write parameters and their hash to a TSV file
# update_mapping_file <- function(params, yaml_hash, mapping_file_path) {
#   # Extract parameters
#   params_values <- c(
#     params$mapp_project,
#     params$mapp_batch,
#     params$polarity,
#     params$actions$run_with_gap_filled,
#     params$actions$scale_data,
#     params$actions$filter_sample_type,
#     params$actions$filter_sample_metadata_one,
#     params$actions$filter_sample_metadata_two,
#     params$actions$filter_variable_metadata_one,
#     params$actions$filter_variable_metadata_two,
#     params$actions$filter_variable_metadata_annotated,
#     params$actions$filter_variable_metadata_num,
#     params$actions$filter_features,
#     # Assuming the following is a correction to avoid repetition and correctly list intended parameters
#     params$actions$short_path,
#     params$actions$run_cytoscape_connector
#   )
  
#   # Create a line with the hash and the parameters separated by tabs
#   line <- paste(yaml_hash, paste(params_values, collapse = "\t"), sep = "\t")
  
#   # Check if the mapping file exists, if not, create it with headers
#   if (!file.exists(mapping_file_path)) {
#     headers <- c("hash", "mapp_project", "mapp_batch", "polarity", "run_with_gap_filled", 
#                  "scale_data", "filter_sample_type", "filter_sample_metadata_one", 
#                  "filter_sample_metadata_two", "filter_variable_metadata_one", 
#                  "filter_variable_metadata_two", "filter_variable_metadata_annotated", 
#                  "filter_variable_metadata_num", "filter_features", "short_path", 
#                  "run_cytoscape_connector")
#     write(paste(headers, collapse = "\t"), mapping_file_path)
#   }
  
#   # Append the new line to the mapping file
#   write(line, mapping_file_path, append = TRUE)
# }


# Revised function to ensure a newline at the end of YAML files
ensure_newline <- function(file_path) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) > 0 && nchar(tail(lines, 1)) != 0) {
    lines <- c(lines, "")
    writeLines(lines, file_path)
  }
}



update_configuration_json <- function(yaml_file_path, json_file_path) {
  # Load the YAML content into an R object, ignoring comments
  params <- yaml.load_file(yaml_file_path)

  # Serialize the R object (YAML content without comments) to a JSON string
  # Using jsonlite::toJSON for consistent ordering and representation
  # auto_unbox and pretty are used to ensure consistency in the string representation
  params_json_str <- toJSON(params, auto_unbox = TRUE, pretty = TRUE)

  # Generate hash from the JSON string representation of the parameters
  config_hash <- digest(params_json_str, algo = "md5")
  
  # Rest of the function implementation...
  configurations <- list()
  
  # Attempt to read the JSON file content and parse it
  if (file.exists(json_file_path) && file.size(json_file_path) > 0) {
    json_content <- readLines(json_file_path, warn = FALSE)
    json_content <- paste(json_content, collapse = "")
    configurations <- fromJSON(txt = json_content)
  }
  
  if (!is.null(configurations[[config_hash]])) {
    message("Configuration already exists. No update needed.")
  } else {
    configurations[[config_hash]] <- params
    write_json(configurations, json_file_path, pretty = TRUE)
    message("New configuration added to JSON.")
  }
  
  return(config_hash)
}