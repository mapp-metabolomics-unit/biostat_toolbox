library(yaml)
library(digest)
library(dplyr)
library(tidyr)
library(purrr)
library(rio)

df = import("params/test/params.yaml")

# Path for yaml file

path_yaml <- "params/test/params.yaml"

# Read yaml file

params <- yaml.load_file(path_yaml)

# Create a hash of the yaml file

hash <- digest(params)

# Save the hash in the yaml file

params$hash <- hash

# Check the type of params

print(class(params))

# Convert the named list to a dataframe

params_df <- as.data.frame(t(params))

# Unesting the dataframe

str(params_df)


unnested_list <- unlist(params_df, recursive = FALSE) 


str(unnested_list)

unnested_list_df = as.data.frame(unnested_list)

str(unnested_list_df)


collapsed_df <- df %>%
  # all rows are listed alphabetically
  arrange(across(everything())) %>%
  summarise(across(everything(), collapse_different)) %>% 
  as.data.frame()

str(collapsed_df)


# Save the collapsed dataframe as a tsv file

# write.table(collapsed_df, file = "params/test/params.tsv", sep = "\t", row.names = FALSE)


fwrite(collapsed_df, file = "params/test/params.csv")



###########################################################################################
###########################################################################################
###########################################################################################

library(yaml)
library(dplyr)

# Revised flattening function with alphabetical sorting
flatten_list <- function(x, prefix = "") {
  flattened <- vector("list")
  
  for (name in names(x)) {
    current_path <- if (prefix == "") name else paste(prefix, name, sep = ".")
    
    if (is.list(x[[name]])) {
      # Recursive call for lists
      flattened <- c(flattened, flatten_list(x[[name]], current_path))
    } else {
      # Handling atomic vectors and single elements
      value <- x[[name]]
      if (is.atomic(value) && length(value) > 1) {
        # Sort and then concatenate vector values into a single string
        value <- paste(sort(value), collapse = "|")
      }
      # Convert logicals to character to ensure consistency
      if (is.logical(value)) {
        value <- as.character(value)
      }
      flattened[[current_path]] <- value
    }
  }
  return(flattened)
}

# Function to convert the YAML content into a single-row dataframe
convert_yaml_to_df <- function(yaml_content) {
  flattened_list <- flatten_list(yaml_content)
  df <- as.data.frame(t(unlist(flattened_list)), stringsAsFactors = FALSE)
  colnames(df) <- names(flattened_list)
  return(df)
}

# Your YAML content
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
options:
  gnps_column_for_boxplots:
    factor_name:
        - 'process'
        - 'sample_type'
        - 'instrument'
        - 'sample_group'
filter_sample_type:
  mode: 'include'
  factor_name: 'sample_type'
  levels:
   - 'BK' 
   - 'BLANK'
   - 'QC'
"

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
  scale_data: TRUE
options:
  gnps_column_for_boxplots:
    factor_name:
        - 'process'
        - 'sample_type'
        - 'instrument'
        - 'sample_group'
filter_sample_type:
  mode: 'include'
  factor_name: 'sample_type'
  levels:
   - 'BK' 
   - 'BLANK'
   - 'QC'
"

# Load and convert YAML content
yaml_content <- yaml.load(your_yaml_string)
df <- convert_yaml_to_df(yaml_content)

# Display the resulting dataframe
print(df)

