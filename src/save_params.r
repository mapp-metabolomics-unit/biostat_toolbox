library(yaml)
library(digest)
library(dplyr)
library(tidyr)
library(purrr)

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


collapsed_df <- unnested_list_df %>%
  # all rows are listed alphabetically
  arrange(across(everything())) %>%
  summarise(across(everything(), collapse_different)) %>% 
  as.data.frame()

str(collapsed_df)


# Save the collapsed dataframe as a tsv file

# write.table(collapsed_df, file = "params/test/params.tsv", sep = "\t", row.names = FALSE)


fwrite(collapsed_df, file = "params/test/params.csv")
