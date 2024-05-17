############################################################################################
############################################################################################
#####################################     PACKAGES     #####################################
############################################################################################
############################################################################################


# install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# install structToolbox and dependencies
# BiocManager::install("structToolbox") v1.10.0

## install additional bioc packages for vignette if needed
# BiocManager::install(c("pmp", "ropls", "BiocFileCache"))

## install additional CRAN packages if needed
# install.packages(c('cowplot', 'openxlsx'))


# We define the following helper function in order to load or install the packages according to the condition

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[, 1])) {
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}

# This one below is to the the default CRAN repo

r <- getOption("repos")
r["CRAN"] <- "http://cran.us.r-project.org"
options(repos = r)
rm(r)

# Package organized alphabetically

usePackage("ape") ### super package
usePackage("base")
usePackage("BiocFileCache")
usePackage("cowplot")
usePackage("crosstalk")
usePackage("data.table")
usePackage("datasets")
usePackage("dbscan")
usePackage("digest")
usePackage("dplyr")
usePackage("DT")
usePackage("emmeans")
usePackage("EnhancedVolcano")
usePackage("fpc")
usePackage("ggdendro")
usePackage("ggh4x")
usePackage("ggplot2")
usePackage("ggraph")
usePackage("ggrepel")
usePackage("ggtree")
usePackage("graphlayouts")
usePackage("gridExtra")
usePackage("gt")
usePackage("heatmaply")
usePackage("here")
usePackage("htmltools")
usePackage("igraph")
usePackage("iheatmapr")
usePackage("janitor")
usePackage("jsonlite")
usePackage("magick")
usePackage("manhattanly")
usePackage("microshades") ### remotes::install_github("KarstensLab/microshades", dependencies=TRUE)
usePackage("modEvA")
usePackage("openxlsx")
usePackage("plotly")
usePackage("pls")
usePackage("pmp")
usePackage("purrr")
usePackage("randomcoloR")
usePackage("randomForest")
usePackage("rcdk")
usePackage("RColorBrewer")
usePackage("readr")
usePackage("reshape2")
usePackage("reticulate")
usePackage("rfPermute")
usePackage("rgl")
usePackage("rockchalk")
usePackage("ropls")
usePackage("this.path")
usePackage("tidyr")
usePackage("tidyverse")
usePackage("tools")
usePackage("treemap")
usePackage("vegan")
usePackage("viridis")
usePackage("webchem")
usePackage("wesanderson")
usePackage("WikidataQueryServiceR")
usePackage("yaml")

# struct(1.10)

# devtools::install_github("jcheng5/d3scatter")

# install.packages("BiocManager")
# BiocManager::install("RCy3")
# remotes::install_gitlab("artemklevtsov/uchardet@devel")


# We use the MAPPstructToolbox package
# Uncomment the lines below to download the MAPPstructToolbox package from github

# library(devtools)
# install_github("mapp-metabolomics-unit/MAPPstructToolbox", force = TRUE)
library(MAPPstructToolbox)

############################################################################################
############################################################################################
################################ LOAD REQUIRED FUNCTIONS  ##################################
############################################################################################
############################################################################################

# We load the required functions from the MAPPstructToolbox package
# these are in the helpers.r file

source("src/helpers.r")


############################################################################################
############################################################################################
################################ LOAD & FORMAT  DATA  ######################################
############################################################################################
############################################################################################


# We set the wd
current_script <- deparse(substitute())
script_path <- file.path(getwd(), current_script)

print(script_path)

if (!exists("params")) {
  my_path_params <- script_path
} ### conserve the path after multiple run
if (exists("params")) {
  setwd(my_path_params)
} ### conserve the path after multiple run

# We call the external params
path_to_params <- "./params/params.yaml"
path_to_params_user <- "./params/params_user.yaml"

# Ensure newlines at the end of the YAML files


ensure_newline(path_to_params)
ensure_newline(path_to_params_user)

# Load the params.yaml file

params <- yaml.load_file(path_to_params)
params_user <- yaml.load_file(path_to_params_user)


# Here we load the user params if they exist

params$paths$docs <- params_user$paths$docs
params$paths$output <- params_user$paths$output
params$operating_system$system <- params_user$operating_system$system
params$operating_system$pandoc <- params_user$operating_system$pandoc


# We generate a hash from the params.yaml file

# yaml_hash <- generate_hash_from_yaml(path_to_params)

# Description of this configuration
# Here we output a fully formatted description of the configuration using the params.yaml file and it's set parameters


# We set the working directory

working_directory <- file.path(params$paths$docs, params$mapp_project, params$mapp_batch)


# # Path to your mapping file
# mapping_file_path <- file.path(params$paths$output, "mapping_file.tsv")

# update_mapping_file(params, yaml_hash, mapping_file_path)



# We set the output directory

if (params$actions$scale_data == "TRUE") {
  scaling_status <- "scaled"
} else {
  scaling_status <- "raw"
}

possible_modes <- c("exclude", "include", "above", "below", "activated", "deactivated")



filter_sample_type_status = formatted_filter_status(params$filter_sample_type)
filter_sample_metadata_one_status = formatted_filter_status(params$filter_sample_metadata_one)
filter_sample_metadata_two_status = formatted_filter_status(params$filter_sample_metadata_two)

filter_variable_metadata_one_status = formatted_filter_status(params$filter_variable_metadata_one)
filter_variable_metadata_two_status = formatted_filter_status(params$filter_variable_metadata_two)
filter_variable_metadata_annotated_status = formatted_filter_status(params$filter_variable_metadata_annotated)
filter_variable_metadata_num_status = formatted_filter_status(params$filter_variable_metadata_num)

filter_sample_metadata_status = paste(filter_sample_type_status, filter_sample_metadata_one_status, filter_sample_metadata_two_status, sep = "_")

filter_variable_metadata_status = paste(filter_variable_metadata_one_status, filter_variable_metadata_two_status, filter_variable_metadata_annotated_status, filter_variable_metadata_num_status, sep = "_")

# We make sure that no multiple _ exists using the sanitize_string function

filter_sample_metadata_status = sanitize_string(filter_sample_metadata_status)
filter_variable_metadata_status = sanitize_string(filter_variable_metadata_status)



#################################################################################################
#################################################################################################
################### Filename and paths establishment ##########################################
#################################################################################################


# The Figures filename is conditionally defined according to the user's choice of filtering the dataset according to CANOPUS NPClassifier classifications or not.


# file_prefix = paste(params$mapp_batch,
#                    params$target$sample_metadata_header,
#                    filter_variable_metadata_status,
#                    filter_sample_metadata_status,
#                    params$polarity,
#                    scaling_status,
#                    sep = "_")

# file_prefix = paste(params$mapp_batch,
#                     params$target$sample_metadata_header,
#                     sep = "_")

file_prefix <- paste("")


filename_box_plots <- paste(file_prefix, "Boxplots.pdf", sep = "")
filename_box_plots_interactive <- paste(file_prefix, "Boxplots_interactive.html", sep = "")
filename_DE <- paste(file_prefix, "DE.rds", sep = "")
filename_DE_original <- paste(file_prefix, "DE_original.rds", sep = "")
filename_DE_description <- paste(file_prefix, "DE_description.txt", sep = "")
filename_DE_original_description <- paste(file_prefix, "DE_original_description.txt", sep = "")
filename_foldchange_pvalues <- paste(file_prefix, "foldchange_pvalues.csv", sep = "")
filename_formatted_peak_table <- paste(file_prefix, "formatted_peak_table.csv", sep = "")
filename_formatted_sample_data_table <- paste(file_prefix, "formatted_sample_data_table.csv", sep = "")
filename_formatted_sample_metadata <- paste(file_prefix, "formatted_sample_metadata.tsv", sep = "")
filename_formatted_variable_metadata <- paste(file_prefix, "formatted_variable_metadata.csv", sep = "")
filename_graphml <- paste(file_prefix, "graphml.graphml", sep = "")
filename_heatmap_pval <- paste(file_prefix, "Heatmap_pval.html", sep = "")
filename_heatmap_rf <- paste(file_prefix, "Heatmap_rf.html", sep = "")
filename_interactive_table <- paste(file_prefix, "interactive_table.html", sep = "")
filename_metaboverse_table <- paste(file_prefix, "metaboverse_table.tsv", sep = "")
filename_params <- paste(file_prefix, "params.yaml", sep = "")
filename_params_user <- paste(file_prefix, "params_user.yaml", sep = "")
filename_PCA <- paste(file_prefix, "PCA.pdf", sep = "")
filename_PCA3D <- paste(file_prefix, "PCA3D.html", sep = "")
filename_PCoA <- paste(file_prefix, "PCoA.pdf", sep = "")
filename_PCoA3D <- paste(file_prefix, "PCoA3D.html", sep = "")
filename_PLSDA <- paste(file_prefix, "PLSDA.pdf", sep = "")
filename_PLSDA_loadings <- paste(file_prefix, "PLSDA_loadings.tsv", sep = "")
filename_PLSDA_VIP_plot <- paste(file_prefix, "PLSDA_VIP.pdf", sep = "")
filename_PLSDA_VIP_table <- paste(file_prefix, "PLSDA_VIP.tsv", sep = "")
filename_R_script <- paste(file_prefix, "R_script_backup.R", sep = "")
filename_random_forest <- paste(file_prefix, "RF_importance.html", sep = "")
filename_random_forest_model <- paste(file_prefix, "RF_model.txt", sep = "")
filename_session_info <- paste(file_prefix, "session_info.txt", sep = "")
filename_summary_stats_table_full <- paste(file_prefix, "summary_stats_table_full.csv", sep = "")
filename_summary_stats_table_selected <- paste(file_prefix, "summary_stats_table_selected.csv", sep = "")
filename_summary_stat_output_selected_cytoscape <- paste(file_prefix, "summary_stats_table_selected_cytoscape.csv", sep = "")
filename_treemap <- paste(file_prefix, "Treemap_interactive.html", sep = "")
filename_volcano <- paste(file_prefix, "Volcano.pdf", sep = "")
filename_volcano_interactive <- paste(file_prefix, "Volcano_interactive.html", sep = "")

# Sirius filenames

sirius_annotations_filename = "compound_identifications.tsv"
canopus_annotations_filename = "canopus_compound_summary.tsv"

################################### load peak table ########################################
############################################################################################

# If params$actions$run_with_gap_filled is set to TRUE, we load the gap filled peak table

if (params$actions$run_with_gap_filled == "TRUE") {
  feature_table <- read_delim(file.path(working_directory, "results", "mzmine", paste0(params$mapp_batch, "_gf_quant.csv")),
    delim = ",", escape_double = FALSE,
    trim_ws = TRUE
  )
} else {
  feature_table <- read_delim(file.path(working_directory, "results", "mzmine", paste0(params$mapp_batch, "_quant.csv")),
    delim = ",", escape_double = FALSE,
    trim_ws = TRUE
  )
}

# The column names are modified using the rename function from the dplyr package

feature_table <- feature_table %>%
  rename(
    "feature_id" = "row ID",
    "feature_mz" = "row m/z",
    "feature_rt" = "row retention time"
  )

# The row m/z and row retention time columns are concatenated to create a new column called `feature_id_full`
feature_table$"feature_id_full" <- paste(feature_table$feature_id,
  round(feature_table$feature_mz, digits = 2),
  round(feature_table$feature_rt, digits = 1),
  sep = "_"
)

# The dataframe is subsetted to keep only columns containing the pattern ` Peak area` and the `feature_id_full` column
# We use dplyr's `select` function and the pipe operator `%>%` to chain the operations.
# We then remove the ` Peak area` pattern from the column names using the rename_with function from the dplyr package
# We then set the `feature_id_full` column as the rownames of the dataframe and transpose it

feature_table_intensities <- feature_table %>%
  select(feature_id, contains(" Peak area")) %>%
  rename_with(~ gsub(" Peak area", "", .x)) %>%
  column_to_rownames(var = "feature_id") %>%
  as.data.frame() %>%
  t()

# We keep the feature_table_intensities dataframe in a separate variable

X <- feature_table_intensities


# We order the X by rownames and by column names

X <- X[order(row.names(X)), ]
X <- X[, order(colnames(X))]

X <- as.data.frame(X)

# Min value imputation (to be checked !!!)

half_min <- min(X[X > 0], na.rm = TRUE) / 2
min <- min(X[X > 0], na.rm = TRUE)


X[X == 0] <- min

# min(X)

# Uncomment for testing purposes
# X <- X[,1:100]


# We keep the feature metadata in a separate dataframe

feature_metadata <- feature_table %>%
  select(feature_id_full, feature_id, feature_mz, feature_rt)

############################### load annotation tables #####################################
############################################################################################

# The Sirius data is loaded
# First we check if a chebied version exists, if not we create it

if (file.exists(file.path(working_directory, "results", "sirius", paste("chebied", sirius_annotations_filename, sep = "_")))) {
  data_sirius <- read_delim(file.path(working_directory, "results", "sirius", paste("chebied", sirius_annotations_filename, sep = "_")),
    delim = "\t", escape_double = FALSE,
    trim_ws = TRUE
  )
} else {
  data_sirius <- read_delim(file.path(working_directory, "results", "sirius", sirius_annotations_filename),
    delim = "\t", escape_double = FALSE,
    trim_ws = TRUE
  )

  # Here we add this step to "standardize" the sirius names to more classical names
  # We first remove duplicates form the Sirius smiles columns

  for_chembiid_smiles <- unique(data_sirius$smiles)

  # We then use the get_chebiid function from the chembiid package to get the ChEBI IDs

  print("Getting ChEBI IDs from smiles ...")

  chebi_ids <- get_chebiid(for_chembiid_smiles, from = "smiles", to = "chebiid", match = "best")

  # And we merge the data_sirius dataframe with the chebi_ids dataframe
  data_sirius <- merge(data_sirius, chebi_ids, by.x = "smiles", by.y = "query")

  # The column names are modified to include the source of the data

  colnames(data_sirius) <- paste("sirius", colnames(data_sirius),  sep = "_")


  # We now build a unique feature_id for each feature in the Sirius data

  data_sirius$feature_id <- sub("^.*_([[:alnum:]]+)$", "\\1", data_sirius$sirius_id)
  data_sirius$feature_id <- as.numeric(data_sirius$feature_id)

  # Since this step takes time we save the output locally

  write.table(data_sirius, file = file.path(working_directory, "results", "sirius", paste("chebied", sirius_annotations_filename, sep = "_")), sep = "\t", row.names = FALSE)
}

# The CANOPUS data is loaded

data_canopus <- read_delim(file.path(working_directory, "results", "sirius", canopus_annotations_filename),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)


# The column names are modified to include the source of the data

colnames(data_canopus) <- paste("canopus", colnames(data_canopus),  sep = "_")

# We now build a unique feature_id for each feature in the Sirius data

data_canopus$feature_id <- sub("^.*_([[:alnum:]]+)$", "\\1", data_canopus$canopus_id)
data_canopus$feature_id <- as.numeric(data_canopus$feature_id)


write.table(data_canopus, file = file.path(working_directory, "results", "sirius", paste("featured", canopus_annotations_filename, sep = "_")), sep = "\t", row.names = FALSE)

# The MetAnnot data is loaded

data_met_annot <- read_delim(file.path(working_directory, "results", "met_annot_enhancer", params$met_annot_enhancer_folder, paste0(params$met_annot_enhancer_folder, "_spectral_match_results_repond.tsv")),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

# The column names are modified to include the source of the data

colnames(data_met_annot) <- paste("met_annot", colnames(data_met_annot),  sep = "_")


# We now build a unique feature_id for each feature in the Metannot data

data_met_annot$feature_id <- data_met_annot$met_annot_feature_id
data_met_annot$feature_id <- as.numeric(data_met_annot$feature_id)



# The GNPS data is loaded. Note that we use the `Sys.glob` function to get the path to the file and expand the wildcard

# At this point we check wether we have to deal with a GNPS2 job or a job from the GNPS legacy interface
# For this we check for the presence of a directory named `nf_output` in file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id)
# If their is such directory then we set the variable gnps2_job to TRUE, else it is FALSE
# Check if the directory exists
gnps2_job <- file.exists(Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id)))

# Print the result
if (gnps2_job) {
  print("This is a GNPS2 job.")
} else {
  print("This is a job from the GNPS legacy interface.")
}


if (gnps2_job) {
  data_gnps_mn <- read_delim(Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "nf_output", "networking", "clustersummary_with_network.tsv")),
    delim = "\t", escape_double = FALSE,
    trim_ws = TRUE
  )
  # The GNPS `Compound_Name` is dropped
  data_gnps_mn <- data_gnps_mn %>%
    select(-contains("Compound_Name"))
  data_gnps_lib <- read_delim(Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "nf_output", "library", "merged_results_with_gnps.tsv")),
    delim = "\t", escape_double = FALSE,
    trim_ws = TRUE
  )
  # Both df are mergeq using the `#Scan#` column
  data_gnps <- merge(data_gnps_mn, data_gnps_lib, by = "#Scan#")
} else {
  data_gnps <- read_delim(Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "clusterinfo_summary", "*.tsv")),
    delim = "\t", escape_double = FALSE,
    trim_ws = TRUE
  )
}


# The column names are modified to include the source of the data

colnames(data_gnps) <- paste("gnps", colnames(data_gnps),  sep = "_")

# We now build a unique feature_id for each feature in the GNPS data

data_gnps$feature_id <- data_gnps$`gnps_cluster index`
data_gnps$feature_id <- as.numeric(data_gnps$feature_id)


# First we check if a chebied version exists, if not we create it

# if (file.exists(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "chebied_DB_result.tsv"))) {
#   data_gnps_annotations = read_delim(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "DB_result", "chebied_DB_result.tsv"),
#   delim = "\t", escape_double = FALSE,
#   trim_ws = TRUE
# )
# } else {
#   data_gnps_annotations = read_delim(Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "DB_result", "*.tsv")),
#     delim = "\t", escape_double = FALSE,
#     trim_ws = TRUE
#   )

#   # Here we add this step to "standardize" the sirius names to more classical names
#   # We first remove duplicates form the Sirius smiles columns

#   for_chembiid_smiles <- unique(data_gnps_annotations$Smiles)

#   # We then use the get_chebiid function from the chembiid package to get the ChEBI IDs

#   print("Getting ChEBI IDs from smiles ...")

#   chebi_ids <- get_chebiid(for_chembiid_smiles, from = "smiles", to = "chebiid", match = "best")
#   str(chebi_ids)
#   # And we merge the data_sirius dataframe with the chebi_ids dataframe
#   data_gnps_annotations <- merge(data_gnps_annotations, chebi_ids, by.x = "Smiles", by.y = "query")

#   # The column names are modified to include the source of the data

#   colnames(data_gnps_annotations) = paste(colnames(data_gnps_annotations), "dbresult_gnps", sep = "_")

#   # We now build a unique feature_id for each feature in the GNPS data

#   data_gnps_annotations$feature_id = data_gnps_annotations$`#Scan#_dbresult_gnps`
#   data_gnps_annotations$feature_id = as.numeric(data_gnps_annotations$feature_id)
#   str(data_gnps_annotations)

#   write.table(data_gnps_annotations, file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "DB_result", "chebied_DB_result.tsv"),
#   sep = "\t", row.names = FALSE)
# }


# The four previous dataframe are merged into one using the common `feature_id` column as key and the tidyverse `reduce` function

list_df <- list(feature_metadata, data_sirius, data_canopus, data_met_annot, data_gnps)
VM <- list_df %>% reduce(full_join, by = "feature_id")

# We take care to convert all N/A values to NA

VM[VM == "N/A"] <- NA

# We add a sanitizing function. first we lowercase all colnames

colnames(VM) <- tolower(colnames(VM))

# We then take care of the # chracter and change it to _

colnames(VM) <- gsub("#", "_", colnames(VM))


VM <- VM %>%
  clean_names(case = "snake")


# The row m/z and row retention time columns are concatenated to create a new column called `feature_id_full_annotated`
VM$"feature_id_full_annotated" <- paste0(
  VM$sirius_chebiasciiname,
  "_[",
  VM$feature_id_full,
  "]",
  sep = ""
)


# Make sure that all column containing score in their name are as.numeric. But we keep all the dataframes columns (we might want to find a more generic way to do this)

VM <- VM %>%
  mutate_at(vars(contains("score")), as.numeric)


# We now convert the VM tibble into a dataframe and set the `feature_id` column as the rownames

VM <- as.data.frame(VM)
row.names(VM) <- VM$feature_id


# We order the VM by rownames and by column names

VM <- VM[order(row.names(VM)), ]
VM <- VM[, order(colnames(VM))]

# Uncomment for testing purposes
# VM = head(VM, 100)


################################ load sample  metadata #####################################
############################################################################################

# Later on ... implement a test stage where we check for the presence of a "species" and "sample_type" column in the metadata file.


# We here load the sample metadata


sample_metadata <- read_delim(file.path(working_directory, "metadata", "treated", paste(params$mapp_batch, "metadata.txt", sep = "_")),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Here we establish a small test which will check if the sample metadata file contains the required columns (filename, sample_id, sample_type and species)

required_columns <- c("filename", "sample_id", "sample_type", "source_taxon")

if (!all(required_columns %in% colnames(sample_metadata))) {
  stop("The sample metadata file does not contain the required columns (filename, sample_id, sample_type and source_taxon). Please check your metadata file and try again.")
}


SM <- data.frame(sample_metadata)



# Sanitize SM colnames
SM <- SM %>%
  clean_names(case = "snake")



# Here we fetch the wikidata QIDs for the source_taxon columns

# Get distinct taxon names (including multiple taxa in a single entry)
distinct_taxa <- SM %>%
  filter(sample_type == "sample") %>%
  mutate(source_taxon = strsplit(source_taxon, ", ")) %>% # Split multiple taxa
  unnest(source_taxon) %>% # Expand multiple taxa into separate rows
  distinct(source_taxon)


taxon_names <- distinct_taxa$source_taxon



# Function to query Wikidata for QIDs based on taxon names
get_taxon_qids <- function(taxon_names) {
  qids <- character(length(taxon_names))
  # i= 1
  for (i in seq_along(taxon_names)) {
    taxon_name <- taxon_names[i]
    query <- paste0('SELECT ?taxon WHERE { ?taxon wdt:P225 "', taxon_name, '". }')
    result <- WikidataQueryServiceR::query_wikidata(query)

    if (!is.null(result$taxon) && length(result$taxon) > 0) {
      qid_full <- result$taxon[1]
      qid_plain <- sub("http://www.wikidata.org/entity/", "", qid_full)
      qids[i] <- qid_plain
    } else {
      qids[i] <- NA
    }
  }

  return(qids)
}


# Get QIDs for distinct taxon names
distinct_qids <- get_taxon_qids(distinct_taxa$source_taxon)

# Combine distinct taxon names with their QIDs
distinct_taxa$source_taxon_qid <- distinct_qids


# Use dplyr to create comma-separated QID column
SM <- SM %>%
  mutate(source_taxon = strsplit(source_taxon, ", ")) %>%
  rowwise() %>%
  mutate(source_taxon_qid = paste(distinct_taxa$source_taxon_qid[distinct_taxa$source_taxon %in% source_taxon], collapse = ", ")) %>%
  ungroup() %>%
  mutate(source_taxon = sapply(source_taxon, paste, collapse = ", ")) %>%
  as.data.frame()



# We take full power over the matrix (sic. Defossez, 2023)
# First we work vertically (within a given SM column)

for (column in names(params$to_combine_vertically)) {
  col_info <- params$to_combine_vertically[[column]]
  col_name <- col_info$factor_name

  # Initialize aggregated groups with original condition variable
  SM[paste(col_info$factor_name, "simplified", sep = "_")] <- as.factor(SM[[col_info$factor_name]])

  # Iterate over each group in params$tocomb
  for (group in names(col_info$groups)) {
    group_info <- col_info$groups[[group]]
    levels <- group_info$levels

    # now make sure to sort the levels
    levels <- sort(levels, decreasing = FALSE)
    # We create a new label for the current group by concatenating the levels value with a "_"
    new_label <- paste(levels, collapse = "_")

    # Combine levels for the current group
    SM[paste(col_info$factor_name, "simplified", sep = "_")] <- combineLevels(SM[[paste(col_info$factor_name, "simplified", sep = "_")]], levs = levels, newLabel = c(new_label))
  }
}

# The function below is used to create metadata combinations
# Then we work horizontally (across SM columns)

if (!is.null(params$to_combine_horizontally$factor_name)) {
  df <- SM %>%
    filter(sample_type == "sample")

  # This line allows us to make sure that the columns will be combined in alphabetical order
  cols <- sort(c(params$to_combine_horizontally$factor_name), decreasing = FALSE)

  for (n in 1:length(cols)) {
    combos <- combn(cols, n, simplify = FALSE)
    for (combo in combos) {
      new_col_name <- paste(combo, collapse = "_")
      df[new_col_name] <- apply(df[combo], 1, paste, collapse = "_")
    }
  }

  # We merge back the resulting df to the original SM dataframe and fill the NA values with "ND"
  SM <- merge(x = SM, y = df, all.x = TRUE)
}

SM[is.na(SM)] <- "ND"

SM <- SM %>%
  remove_rownames() %>%
  column_to_rownames(var = "filename")


# This allows us to both have the filename as rownames and as a unique column
SM$filename <- rownames(SM)
# cuirmoustache

# We order the SM by rownames and by column names.

SM <- SM[order(row.names(SM)), ]
SM <- SM[, order(colnames(SM))]


# glimpse(SM)
# Optional ponderation step.
# To clean

# XSM = merge(x = X, y = SMDF, by = "row.names", all = TRUE)

# X_pond = XSM %>%
#   select(Row.names, grep("peak", colnames(XSM))) %>%
#   column_to_rownames(var = "Row.names")

# X_pond = X_pond[order(row.names(X_pond)), ]
# X_pond = X_pond[, order(colnames(X_pond))]
# SMDF = SMDF[order(row.names(SMDF)), ]
# VM = VM[order(row.names(VM)), ]


if (any(colnames(X) != row.names(VM))) {
  stop("Some columns in X are not present in the rownames of VM. Please check the column names in X and the rownames of VM.")
}

# We repeat for row.names(SMDF) == row.names(X_pond)

if (any(row.names(X) != row.names(SM))) {
  stop("Some rownames in X are not present in the rownames of SM. Please check the rownames in X and the rownames of SM.")
}

# length(unique(row.names(X)))
# length(unique(row.names(SM)))

# # We troubleshoot and find which rownames are not present in both X and SM

# rownames_not_present = setdiff(row.names(SM), row.names(X))


#################################################################################################
#################################################################################################
#################################################################################################

# The DatasetExperiment object is created using the X_pond, SMDF and VM objects.

DE_original <- DatasetExperiment(
  data = X,
  sample_meta = SM,
  variable_meta = VM,
  name = params$dataset_experiment$name,
  description = params$dataset_experiment$description
)


# variable_meta = DE_original$variable_meta

## Filtering steps

if (length(params$feature_to_filter) > 0) {
  filter_by_name_model <- filter_by_name(mode = "exclude", dimension = "variable", names = params$feature_to_filter)

  # apply model sequence
  filter_by_name_result <- model_apply(filter_by_name_model, DE_original)
  DE_filtered_name <- filter_by_name_result@filtered
} else {
  DE_filtered_name <- DE_original
}


## Filtering steps

## Blank filter

if (is.numeric(params$filter_blank$fold_change) == TRUE) {

  filter_blank_model <- blank_filter(
    fold_change = params$filter_blank$fold_change,
    factor_name = params$filter_blank$factor_name,
    blank_label = params$filter_blank$blank_label,
    qc_label = params$filter_blank$qc_label,
    fraction_in_blank = params$filter_blank$fraction_in_blank
  )

  # apply model sequence
  filter_blank_result <- model_apply(filter_blank_model, DE_filtered_name)

  DE_filtered <- filter_blank_result$filtered
} else {
  DE_filtered <- DE_filtered_name
}


if (params$filter_sample_type$mode %in% possible_modes){
  filter_smeta_model <- filter_smeta(
    mode = params$filter_sample_type$mode,
    factor_name = params$filter_sample_type$factor_name,
    levels = params$filter_sample_type$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
} 

if (params$filter_sample_metadata_one$mode %in% possible_modes){
  filter_smeta_model <- filter_smeta(
    mode = params$filter_sample_metadata_one$mode,
    factor_name = params$filter_sample_metadata_one$factor_name,
    levels = params$filter_sample_metadata_one$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
}

if (params$filter_sample_metadata_two$mode %in% possible_modes) {
  filter_smeta_model <- filter_smeta(
    mode = params$filter_sample_metadata_two$mode,
    factor_name = params$filter_sample_metadata_two$factor_name,
    levels = params$filter_sample_metadata_two$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
}


if (params$filter_variable_metadata_one$mode %in% possible_modes) {
  filter_vmeta_model <- filter_vmeta(
    mode = params$filter_variable_metadata_one$mode,
    factor_name = params$filter_variable_metadata_one$factor_name,
    levels = params$filter_variable_metadata_one$levels
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}

if (params$filter_variable_metadata_two$mode %in% possible_modes) {
  filter_vmeta_model <- filter_vmeta(
    mode = params$filter_variable_metadata_two$mode,
    factor_name = params$filter_variable_metadata_two$factor_name,
    levels = params$filter_variable_metadata_two$levels
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}



if (params$filter_variable_metadata_annotated$mode %in% possible_modes) {
  # Convert the "levels" value to NA if it is "NA" as a character string
  if (params$filter_variable_metadata_annotated$levels == "NA") {
    params$filter_variable_metadata_annotated$levels <- NA
  }

  filter_vmeta_model <- filter_vmeta(
    mode = params$filter_variable_metadata_annotated$mode,
    factor_name = params$filter_variable_metadata_annotated$factor_name,
    levels = as.character(params$filter_variable_metadata_annotated$levels)
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}


if (params$filter_variable_metadata_num$mode %in% possible_modes) {
  filter_vmeta_model <- filter_vmeta_num(
    mode = params$filter_variable_metadata_num$mode,
    factor_name = params$filter_variable_metadata_num$factor_name,
    level = params$filter_variable_metadata_num$level
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}


if (params$actions$scale_data == "FALSE") {
  DE <- DE_filtered

} else if (params$actions$scale_data == "TRUE") {
  # Overall Pareto scaling (test)

  M <- pareto_scale()
  M <- model_train(M, DE_filtered)
  M <- model_predict(M, DE_filtered)
  DE <- M$scaled

  # We use the filter_na_count function to filter out features with a number of NAs greater than the threshold

  M <- filter_na_count(threshold = 1, factor_name = "sample_type")
  M <- model_apply(M, DE)

  DE <- M$filtered

  ##### we range all feature from 0 to 1

  # @Manu !!! Why do we have this ?!
  DE$data <- apply(DE$data, 2, modEvA::range01)


  # Min value imputation also after the scaling stage (to be checked !!!)

  # half_min_sec = min(DE$data[DE$data > 0], na.rm = TRUE) / 2

  # min_sec = min(DE$data[DE$data > 0], na.rm = TRUE)

  # DE$data[DE$data == 0] = min_sec


} else {
stop("Please check the value of the 'scale_data' parameter in the params file.")
}

################################################################################################
################################################################################################
######################## structool box formatted data export

message("Outputting X, VM and SM ...")

formatted_peak_table <- DE$data

formatted_variable_metadata <- DE$variable_meta ### need to be filter with only usefull output


# col_filter <- c("feature_id_full", "feature_id", "feature_mz" ,"feature_rt", "sirius_molecularformula","sirius_inchikey2d","InChI_sirius",
# "sirius_name","sirius_smiles", "pubchemids_sirius", "molecularFormula_canopus", "canopus_npc_pathway","NPC.pathway.Probability_canopus",
# "canopus_npc_superclass", "canopus_npc_class","ClassyFire.most.specific.class_canopus","ClassyFire.most.specific.class.Probability_canopus",
# "ClassyFire.level.5_canopus","ClassyFire.subclass_canopus","ClassyFire.class_canopus","ClassyFire.superclass_canopus","ClassyFire.all.classifications_canopus",
# "met_annot_structure_wikidata","met_annot_structure_inchikey","met_annot_structure_inchi","met_annot_structure_smiles","met_annot_structure_molecular_formula",
# "met_annot_short_inchikey","met_annot_structure_taxonomy_npclassifier_01pathway","met_annot_structure_taxonomy_npclassifier_02superclass",
# "met_annot_structure_taxonomy_npclassifier_02superclass","met_annot_organism_wikidata","met_annot_organism_name","met_annot_organism_taxonomy_ottid","met_annot_organism_taxonomy_01domain",
# "met_annot_organism_taxonomy_02kingdom","met_annot_organism_taxonomy_03phylum","met_annot_organism_taxonomy_04class","met_annot_organism_taxonomy_05order",
# "met_annot_organism_taxonomy_06family","met_annot_organism_taxonomy_07tribe","met_annot_organism_taxonomy_08genus","met_annot_organism_taxonomy_09species","met_annot_organism_taxonomy_10varietas",
# "matched_domain_met_annot","matched_kingdom_met_annot","matched_phylum_met_annot","matched_class_met_annot","matched_order_met_annot","matched_family_met_annot","matched_tribe_met_annot",
# "matched_genus_met_annot","matched_species_met_annot","met_annot_score_taxo","score_max_consistency_met_annot","final_score_met_annot","rank_final_met_annot","component_id_met_annot",
# "met_annot_structure_taxonomy_npclassifier_01pathway_consensus","freq_met_annot_structure_taxonomy_npclassifier_01pathway","met_annot_structure_taxonomy_npclassifier_02superclass_consensus",
# "freq_met_annot_structure_taxonomy_npclassifier_02superclass","met_annot_structure_taxonomy_npclassifier_03class_consensus","freq_met_annot_structure_taxonomy_npclassifier_02superclass")


col_filter <- c("feature_id_full", "feature_id", "feature_mz", "feature_rt", "sirius_molecularformula")

formatted_variable_metadata_filtered <- formatted_variable_metadata[col_filter]

formatted_sample_metadata <- DE$sample_meta

formatted_sample_data_table <- merge(DE$sample_meta, DE$data, by = "row.names")


# We work on an export for MetaboAnalyst Pathways analysis

# DE$sample_meta


# DE$data

# # First we merge the sample metadata and the data

# sample_data_table = merge(DE$sample_meta, DE$data, by="row.names")

# # We now drop the useless columns. We use the dplyr syntax

# colnames_to_drop = c("Row.names","condition_detailed","condition_simplified","filename","id","internal_id","sample_type","source_taxon","source_taxon_qid")

# sample_data_table = sample_data_table %>%
#   select(-one_of(colnames_to_drop))  %>%
#   # reorganize the columns
#   select(sample_id, condition, everything())

# # We now filter the variable metadata to keep only the columns and rows we need

# colnames(DE$variable_meta)

# DE$variable_meta$sirius_chebiasciiname

# compound_names = DE$variable_meta  %>%
# select(sirius_chebiasciiname)  %>%
# # NA are dropped by default
# filter(!is.na(sirius_chebiasciiname))

# # We now replace the column names in the sample_data_table with the compound names
# # For this we transpose the sample_data_table and then match

# sample_data_table_transposed = as.data.frame(t(sample_data_table))

# # We now merge the compound names with the sample_data_table_transposed

# merged = merge(sample_data_table_transposed, compound_names, by = "row.names", all = TRUE)

# as.data.frame(merged)

# # We now transpose the merged dataframe defining the rownames column as the first row


# sample_compounds = as.data.frame(t(merged))

# rownames(sample_compounds)

# # Row "sirius_chebiasciiname" is now the first row. We now rename it to "compound_name"



###################################################################################################
######################### rename main folder - short version
# Here we check if the params$paths$out value exist and use it else we use the default output_directory

target_name = paste(as.vector(sort(as.character(unique(DE$sample_meta[[params$target$sample_metadata_header]])), decreasing = FALSE)), collapse = "_vs_")


common_df_path <- file.path(params$paths$output, "params_log.rds")
common_tsv_path <- file.path(params$paths$output, "params_log.tsv")


  # Convert the YAML content to a dataframe row
new_row_df <- convert_yaml_to_single_row_df_with_hash(params)

# Append this row to the common dataframe and save
append_to_common_df_and_save(new_row_df, common_df_path, common_tsv_path)



if (params$paths$output != "") {
  output_directory <- file.path(params$paths$output, new_row_df$hash)
} else {
  output_directory <- file.path(working_directory, "results", "stats", new_row_df$hash)
}


if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
  message("Directory created:", output_directory, "\n")
} else {
  message("Directory already exists:", output_directory, "\n")
}

#################################################################################
#################################################################################
##### write raw data and param

## We save the used params.yaml

message("Writing params.yaml ...")

file.copy(path_to_params, file.path(output_directory, filename_params), overwrite = TRUE)
file.copy(path_to_params_user, file.path(output_directory, filename_params_user), overwrite = TRUE)

# We move to the output directory

setwd(output_directory)

#######################

# The DE and DE_original objects is printed and saved in the output directory

message("Saving DE object ...")

saveRDS(DE, filename_DE)

message("DatasetExperiment object properties: ")

sink(filename_DE_description)

print(DE)

sink()

message("Saving DE_original object ...")

saveRDS(DE_original, filename_DE_original)

message("DatasetExperiment object properties: ")

sink(filename_DE_original_description)

print(DE_original)

sink()


write.table(formatted_peak_table, file = filename_formatted_peak_table, sep = ",", row.names = FALSE)
write.table(formatted_variable_metadata_filtered, file = filename_formatted_variable_metadata, sep = ",", row.names = FALSE)
write.table(formatted_sample_metadata, file = filename_formatted_sample_metadata, sep = "\t", row.names = FALSE)
write.table(formatted_sample_data_table, file = filename_formatted_sample_data_table, sep = ",", row.names = FALSE)
################################################################################################
################################################################################################

# The Figures titles are conditionally defined according to the user's choices and option in the parameters file

# title_PCA = paste("PCA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
# title_PLSDA = paste("PLSDA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
# title_PLSDA_VIP = paste("PLSDA selected Features of Importance", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
# title_PCA3D = paste("PCA3D", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
# title_PCoA = paste("PCoA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
# title_PCoA3D = paste("PCoA3D", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
# title_volcano = paste("Volcano plot", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
# title_treemap = paste("Treemap", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
# title_random_forest = paste("Random Forest results", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
# title_box_plots = paste("Top", params$boxplot$topN, "boxplots", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
# title_heatmap = paste("Heatmap of","top", params$heatmap$topN,"Random Forest filtered features", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")

# title_PCA = paste("PCA", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")

title_PCA <- paste(
  paste("PCA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)


# title_PLSDA = paste("PLSDA", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PLSDA <- paste(
  paste("PLSDA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

title_PLSDA_VIP <- paste("PLSDA selected Features of Importance", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")

# title_PCA3D = paste("PCA3D", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PCA3D <- paste(
  paste("PCA3D", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

# title_PCoA = paste("PCoA", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PCoA <- paste(
  paste("PCoA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

# title_PCoA3D = paste("PCoA3D", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PCoA3D <- paste(
  paste("PCoA3D", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

# title_volcano = paste("Volcano plot", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")
title_treemap <- paste("Treemap", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")
title_random_forest <- paste("Random Forest results", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")
title_box_plots <- paste("Top", params$boxplot$topN, "boxplots", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")
title_heatmap_rf <- paste("Heatmap of", "top", params$heatmap$topN, "Random Forest filtered features", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")

# title_heatmap_pval = paste("Heatmap of significant feature for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_heatmap_pval <- paste(
  paste("Heatmap of significant feature for dataset", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "<br>"
)

title_volcano <- paste(
  paste("Volcano plot of significant feature for dataset", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

#################################################################################################
#################################################################################################
############# Colors definition #################################################################
#################################################################################################
#################################################################################################
wes_palettes_vec <- sample(sort(unique(unlist(wes_palettes[names(wes_palettes)]))))

factor_name_meta <- unlist(unique(DE$sample_meta[params$target$sample_metadata_header]))


# Check that alll members of params$colors$all$key are present in factor_name_meta.
# If not, return the values of params$colors$all$key that are not present in factor_name_meta.

if (!all(params$colors$all$key %in% factor_name_meta)) {
  missing_colors <- params$colors$all$key[!params$colors$all$key %in% factor_name_meta]
  factor_name_meta_str <- paste(unique(factor_name_meta), collapse=", ")
  stop(paste("The following values in params$colors$all$key are not present in the sample metadata:", paste(missing_colors, collapse=", "), "Check the spelling of values in params$colors$all$key, they should match the following available values:", factor_name_meta_str))
}

# We establish a named vector for the whole dataset
if (length(params$colors$all$key) > 0) {
  custom_colors <- setNames(c(params$colors$all$value), c(params$colors$all$key))
} else {
  custom_colors <- wes_palettes_vec[sample(c(1:length(wes_palettes_vec)), length(factor_name_meta))]
  names(custom_colors) <- factor_name_meta
}





#################################################################################################
#################################################################################################
#################################################################################################
##### PCA filtered data #######################################################################

message("Launching PCA calculations ...")


pca_seq_model <- filter_na_count(threshold = 1, factor_name = "sample_type") +
  knn_impute(neighbours = 5) +
  vec_norm() +
  # log_transform(base = 10) +
  mean_centre() +
  PCA(number_components = 3)

# apply model sequence
pca_seq_result <- model_apply(pca_seq_model, DE)

# Fetching the PCA data object
pca_object <- pca_seq_result[length(pca_seq_result)]

# PCA scores plot

pca_scores_plot <- pca_scores_plot(
  factor_name = params$target$sample_metadata_header,
  label_factor = "sample_id",
  ellipse_type = "t",
  ellipse_confidence = 0.9,
  points_to_label = "all"
)

# plot
pca_plot <- chart_plot(pca_scores_plot, pca_object)



fig_PCA <- pca_plot +
  theme_classic() +
  facet_wrap(~ pca_plot$labels$title) +
  ggtitle(title_PCA)


fig_PCA <- fig_PCA +
  scale_colour_manual(name = "Groups", values = custom_colors)



#   theme(plot.title = element_text(hjust = 0.2, vjust = -2)) +


# We merge PCA scores and metadata info in a single df

PCA_meta <- merge(x = pca_object$scores$sample_meta, y = pca_object$scores$data, by = 0, all = TRUE)


fig_PCA3D <- plot_ly(PCA_meta, x = ~PC1, y = ~PC2, z = ~PC3, color = PCA_meta[, params$target$sample_metadata_header], colors = custom_colors)


fig_PCA3D <- fig_PCA3D %>% add_markers()
fig_PCA3D <- fig_PCA3D %>% layout(
  scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")
  ),
  legend = list(title = list(text = params$target$sample_metadata_header)),
  title = title_PCA3D
)


# The files are exported

ggsave(plot = fig_PCA, filename = filename_PCA, width = 10, height = 10)


if (params$operating_system$system == "unix") {
  ### linux version
  fig_PCA3D %>%
    htmlwidgets::saveWidget(file = filename_PCA3D, selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
  ### windows version
  Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
  fig_PCA3D %>%
    htmlwidgets::saveWidget(file = filename_PCA3D, selfcontained = TRUE, libdir = "lib")
  unlink("lib", recursive = FALSE)
}

# #################################################################################################
# #################################################################################################
# #################################################################################################
# ##### PLSDA filtered data #######################################################################


if (params$actions$run_PLSDA == "TRUE") {
  message("Launching PLSDA calculations ...")

  # First we make sure that the sample metadata variable of interest is a factor
  # For now we use DE_original here ... check if this is correct

  DE$sample_meta[, params$target$sample_metadata_header] <- as.factor(DE$sample_meta[, params$target$sample_metadata_header])

  # glimpse(DE_filtered$sample_meta)

  # check the outcome of a Pareto scaling methods


  # # prepare model sequence
  plsda_seq_model <- # autoscale() +
    filter_na_count(threshold = 1, factor_name = params$target$sample_metadata_header) +
    # knn_impute() +
    PLSDA(factor_name = params$target$sample_metadata_header, number_components = 2)

  plsda_seq_result <- model_apply(plsda_seq_model, DE)



  # Fetching the PLSDA data object
  plsda_object <- plsda_seq_result[length(plsda_seq_result)]

  # We merge the plsda_object$vip object with the DE$variable_meta object. We use dplyr syntax and keep a new object called variable_meta_plsda. We keep the rownames of the plsda_object$vip.

  plsda_object_vip <- plsda_object$vip %>%
    rownames_to_column(var = "feature_id") %>%
    mutate(feature_id = as.numeric(feature_id))

  vip_variable_meta <- plsda_object_vip %>%
    left_join(DE$variable_meta, by = "feature_id") %>%
    select(feature_id, feature_id_full_annotated)

  rownames(plsda_object$vip) <- vip_variable_meta$feature_id_full


  C <- pls_scores_plot(factor_name = params$target$sample_metadata_header)

  plsda_plot <- structToolbox::chart_plot(C, plsda_object)




  fig_PLSDA <- plsda_plot + theme_classic() + facet_wrap(~ plsda_plot$labels$title) + ggtitle(title_PLSDA)


  fig_PLSDA <- fig_PLSDA +
    scale_colour_manual(name = "Groups", values = custom_colors)


  # We output the feature importance

  C <- plsda_feature_importance_plot(n_features = 30, metric = "vip")

  vip_plot <- chart_plot(C, plsda_object)



  fig_PLSDA_VIP <- vip_plot + theme_classic() + facet_wrap(~ plsda_plot$labels$title) + ggtitle(title_PLSDA_VIP)

  # We keep the loadings

  loadings <- plsda_object$loadings

  # The rownames of the loadings are the feature names, we keep them as a column
  # We also keep these as integers

  loadings <- loadings %>%
    rownames_to_column(var = "feature_id") %>%
    mutate(feature_id = as.numeric(feature_id))

  # The process is repeated for the vip


  vip <- plsda_object$vip %>%
    # We keep the first column and rename it to VIP
    select(1) %>%
    rename(VIP = 1) %>%
    rownames_to_column(var = "feature")

  # We now merge this vip object with the loading to fetch the correct feature names
  # We assume that the row are in the same order. We use cbind()

  vip <- cbind(vip, loadings)

  # We reorganize the columns to keep feature_id and feature at the beginning. We use dplyr syntax
  # We order by decreasing value of the VIP column

  vip <- vip %>%
    select(feature_id, feature, everything()) %>%
    arrange(desc(VIP))

  # The plots are exported

  ggsave(plot = fig_PLSDA, filename = filename_PLSDA, width = 10, height = 10)
  ggsave(plot = fig_PLSDA_VIP, filename = filename_PLSDA_VIP_plot, width = 20, height = 10)

  # We export the loadings

  # write.table(loadings, file = filename_PLSDA_loadings, sep = "\t", row.names = FALSE)

  # We export the vip

  write.table(vip, file = filename_PLSDA_VIP_table, sep = "\t", row.names = FALSE)

}

#################################################################################################
#################################################################################################
#################################################################################################
##### HClustering ###############################################################################
#################################################################################################

# # prepare model sequence

# MS_hclust = filter_smeta(mode = "exclude", levels = "QC", factor_name = "sample_type") +
#   filter_smeta(mode = "exclude", levels = "BK", factor_name = "sample_type") +
#   # filter_smeta(mode = var_filter_type, levels= var_levels, factor_name = var_filter) +
#  #log_transform(base = 10) +
#   filter_by_name(mode = "include", dimension = "variable", names = names_var)

# # apply model sequence

# DE_MS_hclust = model_apply(MS_hclust, DE)
# DE_MS_hclust = DE_MS_hclust[length(DE_MS_hclust)]

# ######################################################
# ######################################################

# data_RF = DE_MS_hclust@filtered

# data_subset_norm_rf = data_RF$data
# data_subset_norm_rf[sapply(data_subset_norm_rf, is.infinite)] = NA
# data_subset_norm_rf[is.na(data_subset_norm_rf)] = 0

# dist_metabo = vegdist(data_subset_norm_rf, method = "bray") # method="man" # is a bit better
# hclust_metabo = hclust(dist_metabo, method = "complete")

# dendro_metabo_noqc = as.phylo(hclust_metabo)

# #### remove_QC

# plot(dendro_metabo_noqc)

# ############ circular plot

# var_treatment2 = data_RF$sample_meta[var_treatment]
# var_treatment2 = var_treatment2[, 1]
# fact_plot = as.factor(var_treatment2)
# g2 = split(as.factor(row.names(data_RF$sample_meta)), fact_plot)
# tree_plot2 = ggtree::groupOTU(dendro_metabo_noqc, g2, group_name = "grouper")


# # cols = wes_palette("Cavalcanti1") ## few factor
# cols = distinctColorPalette(length(unique(fact_plot))) ## several factor
# # cols = cols = myColorRamp(c("white", "darkgreen", "black", "black", "black"), sort(as.numeric(as.vector(fact_plot)))) ## continiou svalue

# cols = cols[1:length(unique(fact_plot))]

# circ = ggtree(tree_plot2, aes(color = grouper), size = 2) +
#   geom_tiplab(size = 2, offset = 0) +
#   scale_color_manual(values = cols) + theme(legend.position = "none")

# df = data.frame(fact_plot)
# colnames(df) = c("grouper")

# rownames(df) = as.factor(row.names(data_RF$sample_meta))


# if (!(molecular_pathway_target == "all")) {
#   mol_family = molecular_pathway_target
# }
# if (molecular_pathway_target == "all") {
#   mol_family = ""
# }


# title_hclust = paste("hclust", var_treatment, mol_family, sep = " ")


# ggheat_plot = gheatmap(circ, df[, "grouper", drop = F],
#   offset = 0.05, width = 0.1, colnames = FALSE,
#   colnames_angle = 90, colnames_offset_y = 1
# ) + scale_fill_manual(values = cols) + ggtitle(title_hclust)

# ggheat_plot

# # setwd("G:/My Drive/taf/git_repository/andrea-brenna-group/docs/mapp_project_00021/mapp_batch_00037/results/stats")
# ggsave(paste("tree_plot_", var_treatment, "_", mol_family, ".pdf", sep = ""), width = 15, height = 20)


#################################################################################################
#################################################################################################
#################################################################################################
##### PCoA ##########################################################################


message("Launching PCoA calculations ...")

# # prepare model sequence

# MS_PCOA = filter_smeta(mode = "include", levels = params$filters$to_include, factor_name = "sample_type") +
#  #log_transform(base = 10) +
#   filter_by_name(mode = "include", dimension = "variable", names = names_var)

# # apply model sequence
# # Note that for the PCoA we need to use the original data, not the scaled one

# DE_MS_PCOA = model_apply(MS_PCOA, DE_original)
# DE_MS_PCOA = DE_MS_PCOA[length(DE_MS_PCOA)]

######################################################
######################################################

# @Manu explain what is done below filters etc ....


data_RF <- DE # DE_filtered
sample_name <- data_RF$sample_meta$sample_id #### check
data_subset_norm_rf <- data_RF$data
data_subset_norm_rf[sapply(data_subset_norm_rf, is.infinite)] <- NA
data_subset_norm_rf[is.na(data_subset_norm_rf)] <- 0


dist_metabo <- vegdist(data_subset_norm_rf, method = "bray") # method="man" # is a bit better
# Why dont we use it then ???
D3_data_dist <- cmdscale(dist_metabo, k = 3)
D3_data_dist <- data.frame(D3_data_dist)
D3_data_dist$sample_name <- sample_name
D3_data_dist <- D3_data_dist[order(D3_data_dist$sample_name), ]
metadata_merge <- data_RF$sample_meta[order(sample_name), ]

data_PCOA_merge <- data.frame(cbind(D3_data_dist, metadata_merge))

cols <- data_PCOA_merge[params$target$sample_metadata_header]
cols <- cols[, 1]


fig_PCoA <- ggplot(data_PCOA_merge, aes(x = X1, y = X2, color = cols)) +
  geom_point() +
  ggtitle(title_PCoA) +
  theme_classic()



fig_PCoA <- fig_PCoA +
  scale_colour_manual(name = "Groups", values = custom_colors)



#### PCoA 3D

fig_PCoA3D <- plot_ly(
  x = data_PCOA_merge$X1, y = data_PCOA_merge$X2, z = data_PCOA_merge$X3,
  type = "scatter3d", mode = "markers", color = cols, colors = custom_colors,
  hoverinfo = "text",
  text = ~ paste(
    "</br> name: ", data_PCOA_merge$sample_name,
    "</br> num: ", data_PCOA_merge$sample_id
  )
)



fig_PCoA3D <- fig_PCoA3D %>% layout(
  title = title_PCoA3D,
  legend = list(title = list(text = params$target$sample_metadata_header))
)


# The files are exported

ggsave(plot = fig_PCoA, filename = filename_PCoA, width = 10, height = 10)


if (params$operating_system$system == "unix") {
  ### linux version
  fig_PCoA3D %>%
    htmlwidgets::saveWidget(file = filename_PCoA3D, selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
  ### windows version
  Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
  fig_PCoA3D %>%
    htmlwidgets::saveWidget(file = filename_PCoA3D, selfcontained = TRUE, libdir = "lib")
  unlink("lib", recursive = FALSE)
}

#################################################################################################
#################################################################################################
#################################################################################################
##### group detection with DBSCAN 3D ############################################################
#################################################################################################


# df = data_PCOA_merge[, 1:2]

# # Compute DBSCAN using fpc package
# db = fpc::dbscan(df, eps = 0.01, MinPts = 3)
# # Plot DBSCAN results
# plot(db, df, main = "DBSCAN", frame = FALSE)


# db_cluster = db$cluster
# PC1 = data_PCOA_merge$X1
# PC2 = data_PCOA_merge$X2
# PC3 = data_PCOA_merge$X3
# sample_raw_id = row.names(data_PCOA_merge)
# id = data_PCOA_merge$sample_name
# ATTRIBUTE = data_PCOA_merge$age_genotype ### check
# data_dbclust_merge = data.frame(sample_raw_id, ATTRIBUTE, db_cluster, PC1, PC2, PC3)


# write.table(data_dbclust_merge, file = file.path(working_directory, "results", "stats", paste(params$mapp_batch, "_", "data_dbclust_merge", ".csv", sep = "")), sep = ",")

# ellipse = ellipse3d(cov(data_dbclust_merge[c("PC1", "PC2", "PC3")]))


# title_dbclut = paste("Density Cbased Clustering", var_treatment, mol_family, sep = " ")


# DBscan_3D = plot_ly(
#   x = data_dbclust_merge$PC1, y = data_dbclust_merge$PC2, z = data_dbclust_merge$PC3,
#   type = "scatter3d", mode = "markers", color = db_clut,
#   hoverinfo = "text",
#   text = ~ paste(
#     "</br> condition: ", data_dbclust_merge$ATTRIBUTE,
#     "</br> name: ", data_dbclust_merge$sample_raw_id
#   )
# ) %>%
#   layout(title = title_dbclut)


# # setwd("G:/My Drive/taf/git_repository/andrea-brenna-group/docs/mapp_project_00021/mapp_batch_00037/results/stats")
# # Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/bin/quarto/bin/tools") ### for seflcontained

# DBscan_3D %>%
#   htmlwidgets::saveWidget(file = file.path(working_directory, "results", "stats", paste("DBscan_3D", var_treatment, "_", mol_family, ".html", sep = "")), selfcontained = TRUE)

#################################################################################################
#################################################################################################
#################################################################################################
##### Fold Changes and Tukey’s Honest Significant Difference calculations #######################
#################################################################################################
#################################################################################################


message("Launching Fold Changes and Tukey’s Honest Significant Difference calculations ...")


### Here we wil work on outputting pvalues and fc for time series.

# We build a for loop to iterate over the different time points
# This loop generate a set of DE results for each time point

# params = yaml.load_file('/Users/pma/Dropbox/git_repos/mapp-metabolomics-unit/biostat_toolbox/params/params.yaml')


if (params$actions$calculate_multi_series_fc == "TRUE") {
  l <- list()

  for (i in params$multi_series$points) {
    print(i)
    filter_smeta_model <- filter_smeta(
      mode = "include",
      factor_name = params$multi_series$colname,
      levels = i
    )

    # apply model sequence
    filter_smeta_result <- model_apply(filter_smeta_model, DE)

    DE_tp <- filter_smeta_result@filtered
    # assign(paste("DE_filtered", i, sep = "_"), filter_smeta_result@filtered)

    # The formula is defined externally
    formula <- as.formula(paste0(
      "y", "~", params$target$sample_metadata_header, "+",
      "Error(sample_id/",
      params$target$sample_metadata_header,
      ")"
    ))

    # DE$sample_meta

    HSDEM_model <- HSDEM(
      alpha = params$posthoc$p_value,
      formula = formula, mtc = "none"
    )

    HSDEM_result <- model_apply(HSDEM_model, DE_tp)

    HSDEM_result_p_value <- HSDEM_result$p_value

    # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.

    colnames(HSDEM_result_p_value) <- plotrix::pasteCols(sapply(strsplit(colnames(HSDEM_result_p_value), " - "), sort), sep = "_")

    # We now add a specific suffix (`_p_value`) to each of the colnames

    colnames(HSDEM_result_p_value) <- paste0("tp_", i, "_", colnames(HSDEM_result_p_value), "_p_value")

    p_value_column <- colnames(HSDEM_result_p_value)

    # We set the row names as columns row_id to be able to merge the two dataframes

    HSDEM_result_p_value$row_id <- rownames(HSDEM_result_p_value)


    # We build a fold change model

    fold_change_model <- fold_change(
      factor_name = params$target$sample_metadata_header,
      paired = FALSE,
      sample_name = character(0),
      threshold = 0.5,
      control_group = character(0),
      method = "mean",
      conf_level = 0.95
    )


    fold_change_result <- model_apply(fold_change_model, DE_tp)

    # view(DE$data)
    # DE$data[,2]  <- c(-500,-500,-500,-500,500,500,500,500)

    # We suffix the column name of the dataframe with `_fold_change`, using dplyr rename function

    fold_change_result_fold_change <- fold_change_result$fold_change

    # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.
    # n !!!! We need to make sure that the header of metadata variable is not in the colnames of the fold change result

    colnames(fold_change_result_fold_change) <- plotrix::pasteCols(sapply(strsplit(colnames(fold_change_result_fold_change), "/"), sort), sep = "_")


    # We now add a specific suffix (`_p_value`) to each of the colnames

    colnames(fold_change_result_fold_change) <- paste0("tp_", i, "_", colnames(fold_change_result_fold_change), "_fold_change")


    fc_column <- colnames(fold_change_result_fold_change)


    # We set the row names as columns row_id to be able to merge the two dataframes

    fold_change_result_fold_change$row_id <- rownames(fold_change_result_fold_change)

    # # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

    # fold_change_result_fold_change = pivot_longer(fold_change_result_fold_change, cols = -row_id, names_to = "pairs", values_to = "fold_change")


    # We merge the two dataframes according to both the row_id and the pairs columns.

    DE_foldchange_pvalues <- merge(HSDEM_result_p_value, fold_change_result_fold_change, by = "row_id")


    # We add columns corresponding to the Log2 of the fold change column (suffix by fold_change). For this we use mutate_at function from dplyr package. We save the results in new columns with a novel suffix `_log2_FC`.

    message("Calculating logs ...")

    DE_foldchange_pvalues <- DE_foldchange_pvalues %>%
      mutate(across(contains("_fold_change"),
        .fns = list(log2 = ~ log2(.)),
        .names = "{col}_{fn}"
      )) %>%
      mutate(across(contains("_p_value"),
        .fns = list(minus_log10 = ~ -log10(.)),
        .names = "{col}_{fn}"
      ))


    l[[i]] <- DE_foldchange_pvalues
  }

  # We now merge the different dataframes in the list l

  DE_foldchange_pvalues <- Reduce(function(x, y) merge(x, y, by = "row_id"), l)
} else {
  # The formula is defined externally
  formula <- as.formula(paste0(
    "y", "~", params$target$sample_metadata_header, "+",
    "Error(sample_id/",
    params$target$sample_metadata_header,
    ")"
  ))


  # DE$sample_meta

  HSDEM_model <- HSDEM(
    alpha = params$posthoc$p_value,
    formula = formula, mtc = "none"
  )

  HSDEM_result <- model_apply(HSDEM_model, DE)

  HSDEM_result_p_value <- HSDEM_result$p_value


  # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.

  colnames(HSDEM_result_p_value) <- plotrix::pasteCols(sapply(strsplit(colnames(HSDEM_result_p_value), " - "), sort), sep = "_vs_")

  # We now add a specific suffix (`_p_value`) to each of the colnames

  colnames(HSDEM_result_p_value) <- paste0(colnames(HSDEM_result_p_value), "_p_value")

  p_value_column <- colnames(HSDEM_result_p_value)


  # We set the row names as columns row_id to be able to merge the two dataframes

  HSDEM_result_p_value$row_id <- rownames(HSDEM_result_p_value)

  # # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

  # HSDEM_result_p_value_long = pivot_longer(HSDEM_result_p_value, cols = -row_id, names_to = "pairs", values_to = "p_value")



  fold_change_model <- fold_change(
    factor_name = params$target$sample_metadata_header,
    paired = FALSE,
    sample_name = character(0),
    threshold = 0.5,
    control_group = character(0),
    method = "mean",
    conf_level = 0.95
  )

  # Check if this can be important to apply.

  DE_fc <- DE

  DE_fc$data <- DE_fc$data + 1

  fold_change_result <- model_apply(fold_change_model, DE_fc)

  # view(DE$data)
  # DE$data[,2]  <- c(-500,-500,-500,-500,500,500,500,500)

  # We suffix the column name of the dataframe with `_fold_change`, using dplyr rename function

  fold_change_result_fold_change <- fold_change_result$fold_change

  # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.
  # n !!!! We need to make sure that the header of metadata variable is not in the colnames of the fold change result


  colnames(fold_change_result_fold_change) <- plotrix::pasteCols(sapply(strsplit(colnames(fold_change_result_fold_change), "/"), sort), sep = "_vs_")


  # We now add a specific suffix (`_p_value`) to each of the colnames

  colnames(fold_change_result_fold_change) <- paste0(colnames(fold_change_result_fold_change), "_fold_change")


  fc_column <- colnames(fold_change_result_fold_change)


  # We set the row names as columns row_id to be able to merge the two dataframes

  fold_change_result_fold_change$row_id <- rownames(fold_change_result_fold_change)

  # # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

  # fold_change_result_fold_change = pivot_longer(fold_change_result_fold_change, cols = -row_id, names_to = "pairs", values_to = "fold_change")


  # We merge the two dataframes according to both the row_id and the pairs columns.

  DE_foldchange_pvalues <- merge(HSDEM_result_p_value, fold_change_result_fold_change, by = "row_id")


  # We add columns corresponding to the Log2 of the fold change column (suffix by fold_change). For this we use mutate_at function from dplyr package. We save the results in new columns with a novel suffix `_log2_FC`.

  message("Calculating logs ...")

  DE_foldchange_pvalues <- DE_foldchange_pvalues %>%
    mutate(across(contains("_fold_change"),
      .fns = list(log2 = ~ log2(.)),
      .names = "{col}_{fn}"
    )) %>%
    mutate(across(contains("_p_value"),
      .fns = list(minus_log10 = ~ -log10(.)),
      .names = "{col}_{fn}"
    ))
}

# We now merge the DE_foldchange_pvalues with the variable metadata using the row_ID column and the rownames of the variable metadata

DE_foldchange_pvalues <- merge(DE_foldchange_pvalues, DE$variable_meta, by.x = "row_id", by.y = "row.names")


# The file is exported

write.table(DE_foldchange_pvalues, file = filename_foldchange_pvalues, sep = ",")



##############################################################################
##############################################################################
############ Volcano Plots   #################################################
##############################################################################
##############################################################################


###### CrossTalk DT / Plotly - Volcano Plot

message("Launching Volcano Plots calculations ...")


formatted_qids <- paste("wd:", distinct_qids, sep = "") # Add "wd:" prefix
target_taxa <- paste(formatted_qids, collapse = "%0A") # Separate with "%0A" the URLencode equivalent of "\n"

# 'nan' strings in the met_annot_structure_smiles column are replaced by NA
DE_foldchange_pvalues$met_annot_structure_smiles[DE_foldchange_pvalues$met_annot_structure_smiles == "nan"] <- NA

if (gnps2_job) {
  # We first prepare the table for the dt export

  de4dt <- DE_foldchange_pvalues %>%
    select(
      feature_id,
      feature_id_full,
      sirius_chebiasciiname,
      sirius_chebiid,
      sirius_name,
      canopus_npc_pathway,
      canopus_npc_superclass,
      canopus_npc_class,
      canopus_npc_pathway_probability,
      canopus_npc_superclass_probability,
      canopus_npc_class_probability,
      feature_mz,
      feature_rt,
      gnps_component,
      gnps_compound_name,
      sirius_confidencescore,
      sirius_inchi,
      sirius_inchikey2d,
      sirius_molecularformula,
      sirius_adduct,
      sirius_smiles,
      met_annot_structure_inchi,
      met_annot_structure_inchikey,
      met_annot_structure_molecular_formula,
      met_annot_structure_nametraditional,
      met_annot_structure_smiles,
      met_annot_structure_taxonomy_npclassifier_01pathway,
      met_annot_structure_taxonomy_npclassifier_02superclass,
      met_annot_structure_taxonomy_npclassifier_03class,
      met_annot_structure_wikidata,
      met_annot_organism_name,
      met_annot_organism_taxonomy_01domain,
      met_annot_organism_taxonomy_02kingdom,
      met_annot_organism_taxonomy_03phylum,
      met_annot_organism_taxonomy_04class,
      met_annot_organism_taxonomy_05order,
      met_annot_organism_taxonomy_06family,
      met_annot_organism_taxonomy_07tribe,
      met_annot_organism_taxonomy_08genus,
      met_annot_organism_taxonomy_09species,
      met_annot_organism_taxonomy_10varietas,
      met_annot_organism_taxonomy_ottid,
      met_annot_organism_wikidata,
      met_annot_score_taxo,
      contains("p_value_minus_log10"),
      contains("p_value"),
      contains("fold_change_log2"),
      contains("fold_change")
    ) %>%
    # We format the smiles column to be able to display it in the datatable. We make sure this is only applied when sirius_smiles is not NA
    mutate(sirius_chemical_structure = ifelse(!is.na(sirius_smiles),
      sprintf('<img src="https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=%s&zoom=2.0" height="50"></img>', sirius_smiles),
      ""
    )) %>%
    mutate(met_annot_chemical_structure = ifelse(!is.na(met_annot_structure_smiles),
      sprintf('<img src="https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=%s&zoom=2.0" height="50"></img>', met_annot_structure_smiles),
      ""
    )) %>%
    mutate(met_annot_structure_wikidata = ifelse(!is.na(met_annot_structure_wikidata), sprintf('<a href="%s">%s</a>', met_annot_structure_wikidata, met_annot_structure_wikidata), "")) %>%
    mutate(met_annot_organism_wikidata = ifelse(!is.na(met_annot_organism_wikidata), sprintf('<a href="%s">%s</a>', met_annot_organism_wikidata, met_annot_organism_wikidata), "")) %>%
    mutate(sirius_name_url_safe = URLencode(sirius_name)) %>%
    # We then build the link to the PubChem website
    mutate(sirius_name = ifelse(!is.na(sirius_smiles),
      sprintf('<a href="https://pubchem.ncbi.nlm.nih.gov/#query=%s">%s</a>', sirius_name_url_safe, sirius_name),
      ""
    )) %>%
    # We then build the link to the CheBI website
    mutate(sirius_chebiid = ifelse(!is.na(sirius_chebiid),
      sprintf('<a href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s">%s</a>', sirius_chebiid, sirius_chebiid),
      ""
    )) %>%
    # We build a column for WD query
    mutate(wd_occurence_reports = ifelse(!is.na(sirius_inchikey2d), str_glue('<a href="https://query.wikidata.org/embed.html#SELECT%20%20%3Fcompound%20%3FInChIKey%20%3Ftaxon%20%3FtaxonLabel%20%3Fgenus_name%20%3Ffamily_name%20%3Fkingdom_name%20%3Freference%20%3FreferenceLabel%20WITH%20%7B%0A%20%20SELECT%20%3FqueryKey%20%3Fsrsearch%20%3Ffilter%20WHERE%20%7B%0A%20%20%20%20VALUES%20%3FqueryKey%20%7B%0A%20%20%20%20%20%20%22{sirius_inchikey2d}%22%0A%20%20%20%20%7D%0A%20%20%20%20BIND%20%28CONCAT%28substr%28%24queryKey%2C1%2C14%29%2C%20%22%20haswbstatement%3AP235%22%29%20AS%20%3Fsrsearch%29%0A%20%20%20%20BIND%20%28CONCAT%28%22%5E%22%2C%20substr%28%24queryKey%2C1%2C14%29%29%20AS%20%3Ffilter%29%0A%20%20%7D%0A%7D%20AS%20%25comps%20WITH%20%7B%0A%20%20SELECT%20%3Fcompound%20%3FInChIKey%20WHERE%20%7B%0A%20%20%20%20INCLUDE%20%25comps%0A%20%20%20%20%20%20%20%20%20%20%20%20SERVICE%20wikibase%3Amwapi%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Aendpoint%20%22www.wikidata.org%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20wikibase%3Aapi%20%22Search%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrsearch%20%3Fsrsearch%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrlimit%20%22max%22.%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%3Fcompound%20wikibase%3AapiOutputItem%20mwapi%3Atitle.%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%3Fcompound%20wdt%3AP235%20%3FInChIKey%20.%0A%20%20%20%20FILTER%20%28REGEX%28STR%28%3FInChIKey%29%2C%20%3Ffilter%29%29%0A%20%20%7D%0A%7D%20AS%20%25compounds%0AWHERE%20%7B%0A%20%20INCLUDE%20%25compounds%0A%20%20%20VALUES%20%3Ftaxon%20%7B%0A%20%20%20%20%20%20{target_taxa}%0A%20%20%20%20%7D%0A%20%20%7B%0A%20%20%20%20%3Fcompound%20p%3AP703%20%3Fstmt.%0A%20%20%20%20%3Fstmt%20ps%3AP703%20%3Ftaxon.%0A%20%20%20%20%3Fkingdom%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ36732%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fkingdom_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Ffamily%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ35409%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Ffamily_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Fgenus%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ34740%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fgenus_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%7D%0A%20%20OPTIONAL%20%7B%0A%20%20%20%20%3Fstmt%20prov%3AwasDerivedFrom%20%3Fref.%0A%20%20%20%20%3Fref%20pr%3AP248%20%3Freference.%0A%20%20%7D%20%0A%20%20SERVICE%20wikibase%3Alabel%20%7B%20bd%3AserviceParam%20wikibase%3Alanguage%20%22en%22.%20%7D%0A%7D%0ALIMIT%2010000">Biological occurences of this molecule (limited to organism(s) of the current dataset)</a>'), "")) %>%
    # We build a column for WD query
    mutate(wd_occurence_reports_all = ifelse(!is.na(sirius_inchikey2d), str_glue('<a href="https://query.wikidata.org/embed.html#SELECT%20%20%3Fcompound%20%3FInChIKey%20%3Ftaxon%20%3FtaxonLabel%20%3Fgenus_name%20%3Ffamily_name%20%3Fkingdom_name%20%3Freference%20%3FreferenceLabel%20%0AWITH%20%7B%0A%20%20SELECT%20%3FqueryKey%20%3Fsrsearch%20%3Ffilter%20WHERE%20%7B%0A%20%20%20%20VALUES%20%3FqueryKey%20%7B%0A%20%20%20%20%20%20%22{sirius_inchikey2d}%22%0A%20%20%20%20%7D%0A%20%20%20%20BIND%20%28CONCAT%28substr%28%24queryKey%2C1%2C14%29%2C%20%22%20haswbstatement%3AP235%22%29%20AS%20%3Fsrsearch%29%0A%20%20%20%20BIND%20%28CONCAT%28%22%5E%22%2C%20substr%28%24queryKey%2C1%2C14%29%29%20AS%20%3Ffilter%29%0A%20%20%7D%0A%7D%20AS%20%25comps%20WITH%20%7B%0A%20%20SELECT%20%3Fcompound%20%3FInChIKey%20WHERE%20%7B%0A%20%20%20%20INCLUDE%20%25comps%0A%20%20%20%20SERVICE%20wikibase%3Amwapi%20%7B%0A%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Aendpoint%20%22www.wikidata.org%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20wikibase%3Aapi%20%22Search%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrsearch%20%3Fsrsearch%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrlimit%20%22max%22.%0A%20%20%20%20%20%20%3Fcompound%20wikibase%3AapiOutputItem%20mwapi%3Atitle.%0A%20%20%20%20%7D%0A%20%20%20%20%3Fcompound%20wdt%3AP235%20%3FInChIKey%20.%0A%20%20%20%20FILTER%20%28REGEX%28STR%28%3FInChIKey%29%2C%20%3Ffilter%29%29%0A%20%20%7D%0A%7D%20AS%20%25compounds%0AWHERE%20%7B%0A%20%20INCLUDE%20%25compounds%0A%20%20%7B%0A%20%20%20%20%3Fcompound%20p%3AP703%20%3Fstmt.%0A%20%20%20%20%3Fstmt%20ps%3AP703%20%3Ftaxon.%0A%20%20%20%20%3Fkingdom%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ36732%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fkingdom_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Ffamily%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ35409%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Ffamily_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Fgenus%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ34740%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fgenus_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20%0A%20%20%7D%0A%20%20OPTIONAL%20%7B%0A%20%20%20%20%3Fstmt%20prov%3AwasDerivedFrom%20%3Fref.%0A%20%20%20%20%3Fref%20pr%3AP248%20%3Freference.%0A%20%20%7D%20%0A%20%20SERVICE%20wikibase%3Alabel%20%7B%20bd%3AserviceParam%20wikibase%3Alanguage%20%22en%22.%20%7D%0A%7D%0ALIMIT%2010000%0A">All biological occurences of this molecule</a>'), "")) %>%
    # We build a column for the gnps plotter for interactive box plots
    mutate(gnps_plotter_box_plot = str_glue('<a href="http://plotter.gnps2.org/?gnps_quant_table_usi=mzspec%3AGNPS2%3ATASK-{params$gnps_job_id}-nf_output%2Fclustering%2Ffeaturetable_reformated.csv&gnps_metadata_table_usi=mzspec%3AGNPS2%3ATASK-{params$gnps_job_id}-nf_output%2Fmetadata%2Fmerged_metadata.tsv&feature={feature_id}&filter_metadata_column=None&filter_metadata_value=%5B%5D&metadata={params$options$gnps_column_for_boxplots$factor_name}&facet=&groups={params$options$gnps_column_for_boxplots$factor_name}&plot_type=box&color_column={params$options$gnps_column_for_boxplots$factor_name}&color_selection=%5B%5D&points_toggle=False&theme=ggplot2&animation_column=&lat_column=&long_column=&map_animation_column=&map_scope=world">Box plots for {feature_id}</a>')) %>%
    select(
      feature_id,
      feature_id_full,
      sirius_chebiasciiname,
      sirius_chemical_structure,
      sirius_chebiid,
      sirius_name,
      wd_occurence_reports,
      wd_occurence_reports_all,
      canopus_npc_pathway,
      canopus_npc_superclass,
      canopus_npc_class,
      canopus_npc_pathway_probability,
      canopus_npc_superclass_probability,
      canopus_npc_class_probability,
      feature_mz,
      feature_rt,
      gnps_component,
      gnps_compound_name,
      gnps_plotter_box_plot,
      sirius_confidencescore,
      sirius_inchi,
      sirius_inchikey2d,
      sirius_molecularformula,
      sirius_adduct,
      sirius_smiles,
      met_annot_chemical_structure,
      met_annot_structure_inchi,
      met_annot_structure_inchikey,
      met_annot_structure_molecular_formula,
      met_annot_structure_nametraditional,
      met_annot_structure_smiles,
      met_annot_structure_taxonomy_npclassifier_01pathway,
      met_annot_structure_taxonomy_npclassifier_02superclass,
      met_annot_structure_taxonomy_npclassifier_03class,
      met_annot_structure_wikidata,
      met_annot_organism_name,
      met_annot_organism_taxonomy_01domain,
      met_annot_organism_taxonomy_02kingdom,
      met_annot_organism_taxonomy_03phylum,
      met_annot_organism_taxonomy_04class,
      met_annot_organism_taxonomy_05order,
      met_annot_organism_taxonomy_06family,
      met_annot_organism_taxonomy_07tribe,
      met_annot_organism_taxonomy_08genus,
      met_annot_organism_taxonomy_09species,
      met_annot_organism_taxonomy_10varietas,
      met_annot_organism_taxonomy_ottid,
      met_annot_organism_wikidata,
      met_annot_score_taxo,
      contains("p_value_minus_log10"),
      contains("p_value"),
      contains("fold_change_log2"),
      contains("fold_change")
    ) %>%
    # We set the type of the sirius_confidencescore column to numeric
    mutate(sirius_confidencescore = as.numeric(sirius_confidencescore))
} else {
  # We first prepare the table for the dt export

  de4dt <- DE_foldchange_pvalues %>%
    select(
      feature_id,
      feature_id_full,
      sirius_chebiasciiname,
      sirius_chebiid,
      sirius_name,
      canopus_npc_pathway,
      canopus_npc_superclass,
      canopus_npc_class,
      canopus_npc_pathway_probability,
      canopus_npc_superclass_probability,
      canopus_npc_class_probability,
      feature_mz,
      feature_rt,
      gnps_componentindex,
      gnps_gnpslinkout_network,
      gnps_libraryid,
      sirius_confidencescore,
      sirius_inchi,
      sirius_inchikey2d,
      sirius_molecularformula,
      sirius_adduct,
      sirius_smiles,
      met_annot_structure_inchi,
      met_annot_structure_inchikey,
      met_annot_structure_molecular_formula,
      met_annot_structure_nametraditional,
      met_annot_structure_smiles,
      met_annot_structure_taxonomy_npclassifier_01pathway,
      met_annot_structure_taxonomy_npclassifier_02superclass,
      met_annot_structure_taxonomy_npclassifier_03class,
      met_annot_structure_wikidata,
      met_annot_organism_name,
      met_annot_organism_taxonomy_01domain,
      met_annot_organism_taxonomy_02kingdom,
      met_annot_organism_taxonomy_03phylum,
      met_annot_organism_taxonomy_04class,
      met_annot_organism_taxonomy_05order,
      met_annot_organism_taxonomy_06family,
      met_annot_organism_taxonomy_07tribe,
      met_annot_organism_taxonomy_08genus,
      met_annot_organism_taxonomy_09species,
      met_annot_organism_taxonomy_10varietas,
      met_annot_organism_taxonomy_ottid,
      met_annot_organism_wikidata,
      score_taxo,
      contains("p_value_minus_log10"),
      contains("p_value"),
      contains("fold_change_log2"),
      contains("fold_change")
    ) %>%
    # We format the smiles column to be able to display it in the datatable. We make sure this is only applied when sirius_smiles is not NA
    mutate(sirius_chemical_structure = ifelse(!is.na(sirius_smiles),
      sprintf('<img src="https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=%s&zoom=2.0" height="50"></img>', sirius_smiles),
      ""
    )) %>%
    mutate(met_annot_chemical_structure = ifelse(!is.na(met_annot_structure_smiles),
      sprintf('<img src="https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=%s&zoom=2.0" height="50"></img>', met_annot_structure_smiles),
      ""
    )) %>%
    mutate(met_annot_structure_wikidata = ifelse(!is.na(met_annot_structure_wikidata), sprintf('<a href="%s">%s</a>', met_annot_structure_wikidata, met_annot_structure_wikidata), "")) %>%
    mutate(met_annot_organism_wikidata = ifelse(!is.na(met_annot_organism_wikidata), sprintf('<a href="%s">%s</a>', met_annot_organism_wikidata, met_annot_organism_wikidata), "")) %>%
    mutate(cluster_gnps_link = sprintf('<a href="%s">%s</a>', gnps_gnpslinkout_network, gnps_componentindex)) %>%
    # mutate(spectra_gnps_link = sprintf("<a href='%s'>gnps spectrum %s</a>", gnpslinkout_cluster_gnps, feature_id)) %>%
    # We first sanitize the sirius_name column and make it URL safe
    mutate(sirius_name_url_safe = URLencode(sirius_name)) %>%
    # We then build the link to the PubChem website
    mutate(sirius_name = ifelse(!is.na(sirius_smiles),
      sprintf('<a href="https://pubchem.ncbi.nlm.nih.gov/#query=%s">%s</a>', sirius_name_url_safe, sirius_name),
      ""
    )) %>%
    # We then build the link to the CheBI website
    mutate(sirius_chebiid = ifelse(!is.na(sirius_chebiid),
      sprintf('<a href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s">%s</a>', sirius_chebiid, sirius_chebiid),
      ""
    )) %>%
    # We build a column for WD query
    mutate(wd_occurence_reports = ifelse(!is.na(sirius_inchikey2d), str_glue('<a href="https://query.wikidata.org/embed.html#SELECT%20%20%3Fcompound%20%3FInChIKey%20%3Ftaxon%20%3FtaxonLabel%20%3Fgenus_name%20%3Ffamily_name%20%3Fkingdom_name%20%3Freference%20%3FreferenceLabel%20WITH%20%7B%0A%20%20SELECT%20%3FqueryKey%20%3Fsrsearch%20%3Ffilter%20WHERE%20%7B%0A%20%20%20%20VALUES%20%3FqueryKey%20%7B%0A%20%20%20%20%20%20%22{sirius_inchikey2d}%22%0A%20%20%20%20%7D%0A%20%20%20%20BIND%20%28CONCAT%28substr%28%24queryKey%2C1%2C14%29%2C%20%22%20haswbstatement%3AP235%22%29%20AS%20%3Fsrsearch%29%0A%20%20%20%20BIND%20%28CONCAT%28%22%5E%22%2C%20substr%28%24queryKey%2C1%2C14%29%29%20AS%20%3Ffilter%29%0A%20%20%7D%0A%7D%20AS%20%25comps%20WITH%20%7B%0A%20%20SELECT%20%3Fcompound%20%3FInChIKey%20WHERE%20%7B%0A%20%20%20%20INCLUDE%20%25comps%0A%20%20%20%20%20%20%20%20%20%20%20%20SERVICE%20wikibase%3Amwapi%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Aendpoint%20%22www.wikidata.org%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20wikibase%3Aapi%20%22Search%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrsearch%20%3Fsrsearch%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrlimit%20%22max%22.%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%3Fcompound%20wikibase%3AapiOutputItem%20mwapi%3Atitle.%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%3Fcompound%20wdt%3AP235%20%3FInChIKey%20.%0A%20%20%20%20FILTER%20%28REGEX%28STR%28%3FInChIKey%29%2C%20%3Ffilter%29%29%0A%20%20%7D%0A%7D%20AS%20%25compounds%0AWHERE%20%7B%0A%20%20INCLUDE%20%25compounds%0A%20%20%20VALUES%20%3Ftaxon%20%7B%0A%20%20%20%20%20%20{target_taxa}%0A%20%20%20%20%7D%0A%20%20%7B%0A%20%20%20%20%3Fcompound%20p%3AP703%20%3Fstmt.%0A%20%20%20%20%3Fstmt%20ps%3AP703%20%3Ftaxon.%0A%20%20%20%20%3Fkingdom%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ36732%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fkingdom_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Ffamily%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ35409%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Ffamily_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Fgenus%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ34740%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fgenus_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%7D%0A%20%20OPTIONAL%20%7B%0A%20%20%20%20%3Fstmt%20prov%3AwasDerivedFrom%20%3Fref.%0A%20%20%20%20%3Fref%20pr%3AP248%20%3Freference.%0A%20%20%7D%20%0A%20%20SERVICE%20wikibase%3Alabel%20%7B%20bd%3AserviceParam%20wikibase%3Alanguage%20%22en%22.%20%7D%0A%7D%0ALIMIT%2010000">Biological occurences of this molecule (limited to organism(s) of the current dataset)</a>'), "")) %>%
    # We build a column for WD query
    mutate(wd_occurence_reports_all = ifelse(!is.na(sirius_inchikey2d), str_glue('<a href="https://query.wikidata.org/embed.html#SELECT%20%20%3Fcompound%20%3FInChIKey%20%3Ftaxon%20%3FtaxonLabel%20%3Fgenus_name%20%3Ffamily_name%20%3Fkingdom_name%20%3Freference%20%3FreferenceLabel%20%0AWITH%20%7B%0A%20%20SELECT%20%3FqueryKey%20%3Fsrsearch%20%3Ffilter%20WHERE%20%7B%0A%20%20%20%20VALUES%20%3FqueryKey%20%7B%0A%20%20%20%20%20%20%22{sirius_inchikey2d}%22%0A%20%20%20%20%7D%0A%20%20%20%20BIND%20%28CONCAT%28substr%28%24queryKey%2C1%2C14%29%2C%20%22%20haswbstatement%3AP235%22%29%20AS%20%3Fsrsearch%29%0A%20%20%20%20BIND%20%28CONCAT%28%22%5E%22%2C%20substr%28%24queryKey%2C1%2C14%29%29%20AS%20%3Ffilter%29%0A%20%20%7D%0A%7D%20AS%20%25comps%20WITH%20%7B%0A%20%20SELECT%20%3Fcompound%20%3FInChIKey%20WHERE%20%7B%0A%20%20%20%20INCLUDE%20%25comps%0A%20%20%20%20SERVICE%20wikibase%3Amwapi%20%7B%0A%20%20%20%20%20%20bd%3AserviceParam%20wikibase%3Aendpoint%20%22www.wikidata.org%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20wikibase%3Aapi%20%22Search%22%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrsearch%20%3Fsrsearch%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20mwapi%3Asrlimit%20%22max%22.%0A%20%20%20%20%20%20%3Fcompound%20wikibase%3AapiOutputItem%20mwapi%3Atitle.%0A%20%20%20%20%7D%0A%20%20%20%20%3Fcompound%20wdt%3AP235%20%3FInChIKey%20.%0A%20%20%20%20FILTER%20%28REGEX%28STR%28%3FInChIKey%29%2C%20%3Ffilter%29%29%0A%20%20%7D%0A%7D%20AS%20%25compounds%0AWHERE%20%7B%0A%20%20INCLUDE%20%25compounds%0A%20%20%7B%0A%20%20%20%20%3Fcompound%20p%3AP703%20%3Fstmt.%0A%20%20%20%20%3Fstmt%20ps%3AP703%20%3Ftaxon.%0A%20%20%20%20%3Fkingdom%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ36732%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fkingdom_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Ffamily%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ35409%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Ffamily_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20.%0A%20%20%20%20%3Fgenus%20wdt%3AP31%20wd%3AQ16521%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP105%20wd%3AQ34740%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20wdt%3AP225%20%3Fgenus_name%20%3B%0A%20%20%20%20%20%20%20%20%20%20%20%5Ewdt%3AP171%2a%20%3Ftaxon%20%0A%20%20%7D%0A%20%20OPTIONAL%20%7B%0A%20%20%20%20%3Fstmt%20prov%3AwasDerivedFrom%20%3Fref.%0A%20%20%20%20%3Fref%20pr%3AP248%20%3Freference.%0A%20%20%7D%20%0A%20%20SERVICE%20wikibase%3Alabel%20%7B%20bd%3AserviceParam%20wikibase%3Alanguage%20%22en%22.%20%7D%0A%7D%0ALIMIT%2010000%0A">All biological occurences of this molecule</a>'), "")) %>%
    # We build a column for the gnps plotter for interactive box plots
    mutate(gnps_plotter_box_plot = str_glue('<a href="http://plotter.gnps2.org/?gnps_tall_table_usi=mzdata%3AGNPS%3ATASK-{params$gnps_job_id}-feature_statistics%2Fdata_long.csv&gnps_quant_table_usi=&gnps_metadata_table_usi=&feature={feature_id}&filter_metadata_column=None&filter_metadata_value=%5B%5D&metadata={params$options$gnps_column_for_boxplots$factor_name}&facet=&groups=&plot_type=box&color_column={params$options$gnps_column_for_boxplots$factor_name}&color_selection=%5B%5D&points_toggle=True&theme=ggplot2&animation_column=&lat_column=&long_column=&map_animation_column=&map_scope=world">Box plots for {feature_id}</a>')) %>%
    select(
      feature_id,
      feature_id_full,
      sirius_chebiasciiname,
      sirius_chemical_structure,
      sirius_chebiid,
      sirius_name,
      wd_occurence_reports,
      wd_occurence_reports_all,
      canopus_npc_pathway,
      canopus_npc_superclass,
      canopus_npc_class,
      canopus_npc_pathway_probability,
      canopus_npc_superclass_probability,
      canopus_npc_class_probability,
      feature_mz,
      feature_rt,
      cluster_gnps_link,
      gnps_libraryid,
      spectra_gnps_link,
      gnps_plotter_box_plot,
      sirius_confidencescore,
      sirius_inchi,
      sirius_inchikey2d,
      sirius_molecularformula,
      sirius_adduct,
      sirius_smiles,
      met_annot_chemical_structure,
      met_annot_structure_inchi,
      met_annot_structure_inchikey,
      met_annot_structure_molecular_formula,
      met_annot_structure_nametraditional,
      met_annot_structure_smiles,
      met_annot_structure_taxonomy_npclassifier_01pathway,
      met_annot_structure_taxonomy_npclassifier_02superclass,
      met_annot_met_annot_structure_taxonomy_npclassifier_02superclass,
      met_annot_structure_wikidata,
      met_annot_organism_name,
      met_annot_organism_taxonomy_01domain,
      met_annot_organism_taxonomy_02kingdom,
      met_annot_organism_taxonomy_03phylum,
      met_annot_organism_taxonomy_04class,
      met_annot_organism_taxonomy_05order,
      met_annot_organism_taxonomy_06family,
      met_annot_organism_taxonomy_07tribe,
      met_annot_organism_taxonomy_08genus,
      met_annot_organism_taxonomy_09species,
      met_annot_organism_taxonomy_10varietas,
      met_annot_organism_taxonomy_ottid,
      met_annot_organism_wikidata,
      met_annot_score_taxo,
      contains("p_value_minus_log10"),
      contains("p_value"),
      contains("fold_change_log2"),
      contains("fold_change")
    ) %>%
    # We set the type of the sirius_confidencescore column to numeric
    mutate(sirius_confidencescore = as.numeric(sirius_confidencescore))
}
# We outpuzt a generic DT for data exploration of the whole set


### Defining the DT object
DT_volcano <- datatable(de4dt,
  escape = FALSE,
  rownames = FALSE,
  extensions = c("Buttons", "Select"),
  selection = "none",
  filter = "top",
  class = list(stripe = FALSE),
  options =
    list(
      #       initComplete = JS(
      #   "function(settings, json) {",
      #   "$('body').css({'font-family': 'Calibri'});",
      #   "}"
      # ),
      pageLength = 10,
      select = TRUE,
      searching = TRUE,
      scrollX = TRUE,
      scrollY = TRUE,
      dom = "Blfrtip",
      buttons = list(
        list(
          extend = "copy",
          text = "Copy"
          # ,
          # exportOptions = list(modifier = list(selected = TRUE))
        ),
        list(
          extend = "csv",
          text = "CSV"
          # exportOptions = list(modifier = list(selected = TRUE))
        ),
        list(
          extend = "excel",
          text = "Excel"
          # exportOptions = list(modifier = list(selected = TRUE))
        ),
        list(
          extend = "pdf",
          text = "PDF"
          # exportOptions = list(modifier = list(selected = TRUE))
        ),
        list(
          extend = "print",
          text = "Print"
          # exportOptions = list(modifier = list(selected = TRUE))
        )
      ),
      lengthMenu = list(
        c(10, 25, 50, -1),
        c(10, 25, 50, "All")
      )
    )
) %>%
  # formatRound(c("log2_fold_change", "pvalue_minus_log10"), digits = 3) %>%
  # formatSignif(c("log2_fold_change", "pvalue_minus_log10"), digits = 3)  %>%
  formatRound(c("feature_mz", "feature_rt"), digits = 3) %>%
  formatRound(c(
    "canopus_npc_pathway_probability",
    "canopus_npc_superclass_probability",
    "canopus_npc_class_probability",
    "sirius_confidencescore"
  ), digits = 2)



if (params$operating_system$system == "unix") {
  ### linux version
  htmltools::save_html(DT_volcano, file = paste0("DT_full_dataset.html"))
}



if (params$operating_system$system == "windows") {
  ### windows version
  Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
  htmltools::save_html(DT_volcano, file = paste0("DT_full_dataset.html"), libdir = "lib")
  unlink("lib", recursive = FALSE)
}



# Extract prefixes of columns with "_p_value" suffix
conditions <- sub("_p_value$", "", grep("_p_value$", names(DE_foldchange_pvalues), value = TRUE))


# Print message before iterating over conditions
message("Iterating over the following conditions for the Volcano plots generation:\n")


# Iterate over the prefixes
for (condition in conditions) {
  # condition = "Argon_HN_NifH_vs_WT"
  message("Generating Volcano plot for condition: ", condition, "\n")


  # Perform filtering using the prefix as a condition
  de <- de4dt %>%
    filter(!!sym(paste0(condition, "_p_value_minus_log10")) > 0)

  # Print the filtered data
  message("Filtered data for condition: ", condition, "\n")
  # print(head(DE_foldchange_pvalues_signi))
  message("\n")

  condition_parts <- strsplit(condition, "_vs_")[[1]]
  first_part <- condition_parts[1]
  second_part <- condition_parts[2]


  de <- de %>%
    # we rename the day_vs_night_p_value_minus_log10 column to pvalue
    rename(pvalue_minus_log10 = !!sym(paste0(condition, "_p_value_minus_log10"))) %>%
    # we rename the day_vs_night_fold_change_log2 column to log2FoldChange
    rename(log2_fold_change = !!sym(paste0(condition, "_fold_change_log2")))


  # m <- SharedData$new(x, key = ~feature_id)

  m <- SharedData$new(de)

  ### Defining the plotly object


  plotly_volcano <- plot_ly(m, x = ~log2_fold_change, y = ~pvalue_minus_log10) %>%
    add_markers(text = row.names(m)) %>%
    add_markers(text = row.names(m), yaxis = "y2") %>%
    # config(displayModeBar = FALSE) %>%
    layout(
      title = "Hold shift while clicking \n markers for persistent selection",
      margin = list(t = 60)
    ) %>%
    layout(
      title = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters: [", filter_sample_metadata_status, "]"),
      margin = list(
        l = 100, # Left margin in pixels, adjust as needed
        r = 100, # Right margin in pixels, adjust as needed
        t = 100, # Top margin in pixels, adjust as needed
        b = 100 # Bottom margin in pixels, adjust as needed
      )
    ) %>%
    layout(
      title = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters: [", filter_sample_metadata_status, "]"), plot_bgcolor = "#e5ecf6",
      xaxis = list(title = "-log10(pvalue)"),
      yaxis = list(title = paste0("log2(FC) ", first_part), side = "left"),
      yaxis2 = list(title = paste0("log2(FC) ", second_part), side = "right")
    ) %>%
    layout(showlegend = FALSE)

  # gg_plotly_volcano <- ggplotly(plotly_volcano)

  ### Defining the DT object

  DT_volcano <- datatable(m,
    escape = FALSE,
    rownames = FALSE,
    extensions = c("Buttons", "Select"),
    selection = "none",
    filter = "top",
    class = list(stripe = FALSE),
    options =
      list(
        #       initComplete = JS(
        #   "function(settings, json) {",
        #   "$('body').css({'font-family': 'Calibri'});",
        #   "}"
        # ),
        pageLength = 10,
        select = TRUE,
        searching = TRUE,
        scrollX = TRUE,
        scrollY = TRUE,
        dom = "Blfrtip",
        buttons = list(
          list(
            extend = "copy",
            text = "Copy"
            # ,
            # exportOptions = list(modifier = list(selected = TRUE))
          ),
          list(
            extend = "csv",
            text = "CSV"
            # exportOptions = list(modifier = list(selected = TRUE))
          ),
          list(
            extend = "excel",
            text = "Excel"
            # exportOptions = list(modifier = list(selected = TRUE))
          ),
          list(
            extend = "pdf",
            text = "PDF"
            # exportOptions = list(modifier = list(selected = TRUE))
          ),
          list(
            extend = "print",
            text = "Print"
            # exportOptions = list(modifier = list(selected = TRUE))
          )
        ),
        lengthMenu = list(
          c(10, 25, 50, -1),
          c(10, 25, 50, "All")
        )
      )
  ) %>%
    formatRound(c("log2_fold_change", "pvalue_minus_log10"), digits = 3) %>%
    formatSignif(c("log2_fold_change", "pvalue_minus_log10"), digits = 3) %>%
    formatRound(c("feature_mz", "feature_rt"), digits = 3) %>%
    formatRound(c(
      "canopus_npc_pathway_probability",
      "canopus_npc_superclass_probability",
      "canopus_npc_class_probability",
      "sirius_confidencescore"
    ), digits = 2)



  ### Defining the crosstalked object

  # plotly_DT_crosstalked <- bscols(
  #   plotly_volcano %>%
  #     highlight(
  #       color = "green", on = "plotly_selected",
  #       off = "plotly_deselect"
  #     ),
  #   DT_volcano
  # )

  plotly_DT_crosstalked_div <- browsable(div(
    div(
      style = "display: grid; grid-template-columns: 1fr;",
      plotly_volcano %>%
        highlight(
          color = "green", on = "plotly_selected",
          off = "plotly_deselect"
        )
    ), DT_volcano
  ))

  ### Saving the plotly_DT_crosstalked object

  if (params$operating_system$system == "unix") {
    ### linux version
    htmltools::save_html(plotly_DT_crosstalked_div, file = paste0("Volcano_DT_", first_part, "_vs_", second_part, ".html"))
  }



  if (params$operating_system$system == "windows") {
    ### windows version
    Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
    htmltools::save_html(plotly_DT_crosstalked_div, file = paste0("Volcano_DT_", first_part, "_vs_", second_part, ".html"), libdir = "lib")
    unlink("lib", recursive = FALSE)
  }

  # We now generate the associated ggplots

  log2_fold_change_threshold <- 0.25
  pvalue_minus_log10_threshold <- 1.3

  # add a column of NAs
  de$diffexpressed <- "NO"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2_fold_change < -log2_fold_change_threshold & de$pvalue_minus_log10 > pvalue_minus_log10_threshold] <- first_part
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  de$diffexpressed[de$log2_fold_change > log2_fold_change_threshold & de$pvalue_minus_log10 > pvalue_minus_log10_threshold] <- second_part

  # We define a vector of columns to use as labels

  label_columns <- c("sirius_chebiasciiname", "canopus_npc_class", "canopus_npc_superclass", "canopus_npc_pathway")


  # Now we iterate over the vectors

  for (i in 1:length(label_columns)) {
    # We define the label column

    label_column <- label_columns[i]

    # We define the ggplot object

    de$delabel <- NA
    de$delabel[de$diffexpressed != "NO"] <- de[[label_column]][de$diffexpressed != "NO"]

    # cols <- setNames(c(params$colors$volcano), c(first_part, second_part))

    cols <- custom_colors[c(first_part, second_part)]

    # Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function

    # plot adding up all layers we have seen so far
    gg_volcano <- ggplot(data = de, aes(x = log2_fold_change, y = pvalue_minus_log10, col = diffexpressed, label = delabel)) +
      geom_point() +
      theme_minimal() +
      geom_label_repel() +
      scale_colour_manual(name = "Differentially\nExpressed", values = cols) +
      geom_vline(xintercept = c(-log2_fold_change_threshold, log2_fold_change_threshold), col = "grey", linewidth = 0.2) +
      geom_hline(yintercept = pvalue_minus_log10_threshold, col = "grey", linewidth = 0.2) +
      labs(title = title_volcano) # Set the ggplot title

    # We save the plot in a pdf file
    tryCatch(
      {
        ggsave(plot = gg_volcano, filename = paste0("Volcano_", first_part, "_vs_", second_part, "_", label_column, ".pdf"), width = 10, height = 10)
      },
      error = function(e) {}
    )
  }
}



##############################################################################
##############################################################################
############ Treemaps fold change ############################################
##############################################################################
##############################################################################



treat_npclassifier_json <- function(taxonomy) {
  taxonomy_classes <- taxonomy$Class %>%
    rbind()
  rownames(taxonomy_classes) <- "id_class"
  taxonomy_classes <- taxonomy_classes %>%
    t() %>%
    data.frame() %>%
    mutate(
      class = rownames(.),
      id_class = as.numeric(id_class)
    )

  taxonomy_superclasses <- taxonomy$Superclass %>%
    rbind()
  rownames(taxonomy_superclasses) <- "id_superclass"
  taxonomy_superclasses <- taxonomy_superclasses %>%
    t() %>%
    data.frame() %>%
    mutate(
      superclass = rownames(.),
      id_superclass = as.numeric(id_superclass)
    )

  taxonomy_pathways <- taxonomy$Pathway %>%
    rbind()
  rownames(taxonomy_pathways) <- "id_pathway"
  taxonomy_pathways <- taxonomy_pathways %>%
    t() %>%
    data.frame() %>%
    mutate(
      pathway = rownames(.),
      id_pathway = as.numeric(id_pathway)
    )

  taxonomy_hierarchy_class <- taxonomy$Class_hierarchy

  id_pathway <- list()
  id_superclass <- list()
  id_class <- list()

  for (i in seq_len(length(taxonomy_hierarchy_class))) {
    id_pathway[[i]] <- taxonomy_hierarchy_class[[i]]$Pathway
    id_superclass[[i]] <- taxonomy_hierarchy_class[[i]]$Superclass
    id_class[[i]] <- names(taxonomy_hierarchy_class[i])
  }

  zu <- cbind(id_pathway, id_superclass, id_class) %>%
    data.frame() %>%
    mutate(id_class = as.numeric(id_class)) %>%
    unnest(id_superclass) %>%
    unnest(id_pathway)

  ## No idea why would this be needed... class already has everything?

  id_pathway_2 <- list()
  id_superclass <- list()

  taxonomy_hierarchy_superclass <- taxonomy$Super_hierarchy

  for (i in seq_len(length(taxonomy_hierarchy_superclass))) {
    id_pathway_2[[i]] <- taxonomy_hierarchy_superclass[[i]]$Pathway
    id_superclass[[i]] <- names(taxonomy_hierarchy_superclass[i])
  }

  zu_2 <- cbind(id_pathway_2, id_superclass) %>%
    data.frame() %>%
    mutate(id_superclass = as.numeric(id_superclass)) %>%
    unnest(id_pathway_2)

  taxonomy_semicleaned <- full_join(zu, taxonomy_classes) %>%
    full_join(., taxonomy_superclasses) %>%
    full_join(., taxonomy_pathways) %>%
    distinct(class, superclass, pathway)
  return(taxonomy_semicleaned)
}



# ################################### function
# ################################# treemap shaper

dt_for_treemap <- function(datatable, parent_value, value, count) {
  parent_value <- enquo(parent_value)
  value <- enquo(value)
  count <- enquo(count)

  datatable <- data.frame(datatable %>%
    group_by(!!parent_value, !!value, ) %>%
    summarise(count = sum(as.numeric(!!count), na.rm = T)))

  datatable <- datatable %>%
    select(!!parent_value, !!value, count) %>% # create id labels for each row # Notre the !! to pass aruguments to a dplyr function
    rename(
      parent.value = !!parent_value,
      value = !!value
    ) %>%
    mutate(ids = ifelse(parent.value == "", value,
      paste0(value, "-", parent.value) # Notre that here we are passing argument to a non dplyr function call
    )) %>%
    select(ids, everything())

  par_info <- datatable %>% dplyr::group_by(parent.value) %>% # group by parent
    dplyr::summarise(count = sum(as.numeric(count), na.rm = T)) %>% # parent total
    rename(value = parent.value) %>% # parent labels for the item field
    mutate(parent.value = "", ids = value) %>% # add missing fields for my_data
    select(names(datatable)) # put cols in same order as my_data

  data_for_plot <- rbind(datatable, par_info)

  return(data_for_plot)
}
# ###################################################################################
# ###################################################################################

dt_for_treemap_mean <- function(datatable, parent_value, value, count) {
  parent_value <- enquo(parent_value)
  value <- enquo(value)
  count <- enquo(count)

  datatable <- data.frame(datatable %>%
    group_by(!!parent_value, !!value, ) %>%
    summarise(count = mean(as.numeric(!!count), na.rm = T)))

  datatable <- datatable %>%
    select(!!parent_value, !!value, count) %>% # create id labels for each row # Notre the !! to pass aruguments to a dplyr function
    rename(
      parent.value = !!parent_value,
      value = !!value
    ) %>%
    mutate(ids = ifelse(parent.value == "", value,
      paste0(value, "-", parent.value) # Notre that here we are passing argument to a non dplyr function call
    )) %>%
    select(ids, everything())

  par_info <- datatable %>% dplyr::group_by(parent.value) %>% # group by parent
    dplyr::summarise(count = mean(as.numeric(count), na.rm = T)) %>% # parent total
    rename(value = parent.value) %>% # parent labels for the item field
    mutate(parent.value = "", ids = value) %>% # add missing fields for my_data
    select(names(datatable)) # put cols in same order as my_data

  data_for_plot <- rbind(datatable, par_info)

  return(data_for_plot)
}

if (params$actions$run_fc_treemaps == "TRUE") {
  message("Great ! You decided to launch the fc treemaps calculations :) :\n")
  ############################ version 2
  # Create a data frame
  library(jsonlite)

  # Specify the URL of the JSON file
  url <- "https://raw.githubusercontent.com/mwang87/NP-Classifier/master/Classifier/dict/index_v1.json"

  # Load the JSON file
  json_data <- fromJSON(url)

  npclassifier_origin <- treat_npclassifier_json(json_data)


  # Aggregate rows by concatenating values in superclass and path columns
  npclassifier_newpath <- aggregate(cbind(superclass, pathway) ~ class, data = npclassifier_origin, FUN = function(x) paste(unique(unlist(strsplit(x, " x "))), collapse = " x "))
  colnames(npclassifier_newpath) <- c("canopus_npc_class", "canopus_npc_superclass", "canopus_npc_pathway")
  npclassifier_newpath$canopus_npc_superclass[grep(" x ", npclassifier_newpath$canopus_npc_pathway)] <- paste(npclassifier_newpath$canopus_npc_superclass[grep(" x ", npclassifier_newpath$canopus_npc_pathway)], "x")

  # Alternatively we generate a new df where all class-superclass pairs are distinct and we add a column with the corresponding pathway (we keep the first occurence). We rename the columns (canopus_npc_class = class, canopus_npc_superclass = superclass, canopus_npc_pathway = pathway).We return a data.frame.

  npclassifier_newpath_simple <- npclassifier_origin %>%
    distinct(class, .keep_all = TRUE) %>%
    rename(canopus_npc_class = class, canopus_npc_superclass = superclass, canopus_npc_pathway = pathway) %>%
    na.omit() %>%
    data.frame()



  # Here we list the distinct values in the npclassifier_newpath$canopus_npc_pathway and order them alphabetically

  # npclassifier_newpath  %>%
  # distinct(canopus_npc_pathway)  %>%
  # arrange(canopus_npc_pathway)

  # Check wether this line is used or not ?
  index <- sort(unique(paste(npclassifier_newpath$canopus_npc_superclass, npclassifier_newpath$canopus_npc_pathway)))


  # DE_foldchange_pvalues_signi <- DE_foldchange_pvalues[DE_foldchange_pvalues$C_WT_p_value < 0.05,]


  # Extract prefixes of columns with "_p_value" suffix
  conditions <- sub("_p_value$", "", grep("_p_value$", names(DE_foldchange_pvalues), value = TRUE))

  # Print message before iterating over conditions
  message("Iterating over the following conditions for the treemaps generation:\n")


  # Iterate over the prefixes
  for (condition in conditions) {
    message("Generating treemaps for condition: ", condition)

    # Perform filtering using the prefix as a condition
    DE_foldchange_pvalues_signi <- DE_foldchange_pvalues %>%
      filter(!!sym(paste0(condition, "_p_value")) < 0.05)

    # Print the filtered data
    message("Filtered data for condition: ", condition, ":\n")
    # print(head(DE_foldchange_pvalues_signi))
    message("\n")

    condition_parts <- strsplit(condition, "_vs_")[[1]]
    first_part <- condition_parts[1]
    second_part <- condition_parts[2]

    if (gnps2_job) {
    mydata_meta <- select(
      DE_foldchange_pvalues_signi, "sirius_inchikey2d", "row_id", "sirius_name", "sirius_smiles",
      "gnps_cluster_index", "feature_rt", "feature_mz", "sirius_adduct", "sirius_chebiasciiname", "sirius_chebiid", "sirius_molecularformula", "gnps_component"
    )    } else {
    mydata_meta <- select(
      DE_foldchange_pvalues_signi, "sirius_inchikey2d", "row_id", "sirius_name", "sirius_smiles",
      "gnps_cluster_index", "feature_rt", "feature_mz", "sirius_adduct", "sirius_chebiasciiname", "sirius_chebiid", "sirius_molecularformula", "gnps_componentindex", "gnpslinkout_cluster_gnps", "gnps_libraryid"
    )
    }

    mydata_meta$name_comp <- "unknown"
    mydata_meta$name_comp[!is.na(mydata_meta$sirius_inchikey2d)] <- mydata_meta$sirius_name[!is.na(mydata_meta$sirius_inchikey2d)]



    mydata1 <- select(
      DE_foldchange_pvalues_signi,
      !!sym(paste0(condition, "_fold_change_log2")), "sirius_name", "row_id",
      "canopus_npc_class"
    ) %>%
      # this line remove rows with NA in the canopus_npc_class column using the filter function
      filter(!is.na(canopus_npc_class))

    mydata1 <- merge(mydata1, npclassifier_newpath, by = "canopus_npc_class")


    # mydata1 <- mydata1 %>%
    # mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>%
    # mutate_if(is.numeric, function(x) ifelse(is.nan(x), 0, x))


    mydata1_neg <- mydata1 %>%
      filter(!!sym(paste0(condition, "_fold_change_log2")) < 0)

    mydata1_pos <- mydata1 %>%
      filter(!!sym(paste0(condition, "_fold_change_log2")) >= 0)


    # Check if the data frame has zero rows
    if (nrow(mydata1_pos) == 0) {
      # Recycle the original column names and create a new data frame with zeros
      mydata1_pos <- tibble(
        !!!setNames(rep(0, length(names(mydata1_pos))), names(mydata1_pos))
      )
    } else {
      # Data frame already has rows, no need to fill with zeros
      # You can add additional code here to perform operations on the existing data
    }
    # Check if the data frame has zero rows
    if (nrow(mydata1_neg) == 0) {
      # Recycle the original column names and create a new data frame with zeros
      mydata1_neg <- tibble(
        !!!setNames(rep(0, length(names(mydata1_neg))), names(mydata1_neg))
      )
    } else {
      # Data frame already has rows, no need to fill with zeros
      # You can add additional code here to perform operations on the existing data
    }

    # Aggregate the data
    ####
    # We protect the code with a tryCatch to avoid errors if the data is empty. This can hapen when no classified features are returned fopr a specific condition. This should return an empty treemap. Beware !!!!

    mydata1 <- mydata1[!is.na(mydata1$canopus_npc_pathway), ]
    mydata1$counter <- 1

    # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
    mydata1_neg <- mydata1_neg[!is.na(mydata1_neg$canopus_npc_pathway), ]
    mydata1_neg$counter <- 1
    mydata1_neg$fold_dir <- paste("neg", mydata1_neg$canopus_npc_superclass, sep = "_")
    # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
    mydata1_pos <- mydata1_pos[!is.na(mydata1_pos$canopus_npc_superclass), ]
    mydata1_pos$counter <- 1
    mydata1_pos$fold_dir <- paste("pos", mydata1_pos$canopus_npc_superclass, sep = "_")

    #####################################################################
    #####################################################################


    dt_se_prop_prep_count_tot <- dt_for_treemap(
      datatable = mydata1,
      parent_value = canopus_npc_pathway,
      value = canopus_npc_superclass,
      count = counter
    )


    dt_se_prop_prep_fold_tot <- dt_for_treemap_mean(
      datatable = mydata1,
      parent_value = canopus_npc_pathway,
      value = canopus_npc_superclass,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_tot <- dt_se_prop_prep_fold_tot %>%
      select(-c("value", "parent.value"))
    matt_class_fig_tot <- merge(dt_se_prop_prep_count_tot, dt_se_prop_prep_fold_tot, by = "ids")

    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_pos <- dt_for_treemap(
      datatable = mydata1_pos,
      parent_value = canopus_npc_superclass,
      value = fold_dir,
      count = counter
    )

    dt_se_prop_prep_fold_pos <- dt_for_treemap_mean(
      datatable = mydata1_pos,
      parent_value = canopus_npc_superclass,
      value = fold_dir,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_pos <- dt_se_prop_prep_fold_pos %>%
      select(-c("value", "parent.value"))
    matt_class_fig_pos_dir <- merge(dt_se_prop_prep_count_pos, dt_se_prop_prep_fold_pos, by = "ids")


    matt_class_fig_pos_dir <- matt_class_fig_pos_dir[!(matt_class_fig_pos_dir$parent.value == ""), ]
    matt_class_fig_pos_dir <- na.omit(matt_class_fig_pos_dir)

    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_neg <- dt_for_treemap(
      datatable = mydata1_neg,
      parent_value = canopus_npc_superclass,
      value = fold_dir,
      count = counter
    )

    dt_se_prop_prep_fold_neg <- dt_for_treemap_mean(
      datatable = mydata1_neg,
      parent_value = canopus_npc_superclass,
      value = fold_dir,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_neg <- dt_se_prop_prep_fold_neg %>%
      select(-c("value", "parent.value"))
    matt_class_fig_neg_dir <- merge(dt_se_prop_prep_count_neg, dt_se_prop_prep_fold_neg, by = "ids")

    matt_class_fig_neg_dir <- matt_class_fig_neg_dir[!(matt_class_fig_neg_dir$parent.value == ""), ]
    matt_class_fig_neg_dir <- na.omit(matt_class_fig_neg_dir)

    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_pos_sirius <- dt_for_treemap(
      datatable = mydata1_pos,
      parent_value = fold_dir,
      value = row_id,
      count = counter
    )

    dt_se_prop_prep_fold_pos_sirius <- dt_for_treemap_mean(
      datatable = mydata1_pos,
      parent_value = fold_dir,
      value = row_id,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_pos_sirius <- dt_se_prop_prep_fold_pos_sirius %>%
      select(-c("value", "parent.value"))
    matt_class_fig_pos_dir_sirius <- merge(dt_se_prop_prep_count_pos_sirius, dt_se_prop_prep_fold_pos_sirius, by = "ids")

    matt_class_fig_pos_dir_sirius <- matt_class_fig_pos_dir_sirius[!(matt_class_fig_pos_dir_sirius$parent.value == ""), ]
    matt_class_fig_pos_dir_sirius <- na.omit(matt_class_fig_pos_dir_sirius)


    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_neg_sirius <- dt_for_treemap(
      datatable = mydata1_neg,
      parent_value = fold_dir,
      value = row_id,
      count = counter
    )

    dt_se_prop_prep_fold_neg_sirius <- dt_for_treemap_mean(
      datatable = mydata1_neg,
      parent_value = fold_dir,
      value = row_id,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_neg_sirius <- dt_se_prop_prep_fold_neg_sirius %>%
      select(-c("value", "parent.value"))
    matt_class_fig_neg_dir_sirius <- merge(dt_se_prop_prep_count_neg_sirius, dt_se_prop_prep_fold_neg_sirius, by = "ids")

    matt_class_fig_neg_dir_sirius <- matt_class_fig_neg_dir_sirius[!(matt_class_fig_neg_dir_sirius$parent.value == ""), ]
    matt_class_fig_neg_dir_sirius <- na.omit(matt_class_fig_neg_dir_sirius)



    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################

    matttree <- rbind(matt_class_fig_tot, matt_class_fig_pos_dir, matt_class_fig_neg_dir, matt_class_fig_pos_dir_sirius, matt_class_fig_neg_dir_sirius)
    matttree$labels_adjusted <- matttree$value
    matttree$labels_adjusted[grep("pos_", matttree$labels_adjusted)] <- ""
    matttree$labels_adjusted[grep("neg_", matttree$labels_adjusted)] <- ""
    matttree$labels_adjusted <- gsub(" x", " ", matttree$labels_adjusted)

    # We rename the count.x column as count and the count.y column as foldchange_log2
    matttree <- matttree %>%
      rename(count = count.x) %>%
      rename(foldchange_log2 = count.y)



    matttree <- merge(matttree, mydata_meta, by.x = "labels_adjusted", by.y = "row_id", all.x = T)

    matttree$labels_adjusted[!is.na(matttree$name_comp)] <- matttree$name_comp[!is.na(matttree$name_comp)]
    matttree$value[matttree$labels_adjusted == "unknown"] <- ""
    matttree$value[matttree$labels_adjusted == "unknown"] <- ""


    #####################################################################

    # The follow function creates a new hyperlink column based on the labels_adjusted columns

    # matttree$hl <- paste0("https://en.wikipedia.org/wiki/", matttree$labels_adjusted)

    # # <a href='https://example.com/box1' target='_blank'>Box 1</a>
    # matttree$full_hl <- paste0("<a href='", matttree$hl, "' target='_blank'>", matttree$labels_adjusted, "</a>")
    # matttree$full_hl <- paste0(
    #   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
    # )

    # matttree$hl <- paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", matttree$sirius_inchikey2d, "&sort=annothitcnt")

    # # <a href='https://example.com/box1' target='_blank'>Box 1</a>
    # matttree$full_hl <- paste0(
    #   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
    # )

    # <a href='https://example.com/box1' target='_blank'>Box 1</a>
    matttree$smiles_url <- paste0(
      "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=", matttree$sirius_smiles, "&zoom=2.0&annotate=cip"
    )

    # Generate hl URL only if sirius_inchikey2d is not NA
    matttree$hl <- ifelse(!is.na(matttree$sirius_inchikey2d),
      paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", matttree$sirius_inchikey2d, "&sort=annothitcnt"),
      NA
    )

    # Generate full_hl hyperlink only if hl is not NA
    matttree$full_hl <- paste0(
      "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
    )

    # Generate hl URL only if sirius_inchikey2d is not NA
    matttree$chebi_hl <- ifelse(!is.na(matttree$sirius_chebiid),
      paste0("https://www.ebi.ac.uk/chebi/searchId.do?chebiId=", matttree$sirius_chebiid),
      NA
    )

    # Generate full_hl hyperlink only if hl is not NA
    matttree$chebi_hl_formatted <- ifelse(!is.na(matttree$sirius_chebiid),
      paste0(
        "<a href='", matttree$chebi_hl, "' target='_blank' style='color: black;'>", matttree$sirius_chebiid, "</a>"
      ), ""
    )

    # # Generate full_hl hyperlink only if hl is not NA
    # matttree$gnps_hl_formatted <- ifelse(!is.na(matttree$gnps_cluster_index),
    #   paste0(
    #     "<a href='", matttree$gnpslinkout_cluster_gnps, "' target='_blank' style='color: black;'>", matttree$gnps_cluster_index, "</a>"
    #   ), ""
    # )

    # Generate smiles_url only if sirius_smiles is not NA
    matttree$smiles_url <- ifelse(!is.na(matttree$sirius_smiles),
      paste0("https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=", matttree$sirius_smiles, "&zoom=2.0&annotate=cip"),
      NA
    )
    # Generate clickable smiles_url only if smiles_url is not NA
    matttree$smiles_clickable_url <- ifelse(!is.na(matttree$smiles_url),
      paste0("<a href='", matttree$smiles_url, "' target='_blank' style='color: black;'>", matttree$sirius_smiles, "</a>"),
      NA
    )


    # "sirius_molecularformula", "gnps_componentindex", "gnpslinkout_cluster_gnps", "LibraryID_GNPS"
    # mattree$smiles_clickable_url <- paste0("<a href=", matttree$smiles_url, " target='_blank' rel='noopener noreferrer'>", matttree$sirius_smiles, "</a>")


    # Here we replace all NA by empty cells in the matttree$smiles_clickable_url column

    matttree$smiles_clickable_url[is.na(matttree$smiles_clickable_url)] <- ""
    matttree$sirius_chebiid[is.na(matttree$sirius_chebiid)] <- ""
    matttree$sirius_chebiasciiname[is.na(matttree$sirius_chebiasciiname)] <- ""
    # matttree$gnps_libraryid[is.na(matttree$gnps_libraryid)] <- ""


    # Create a new column in the data frame to store the colors for each value
    matttree$colors <- NA

    # Assign specific colors to the classes
    matttree$colors[matttree$parent.value == "Alkaloids" | matttree$value == "Alkaloids"] <- "#514300"
    matttree$colors[matttree$parent.value == "Alkaloids x Amino acids and Peptides" | matttree$value == "Alkaloids x Amino acids and Peptides"] <- "#715e00"
    matttree$colors[matttree$parent.value == "Alkaloids x Terpenoids" | matttree$value == "Alkaloids x Terpenoids"] <- "#756101"
    matttree$colors[matttree$parent.value == "Amino acids and Peptides" | matttree$value == "Amino acids and Peptides"] <- "#ca5a04"
    matttree$colors[matttree$parent.value == "Amino acids and Peptides x Polyketides" | matttree$value == "Amino acids and Peptides x Polyketides"] <- "#d37f3e"
    matttree$colors[matttree$parent.value == "Amino acids and Peptides x Shikimates and Phenylpropanoids" | matttree$value == "Amino acids and Peptides x Shikimates and Phenylpropanoids"] <- "#ca9f04"
    matttree$colors[matttree$parent.value == "Carbohydrates" | matttree$value == "Carbohydrates"] <- "#485f2f"
    matttree$colors[matttree$parent.value == "Fatty acids" | matttree$value == "Fatty acids"] <- "#612ece"
    matttree$colors[matttree$parent.value == "Polyketides" | matttree$value == "Polyketides"] <- "#865993"
    matttree$colors[matttree$parent.value == "Polyketides x Terpenoids" | matttree$value == "Polyketides x Terpenoids"] <- "#6a5c8a"
    matttree$colors[matttree$parent.value == "Shikimates and Phenylpropanoids" | matttree$value == "Shikimates and Phenylpropanoids"] <- "#6ba148"
    matttree$colors[matttree$parent.value == "Terpenoids" | matttree$value == "Terpenoids"] <- "#63acf5"


    # To check what this is doing
    matttree <- matttree[order(matttree$value), ]


    #########################################################
    #########################################################
    
    if (gnps2_job) {
    txt <- as.character(paste0
    (
      "feature id: ", matttree$gnps_cluster_index, "<br>",
      "component id: ", matttree$gnps_component, "<br>",
      "name: ", matttree$labels_adjusted, "<br>",
      "m/z: ", round(matttree$feature_mz, 4), "<br>",
      "RT: ", round(matttree$feature_rt, 2), "<br>",
      "MF: ", matttree$sirius_molecularformula, "<br>",
      "adduct: ", matttree$sirius_adduct, "<br>",
      "FC (log 2): ", round(matttree$foldchange_log2, 2),
      "<extra></extra>"
    ))
    } else {
    txt <- as.character(paste0
    (
      "feature id: ", matttree$gnps_cluster_index, "<br>",
      "component id: ", matttree$gnps_componentindex, "<br>",
      "name: ", matttree$labels_adjusted, "<br>",
      "m/z: ", round(matttree$feature_mz, 4), "<br>",
      "RT: ", round(matttree$feature_rt, 2), "<br>",
      "MF: ", matttree$sirius_molecularformula, "<br>",
      "adduct: ", matttree$sirius_adduct, "<br>",
      "FC (log 2): ", round(matttree$foldchange_log2, 2),
      "<extra></extra>"
    ))
    }



    matttree$txt <- txt



    fig_treemap_qual <- plot_ly(
      data = matttree,
      type = "treemap",
      ids = ~value,
      labels = ~ paste0("<b>", matttree$full_hl, "</b><br>", matttree$smiles_clickable_url, "<br><b>", matttree$sirius_chebiasciiname, "</b><br>", matttree$chebi_hl_formatted, "<br>", "</a>"),
      parents = ~parent.value,
      values = ~count,
      branchvalues = "total",
      maxdepth = 3,
      hovertemplate = ~txt,
      marker = list(
        colors = matttree$colors # Use the colors column from the data frame
      )
    ) %>%
      layout(
        title = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters: [", filter_sample_metadata_status, "]"),
        margin = list(
          l = 100, # Left margin in pixels, adjust as needed
          r = 100, # Right margin in pixels, adjust as needed
          t = 100, # Top margin in pixels, adjust as needed
          b = 100 # Bottom margin in pixels, adjust as needed
        )
      )

    fig_treemap_quan <- plot_ly(
      data = matttree,
      type = "treemap",
      ids = ~value,
      labels = ~ paste0("<b>", matttree$full_hl, "</b><br>", matttree$smiles_clickable_url, "<br><b>", matttree$sirius_chebiasciiname, "</b><br>", matttree$chebi_hl_formatted, "<br>", "</a>"),
      parents = ~parent.value,
      values = ~count,
      branchvalues = "total",
      maxdepth = 4,
      hovertemplate = ~txt,
      marker = list(
        colors = matttree$foldchange_log2,
        colorscale = list(
          c(0, 0.5, 1),
          c(custom_colors[first_part], "#FFFFFF", custom_colors[second_part])
        ),
        cmin = max(abs(matttree$foldchange_log2)) * (-1),
        cmax = max(abs(matttree$foldchange_log2)),
        showscale = TRUE,
        colorbar = list(
          # the title html is set to add a line return
          title = "",
          tickmode = "array",
          tickvals = c((quantile(abs(matttree$foldchange_log2), probs = 0.75) * (-1)), 0, (quantile(abs(matttree$foldchange_log2), probs = 0.75))),
          ticktext = c(
            paste0("<b>", first_part, "</b>"),
            "",
            paste0("<b>", second_part, "</b>")
          ),
          len = 0.5,
          thickness = 30,
          outlinewidth = 1,
          tickangle = 270
        ),
        reversescale = FALSE # Set to FALSE to maintain the color gradient order
      )
    ) %>%
      layout(
        title = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters: [", filter_sample_metadata_status, "]"),
        margin = list(
          l = 100, # Left margin in pixels, adjust as needed
          r = 100, # Right margin in pixels, adjust as needed
          t = 100, # Top margin in pixels, adjust as needed
          b = 100 # Bottom margin in pixels, adjust as needed
        )
      )

    # We now save the treemap as a html file locally

    if (params$operating_system$system == "unix") {
      ### linux version
      htmlwidgets::saveWidget(fig_treemap_qual, file = paste0("Treemap_", first_part, "_vs_", second_part, "_qual.html"), selfcontained = TRUE) # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_qual.html")

      htmlwidgets::saveWidget(fig_treemap_quan, file = paste0("Treemap_", first_part, "_vs_", second_part, "_quan.html"), selfcontained = TRUE) # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_quan.html")
    }


    if (params$operating_system$system == "windows") {
      ### windows version
      Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
      htmlwidgets::saveWidget(fig_treemap_qual, file = paste0("Treemap_", first_part, "_vs_", second_part, "_qual.html"), selfcontained = TRUE, libdir = "lib") # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_qual.html")
      unlink("lib", recursive = FALSE)

      htmlwidgets::saveWidget(fig_treemap_quan, file = paste0("Treemap_", first_part, "_vs_", second_part, "_quan.html"), selfcontained = TRUE, libdir = "lib") # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_qual.html")
      unlink("lib", recursive = FALSE)
    }
  }
}


# mydata_meta <- select(DE_foldchange_pvalues, "sirius_inchikey2d", "row_id","sirius_name","sirius_smiles",
# "gnps_cluster_index","feature_rt","feature_mz")
# mydata_meta$name_comp<- "unknown"
# mydata_meta$name_comp[!is.na(mydata_meta$sirius_inchikey2d)] <- mydata_meta$sirius_name[!is.na(mydata_meta$sirius_inchikey2d)]




# mydata1 <- select(DE_foldchange_pvalues,
# "C_WT_fold_change_log2","sirius_name", "row_id",
# "canopus_npc_class")  %>%
# # this line remove rows with NA in the canopus_npc_class column using the filter function
# filter(!is.na(canopus_npc_class))  %>%
# filter(!is.na(C_WT_fold_change_log2))


# mydata1 <- merge(mydata1,npclassifier_newpath,by="canopus_npc_class")



# mydata1 <- mydata1 %>%
# mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>%
# mutate_if(is.numeric, function(x) ifelse(is.nan(x), 0, x))


# mydata1_neg <- mydata1[mydata1$C_WT_fold_change_log2 < 0,]
# #mydata1_neg$C_WT_fold_change_log2 <- abs(mydata1_neg$C_WT_fold_change_log2)
# mydata1_pos <- mydata1[mydata1$C_WT_fold_change_log2 >= 0,]

# # Aggregate the data
# ####
# mydata1 = mydata1[!is.na(mydata1$canopus_npc_pathway), ]
# mydata1$counter = 1

# # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
# mydata1_neg = mydata1_neg[!is.na(mydata1_neg$canopus_npc_pathway), ]
# mydata1_neg$counter = 1
# mydata1_neg$fold_dir <- paste("neg",mydata1_neg$canopus_npc_superclass,sep="_")
# # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
# mydata1_pos = mydata1_pos[!is.na(mydata1_pos$canopus_npc_superclass), ]
# mydata1_pos$counter = 1
# mydata1_pos$fold_dir <- paste("pos",mydata1_pos$canopus_npc_superclass,sep="_")

# #####################################################################
# #####################################################################

# dt_se_prop_prep_count_tot = dt_for_treemap(
#   datatable = mydata1,
#   parent_value = canopus_npc_pathway,
#   value = canopus_npc_superclass,
#   count = counter
# )


# dt_se_prop_prep_fold_tot = dt_for_treemap_mean(
#   datatable = mydata1,
#   parent_value = canopus_npc_pathway,
#   value = canopus_npc_superclass,
#   count = C_WT_fold_change_log2
# )

# dt_se_prop_prep_fold_tot <- dt_se_prop_prep_fold_tot %>%
# select(-c("value","parent.value"))
# matt_class_fig_tot <- merge(dt_se_prop_prep_count_tot,dt_se_prop_prep_fold_tot,by="ids")

# #####################################################################
# #####################################################################
# #####################################################################
# #####################################################################

# dt_se_prop_prep_count_pos = dt_for_treemap(
#   datatable = mydata1_pos,
#   parent_value = canopus_npc_superclass,
#   value = fold_dir,
#   count = counter
# )

# dt_se_prop_prep_fold_pos = dt_for_treemap_mean(
#   datatable = mydata1_pos,
#   parent_value = canopus_npc_superclass,
#   value = fold_dir,
#   count = C_WT_fold_change_log2
# )

# dt_se_prop_prep_fold_pos <- dt_se_prop_prep_fold_pos %>%
# select(-c("value","parent.value"))
# matt_class_fig_pos_dir <- merge(dt_se_prop_prep_count_pos,dt_se_prop_prep_fold_pos,by="ids")

# matt_class_fig_pos_dir <- matt_class_fig_pos_dir[!(matt_class_fig_pos_dir$parent.value == ""),]
# matt_class_fig_pos_dir <- na.omit(matt_class_fig_pos_dir)

# #####################################################################
# #####################################################################

# dt_se_prop_prep_count_neg = dt_for_treemap(
#   datatable = mydata1_neg,
#   parent_value = canopus_npc_superclass,
#   value = fold_dir,
#   count = counter
# )

# dt_se_prop_prep_fold_neg = dt_for_treemap_mean(
#   datatable = mydata1_neg,
#   parent_value = canopus_npc_superclass,
#   value = fold_dir,
#   count = C_WT_fold_change_log2
# )

# dt_se_prop_prep_fold_neg <- dt_se_prop_prep_fold_neg %>%
# select(-c("value","parent.value"))
# matt_class_fig_neg_dir <- merge(dt_se_prop_prep_count_neg,dt_se_prop_prep_fold_neg,by="ids")

# matt_class_fig_neg_dir <- matt_class_fig_neg_dir[!(matt_class_fig_neg_dir$parent.value == ""),]
# matt_class_fig_neg_dir <- na.omit(matt_class_fig_neg_dir)

# #####################################################################
# #####################################################################
# #####################################################################
# #####################################################################

# dt_se_prop_prep_count_pos_sirius = dt_for_treemap(
#   datatable = mydata1_pos,
#   parent_value = fold_dir,
#   value = row_id,
#   count = counter
# )

# dt_se_prop_prep_fold_pos_sirius = dt_for_treemap_mean(
#   datatable = mydata1_pos,
#   parent_value = fold_dir,
#   value = row_id,
#   count = C_WT_fold_change_log2
# )

# dt_se_prop_prep_fold_pos_sirius <- dt_se_prop_prep_fold_pos_sirius %>%
# select(-c("value","parent.value"))
# matt_class_fig_pos_dir_sirius <- merge(dt_se_prop_prep_count_pos_sirius,dt_se_prop_prep_fold_pos_sirius,by="ids")

# matt_class_fig_pos_dir_sirius <- matt_class_fig_pos_dir_sirius[!(matt_class_fig_pos_dir_sirius$parent.value == ""),]
# matt_class_fig_pos_dir_sirius <- na.omit(matt_class_fig_pos_dir_sirius)


# #####################################################################
# #####################################################################

# dt_se_prop_prep_count_neg_sirius = dt_for_treemap(
#   datatable = mydata1_neg,
#   parent_value = fold_dir,
#   value = row_id,
#   count = counter
# )

# dt_se_prop_prep_fold_neg_sirius = dt_for_treemap_mean(
#   datatable = mydata1_neg,
#   parent_value = fold_dir,
#   value = row_id,
#   count = C_WT_fold_change_log2
# )

# dt_se_prop_prep_fold_neg_sirius <- dt_se_prop_prep_fold_neg_sirius %>%
# select(-c("value","parent.value"))
# matt_class_fig_neg_dir_sirius <- merge(dt_se_prop_prep_count_neg_sirius,dt_se_prop_prep_fold_neg_sirius,by="ids")

# matt_class_fig_neg_dir_sirius <- matt_class_fig_neg_dir_sirius[!(matt_class_fig_neg_dir_sirius$parent.value == ""),]
# matt_class_fig_neg_dir_sirius <- na.omit(matt_class_fig_neg_dir_sirius)



# #####################################################################
# #####################################################################
# #####################################################################
# #####################################################################

# matttree <- rbind(matt_class_fig_tot,matt_class_fig_pos_dir,matt_class_fig_neg_dir,matt_class_fig_pos_dir_sirius,matt_class_fig_neg_dir_sirius)
# matttree$labels_adjusted <- matttree$value
# matttree$labels_adjusted[grep("pos_",matttree$labels_adjusted)] <- "+"
# matttree$labels_adjusted[grep("neg_",matttree$labels_adjusted)] <- "-"
# matttree$labels_adjusted <- gsub(" x"," ",matttree$labels_adjusted)

# matttree <- merge(matttree,mydata_meta,by.x="labels_adjusted",by.y="row_id",all.x =T)

# matttree$labels_adjusted[!is.na(matttree$name_comp)] <- matttree$name_comp[!is.na(matttree$name_comp)]
# matttree$value[matttree$labels_adjusted== "unknown"] <- ""
# matttree$value[matttree$labels_adjusted== "unknown"] <- ""
# #####################################################################

# # The follow function creates a new hyperlink column based on the labels_adjusted columns

# # matttree$hl <- paste0("https://en.wikipedia.org/wiki/", matttree$labels_adjusted)

# # # <a href='https://example.com/box1' target='_blank'>Box 1</a>
# # matttree$full_hl <- paste0("<a href='", matttree$hl, "' target='_blank'>", matttree$labels_adjusted, "</a>")
# # matttree$full_hl <- paste0(
# #   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
# # )

# matttree$hl <- paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", matttree$sirius_inchikey2d, "&sort=annothitcnt")

# # <a href='https://example.com/box1' target='_blank'>Box 1</a>
# matttree$full_hl <- paste0(
#   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
# )

# # <a href='https://example.com/box1' target='_blank'>Box 1</a>
# matttree$smiles_url <- paste0(
#   "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=", matttree$sirius_smiles, "&zoom=2.0&annotate=cip"
# )


# matttree <- matttree[order(matttree$value),]
# ####################### annoation structure
# #########################################################

# ############# ad mol
# #mol_list_2d <- list()
# #for (i in c(1:length(matttree$sirius_smiles))) {
# #  psml <- parse.smiles(matttree$sirius_smiles[i], omit.nulls = TRUE)
# #  if (length(psml) == 0) {
# #    psml <- parse.smiles("C")
# #  }

# #  img <- view.image.2d(psml[1][[1]])

# #  r <- as.raster(img)
# #  mol_list_2d[[i]] <- r
# #}

# #########################################################
# #########################################################

# txt <- as.character(paste0
# ("feature id: ",matttree$gnps_cluster_index,"<br>",
#  "RT: ", round(matttree$feature_rt,2),"<br>",
#  "m/z: ", round(matttree$feature_mz,4),
#  "<extra></extra>"
# ))
# matttree$txt <- txt

# fig_treemap = plot_ly(
#   data = matttree,
#   type = "treemap",
#   ids = ~value,
#   labels = ~matttree$full_hl,
#   parents = ~parent.value,
#   values = ~count.x,
#   branchvalues = "total",
#   maxdepth=3,
#   hovertemplate = ~txt
# )



# # d3 <- htmltools::htmlDependency(
# #   "d3", "7.3",
# #   src = c(href = "https://cdnjs.cloudflare.com/ajax/libs/d3/7.3.0/"),
# #   script = "d3.min.js"
# # )

# # p = plot_ly(
# #   data = matttree,
# #   type = "treemap",
# #   ids = ~value,
# #   labels = ~matttree$full_hl,
# #   parents = ~parent.value,
# #   values = ~count.x,
# #   branchvalues = "total",
# #   maxdepth=3
# # ) %>%
# # add_text(x = matttree$value, y = matttree$value, customdata = ~matttree$smiles_url, text = ~matttree$value) %>%
# # htmlwidgets::onRender(readLines("/Users/pma/Dropbox/git_repos/mapp-metabolomics-unit/biostat_toolbox/hover_tooltip.js"))

# # p$dependencies <- c(p$dependencies, list(d3))
# # p


# # library(htmlwidgets)
# # library(magrittr)
# # library(plotly)

# # x <- 1:3
# # y <- 1:3

# # artists <- c("Bethoven", "Mozart", "Bach")

# # image_links <- c(
# #   "https://upload.wikimedia.org/wikipedia/commons/6/6f/Beethoven.jpg",
# #   "https://upload.wikimedia.org/wikipedia/commons/4/47/Croce-Mozart-Detail.jpg",
# #   "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=CC1C(C(C(C(%3DO)C(CC(C(C(C(C(C(%3DO)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O&zoom=2.0&annotate=cip"
# # )


# # x <- 1:3
# # y <- 1:3

# # # hoverinfo = "none" will hide the plotly.js tooltip, but the
# # # plotly_hover event will still fire
# # p <- plot_ly(hoverinfo = "none") %>%
# #   add_text(x = x, y = y, customdata = image_links, text = artists) %>%
# #   htmlwidgets::onRender(readLines("/Users/pma/Dropbox/git_repos/mapp-metabolomics-unit/biostat_toolbox/hover_tooltip.js"))

# # p$dependencies <- c(p$dependencies, list(d3))
# # p

# ## ChatGPT adapted prompts for treemap
# ######################################
# ######################################

# # usePackage('highcharter')


# # library(highcharter)

# # # Sample data
# # labels <- c("Beethoven", "Mozart", "Chemical Structure")
# # sizes <- c(20, 30, 50)
# # urls <- c(
# #   "https://upload.wikimedia.org/wikipedia/commons/6/6f/Beethoven.jpg",
# #   "https://upload.wikimedia.org/wikipedia/commons/4/47/Croce-Mozart-Detail.jpg",
# #   "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=CC1C(C(C(C(%3DO)C(CC(C(C(C(C(C(%3DO)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O&zoom=2.0&annotate=cip"
# # )

# # # Create the data frame
# # data <- data.frame(
# #   name = labels,
# #   value = sizes,
# #   link = urls,
# #   stringsAsFactors = FALSE
# # )

# # # Create the treemap using highchart
# # treemap <- highchart() %>%
# #   hc_chart(type = "treemap") %>%
# #   hc_title(text = "Treemap with Hover Images") %>%
# #   hc_tooltip(
# #     useHTML = TRUE,
# #     pointFormat = "<b>{point.name}</b><br/><img src='{point.link}' width='100' height='100'>"
# #   ) %>%
# #   hc_plotOptions(
# #     series = list(
# #       borderWidth = 0,
# #       dataLabels = list(enabled = FALSE),
# #       cursor = "pointer"
# #     )
# #   ) %>%
# #   hc_series(
# #     data = list_parse(data),
# #     levels = list(
# #       list(level = 1, layoutAlgorithm = "stripes"),
# #       list(level = 2)
# #     )
# #   )

# # # Display the treemap
# # treemap


# # htmlwidgets::saveWidget(widget, file = "fig_treemap.html", selfcontained = TRUE)





# # library(treemap)
# # library(gridSVG)
# # library(XML
# # df <- data.frame(
# #     character=c("Homer","Marge", "Bart", "Lisa","Maggie", "Moe", "Blinky","Bumblebee Man","Duffman","Maude Flanders","Ned Flanders","Rod Flanders","Todd","Jimbo","Otto Mann","Snowball")
# #     ,score=c(268,267,495, 432, 219, 373, 152, 356, 461, 116,107,165, 305,228, 461, 608)

# # tm <- treemap( df, index = "character", vSize = "score"
# # svg <- grid.export()$sv
# # # see http://stackoverflow.com/questions/10688516/fill-svg-path-with-a-background-image-without-knowing-heightwidth?rq=1
# # pattern <- newXMLNode(
# #     "defs"
# #     ,.children = list(
# #         newXMLNode(
# #             "pattern"
# #             , attrs = c(
# #                 id = "img_homer"
# #                 ,patternUnits="userSpaceOnUse"
# #                 ,patternTransform="translate(0, 0) scale(1, -1) rotate(0)"
# #                 ,width="106"
# #                 ,height="98"
# #             )
# #             , .children = newXMLNode(
# #                 "image"
# #                 , attrs = c(
# #                     "xlink:href" = "http://i.imgur.com/JP4s21O.jpg"
# #                     ,width = 106
# #                     ,height = 80
# #                 )
# #             )
# #         )
# #     )

# # addChildren( svg, pattern
# # homer <- getNodeSet(
# #     getNodeSet( svg, "//*[contains(@id,'data.2')]")[[1]]
# #     ,"//*[local-name()='rect']"
# # )[[5]
# # homer_attrs <- xmlAttrs(homer)
# # homer_attrs[["fill"]] <- "url(#img_homer)"
# # xmlAttrs(homer) <- homer_attr
# # library(htmltools)
# # browsable(HTML(saveXML(svg)))





# # g <- ggplot(iris, aes(x = Sepal.Length,
# #                       y = Petal.Length,
# #                       color = Species,
# #                       text = Species)) + geom_point()
# # p <- ggplotly(g, tooltip = "text") %>% partial_bundle()


# # p %>% htmlwidgets::onRender(readLines("/Users/pma/Dropbox/git_repos/mapp-metabolomics-unit/biostat_toolbox/hover_tooltip_adapted.js"))

# fig_treemap <- plot_ly(
#   data = matttree,
#   type = "treemap",
#   ids = ~value,
#   labels = ~matttree$full_hl,
#   parents = ~parent.value,
#   values = ~count.x,
#   branchvalues = "total",
#   maxdepth = 3,
#   marker = list(
#     colors = matttree$count.y,
#     colorscale = list(
#       c(0, 0.5, 1),
#       c("#A89639", "#FFFFFF", "#337AB7")),
#     cmin = max(abs(matttree$count.y)) * (-1),
#     cmax = max(abs(matttree$count.y)),
#     showscale = TRUE,
#     colorbar = list(
#     # the title html is set to add a line return
#      title = '<b>Increased in </b> <br> <br> <br>',
#       tickmode = "array",
#       tickvals = c((quantile(abs(matttree$count.y),probs=0.5)*(-1)), 0, (quantile(abs(matttree$count.y),probs=0.5))),
#       ticktext = c("wild-type","", "Control"),
#       len = 0.5,
#       thickness = 30,
#       outlinewidth = 1,
#       tickangle = 270
#     ),
#     reversescale = FALSE  # Set to FALSE to maintain the color gradient order
#   )
# )%>%
#   layout(
#     title = "<b>Metabolomic Variations Across Treatments</b>"
#   )

# fig_treemap

# params$posthoc$p_value


# # We now save the treempa as a html file locally

# htmlwidgets::saveWidget(p, file = "fig_treemap.html", selfcontained = TRUE)


#############################################################################
#############################################################################
############## Tree Map #####################################################
#############################################################################
#############################################################################

message("Preparing Tree Map ...")

# glimpse(DE_foldchange_pvalues)

# Here we select the features that are significant
# for this we filter for values above the p_value threshold in the column selected using the `p_value_column` variable
# We use the dplyr and pipes syntax to do this
# Note the as.symbol() function to convert the string to a symbol As per https://stackoverflow.com/a/48219802/4908629

matt_donust <- DE_foldchange_pvalues %>%
  filter(if_any(ends_with("_p_value"), ~ .x < params$posthoc$p_value))

# matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
matt_donust2 <- matt_donust[!is.na(matt_donust$canopus_npc_superclass), ]
matt_donust2$counter <- 1



dt_for_treemap <- function(datatable, parent_value, value, count) {
  parent_value <- enquo(parent_value)
  value <- enquo(value)
  count <- enquo(count)

  datatable <- data.frame(datatable %>%
    group_by(!!parent_value, !!value, ) %>%
    summarise(count = sum(as.numeric(!!count))))

  datatable <- datatable %>%
    select(!!parent_value, !!value, count) %>% # create id labels for each row # Notre the !! to pass aruguments to a dplyr function
    rename(
      parent.value = !!parent_value,
      value = !!value
    ) %>%
    mutate(ids = ifelse(parent.value == "", value,
      paste0(value, "-", parent.value) # Notre that here we are passing argument to a non dplyr function call
    )) %>%
    select(ids, everything())

  par_info <- datatable %>% dplyr::group_by(parent.value) %>% # group by parent
    dplyr::summarise(count = sum(as.numeric(count))) %>% # parent total
    rename(value = parent.value) %>% # parent labels for the item field
    mutate(parent.value = "", ids = value) %>% # add missing fields for my_data
    select(names(datatable)) # put cols in same order as my_data

  data_for_plot <- rbind(datatable, par_info)

  return(data_for_plot)
}

dt_se_prop_prep_tm <- dt_for_treemap(
  datatable = matt_donust2,
  parent_value = canopus_npc_superclass,
  value = canopus_npc_class,
  count = counter
)


fig_treemap <- plot_ly(
  data = dt_se_prop_prep_tm,
  type = "treemap",
  labels = ~value,
  parents = ~parent.value,
  values = ~count,
  branchvalues = "total"
)

# Why "significant ? According to what ?

fig_treemap <- fig_treemap %>%
  layout(title = list(text = title_treemap, y = 0.02))


# The files is exported
# The title should be updated !!!



if (params$operating_system$system == "unix") {
  ### linux version
  fig_treemap %>%
    htmlwidgets::saveWidget(file = filename_treemap, selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
  ### windows version
  Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
  fig_treemap %>%
    htmlwidgets::saveWidget(file = filename_treemap, selfcontained = TRUE, libdir = "lib")
  unlink("lib", recursive = FALSE)
}


#############################################################################
#############################################################################
############## Random Forest ################################################
#############################################################################
#############################################################################

message("Launching Random Forest calculations ...")



# Here we traduce to fit Manu's inputs ... to be updated later

features_of_importance <- DE_foldchange_pvalues %>%
  filter((!!as.symbol(p_value_column)) < params$posthoc$p_value) %>%
  select(feature_id) %>%
  # we output the data as a vector
  pull()


#  We select all columns except the params$target$sample_metadata_header columns in
# data_subset_norm_rf_filter and we prefix the column names with an X.
# We use the dplyr syntax to do this and the rename function to rename the columns
# We then subset the data to keep only the columns that are in the imp_filter1 variable

data_subset_for_RF <- DE$data %>%
  select(all_of(as.character(features_of_importance))) %>%
  rename_all(~ paste0("X", .)) %>%
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names") %>%
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names) %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X")) %>%
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header)))

# We define the formula externally to inject the external variable # params$target$sample_metadata_header

formula <- as.formula(paste0(params$target$sample_metadata_header, " ~ ."))

# We launch the rfPermute function

data.rp <- rfPermute(formula, data = data_subset_for_RF, na.action = na.omit, ntree = 500, num.rep = 500)

imp_table_rf <- data.frame(data.rp$pval)
imp_table_rf <- importance(data.rp)
imp_table_rf <- data.frame(imp_table_rf)


sink(filename_random_forest_model)
summary(data.rp)
# f = plotImportance(data.rp, plot.type = "bar", plot = FALSE)

sink()

########### plot importance
#
sorted_indices <- order(-imp_table_rf$MeanDecreaseGini)
# Load the required libraries
# Sort the data based on MeanDecreaseGini
imp_table_rf <- imp_table_rf[sorted_indices, ]

# Create the plotly bar plot
fig_rf <- plot_ly(
  data = imp_table_rf,
  x = ~MeanDecreaseGini,
  y = ~ reorder(row.names(imp_table_rf), -MeanDecreaseGini, decreasing = TRUE), # Use reorder to maintain sorting order
  type = "bar",
  orientation = "h"
) %>%
  layout(
    title = title_random_forest,
    xaxis = list(title = "Importance", tickfont = list(size = 12)), # Adjust the label size here (e.g., size = 12)
    yaxis = list(title = "Features", tickfont = list(size = 5)), # Adjust the label size here (e.g., size = 10)
    margin = list(l = 100, r = 20, t = 50, b = 70),
    showlegend = FALSE
  )
fig_rf
# The file is exported
# The title should be updated !!!


if (params$operating_system$system == "unix") {
  ### linux version
  fig_rf %>%
    htmlwidgets::saveWidget(file = filename_random_forest, selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
  ### windows version
  Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
  fig_rf %>%
    htmlwidgets::saveWidget(file = filename_random_forest, selfcontained = TRUE, libdir = "lib")
  unlink("lib", recursive = FALSE)
}



#############################################################################
#############################################################################
############## Random Forest selected Box Plots #############################
#############################################################################
#############################################################################

# message("Preparing RF-selected Box plots ...")


# imp.scaled = rfPermute::importance(data.rp, scale = TRUE)
# imp.scaled = data.frame(imp.scaled)
# imp.scaled = imp.scaled[order(imp.scaled$MeanDecreaseGini, decreasing = TRUE), ]
# # boxplot_top_N = row.names(imp.scaled)[1:params$boxplot$topN]
# # This below using head is safer in case the number of features is smaller than the number of features to plot
# boxplot_top_N = row.names(head(imp.scaled, n = params$boxplot$topN))

# data_subset_boxplot = data_subset_for_RF %>%
#   select(all_of(boxplot_top_N), params$target$sample_metadata_header)

# # We now establish a side by side box plot for each columns of the data_subset_norm_boxplot
# # We use the melt function to reshape the data to a long format
# # We then use the ggplot2 syntax to plot the data and the facet_wrap function to plot the data side by side

# # Gather value columns into key-value pairs
# df_long <- tidyr::gather(data_subset_boxplot, key = "variable", value = "value", -params$target$sample_metadata_header)


# # Create boxplots faceted by variable and colored by age
# # Note how we use the get function to access the variable name
# # See here for details  https://stackoverflow.com/a/22309328/4908629
# p = ggplot(df_long, aes(x = get(params$target$sample_metadata_header), y = value, fill = get(params$target$sample_metadata_header))) +
#   geom_boxplot() +
#   facet_wrap(~ variable, ncol = 4) +
#   theme_minimal()+
#   ggtitle(title_box_plots)


# fig_boxplot = p + facet_wrap(~variable, scales = "free", dir = "v") + theme(
#   legend.position = "top",
#   legend.title = element_blank()
# )

# # fig_boxplotly = data_subset_boxplot %>%
# #   split(data_subset_boxplot$variable) %>%
# #   map(~ {
# #     plot_ly(data = .x, x = .x$treatment, y = .x$value, color = .x$treatment, type = "box") %>%
# #       layout(showlegend = F, xaxis = list(title = .x$variable[1])) %>%
# #       layout(xaxis = list(titlefont = list(size = 10), tickfont = list(size = 10))) %>%
# #       layout(yaxis = list(titlefont = list(size = 12), tickfont = list(size = 12)))
# #   }) %>%
# #   subplot(margin = 0.02, nrows = 4, titleX = TRUE)  %>%
# #   layout(title = title_box_plots,
# #   legend = list(title = list(text = paste("<b>class </b>")))) # Find a way to define the legend title

# # The files are exported

# ggsave(plot = fig_boxplot, filename = filename_box_plots, width = 10, height = 10)
# # fig_boxplotly %>%
# #     htmlwidgets::saveWidget(file = filename_box_plots_interactive, selfcontained = TRUE)



#############################################################################
#############################################################################
############## p-Value selected Box Plots #############################
#############################################################################
#############################################################################

message("Preparing p-value selected Box plots ...")

features_of_importance_boxplots <- DE_foldchange_pvalues %>%
  # we keep only the features that have a p-value lower than the threshold
  # filter((!!as.symbol(p_value_column)) < params$posthoc$p_value)  %>%
  # We keep only the lowest top n = params$boxplot$topN in the p_value_column
  top_n(-params$boxplot$topN, !!as.symbol(p_value_column)) %>%
  # we order the features by increasing p-value
  arrange(!!as.symbol(p_value_column)) %>%
  select(feature_id) %>%
  # we output the data as a vector
  pull()


data_subset_for_boxplots <- DE$data %>%
  select(all_of(as.character(features_of_importance_boxplots))) %>%
  rename_all(~ paste0("X", .)) %>%
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names") %>%
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names) %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X")) %>%
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header))) %>%
  # Finally we remove the X from the columns names
  rename_all(~ gsub("X", "", .))

# We now establish a side by side box plot for each columns of the data_subset_norm_boxplot
# We use the melt function to reshape the data to a long format
# We then use the ggplot2 syntax to plot the data and the facet_wrap function to plot the data side by side

# Gather value columns into key-value pairs
df_long <- tidyr::gather(data_subset_for_boxplots, key = "variable", value = "value", -params$target$sample_metadata_header)

# Here we merge the df_long with the DE$variable_meta data frame to get the variable type

df_long_informed <- merge(df_long, DE_foldchange_pvalues, by.x = "variable", by.y = "feature_id")


p <- ggplot(df_long_informed, aes(x = !!sym(params$target$sample_metadata_header), y = value, fill = !!sym(params$target$sample_metadata_header))) +
  geom_boxplot() +
  facet_wrap(~feature_id_full_annotated, ncol = 4) +
  # theme_minimal()+
  ggtitle(title_box_plots) +
  geom_point(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.3) # Add data points with jitter for better visibility


ridiculous_strips <- strip_themed(
  # Horizontal strips
  background_x = elem_list_rect(),
  text_x = elem_list_text(face = c("bold", "italic")),
  by_layer_x = TRUE,
  # Vertical strips
  background_y = elem_list_rect(
    fill = c("gold", "tomato", "deepskyblue")
  ),
  text_y = elem_list_text(angle = c(0, 90)),
  by_layer_y = FALSE
)

fig_boxplot <- p + facet_wrap2(~ sirius_chebiasciiname + feature_id_full, labeller = label_value, strip = ridiculous_strips) + theme(
  legend.position = "top",
  legend.title = element_blank()
)


fig_boxplot <- fig_boxplot +
  scale_fill_manual(name = "Groups", values = custom_colors)


# Display the modified plot
print(fig_boxplot)

# The files are exported

ggsave(plot = fig_boxplot, filename = filename_box_plots, width = 10, height = 10)


####
# We now create individual box plots for each selected variable


output_directory_bp <- "./selected_boxplots/"

# Create the directory if it doesn't exist
if (!dir.exists(output_directory_bp)) {
  dir.create(output_directory_bp, recursive = TRUE)
}

# Create and save individual box plots for each selected variable
for (var in features_of_importance_boxplots) {
  # Filter data for the current variable using dplyr
  data_for_plot <- df_long_informed %>%
    filter(variable == var)


  # Round the p-value to 5 digits
  rounded_p_value <- round(pull(data_for_plot, !!as.name(p_value_column)), 5)

  # Create the plot for the current variable (simple box plot)
  p <- ggplot(data_for_plot, aes(x = !!sym(params$target$sample_metadata_header), y = value, fill = !!sym(params$target$sample_metadata_header))) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) + # Add data points with jitter for better visibility
    # ggtitle(paste("Box Plot for", "\n",
    # "Compound name: ", data_for_plot$sirius_chebiasciiname[1], "\n",
    # "Feature details: ", data_for_plot$feature_id_full[1]))
    labs(
      x = params$target$sample_metadata_header,
      y = "Normalized Intensity",
      title = paste("Compared intensities for feature:", var),
      subtitle = paste(
        "\n",
        "Compound name: ", data_for_plot$sirius_chebiasciiname[1], "\n",
        "Feature details: ", data_for_plot$feature_id_full[1]
      ),
      caption = paste("Calculated p-value is ~ ", rounded_p_value)
    ) +
    theme(
      plot.caption = element_text(hjust = 0, face = "italic"), # Default is hjust=1
      plot.title.position = "plot", # NEW parameter. Apply for subtitle too.
      plot.caption.position = "plot"
    ) # NEW parameter


  p <- p +
    scale_fill_manual(name = "Groups", values = custom_colors)


  # Save the plot to a file with a unique filename for each variable
  filename <- paste(output_directory_bp, "boxplot_", gsub(" ", "_", var), ".png", sep = "")
  ggsave(plot = p, filename = filename, width = 8, height = 8)
}

#############################################################################
############## Pvalue filtered Heat Map  #############################
#############################################################################
#############################################################################

message("Preparing p-value filtered Heatmap ...")

features_of_importance <- DE_foldchange_pvalues %>%
  filter((!!as.symbol(p_value_column)) < params$posthoc$p_value) %>%
  select(feature_id) %>%
  # we output the data as a vector
  pull()

data_subset_for_pval_hm <- DE$data %>%
  select(all_of(as.character(features_of_importance))) %>%
  rename_all(~ paste0("X", .)) %>%
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names") %>%
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names) %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X")) %>%
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header))) %>%
  # Finally we remove the X from the columns names
  rename_all(~ gsub("X", "", .))

data_subset_for_pval_hm_peak_height <- DE_original$data %>%
  select(all_of(as.character(features_of_importance))) %>%
  rename_all(~ paste0("X", .)) %>%
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names") %>%
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names) %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X")) %>%
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header))) %>%
  # Finally we remove the X from the columns names
  rename_all(~ gsub("X", "", .))


# data_subset_for_pval_hm_sel = data_subset_for_pval_hm %>%
#   select(params$target$sample_metadata_header)

data_subset_for_pval_hm <- data_subset_for_pval_hm[, colnames(data_subset_for_pval_hm) %in% features_of_importance]

data_subset_for_pval_hm_peak_height <- data_subset_for_pval_hm_peak_height[, colnames(data_subset_for_pval_hm_peak_height) %in% features_of_importance]

# my_sample_col = DE$sample_meta$sample_id

# data_subset_for_Pval = data_subset_for_Pval[, colnames(data_subset_for_Pval) %in% imp_filter2X]
# # my_sample_col = DE$sample_meta$sample_id

my_sample_col <- paste(DE$sample_meta$sample_id, DE$sample_meta[[params$target$sample_metadata_header]], sep = "_")


# We filter the annotation table (DE$variable_meta) to keep only the features of interest identified in the (features_of_importance). We use dplyr

selected_variable_meta <- DE$variable_meta %>%
  filter(feature_id %in% features_of_importance)
# %>%
# select(feature_id, canopus_npc_pathway, canopus_npc_superclass) %>%
# mutate(canopus_npc_pathway = paste(canopus_npc_pathway, canopus_npc_superclass, sep = "_")) %>%
# select(feature_id, canopus_npc_pathway) %>%
# column_to_rownames("feature_id")

selected_variable_meta_NPC <- DE$variable_meta %>%
  filter(feature_id %in% features_of_importance) %>%
  select(feature_id, canopus_npc_superclass, canopus_npc_pathway, canopus_npc_class) %>%
  mutate(NPC.superclass_merged_canopus = paste(canopus_npc_pathway, canopus_npc_superclass, sep = "_")) %>%
  mutate(NPC.class_merged_canopus = paste(NPC.superclass_merged_canopus, canopus_npc_class, sep = "_")) %>%
  select(NPC.class_merged_canopus, NPC.superclass_merged_canopus, canopus_npc_pathway)

selected_variable_meta_NPC_simple <- DE$variable_meta %>%
  filter(feature_id %in% features_of_importance) %>%
  select(canopus_npc_class, canopus_npc_superclass, canopus_npc_pathway)

selected_variable_meta_NPC_simple_ordered <- DE$variable_meta %>%
  filter(feature_id %in% features_of_importance) %>%
  select(canopus_npc_pathway, canopus_npc_superclass, canopus_npc_class)

selected_variable_meta_NPC_simple <- DE$variable_meta %>%
  filter(feature_id %in% features_of_importance) %>%
  select(canopus_npc_class, canopus_npc_superclass, canopus_npc_pathway)



npclassifier_origin_ordered <- npclassifier_origin %>%
  select(pathway, superclass, class)

# ByPal = colorRampPalette(c(wes_palette("Zissou1")))

# data_subset_for_Pval = apply(data_subset_for_Pval, 2, as.numeric)
# # heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")


data_subset_for_pval_hm_mat <- apply(data_subset_for_pval_hm, 2, as.numeric)
# heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")
data_subset_for_pval_hm_mat <- data_subset_for_pval_hm
data_subset_for_pval_hm_mat[] <- lapply(data_subset_for_pval_hm_mat, as.numeric)


# heatmap_filtered_pval = heatmaply(
#   percentize(data_subset_for_pval_hm),
#   seriate = "none", # none , GW , mean, OLO
#   col_side_colors = data.frame(selected_variable_meta_NPC, check.names = FALSE),
#   col_side_palette = ByPal,
#   row_side_colors = data_subset_for_pval_hm_sel,
#   # row_side_palette = ByPal,
#   labRow = as.vector(as.character(my_sample_col)), # [vec_plot]
#   labCol = selected_variable_meta$feature_id_full_annotated,
#   subplot_margin = 0.01,
#   scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
#     low = "lightsteelblue2",
#     # mid = "goldenrod1",
#     high = "firebrick3",
#     midpoint = 0.5,
#     limits = c(0, 1),
#     position = "left"
#   ),
#   hide_colorbar = TRUE,
#   branches_lwd = 0.3,
#   k_row = 4,
#   distfun_row = "pearson",
#   distfun_col = "pearson",
#   fontsize_row = 9,
#   fontsize_col = 9,
#   Rowv = FALSE,
#   side_color_colorbar_len = 1
#   # # # # # Colv = NULL,
#   # plot_method = "plotly"
# )  %>% layout(
#   title = list(text = title_heatmap_pval, font = list(size = 14), x = 0.1),
#   margin = list(t = 150, b = 20) # Adjust the top margin value (e.g., 80) to move the title to the top
# )

# heatmap_filtered_pval



# # The file is exported


# if (params$operating_system$system == "unix") {
# ### linux version

# heatmap_filtered_pval %>%
#     htmlwidgets::saveWidget(file = filename_heatmap_pval, selfcontained = TRUE)
# }

# if (params$operating_system$system == "windows") {
# ### windows version
# Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
# heatmap_filtered_pval %>%
#     htmlwidgets::saveWidget(file = filename_heatmap_pval, selfcontained = TRUE,libdir = "lib")
# unlink("lib", recursive = FALSE)

# }

#### Iheatmapr


target_metadata <- as.factor(DE$sample_meta[[params$target$sample_metadata_header]])


##########################


custom_colors_heatmap <- custom_colors


# Define the vector of colors
micro_cvd_gray <- rev(c(microshades_palette("micro_cvd_gray")))
micro_cvd_purple <- rev(c(microshades_palette("micro_cvd_purple")))
micro_cvd_blue <- rev(c(microshades_palette("micro_cvd_blue")))
micro_cvd_orange <- rev(c(microshades_palette("micro_cvd_orange")))
micro_cvd_green <- rev(c(microshades_palette("micro_cvd_green")))
micro_cvd_turquoise <- rev(c(microshades_palette("micro_cvd_turquoise")))
micro_orange <- rev(c(microshades_palette("micro_orange")))
micro_purple <- rev(c(microshades_palette("micro_purple")))

# Choose the column to which you want to assign the vector of colors (e.g., "column4")

hex_custom <- data.frame(
  micro_cvd_gray = micro_cvd_gray,
  micro_cvd_purple = micro_cvd_purple,
  micro_cvd_blue = micro_cvd_blue,
  micro_cvd_orange = micro_cvd_orange,
  micro_cvd_green = micro_cvd_green,
  micro_cvd_turquoise = micro_cvd_turquoise,
  micro_orange = micro_orange,
  micro_purple = micro_purple
)

# Custom function adapted from https://github.com/KarstensLab/microshades

# We create a fixed color scale function


fixed_custom_create_color_dfs <- function(mdf,
                                          selected_groups = c(
                                            "Proteobacteria",
                                            "Actinobacteria",
                                            "Bacteroidetes",
                                            "Firmicutes"
                                          ),
                                          top_n_subgroups = 4,
                                          group_level = "Phylum",
                                          subgroup_level = "Genus",
                                          cvd = FALSE,
                                          top_orientation = FALSE) {
  # Throws error if too many subgroups
  if (top_n_subgroups > 4) {
    stop("'top_n_subgroups' exceeds MAX value 4")
  }

  if (class(mdf) != "data.frame") {
    stop("mdf argument must be a data frame")
  }
  if (!is.null(mdf$group)) {
    stop("'group' column name already exists; consider renaming or removing")
  }

  if (is.null(mdf[[group_level]])) {
    stop("'group_level' does not exist")
  }

  if (is.null(mdf[[subgroup_level]])) {
    stop("'subgroup_level' does not exist")
  }

  # Here we add a security check to make sure that the Others is present in mdf[[group_level]]. Else we add it directly to the mdf dataframe

  if ("Other" %in% mdf[[group_level]]) {
    print("Other is present in the dataframe")
  } else {
    print("Other is not present in the dataframe. We add it directly to the dataframe")
    row_to_insert <- data.frame(
      canopus_npc_pathway = "Other",
      canopus_npc_superclass = "Other",
      canopus_npc_class = "Other",
      Abundance = 1
    )
    mdf <- mdf %>%
      rows_insert(row_to_insert)
  }

  # Create new column for group level -----
  # Add "Other" category immediately
  col_name_group <- paste0("Top_", group_level)
  mdf[[col_name_group]] <- "Other"

  # Index and find rows containing the selected groups
  rows_to_change <- mdf[[group_level]] %in% selected_groups
  taxa_names_mdf <- row.names(mdf[rows_to_change, ])
  mdf[taxa_names_mdf, col_name_group] <-
    as.character(mdf[taxa_names_mdf, group_level])

  if ("Other" %in% selected_groups) {
    # Create factor for the group level column
    mdf[[col_name_group]] <- factor(mdf[[col_name_group]],
      levels = c(selected_groups)
    )
  } else {
    # Create factor for the group level column
    mdf[[col_name_group]] <- factor(mdf[[col_name_group]],
      levels = c("Other", selected_groups)
    )
  }

  # Check to make sure the selected_groups specified all exist in the dataset
  if (sum(selected_groups %in% as.character(unique(mdf[[col_name_group]]))) != length(selected_groups)) {
    stop("some 'selected_groups' do not exist in the dataset. Consider SILVA 138 c('Proteobacteria', 'Actinobacteriota', 'Bacteroidota', 'Firmicutes')")
  }

  # Rename missing genera
  mdf_unknown_subgroup <- mdf %>%
    mutate(!!sym(subgroup_level) := fct_na_value_to_level(!!sym(subgroup_level), "Unknown")) ## fct_na_value_to_level

  # Rank group-subgroup categories by ranked abundance and add order
  # Ranked abundance aggregated using sum() function
  col_name_subgroup <- paste0("Top_", subgroup_level)
  subgroup_ranks <- mdf_unknown_subgroup %>%
    group_by_at(c(paste(subgroup_level), paste(col_name_group))) %>%
    summarise(rank_abundance = sum(Abundance)) %>%
    arrange(desc(rank_abundance)) %>%
    group_by_at(c(paste(col_name_group))) %>%
    mutate(order = row_number()) %>%
    ungroup()

  # Correctly keep "Other" for lower abundant genera
  # Pseudocode:
  # - set all (top) subgroups to "Other"
  # - change subgroups back to actual subgroups (e.g., Genus) if it is in the
  #   top N number of subgroups passed into `top_n_subgroups` (e.g., 4)
  subgroup_ranks[[col_name_subgroup]] <- "Other"
  rows_to_change <- subgroup_ranks$order <= top_n_subgroups
  subgroup_ranks[rows_to_change, col_name_subgroup] <-
    as.vector(subgroup_ranks[rows_to_change, subgroup_level])

  # Generate group-subgroup categories -----
  # There are `top_n_subgroups` additional groups because each group level has
  # an additional subgroup of "Other"
  # E.g., 4 selected_groups + 1 Other, 4 top_n_groups + 1 Other => 25 groups
  group_info <- subgroup_ranks %>%
    mutate(group = paste(!!sym(col_name_group),
      !!sym(col_name_subgroup),
      sep = "-"
    ))

  # Ensure that the "Other" subgroup is always the lightest shade
  group_info$order[group_info[[col_name_subgroup]] == "Other"] <- top_n_subgroups + 1

  # Merge group info back to df -----
  # Get relevant columns from data frame with group info
  group_info_to_merge <-
    group_info[, c(
      col_name_group, subgroup_level,
      col_name_subgroup, "group"
    )]
  mdf_group <- mdf_unknown_subgroup %>%
    left_join(group_info_to_merge, by = c(col_name_group, subgroup_level))

  # Get beginning of color data frame with top groups/subgroups
  # E.g., 4 selected_groups + 1 Other, 4 top_n_groups + 1 Other => 25 groups
  prep_cdf <- group_info %>%
    select(all_of(c("group", "order", col_name_group, col_name_subgroup))) %>%
    filter(order <= top_n_subgroups + 1) %>% # "+ 1" for other subgroup
    arrange(!!sym(col_name_group), order)

  # Prepare hex colors -----

  # Generates default 5 row x 6 cols of 5 colors for 6 phylum categories
  # Parameter for number of selected phylum
  # "+ 1" is for "Other" group
  num_group_colors <- length(selected_groups) + 1

  # hex_df <- default_hex(num_group_colors, cvd)

  hex_df <- hex_custom %>%
    rownames_to_column("order") %>%
    mutate(order = as.numeric(order))


  # Add hex codes in ranked way
  # creates nested data frame
  # https://tidyr.tidyverse.org/articles/nest.html
  # https://tidyr.tidyverse.org/reference/nest.html
  cdf <- prep_cdf %>%
    group_by_at(c(paste(col_name_group))) %>%
    tidyr::nest() %>%
    arrange(!!sym(col_name_group))

  # Loop through top group and add colors by nested data frame
  # Higher row number = less abundant = lighter color

  # if ("Other" %in% mdf[[col_name_group]])
  # {
  #   start <- 1
  # } else
  # {
  #   start <- 2
  #   num_group_colors <- num_group_colors -1
  # }

  # for (i in 1:num_group_colors) {
  #   cdf$data[[i]]$hex <- hex_df[1:length(cdf$data[[i]]$group),start]
  #   start = start + 1
  # }
  # 1 is micro_cvd_gray

  # Define a function to create a pathway tibble
  create_pathway_tibble <- function(pathway_name, hex_column_name) {
    cdf %>%
      filter(Top_canopus_npc_pathway == pathway_name) %>%
      pull(data) %>%
      as.data.frame() %>%
      left_join(hex_df, by = "order") %>%
      select(group, order, Top_canopus_npc_superclass, !!hex_column_name := !!sym(hex_column_name)) %>%
      rename(hex = !!sym(hex_column_name)) %>%
      as_tibble() %>%
      distinct()
  }


  # List of pathway names and corresponding hex column names
  pathway_data <- list(
    "Terpenoids" = "micro_cvd_purple",
    "Fatty acids" = "micro_cvd_blue",
    "Polyketides" = "micro_cvd_orange",
    "Alkaloids" = "micro_cvd_green",
    "Shikimates and Phenylpropanoids" = "micro_cvd_turquoise",
    "Amino acids and Peptides" = "micro_orange",
    "Carbohydrates" = "micro_purple",
    "Other" = "micro_cvd_gray"
  )
  # Subset pathway_data to match the levels in cdf$Top_canopus_npc_pathway
  valid_pathway_names <- levels(cdf$Top_canopus_npc_pathway)
  valid_pathway_data <- pathway_data[names(pathway_data) %in% valid_pathway_names]


  # Loop through the pathway_data and create tibbles
  pathway_tibbles <- list()
  for (pathway_name in names(valid_pathway_data)) {
    hex_column_name <- valid_pathway_data[[pathway_name]]
    pathway_tibbles[[pathway_name]] <- create_pathway_tibble(pathway_name, hex_column_name)
  }

  # Convert the list of tibbles to a tibble
  cdf <- tibble(
    Top_canopus_npc_pathway = names(pathway_tibbles),
    data = map(pathway_tibbles, as_tibble)
  )

  # Unnest colors and groups and polish for output
  cdf <- cdf %>%
    ungroup() %>%
    arrange(desc(row_number())) %>%
    tidyr::unnest(data) %>%
    select(
      !!sym(col_name_group),
      !!sym(col_name_subgroup),
      group, hex, order
    ) %>%
    mutate_all(as.character) # Remove factor from hex codes

  cdf <- cdf %>% filter(!is.na(hex))

  if (top_orientation) {
    level_assign <- unique(cdf$group)
  } else {
    level_assign <- unique(rev(cdf$group))
  }

  mdf_group$group <- factor(mdf_group$group, levels = level_assign)


  # Return final objects -----
  list(
    mdf = mdf_group,
    cdf = cdf
  )
}

# mdf = selected_variable_meta_NPC_simple_resolved
# # selected_groups = c("Terpenoids", "Fatty acids", "Polyketides", "Alkaloids", "Shikimates and Phenylpropanoids", "Amino acids and Peptides", "Carbohydrates")
# selected_groups = selected_variable_meta_NPC_simple_resolved_count
# group_level = "canopus_npc_pathway"
# subgroup_level = "canopus_npc_superclass"
# cvd = TRUE
# top_n_subgroups = 4
# top_orientation = FALSE



selected_variable_meta_NPC_simple_resolved <- DE$variable_meta %>%
  filter(feature_id %in% features_of_importance) %>%
  select(canopus_npc_class) %>%
  rownames_to_column("index") %>%
  left_join(npclassifier_newpath_simple, by = "canopus_npc_class") %>%
  column_to_rownames("index") %>%
  select(canopus_npc_pathway, canopus_npc_superclass, canopus_npc_class) %>%
  # We convert NA in the canopus_npc_pathway column to "Other"
  mutate(canopus_npc_pathway = ifelse(is.na(canopus_npc_pathway), "Other", canopus_npc_pathway))


selected_variable_meta_NPC_simple_resolved$Abundance <- 1


# show_col(cdf_variable_meta_NPC_simple_resolved_colored$hex, cex_label = 0.5)

### Selected dataset

selected_variable_meta_NPC_simple_resolved_count <- selected_variable_meta_NPC_simple_resolved %>%
  group_by(canopus_npc_pathway) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  select(canopus_npc_pathway) %>%
  pull()


selected_variable_meta_NPC_simple_resolved_colored <- fixed_custom_create_color_dfs(selected_variable_meta_NPC_simple_resolved, selected_groups = selected_variable_meta_NPC_simple_resolved_count, group_level = "canopus_npc_pathway", subgroup_level = "canopus_npc_superclass", cvd = TRUE)

# Extract
mdf_selected_variable_meta_NPC_simple_resolved_colored <- selected_variable_meta_NPC_simple_resolved_colored$mdf
cdf_selected_variable_meta_NPC_simple_resolved_colored <- selected_variable_meta_NPC_simple_resolved_colored$cdf


col_order_np_pathway <- mdf_selected_variable_meta_NPC_simple_resolved_colored %>%
  distinct(group, .keep_all = TRUE) %>%
  # we now merge the df with the cdf_selected_variable_meta_NPC_simple_resolved_colored_plus df to get the hex color code
  # with the Top_canopus_npc_pathway column on the left and the canopus_npc_pathway column on the right
  left_join(cdf_selected_variable_meta_NPC_simple_resolved_colored, by = "group") %>%
  # arrange(canopus_npc_pathway, order)  %>%
  arrange(ifelse(canopus_npc_pathway == "Other", 2, 1), canopus_npc_pathway, order) %>%
  # left_join(df_col_np_pathway, by.x = "Top_canopus_npc_pathway", by.x = "canopus_npc_pathway") %>%
  select(hex, group)


col_np_pathway <- col_order_np_pathway %>%
  select(hex) %>%
  as.vector() %>%
  unlist() %>%
  rev()

order_np_pathway <- col_order_np_pathway %>%
  select(group) %>%
  as.vector() %>%
  unlist() %>%
  rev()


mdf_selected_variable_meta_NPC_simple_resolved_colored$group <- factor(mdf_selected_variable_meta_NPC_simple_resolved_colored$group, levels = order_np_pathway)

# The grid parameters are defined

grid_params <- setup_colorbar_grid(
  nrows = 2,
  y_length = 0.4,
  x_spacing = 0.3,
  y_spacing = 0.5,
  x_start = 1.1,
  y_start = 0.8
)

# We create the hover text for the heatmap

dt <- as.data.frame(t(data_subset_for_pval_hm_mat))

if (gnps2_job) {
  values_mat <- dt %>%
  rownames_to_column("feature_id") %>%
  select(feature_id) %>%
  mutate(feature_id = as.numeric(feature_id)) %>%
  left_join(DE$variable_meta, by = "feature_id") %>%
  select(feature_id, gnps_component, sirius_adduct, sirius_chebiasciiname, sirius_name, canopus_npc_pathway, canopus_npc_superclass, canopus_npc_class) %>%
  mutate(hover_text = paste0(
    "<br>", "feature_id: ", feature_id,
    "<br>", "gnps_component: ", gnps_component,
    "<br>", "Adduct sirius: ", sirius_adduct,
    "<br>", "CheBI name: ", sirius_chebiasciiname,
    "<br>", "Sirius name: ", sirius_name,
    "<br>", "Pathway: ", canopus_npc_pathway,
    "<br>", "Superclass: ", canopus_npc_superclass,
    "<br>", "Class: ", canopus_npc_class
  )) %>%
  select(feature_id, hover_text) %>%
  pivot_wider(names_from = feature_id, values_from = hover_text) %>%
  # we now repeat the hover_text for each row of the matrix. We use dplyr to do that
  mutate(count = nrow(data_subset_for_pval_hm_mat)) %>%
  uncount(count) %>%
  as.matrix()
} else {
  values_mat <- dt %>%
  rownames_to_column("feature_id") %>%
  select(feature_id) %>%
  mutate(feature_id = as.numeric(feature_id)) %>%
  left_join(DE$variable_meta, by = "feature_id") %>%
  select(feature_id, gnps_componentindex, sirius_adduct, sirius_chebiasciiname, sirius_name, canopus_npc_pathway, canopus_npc_superclass, canopus_npc_class) %>%
  mutate(hover_text = paste0(
    "<br>", "feature_id: ", feature_id,
    "<br>", "gnps_componentindex: ", gnps_componentindex,
    "<br>", "Adduct sirius: ", sirius_adduct,
    "<br>", "CheBI name: ", sirius_chebiasciiname,
    "<br>", "Sirius name: ", sirius_name,
    "<br>", "Pathway: ", canopus_npc_pathway,
    "<br>", "Superclass: ", canopus_npc_superclass,
    "<br>", "Class: ", canopus_npc_class
  )) %>%
  select(feature_id, hover_text) %>%
  pivot_wider(names_from = feature_id, values_from = hover_text) %>%
  # we now repeat the hover_text for each row of the matrix. We use dplyr to do that
  mutate(count = nrow(data_subset_for_pval_hm_mat)) %>%
  uncount(count) %>%
  as.matrix()
}



# We change the numeric into E-notation and we round the values to 2 decimals.


data_subset_for_pval_hm_peak_height <- format(data_subset_for_pval_hm_peak_height, digits = 2, scientific = TRUE)


combined_matrix <- matrix(paste(as.matrix(data_subset_for_pval_hm_peak_height), values_mat, sep = "<br>"), nrow = nrow(data_subset_for_pval_hm_peak_height), ncol = ncol(data_subset_for_pval_hm_peak_height))



##########################



iheatmap <- iheatmapr::main_heatmap(as.matrix(t(data_subset_for_pval_hm_mat)), ### add heat map top 100
  name = "Intensity",
  # layout = list(margin = list(b = 80)),
  colorbar_grid = grid_params,
  colors = "GnBu",
  show_colorbar = TRUE,
  text = t(combined_matrix),
  layout = list(
    title = list(text = title_heatmap_pval, font = list(size = 14), x = 0.1),
    margin = list(t = 160, r = 80, b = 80, l = 80)
  )
) %>%
  add_row_labels(
    tickvals = NULL,
    ticktext = selected_variable_meta$feature_id_full_annotated,
    side = "left",
    buffer = 0.01,
    textangle = 0,
    size = 0.45,
    font = list(size = 9)
  ) %>%
  # add_row_annotation(selected_variable_meta_NPC_simple,
  #   side = "right",
  #   buffer = 0.005
  # colors = list(ByPal, ByPal, ByPal) %>%
  # add_row_annotation(data.frame("NPC_Class" = selected_variable_meta_NPC_simple_resolved$canopus_npc_class),
  #   side = "right",
  #   buffer = 0.005,
  #   colors = list("NPC_Class" = col_np_class)
  # ) %>%
  # add_row_annotation(data.frame("NPC_SuperClass" = selected_variable_meta_NPC_simple_resolved$canopus_npc_superclass),
  #   side = "right",
  #   buffer = 0.005,
  #   colors = list("NPC_SuperClass" = col_np_superclass)
  # ) %>%
  add_row_annotation(data.frame("Classification" = mdf_selected_variable_meta_NPC_simple_resolved_colored$group),
    side = "right",
    buffer = 0.05,
    colors = list("Classification" = col_np_pathway)
  ) %>%
  # add_row_annotation(selected_variable_meta_NPC_simple$canopus_npc_superclass,
  #   side = "right",
  #   buffer = 0.005
  # colors = list(ByPal, ByPal, ByPal) %>%
  # add_row_annotation(selected_variable_meta_NPC_simple$canopus_npc_pathway,
  #   side = "right",
  #   buffer = 0.005
  # colors = list(ByPal, ByPal, ByPal) %>%
  add_row_clustering(side = "right") %>%
  add_col_annotation(data.frame("Condition" = target_metadata),
    colors = list("Condition" = custom_colors_heatmap),
    buffer = 0.01
  ) %>%
  add_col_clustering() %>%
  add_col_labels(
    tickvals = NULL,
    ticktext = my_sample_col,
    textangle = -90,
    size = 0.2,
    font = list(size = 10)
  )

# iheatmap

# library(iheatmapr)
# data(measles, package = "iheatmapr")

# main_heatmap(measles, name = "Measles<br>Cases", x_categorical = FALSE,
#              layout = list(font = list(size = 8))) %>%
#   add_col_groups(ifelse(1930:2001 < 1961,"No","Yes"),
#                   side = "bottom", name = "Vaccine<br>Introduced?",
#                   title = "Vaccine?",
#                   colors = c("lightgray","blue")) %>%
#   add_col_labels(ticktext = seq(1930,2000,10),font = list(size = 8)) %>%
#   add_row_labels(size = 0.3,font = list(size = 6)) %>%
#   add_col_summary(layout = list(title = "Average<br>across<br>states"),
#                   yname = "summary")  %>%
#   add_col_title("Measles Cases from 1930 to 2001", side= "top") %>%
#   add_row_summary(groups = TRUE,
#                   type = "bar",
#                   layout = list(title = "Average<br>per<br>year",
#                                 font = list(size = 8)))



# The file is exported

iheatmap %>% save_iheatmap(file = filename_heatmap_pval) # Save interactive HTML



# unique(selected_variable_meta_NPC_simple_resolved$canopus_npc_superclass)

# col_np_superclass


# hex_values <-c(microshades_palette("micro_green", 5, lightest = FALSE),
#                microshades_palette("micro_blue", 3, lightest = FALSE),
#                microshades_palette("micro_purple", 3, lightest = FALSE))


# library(RColorBrewer)

# col.g <- c(brewer.pal(9,"Greens"))[1:102] # select 5 colors from class Greens
# col.r <- c(brewer.pal(9,"Reds"))[6:9] # select 4 colors from class Reds
# col.b <- c(brewer.pal(9,"Blues"))[6:9] # select 4 colors from class Blues
# my.cols <- c(col.g,col.r,col.b)


# # df_col <- treepalette(
# #   selected_variable_meta_NPC_simple_ordered,
# #   index = names(selected_variable_meta_NPC_simple_ordered),
# #   method = "HCL",
# #   palette = NULL,
# #   return.parameters = TRUE,
# #   prepare.dat = TRUE
# # )

# df_col_np <- treepalette(
#   npclassifier_origin_ordered,
#   index = names(npclassifier_origin_ordered),
#   method = "HCL",
#   palette = NULL,
#   palette.HCL.options = list(hue_perm = TRUE),
#   return.parameters = TRUE,
#   prepare.dat = TRUE
# )

# show_col(df_col_np$HCL.color, cex_label = 0.5)


# npclassifier_origin_ordered

# # we drop the row with the NA value in all columns

# npclassifier_origin_ordered = npclassifier_origin_ordered[complete.cases(npclassifier_origin_ordered), ]



# library(scales)
# pal <- rgb(ddf$r, ddf$g, ddf$b)



# show_col(paired, cex_label = 0.5)



# selected_variable_meta_NPC

# colnames(data_subset_for_pval_hm_mat) = selected_variable_meta$feature_id_full_annotated

# typeof(measles)
# DE$sample_meta[[params$target$sample_metadata_header]]


# Indometh_matrix <- acast(Indometh, Subject ~ time, value.var = "conc")
# Indometh_matrix <- Indometh_matrix[as.character(1:6),]
# rownames(Indometh_matrix) <- paste("Patient",rownames(Indometh_matrix))
# Indometh_patient_cor <- cor(t(Indometh_matrix))

# patient_max_conc <- apply(Indometh_matrix,1,max)
# patient_min_conc <- apply(Indometh_matrix,1,min)
# patient_groups <- c("A","A","B","A","B","A") # Arbitrary groups

# example_heatmap <- main_heatmap(Indometh_patient_cor, name = "Correlation")

# example_heatmap
# typeof(Indometh_patient_cor)

# data_subset_for_pval_hm

# rownames(data_subset_for_pval_hm) <- my_sample_col

#############################################################################
#############################################################################
############## Random Forest filtered Heat Map  #############################
#############################################################################
#############################################################################

# message("Preparing Random Forest filtered Heatmap ...")

# imp_table_rf_order = imp_table_rf[order(imp_table_rf$MeanDecreaseGini, decreasing = TRUE), ] # imp_table_rf[order(imp_table_rf$MeanDecreaseGini,decreasing=TRUE),]

# imp_filter2X = row.names(imp_table_rf_order)[1:params$heatmap$topN]
# imp_filter2 = gsub("X", "", imp_filter2X)

# data_subset_for_RF = data_subset_for_RF[, colnames(data_subset_for_RF) %in% imp_filter2X]
# # my_sample_col = DE$sample_meta$sample_id

# my_sample_col = paste(DE$sample_meta$sample_id, DE$sample_meta[[params$target$sample_metadata_header]], sep = "_")

# annot_col = data.frame(paste(DE$variable_meta$canopus_npc_pathway, DE$variable_meta$canopus_npc_superclass, sep = "_"), DE$variable_meta$canopus_npc_pathway)

# colnames(annot_col) = c("Superclass", "Pathway")

# rownames(annot_col) = DE$variable_meta$feature_id
# annot_col_filter = annot_col[rownames(annot_col) %in% imp_filter2, ]


# ByPal = colorRampPalette(c(wes_palette("Zissou1")))

# data_subset_for_RF = apply(data_subset_for_RF, 2, as.numeric)
# # heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")


# heatmap_filtered_rf = heatmaply(
#   percentize(data_subset_for_RF),
#   seriate = "mean", # none , GW , mean, OLO
#   col_side_colors = data.frame(annot_col_filter, check.names = FALSE),
#   col_side_palette = ByPal,
#   labRow = as.vector(as.character(my_sample_col)), # [vec_plot]
#   subplot_margin = 0.01,
#   scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
#     low = "lightsteelblue2",
#     # mid = "goldenrod1",
#     high = "firebrick3",
#     midpoint = 0.5,
#     limits = c(0, 1)
#   ),
#   fontsize_col = 5,
#   branches_lwd = 0.3,
#   k_row = 4,
#   distfun_row = "pearson",
#   distfun_col = "pearson"
# )

# heatmap_filtered_rf = heatmap_filtered_rf %>% layout(title = list(text = title_heatmap_rf, y = 0.05))

# # The file is exported


# if (params$operating_system$system == "unix") {
# ### linux version

# heatmap_filtered_rf %>%
#     htmlwidgets::saveWidget(file = filename_heatmap_rf, selfcontained = TRUE)
# }

# if (params$operating_system$system == "windows") {
# ### windows version
# Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
# heatmap_filtered_rf %>%
#     htmlwidgets::saveWidget(file = filename_heatmap_rf, selfcontained = TRUE,libdir = "lib")
# unlink("lib", recursive = FALSE)

# }

#############################################################################
#############################################################################
############## Summary Table ################################################
#############################################################################
#############################################################################

message("Outputing Summary Table ...")

# Output is not clean. Feature id are repeated x times.
# To tidy ---

# feature_id = DE$variable_meta$feature_id_full
# sample_raw_id = paste("X", DE$variable_meta$feature_id_full, sep = "")
# sirius_name = DE$variable_meta$sirius_name
# sirius_smiles = DE$variable_meta$sirius_smiles
# InChI_sirius = DE$variable_meta$InChI_sirius
# canopus_npc_pathway = DE$variable_meta$canopus_npc_pathway
# canopus_npc_superclass = DE$variable_meta$canopus_npc_superclass
# Annotation_merge = data.frame(feature_id, sample_raw_id, sirius_name, sirius_smiles, InChI_sirius, canopus_npc_pathway, canopus_npc_superclass)

# RF_importance = imp_table_rf$MeanDecreaseGini
# sample_raw_id = row.names(imp_table_rf)
# RF_importance_merge = data.frame(sample_raw_id, RF_importance)

# summary_matt1 = merge(Annotation_merge, RF_importance_merge, by = "sample_raw_id", all = TRUE)

# sample_raw_id = paste("X", matt_volcano_tot$mol, sep = "")
# group = matt_volcano_tot$X1
# posthoc_Pvalue = matt_volcano_tot$p.value
# log10P = matt_volcano_tot$log10P
# posthoc_result_merge = data.frame(sample_raw_id, group, posthoc_Pvalue, log10P)
# split_group = split(posthoc_result_merge, group)
# matt_split = do.call("cbind", split_group)
# matt_split$sample_raw_id = matt_split[, 1]

# summary_stat_output = merge(summary_matt1, matt_split, by = "sample_raw_id", all = TRUE)

summary_stat_output_full <- DE_foldchange_pvalues

# We filter the DE_foldchange_pvalues table to only keep the top N features (any column ending with _p_value string should have a value < 0.05)
# We use the dplyr synthax to filter the table
# We need to make sure to remove the rownames() before exporting


summary_stat_output_selected <- DE_foldchange_pvalues %>%
  filter(if_any(ends_with("_p_value"), ~ . < params$posthoc$p_value)) %>%
  arrange(across(ends_with("_p_value"))) %>%
  select(
    feature_id_full,
    feature_id,
    feature_mz,
    feature_rt,
    contains("p_value"),
    contains("fold"),
    canopus_npc_pathway,
    canopus_npc_superclass,
    canopus_npc_class,
    sirius_name,
    # gnps_libraryid,
    contains("smiles", ignore.case = TRUE),
    contains("inchi", ignore.case = TRUE),
    contains("inchikey", ignore.case = TRUE)
  )

summary_stat_output_selected_cytoscape <- DE_foldchange_pvalues %>%
  select(
    feature_id,
    feature_id_full,
    feature_id_full_annotated,
    feature_mz,
    feature_rt,
    contains("p_value"),
    contains("fold"),
    sirius_name,
    sirius_chebiasciiname,
    sirius_chebiid,
    sirius_confidencescore,
    sirius_csi_fingeridscore,
    sirius_siriusscore,
    sirius_zodiacscore,
    sirius_inchi,
    sirius_inchikey2d,
    sirius_molecularformula,
    sirius_adduct,
    sirius_smiles,
    canopus_npc_pathway,
    canopus_npc_class,
    canopus_npc_superclass,
    met_annot_structure_exact_mass,
    met_annot_structure_inchi,
    met_annot_structure_inchikey,
    met_annot_short_inchikey,
    met_annot_structure_smiles,
    met_annot_structure_molecular_formula,
    met_annot_structure_nametraditional,
    met_annot_structure_wikidata,
    met_annot_structure_taxonomy_npclassifier_01pathway,
    met_annot_structure_taxonomy_npclassifier_02superclass,
    met_annot_structure_taxonomy_npclassifier_02superclass,
    met_annot_structure_taxonomy_npclassifier_01pathway_consensus,
    met_annot_structure_taxonomy_npclassifier_02superclass_consensus,
    met_annot_structure_taxonomy_npclassifier_03class_consensus,
    met_annot_organism_name,
    met_annot_organism_taxonomy_01domain,
    met_annot_organism_taxonomy_02kingdom,
    met_annot_organism_taxonomy_03phylum,
    met_annot_organism_taxonomy_04class,
    met_annot_organism_taxonomy_05order,
    met_annot_organism_taxonomy_06family,
    met_annot_organism_taxonomy_07tribe,
    met_annot_organism_taxonomy_08genus,
    met_annot_organism_taxonomy_09species,
    met_annot_organism_taxonomy_10varietas,
    met_annot_organism_taxonomy_ottid,
    met_annot_organism_wikidata,
    met_annot_score_taxo
  )

#### ad rf importance

imp_table_rf$feature_id <- gsub("X", "", row.names(imp_table_rf))
summary_stat_output_selected_cytoscape <- merge(summary_stat_output_selected_cytoscape, imp_table_rf, by = "feature_id")
# glimpse(summary_stat_output_selected)


# We also prepare Metaboverse outputs from the fc and pvalues tables


metaboverse_table <- DE_foldchange_pvalues

# We then keep the keep the first occurence of the sirius_chebiasciiname

metaboverse_table <- metaboverse_table %>%
  distinct(sirius_chebiasciiname, .keep_all = TRUE)

# We now format the table for Metaboverse
# For this we apply the foillowing steps:
# 1. We select the columns we want to keep (sirius_chebiasciiname, Co_KO_p_value, Co_KO_fold_change_log2)
# 2. We rename the columns to the names Metaboverse expects. Using the rename_with and gsub we replace the the _p_value and _fold_change_log2 suffixes to _stat and _fc
# 3. We reorganize the columns to the order Metaboverse expects (sirius_chebiasciiname, _stat, _fc)
# 4. We remove any rows containing NA values in the dataframe
# 5. We replace the name of the `sirius_chebiasciiname` column by an empty string


metaboverse_table <- metaboverse_table %>%
  select(sirius_chebiasciiname, ends_with("_fold_change_log2"), ends_with("_p_value")) %>%
  # rename_with(~gsub("_p_value", "_stat", .)) %>%
  # rename_with(~gsub("_fold_change_log2", "_fc", .)) %>%
  # select(sirius_chebiasciiname, Co_KO_fc, Co_KO_stat) %>%
  # We remove row containing the `Inf` value
  # filter(!grepl('Inf', ends_with('_fold_change_log2')))  %>%
  # We remove any row containing the `Inf` value across all columns of the dataframe
  filter(if_any(everything(), ~ !str_detect(., "Inf"))) %>%
  # filter(!grepl('Inf', ends_with('_fold_change_log2')))  %>%
  na.omit()

colnames(metaboverse_table)[1] <- ""

# We now sort columns alphabetically

metaboverse_table <- metaboverse_table[, order(colnames(metaboverse_table))]




# The file is exported

write.table(summary_stat_output_full, file = filename_summary_stats_table_full, sep = ",", row.names = FALSE)
write.table(summary_stat_output_selected, file = filename_summary_stats_table_selected, sep = ",", row.names = FALSE)
write.table(summary_stat_output_selected_cytoscape, file = filename_summary_stat_output_selected_cytoscape, sep = ",", row.names = FALSE)
write.table(metaboverse_table, file = filename_metaboverse_table, sep = "\t", row.names = FALSE, quote = FALSE)


#############################################################################
#############################################################################
######################################## summmary table with structure


summary_stat_output_selected_simple <- DE_foldchange_pvalues %>%
  filter(if_any(ends_with("_p_value"), ~ . < params$posthoc$p_value)) %>%
  arrange(across(ends_with("_p_value"))) %>%
  select(
    feature_id,
    feature_id_full,
    contains("p_value"),
    canopus_npc_pathway,
    canopus_npc_superclass,
    canopus_npc_class,
    sirius_name,
    sirius_smiles
  )


# mol_list_2d <- list()

# for (i in c(1:nrow(summary_stat_output_selected_simple))) {
#   psml <- parse.smiles(summary_stat_output_selected_simple$sirius_smiles[i], omit.nulls = TRUE)
#   if (length(psml) == 0) {
#     psml <- parse.smiles("C")
#   }
# img <- view.image.2d(psml[1][[1]])
# test <-  as.matrix(as.raster(img))
# matrice_df <- melt(test)
# # Renommer les colonnes
# colnames(matrice_df) <- c("y", "x", "couleur")
# mol_list_2d[[i]] = ggplot(matrice_df, aes(x = x, y = rev(y), fill = couleur)) +
#                    geom_tile() +
#                    scale_fill_identity() +
#                    labs(x = "Axe X", y = "Axe Y") +
#                   theme(legend.position="none") + theme_void()
# }




# summary_stat_output_selected_simple$plots  <- mol_list_2d
# summary_stat_output_selected_simple$ggplot  <- rep(NA,length(mol_list_2d))

# tab_1 <- summary_stat_output_selected_simple %>%
#     select(-plots) %>%
#     gt() %>%
#   text_transform(locations = cells_body(c(ggplot)),
#                  fn = function(x) {
#                   map(summary_stat_output_selected_simple$plots, ggplot_image, height = px(100))
#                  }
#                  )

# tab_1 %>% gtsave(filename = filename_interactive_table, inline_css = TRUE) ### add path to save

#############################################################################
#############################################################################
############## GraphML output ################################################
#############################################################################
#############################################################################

message("Generating GraphML output ...")


# We first load the GNPS graphml file

if (gnps2_job) {
  graphml_dir <- file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "nf_output", "networking")
  graphml_file <- list.files(path = graphml_dir, pattern = "network.graphml$")

} else {
  graphml_dir <- file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "gnps_molecular_network_graphml")
  graphml_file <- list.files(path = graphml_dir, pattern = "\\.graphml$")

}


graphml_file_path <- file.path(graphml_dir, graphml_file)



g <- read.graph(file = graphml_file_path, format = "graphml")
# net_gnps = igraph::simplify(g, remove.multiple = FALSE, edge.attr.comb = "ignore")

df_from_graph_edges_original <- igraph::as_data_frame(g, what = c("edges"))
df_from_graph_vertices_original <- igraph::as_data_frame(g, what = c("vertices"))


# We define drop the from and to columns from the edges dataframe
# And then rename the node 1 and node 2 columns to from and to, respectively
# These columns are placed at the beginning of the dataframe
# and converted to numerics
if (gnps2_job) {
  df_from_graph_edges <- df_from_graph_edges_original %>%
    select(-from, -to) %>%
    rename(from = scan1, to = scan2) %>%
    select(from, to, everything()) %>%
    mutate_at(vars(from, to), as.numeric)
} else {
  df_from_graph_edges <- df_from_graph_edges_original %>%
    select(-from, -to) %>%
    rename(from = node1, to = node2) %>%
    select(from, to, everything()) %>%
    mutate_at(vars(from, to), as.numeric)
}

# the id column of the vertices dataframe is converted to numerics

df_from_graph_vertices <- df_from_graph_vertices_original %>%
  mutate_at(vars(id), as.numeric)


# We then add the attributes to the vertices dataframe
# For this we merge the vertices dataframe with the VM output using the id column and the feature_id column, respectively

# vm_minus_gnps = DE_original$variable_meta  %>%
# select(-contains("_gnps"))


# df_from_graph_vertices_plus = merge(df_from_graph_vertices, vm_minus_gnps, by.x = "id", by.y = "feature_id", all.x = T)

# glimpse(df_from_graph_vertices_plus)


# Now we will add the results of the statistical outputs to the vertices dataframe
# For this we merge the vertices dataframe with the summary_stat_output using the id column and the feature_id column, respectively

# First we clean the summary_stat_output dataframe
# For this we remove columns that are not needed. The one containing the sirius and canopus pattern in the column names. Indeed they arr already present in the VM dataframe

# summary_stat_output_red = summary_stat_output_full %>%
#   select(-contains("_sirius")) %>%
#   select(-contains("_canopus")) %>%
#   select(-contains("_met_annot")) %>%
#   select(-contains("_gnps"))  %>%
#   #select(-ends_with("_id"))  %>%
#   select(-ends_with("_mz"))  %>%
#   select(-ends_with("_rt"))

DE_original_features <- DE_original$variable_meta 
# %>%
#   select(feature_id)

# We merge DE_original_features and summary_stat_output_selected_cytoscape but make sure to drop duplicated columns from the summary_stat_output_selected_cytoscape dataframe

# Identify the names of columns that are duplicated between the two data frames, excluding the merging key column.
common_cols <- setdiff(intersect(names(DE_original_features), names(summary_stat_output_selected_cytoscape)), "feature_id")

# Remove Duplicated Columns from One DataFrame
summary_stat_output_selected_cytoscape <- summary_stat_output_selected_cytoscape %>%
  select(-all_of(common_cols))


df_from_graph_vertices_plus <- DE_original_features %>%
  left_join(summary_stat_output_selected_cytoscape, by = "feature_id")



# We merge the data from the DE$data dataframe with the DE$sample_meta dataframe using rownames as the key

merged_D_SM <- merge(DE_original$sample_meta, DE_original$data, by = "row.names", all = TRUE)

# We replace NA values with 0 in the merged dataframe
# Check why we dont do this before ??
merged_D_SM[is.na(merged_D_SM)] <- 0


# The function below allows to group data by multiple factors and return a dataframe with the mean of each group


dfList <- list()

for (i in params$to_mean$factor_name) {
  dfList[[i]] <- merged_D_SM %>%
    group_by(!!as.symbol(i)) %>%
    summarise(across(colnames(DE_original$data), mean),
      .groups = "drop"
    ) %>%
    select(!!all_of(i), colnames(DE_original$data)) %>%
    pivot_longer(-!!i) %>%
    pivot_wider(names_from = all_of(i), values_from = value) %>%
    # We prefix all columns with the factor name
    rename_with(.cols = -name, ~ paste0("mean_int", "_", i, "_", .x))
}


flat_dfList <- reduce(dfList, full_join, by = "name")

# We now add the raw feature list to the dataframe

full_flat_dfList <- merge(flat_dfList, t(DE_original$data), by.x = "name", by.y = "row.names", all = TRUE)


#  We add the raw feature list

df_from_graph_vertices_plus_plus <- merge(df_from_graph_vertices_plus, full_flat_dfList, by.x = "feature_id", by.y = "name", all.x = T)


# We set back the id column as the first column of the dataframe

df_from_graph_vertices_plus_plus <- df_from_graph_vertices_plus_plus %>%
  select(feature_id, everything())

node_size <- df_from_graph_vertices_plus_plus$MeanDecreaseGini
node_size[is.na(node_size)] <- min(node_size, na.rm = T)

df_from_graph_vertices_plus_plus$node_size <- node_size
# We then add the attributes to the edges dataframe and generate the igraph object

# In the case when we have been filtering the X data we will add the filtered X data to the vertices dataframe prior to merging.


# We then make sure to have the df_from_graph_vertices_plus_plus dataframe ordered by decreasing value of the MeanDecreaseGini column

df_from_graph_vertices_plus_plus <- df_from_graph_vertices_plus_plus %>%
  arrange(desc(MeanDecreaseGini))


generated_g <- graph_from_data_frame(df_from_graph_edges, directed = FALSE, vertices = df_from_graph_vertices_plus_plus)


################################################################################
################################################################################
##### add annotations to igraph


# The file is exported


write_graph(generated_g, file = filename_graphml, format = "graphml")


message("... the R session info file ...")

sink(filename_session_info)
sessionInfo()
sink()

message("... and the R script file !")

print(getwd())

setwd(script_path)

# Print the current wd
print(getwd())

get_filename <- function() {
  c_args <- commandArgs()
  r_file <- c_args[grepl("\\.R$", c_args, ignore.case = TRUE)]
  r_file <- gsub("--file=", "", r_file)
  r_file <- normalizePath(r_file)
  return(r_file)
}

script_name <- get_filename()

file.copy(script_name, file.path(output_directory, filename_R_script), overwrite = TRUE)


message("Done !")




######### run cytoscape
######### Mind that here you need to have a Cytoscape instance up and running in parralel !

# We launch the Cytoscape connector conditionally 

if (params$actions$run_cytoscape_connector == TRUE) {
  message("Running Cytoscape connector ...")
  # We launch the Cytoscape connector

library(RCy3) # Should be launched at the end as else it will cause conflicts with the previously loaded packages.
library(pdftools)
setwd(output_directory)

cytoscapePing()

createNetworkFromIgraph(generated_g, "myIgraph")



setNodeSizeMapping("node_size", sizes = c(20, 150), style.name = "Solid")


colnames_piechart <- paste("mean_int", params$target$sample_metadata_header, names(custom_colors), sep = "_")
colnames_piechart <- gsub(" ", ".", colnames_piechart)

# Note that since we have a named vectors. The vector needs to be unnamed when passed to the colors argument

setNodeCustomPieChart(colnames_piechart,
  colors = unname(custom_colors)
, style.name = "Solid"
)

# The edge label is set 

setEdgeLabelMapping('EdgeAnnotation', style.name = "Solid")

# The Node label is set

setNodeLabelMapping('parent_mass_gnps', style.name = "Solid")



saveSession("cytoscape_piechart")
exportPDF(filename = "cytoscape_piechart", pageSize = "Auto")


size_gradient <- c(1:4)
pdf("color_legend.pdf")
# Create a legend
plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
legend("topleft",
  legend = names(custom_colors), pch = 15, pt.cex = 2, cex = 1, bty = "n",
  col = custom_colors
)
mtext(params$target$sample_metadata_header, at = 0, cex = 1)
legend("topright",
  legend = size_gradient, pch = 16, pt.cex = size_gradient * 0.7, cex = 1, bty = "n",
  col = "black"
)
mtext("Importance", at = 0.95, cex = 1)
dev.off()


files <- c("cytoscape_piechart.pdf", "color_legend.pdf") # get pdf filenames
pdfs <- Reduce(c, lapply(files, image_read_pdf)) # read in and combine
montage <- image_montage(pdfs, tile = "1x2", geometry = "x1200") # create pages of 4
image_write(montage, format = "pdf", "cytoscape_piechart_legend.pdf") # save as single pdf

file.remove("color_legend.pdf")



# Read the PDF files individually
cytoscape_piechart <- image_read_pdf("cytoscape_piechart.pdf", density = 300)
color_legend <- image_read_pdf("color_legend.pdf", density = 300)

# Check the number of images/pages in each file
print(image_info(cytoscape_piechart))
print(image_info(color_legend))

montage <- image_montage(c(cytoscape_piechart, color_legend), tile = "1x2")
# Example with a simple border specified, adjust values as needed
montage <- image_montage(c(cytoscape_piechart, color_legend))

# Save the combined content as a single PDF
image_write(montage, path = "cytoscape_piechart_legend4.pdf", format = "pdf")


}