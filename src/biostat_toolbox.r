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
# BiocManager::install("structToolbox")

## install additional bioc packages for vignette if needed
# BiocManager::install(c("pmp", "ropls", "BiocFileCache"))

## install additional CRAN packages if needed
# install.packages(c('cowplot', 'openxlsx'))


# We define the following helper function in order to load or install the packages according to the condition

usePackage = function(p) {
  if (!is.element(p, installed.packages()[, 1])) {
    install.packages(p, dep = TRUE)
  }
  require(p, character.only = TRUE)
}
 
# This one below is to the the default CRAN repo

r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
rm(r)

# To organize alphabetically

usePackage("ape")
usePackage("base")
usePackage("BiocFileCache")
usePackage("cowplot")
usePackage("data.table")
usePackage("dbscan")
usePackage("tools")
usePackage("dplyr")
usePackage("emmeans")
usePackage("EnhancedVolcano")
usePackage("fpc")
usePackage("funModeling")
usePackage("ggdendro")
usePackage("ggplot2")
usePackage("ggraph")
usePackage("ggrepel")
usePackage("ggtree")
usePackage("graphlayouts")
usePackage("gridExtra")
usePackage("gt")
usePackage("heatmaply")
usePackage("here")
usePackage("igraph")
usePackage("manhattanly")
usePackage("openxlsx")
usePackage("plotly")
usePackage("pls")
usePackage("pmp")
usePackage("purrr")
usePackage("randomcoloR")
usePackage("randomForest")
usePackage("rcdk")
usePackage("readr")
usePackage("reticulate")
usePackage("rfPermute")
usePackage("rgl")
usePackage("rockchalk")
usePackage("ropls")
usePackage("structToolbox")
usePackage("this.path")
usePackage("tidyr")
usePackage("tidyverse")
usePackage("vegan")
usePackage("viridis")
usePackage("webchem")
usePackage("wesanderson")
usePackage("yaml")
usePackage("ggh4x")
usePackage("iheatmapr")


# We use the MAPPstructToolbox package 
# Uncomment the lines below to download the MAPPstructToolbox package from github

# library(devtools)
# install_github("mapp-metabolomics-unit/MAPPstructToolbox", force = TRUE)
library(MAPPstructToolbox)


############################################################################################
############################################################################################
################################ LOAD & FORMAT  DATA  ######################################
############################################################################################
############################################################################################


# We set the wd
current_script <- deparse(substitute())
script_path <- file.path(getwd(), current_script)

print(script_path)

if(!exists("params"))   {my_path_params <- script_path} ### conserve the path after multiple run
if(exists("params"))   {setwd(my_path_params)} ### conserve the path after multiple run

# We call the external params
path_to_params = "./params/params.yaml"
path_to_params_user = "./params/params_user.yaml"

params = yaml.load_file(path_to_params)
params_user = yaml.load_file(path_to_params_user)

# Here we load the user params if they exist

params$paths$docs = params_user$paths$docs
params$paths$output = params_user$paths$output
params$operating_system$system = params_user$operating_system$system
params$operating_system$pandoc = params_user$operating_system$pandoc

# We set the working directory

working_directory = file.path(params$paths$docs, params$mapp_project, params$mapp_batch)

# We set the output directory 

if (params$actions$scale_data == "TRUE") {
scaling_status = "scaled"
} else { scaling_status = "raw" }

if (params$actions$filter_sample_metadata_one == "TRUE" & params$actions$filter_sample_metadata_two == "TRUE") {
filter_sample_metadata_status = paste(params$filter_sample_metadata_one$mode,
params$filter_sample_metadata_one$factor_name,
paste(params$filter_sample_metadata_one$levels, collapse = "_"),
params$filter_sample_metadata_two$mode,
params$filter_sample_metadata_two$factor_name,
paste(params$filter_sample_metadata_two$levels, collapse = "_"),
sep = "_")
} else if (params$actions$filter_sample_metadata_one == "TRUE") {
filter_sample_metadata_status = paste(params$filter_sample_metadata_one$mode,
params$filter_sample_metadata_one$factor_name,
paste(params$filter_sample_metadata_one$levels, collapse = "_"),
sep = "_") 
} else { filter_sample_metadata_status = "" }


if (params$actions$filter_variable_metadata_one == "TRUE" & params$actions$filter_variable_metadata_two == "TRUE") {
filter_variable_metadata_status = paste(params$filter_variable_metadata_one$mode,
params$filter_variable_metadata_one$factor_name,
paste(params$filter_variable_metadata_one$levels, collapse = "_"),
params$filter_variable_metadata_two$mode,
params$filter_variable_metadata_two$factor_name,
paste(params$filter_variable_metadata_two$levels, collapse = "_"),
sep = "_")
} else if (params$actions$filter_variable_metadata_one == "TRUE") {
filter_variable_metadata_status = paste(params$filter_variable_metadata_one$mode,
params$filter_variable_metadata_one$factor_name,
paste(params$filter_variable_metadata_one$levels, collapse = "_"),
sep = "_") 
} else { filter_variable_metadata_status = "" }


if (params$actions$filter_variable_metadata_annotated == "TRUE") {
filter_variable_metadata_status = paste0(filter_variable_metadata_status, "_only_chebi_annotated")
}



#################################################################################################
#################################################################################################
################### Filename and paths establishment ##########################################
#################################################################################################


# The Figures filename is conditionally defined according to the user's choice of filtering the dataset according to CANOPUS NPClassifier classifications or not.


#file_prefix = paste(params$mapp_batch, 
#                    params$target$sample_metadata_header, 
#                    filter_variable_metadata_status, 
#                    filter_sample_metadata_status, 
#                    params$polarity, 
#                    scaling_status, 
#                    sep = "_")

# file_prefix = paste(params$mapp_batch, 
#                     params$target$sample_metadata_header,
#                     sep = "_")

file_prefix = paste("")


filename_box_plots <- paste(file_prefix, "Boxplots.pdf", sep = "")
filename_box_plots_interactive <- paste(file_prefix, "Boxplots_interactive.html", sep = "")
filename_DE_model <- paste(file_prefix, "DE_description.txt", sep = "")
filename_foldchange_pvalues <- paste(file_prefix, "foldchange_pvalues.csv", sep = "")
filename_formatted_peak_table <- paste(file_prefix, "formatted_peak_table.csv", sep = "")
filename_formatted_sample_data_table <- paste(file_prefix, "formatted_sample_data_table.csv", sep = "")
filename_formatted_sample_metadata <- paste(file_prefix, "formatted_sample_metadata.csv", sep = "")
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
filename_PLSDA_VIP <- paste(file_prefix, "PLSDA_VIP.pdf", sep = "")
filename_R_script <- paste(file_prefix, "R_script_backup.R", sep = "")
filename_random_forest <- paste(file_prefix, "RF_importance.html", sep = "")
filename_random_forest_model <- paste(file_prefix, "RF_model.txt", sep = "")
filename_session_info <- paste(file_prefix, "session_info.txt", sep = "")
filename_summary_stats_table_full <- paste(file_prefix, "summary_stats_table_full.csv", sep = "")
filename_summary_stats_table_selected <- paste(file_prefix, "summary_stats_table_selected.csv", sep = "")
filename_treemap <- paste(file_prefix, "Treemap_interactive.html", sep = "")
filename_volcano <- paste(file_prefix, "Volcano.pdf", sep = "")
filename_volcano_interactive <- paste(file_prefix, "Volcano_interactive.html", sep = "")

################################### load peak table ########################################
############################################################################################



feature_table = read_delim(file.path(working_directory,  "results", "mzmine", paste0(params$mapp_batch, "_quant.csv")),
  delim = ",", escape_double = FALSE,
  trim_ws = TRUE
)

# The column names are modified using the rename function from the dplyr package

feature_table = feature_table %>%
  rename("feature_id" = "row ID",
    "feature_mz" = "row m/z",
    "feature_rt" = "row retention time"
  )

# The row m/z and row retention time columns are concatenated to create a new column called `feature_id_full`
feature_table$"feature_id_full" = paste(feature_table$feature_id,
  round(feature_table$feature_mz, digits = 2),
  round(feature_table$feature_rt, digits = 1),
  sep = "_"
)

# The dataframe is subsetted to keep only columns containing the pattern ` Peak area` and the `feature_id_full` column
# We use dplyr's `select` function and the pipe operator `%>%` to chain the operations.
# We then remove the ` Peak area` pattern from the column names using the rename_with function from the dplyr package
# We then set the `feature_id_full` column as the rownames of the dataframe and transpose it

feature_table_intensities = feature_table %>%
  select(feature_id, contains(" Peak area")) %>%
  rename_with(~gsub(" Peak area", "", .x)) %>%
  column_to_rownames(var = "feature_id") %>%
  as.data.frame() %>%
  t()

# We keep the feature_table_intensities dataframe in a separate variable

X = feature_table_intensities


# We order the X by rownames and by column names

X = X[order(row.names(X)), ]
X = X[, order(colnames(X))]

X = as.data.frame(X)

# Min value imputation (to be checked !!!)

half_min = min(X[X > 0], na.rm = TRUE) / 2
min = min(X[X > 0], na.rm = TRUE)


X[X == 0] = min


# Uncomment for testing purposes
# X <- X[,1:100]


# We keep the feature metadata in a separate dataframe

feature_metadata = feature_table %>%
  select(feature_id_full, feature_id, feature_mz, feature_rt)

############################### load annotation tables #####################################
############################################################################################

# The Sirius data is loaded
# First we check if a chebied version exists, if not we create it

if (file.exists(file.path(working_directory, "results", "sirius", paste("chebied", params$filenames$sirius_annotations, sep = "_")))) {
  data_sirius <- read_delim(file.path(working_directory, "results", "sirius", paste("chebied", params$filenames$sirius_annotations, sep = "_")),
    delim = "\t", escape_double = FALSE,
    trim_ws = TRUE
  )
} else {
  data_sirius <- read_delim(file.path(working_directory, "results", "sirius", params$filenames$sirius_annotations),
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

  colnames(data_sirius) <- paste(colnames(data_sirius), "sirius", sep = "_")


  # We now build a unique feature_id for each feature in the Sirius data

  data_sirius$feature_id <- sub("^.*_([[:alnum:]]+)$", "\\1", data_sirius$id_sirius)
  data_sirius$feature_id <- as.numeric(data_sirius$feature_id)

  # Since this step takes time we save the output locally

  write.table(data_sirius, file = file.path(working_directory, "results", "sirius", paste("chebied", params$filenames$sirius_annotations, sep = "_")), sep = "\t", row.names = FALSE)
}

# The CANOPUS data is loaded

data_canopus = read_delim(file.path(working_directory, "results", "sirius", params$filenames$canopus_annotations),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)


# The column names are modified to include the source of the data

colnames(data_canopus) = paste(colnames(data_canopus), "canopus", sep = "_")

# We now build a unique feature_id for each feature in the Sirius data

data_canopus$feature_id = sub("^.*_([[:alnum:]]+)$", "\\1", data_canopus$id_canopus)
data_canopus$feature_id = as.numeric(data_canopus$feature_id)


write.table(data_canopus, file = file.path(working_directory, "results", "sirius", paste("featured", params$filenames$canopus_annotations, sep = "_")), sep = "\t", row.names = FALSE)

# The MetAnnot data is loaded

data_metannot = read_delim(file.path(working_directory, "results", "met_annot_enhancer", params$met_annot_enhancer_folder, paste0(params$met_annot_enhancer_folder, "_spectral_match_results_repond.tsv")),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

# The column names are modified to include the source of the data

colnames(data_metannot) = paste(colnames(data_metannot), "metannot", sep = "_")


# We now build a unique feature_id for each feature in the Metannot data

data_metannot$feature_id = data_metannot$feature_id_metannot
data_metannot$feature_id = as.numeric(data_metannot$feature_id)



# The GNPS data is loaded. Note that we use the `Sys.glob` function to get the path to the file and expand the wildcard

data_gnps = read_delim(Sys.glob(file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "clusterinfo_summary", "*.tsv")),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

# The column names are modified to include the source of the data

colnames(data_gnps) = paste(colnames(data_gnps), "gnps", sep = "_")

# We now build a unique feature_id for each feature in the GNPS data

data_gnps$feature_id = data_gnps$`cluster index_gnps`
data_gnps$feature_id = as.numeric(data_gnps$feature_id)


# The four previous dataframe are merged into one using the common `feature_id` column as key and the tidyverse `reduce` function

list_df = list(feature_metadata, data_sirius, data_canopus, data_metannot, data_gnps)
VM = list_df %>% reduce(full_join, by='feature_id')


# The row m/z and row retention time columns are concatenated to create a new column called `feature_id_full_annotated`
VM$"feature_id_full_annotated" = paste0(
  VM$chebiasciiname_sirius,
  "_[",
  VM$feature_id_full,
  "]",
  sep = ""
)

# We now convert the VM tibble into a dataframe and set the `feature_id_full` column as the rownames

VM = as.data.frame(VM)
row.names(VM) = VM$feature_id


# We order the VM by rownames and by column names

VM = VM[order(row.names(VM)), ]
VM = VM[, order(colnames(VM))]

# Uncomment for testing purposes
# VM = head(VM, 100)


################################ load sample  metadata #####################################
############################################################################################

# Later on ... implement a test stage where we check for the presence of a "species" and "sample_type" column in the metadata file.


# We here load the sample metadata


sample_metadata = read_delim(file.path(working_directory, "metadata", "treated", paste(params$mapp_batch,  "metadata.txt", sep = "_")),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Here we establish a small test which will check if the sample metadata file contains the required columns (filename, sample_id, sample_type and species)

if (!all(c("filename", "sample_id", "sample_type", "species") %in% colnames(sample_metadata))) {
  stop("The sample metadata file does not contain the required columns (filename, sample_id, sample_type and species). Please check your metadata file and try again.")
}

SM = data.frame(sample_metadata)


# We take full power over the matrix (sic. Defossez, 2023)
# First we work horizontally

for (column in names(params$to_combine_horizontally)) {

  # column = 'column1'
  col_info <- params$to_combine_horizontally[[column]]
  col_name <- col_info$name

  # Initialize aggregated groups with original condition variable
  SM[paste(col_info$name, "simplified", sep = "_")] <- as.factor(SM[[col_info$name]])

  # Iterate over each group in params$tocomb
  for (group in names(col_info$groups)) {
    # group = 'group1'
    group_info <- col_info$groups[[group]]
    cols <- group_info$cols
    new_label <- group_info$name
    
    # Combine levels for the current group
    SM[paste(col_info$name, "simplified", sep = "_")] <- combineLevels(SM[[paste(col_info$name, "simplified", sep = "_")]], levs = cols, newLabel = c(new_label))
}
}

# Then we work vertically
# The function below is used to create metadata combinations

df = SM %>% 
  filter(sample_type == "sample")

# This line allows us to make sure that the columns will be combined in alphabetical order
cols = sort(c(params$colnames$to_combine), decreasing = FALSE)

for (n in 1:length(cols)) {
  combos = combn(cols, n, simplify = FALSE)
  for (combo in combos) {
    new_col_name = paste(combo, collapse = "_")
    df[new_col_name] = apply(df[combo], 1, paste, collapse = "_")
  }
}

# We merge back the resulting df to the original SM dataframe and fill the NA values with "ND"
SM = merge(x = SM, y = df,  all.x = TRUE)
SM[is.na(SM)] = "ND"


SM = SM %>%
  remove_rownames() %>%
  column_to_rownames(var = "filename")


# This allows us to both have the filename as rownames and as a unique column
SM$filename = rownames(SM)
#cuirmoustache

# We order the SM by rownames and by column names.

SM = SM[order(row.names(SM)), ]
SM = SM[, order(colnames(SM))]


# glimpse(SM)
# Optional ponderation step.
#To clean

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

DE_original = DatasetExperiment(
  data = X,
  sample_meta = SM,
  variable_meta = VM,
  name = params$dataset_experiment$name,
  description = params$dataset_experiment$description
)


## Filtering steps

if (params$actions$filter_features == "TRUE") {

filter_by_name_model = filter_by_name(mode='exclude',dimension='variable',names=params$feature_to_filter)

# apply model sequence
filter_by_name_result = model_apply(filter_by_name_model,DE_original)
DE_filtered_name = filter_by_name_result@filtered

} else {

DE_filtered_name = DE_original

}

## Filtering steps

if (params$actions$filter_sample_type == "TRUE") {

filter_smeta_model <- filter_smeta(mode = params$filter_sample_type$mode,
                          factor_name = params$filter_sample_type$factor_name,
                          levels = params$filter_sample_type$levels)

# apply model sequence
filter_smeta_result = model_apply(filter_smeta_model, DE_filtered_name)

DE_filtered = filter_smeta_result@filtered

}

if (params$actions$filter_sample_metadata_one == "TRUE") {

filter_smeta_model <- filter_smeta(mode = params$filter_sample_metadata_one$mode,
                          factor_name = params$filter_sample_metadata_one$factor_name,
                          levels = params$filter_sample_metadata_one$levels)

# apply model sequence
filter_smeta_result = model_apply(filter_smeta_model, DE_filtered)

DE_filtered = filter_smeta_result@filtered

}

if (params$actions$filter_sample_metadata_two == "TRUE") {

filter_smeta_model <- filter_smeta(mode = params$filter_sample_metadata_two$mode,
                          factor_name = params$filter_sample_metadata_two$factor_name,
                          levels = params$filter_sample_metadata_two$levels)

# apply model sequence
filter_smeta_result = model_apply(filter_smeta_model, DE_filtered)

DE_filtered = filter_smeta_result@filtered

}



if (params$actions$filter_variable_metadata_one == "TRUE") {

filter_vmeta_model <- filter_vmeta(mode = params$filter_variable_metadata_one$mode,
                          factor_name = params$filter_variable_metadata_one$factor_name,
                          levels = params$filter_variable_metadata_one$levels)

# apply model sequence
filter_vmeta_result = model_apply(filter_vmeta_model, DE_filtered)

DE_filtered = filter_vmeta_result@filtered

}

if (params$actions$filter_variable_metadata_two == "TRUE") {

filter_vmeta_model <- filter_vmeta(mode = params$filter_variable_metadata_two$mode,
                          factor_name = params$filter_variable_metadata_two$factor_name,
                          levels = params$filter_variable_metadata_two$levels)

# apply model sequence
filter_vmeta_result = model_apply(filter_vmeta_model, DE_filtered)

DE_filtered = filter_vmeta_result@filtered

}



if (params$actions$filter_variable_metadata_annotated == "TRUE") {

# Convert the "levels" value to NA if it is "NA" as a character string
if (params$filter_variable_metadata_annotated$levels == "NA") {
  params$filter_variable_metadata_annotated$levels <- NA
}

filter_vmeta_model <- filter_vmeta(mode = params$filter_variable_metadata_annotated$mode,
                          factor_name = params$filter_variable_metadata_annotated$factor_name,
                          levels = as.character(params$filter_variable_metadata_annotated$levels)
                          )

# apply model sequence
filter_vmeta_result = model_apply(filter_vmeta_model, DE_filtered)

DE_filtered = filter_vmeta_result@filtered

}


if (params$actions$scale_data == "FALSE") {

DE = DE_filtered

# We display the properties of the DatasetExperiment object to the user.
message("DatasetExperiment object properties: ")

sink(filename_DE_model)

print(DE)

sink() } else if (params$actions$scale_data == "TRUE") {

# Overall Pareto scaling (test)

M = pareto_scale()
M = model_train(M,DE_filtered)
M = model_predict(M,DE_filtered)
DE = M$scaled

##### we range  all feature to o to 1

DE$data <- apply(DE$data,2,funModeling::range01)

# Min value imputation also after the scaling stage (to be checked !!!)

half_min_sec = min(DE$data[DE$data > 0], na.rm = TRUE) / 2

min_sec = min(DE$data[DE$data > 0], na.rm = TRUE)

DE$data[DE$data == 0] = min_sec

################################################################################################
################################################################################################
######################## structool box formatted data export 

message("Outputting X, VM and SM ...")

formatted_peak_table <- DE_filtered$data

formatted_variable_metadata <- DE_filtered$variable_meta ### need to be filter with only usefull output


# col_filter <- c("feature_id_full", "feature_id", "feature_mz" ,"feature_rt", "molecularFormula_sirius","InChIkey2D_sirius","InChI_sirius",
# "name_sirius","smiles_sirius", "pubchemids_sirius", "molecularFormula_canopus", "NPC.pathway_canopus","NPC.pathway.Probability_canopus",
# "NPC.superclass_canopus", "NPC.class_canopus","ClassyFire.most.specific.class_canopus","ClassyFire.most.specific.class.Probability_canopus",
# "ClassyFire.level.5_canopus","ClassyFire.subclass_canopus","ClassyFire.class_canopus","ClassyFire.superclass_canopus","ClassyFire.all.classifications_canopus",
# "structure_wikidata_metannot","structure_inchikey_metannot","structure_inchi_metannot","structure_smiles_metannot","structure_molecular_formula_metannot",
# "short_inchikey_metannot","structure_taxonomy_npclassifier_01pathway_metannot","structure_taxonomy_npclassifier_02superclass_metannot",
# "structure_taxonomy_npclassifier_03class_metannot","organism_wikidata_metannot","organism_name_metannot","organism_taxonomy_ottid_metannot","organism_taxonomy_01domain_metannot",
# "organism_taxonomy_02kingdom_metannot","organism_taxonomy_03phylum_metannot","organism_taxonomy_04class_metannot","organism_taxonomy_05order_metannot",
# "organism_taxonomy_06family_metannot","organism_taxonomy_07tribe_metannot","organism_taxonomy_08genus_metannot","organism_taxonomy_09species_metannot","organism_taxonomy_10varietas_metannot",
# "matched_domain_metannot","matched_kingdom_metannot","matched_phylum_metannot","matched_class_metannot","matched_order_metannot","matched_family_metannot","matched_tribe_metannot",
# "matched_genus_metannot","matched_species_metannot","score_taxo_metannot","score_max_consistency_metannot","final_score_metannot","rank_final_metannot","component_id_metannot",
# "structure_taxonomy_npclassifier_01pathway_consensus_metannot","freq_structure_taxonomy_npclassifier_01pathway_metannot","structure_taxonomy_npclassifier_02superclass_consensus_metannot",
# "freq_structure_taxonomy_npclassifier_02superclass_metannot","structure_taxonomy_npclassifier_03class_consensus_metannot","freq_structure_taxonomy_npclassifier_03class_metannot")


col_filter <- c("feature_id_full", "feature_id", "feature_mz" ,"feature_rt", "molecularFormula_sirius","freq_structure_taxonomy_npclassifier_03class_metannot")

formatted_variable_metadata_filtered = formatted_variable_metadata[col_filter]

formatted_sample_metadata = DE_filtered$sample_meta

formatted_sample_data_table = merge(DE_filtered$sample_meta, DE_filtered$data, by="row.names")



###################################################################################################
######################### rename main folder - short version
# Here we check if the params$paths$out value exist and use it else we use the default output_directory


target_name <- paste(as.vector(unique(DE$sample_meta[[params$target$sample_metadata_header]])), collapse = "_vs_")


if (params$paths$output != "") {
  output_directory <- file.path(params$paths$output, paste(params$target$sample_metadata_header, target_name, filter_variable_metadata_status, filter_sample_metadata_status, scaling_status, sep = "_"), sep = "")
} else {
  output_directory <- file.path(working_directory, "results", "stats", paste(params$target$sample_metadata_header, target_name, filter_variable_metadata_status,filter_sample_metadata_status, scaling_status, sep = "_"), sep = "")
}

output_directory <- gsub("__","_",output_directory)

dir.create(output_directory)


#################################################################################
#################################################################################
##### write raw data and param 

## We save the used params.yaml

message("Writing params.yaml ...")


file.copy(path_to_params, file.path(output_directory,filename_params), overwrite = TRUE)
file.copy(path_to_params_user, file.path(output_directory,filename_params_user), overwrite = TRUE)

# We move to the output directory

setwd(output_directory)

# We display the properties of the DatasetExperiment object to the user.
message("DatasetExperiment object properties: ")

sink(filename_DE_model)

print(DE)

sink() 
} else {
  stop("Please check the value of the 'scale_data' parameter in the params file.")
}


#######################


write.table(formatted_peak_table, file = filename_formatted_peak_table, sep = ",", row.names = FALSE)
write.table(formatted_variable_metadata_filtered, file = filename_formatted_variable_metadata, sep = ",", row.names = FALSE)
write.table(formatted_sample_metadata, file = filename_formatted_sample_metadata, sep = ",", row.names = FALSE)
write.table(formatted_sample_data_table, file = filename_formatted_sample_data_table, sep = ",", row.names = FALSE)
################################################################################################
################################################################################################

# The Figures titles are conditionally defined according to the user's choices and option in the parameters file

#title_PCA = paste("PCA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ") 
#title_PLSDA = paste("PLSDA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
#title_PLSDA_VIP = paste("PLSDA selected Features of Importance", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ") 
#title_PCA3D = paste("PCA3D", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
#title_PCoA = paste("PCoA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ") 
#title_PCoA3D = paste("PCoA3D", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$target$sample_metadata_header, sep = " ")
#title_volcano = paste("Volcano plot", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
#title_treemap = paste("Treemap", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
#title_random_forest = paste("Random Forest results", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
#title_box_plots = paste("Top", params$boxplot$topN, "boxplots", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
#title_heatmap = paste("Heatmap of","top", params$heatmap$topN,"Random Forest filtered features", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")

# title_PCA = paste("PCA", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")

title_PCA = paste(
  paste("PCA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)


# title_PLSDA = paste("PLSDA", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PLSDA = paste(
  paste("PLSDA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

title_PLSDA_VIP = paste("PLSDA selected Features of Importance", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

# title_PCA3D = paste("PCA3D", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PCA3D = paste(
  paste("PCA3D", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

# title_PCoA = paste("PCoA", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PCoA = paste(
  paste("PCoA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

# title_PCoA3D = paste("PCoA3D", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_PCoA3D = paste(
  paste("PCoA3D", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

title_volcano = paste("Volcano plot", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")
title_treemap = paste("Treemap", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")
title_random_forest = paste("Random Forest results", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")
title_box_plots = paste("Top", params$boxplot$topN, "boxplots", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")
title_heatmap_rf = paste("Heatmap of","top", params$heatmap$topN,"Random Forest filtered features", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

# title_heatmap_pval = paste("Heatmap of significant feature for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status,scaling_status, sep = " ")

title_heatmap_pval = paste(
  paste("Heatmap of significant feature for dataset", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)


#################################################################################################
#################################################################################################
#################################################################################################
##### PCA filtered data #######################################################################

message("Launching PCA calculations ...")


pca_seq_model <- filter_na_count(threshold = 1, factor_name = "sample_type") +
              knn_impute(neighbours = 5) +
              vec_norm() +
            #log_transform(base = 10) +
              mean_centre() +
              PCA(number_components = 3)

# apply model sequence
pca_seq_result = model_apply(pca_seq_model, DE)

# Fetching the PCA data object
pca_object = pca_seq_result[length(pca_seq_result)]

# PCA scores plot

pca_scores_plot = pca_scores_plot(
  factor_name = params$target$sample_metadata_header,
  label_factor = "sample_id",
  ellipse_type = "t",
  ellipse_confidence = 0.9,
  points_to_label = "all"
)

# plot
pca_plot = chart_plot(pca_scores_plot, pca_object)


# pca_plot$labels$title <- title_PCA

fig_PCA = pca_plot + 
theme_classic() + 
facet_wrap(~ pca_plot$labels$title) +
ggtitle(title_PCA)
#   theme(plot.title = element_text(hjust = 0.2, vjust = -2)) +


# We merge PCA scores and metadata info in a single df

PCA_meta = merge(x = pca_object$scores$sample_meta, y = pca_object$scores$data, by = 0, all = TRUE)


fig_PCA3D = plot_ly(PCA_meta, x = ~PC1, y = ~PC2, z = ~PC3, color = PCA_meta[,params$target$sample_metadata_header])
fig_PCA3D = fig_PCA3D %>% add_markers()
fig_PCA3D = fig_PCA3D %>% layout(scene = list(
  xaxis = list(title = "PC1"),
  yaxis = list(title = "PC2"),
  zaxis = list(title = "PC3")
),
legend = list(title=list(text=params$target$sample_metadata_header)),
title = title_PCA3D
)


# The files are exported

ggsave(plot = fig_PCA, filename = filename_PCA , width = 10, height = 10)


if (params$operating_system$system == "unix") {
### linux version
fig_PCA3D %>%
    htmlwidgets::saveWidget(file = filename_PCA3D, selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
### windows version
Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
fig_PCA3D %>%
    htmlwidgets::saveWidget(file = filename_PCA3D, selfcontained = TRUE,libdir = "lib")
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

DE$sample_meta[,params$target$sample_metadata_header] = as.factor(DE$sample_meta[,params$target$sample_metadata_header])

# glimpse(DE_filtered$sample_meta)

# check the outcome of a Pareto scaling methods


# # prepare model sequence
plsda_seq_model = # autoscale() +
                  filter_na_count(threshold=3,factor_name=params$target$sample_metadata_header) +
                  # knn_impute() +
                  PLSDA(factor_name=params$target$sample_metadata_header, number_components=2)

plsda_seq_result = model_apply(plsda_seq_model,DE)

# Fetching the PLSDA data object
plsda_object = plsda_seq_result[length(plsda_seq_result)]

C = pls_scores_plot(factor_name = params$target$sample_metadata_header)

plsda_plot = chart_plot(C,plsda_object)


fig_PLSDA = plsda_plot + theme_classic() + facet_wrap(~ plsda_plot$labels$title) + ggtitle(title_PLSDA)

# We output the feature importance

C = plsda_feature_importance_plot(n_features=30, metric='vip')

vip_plot <- chart_plot(C,plsda_object)

fig_PLSDA_VIP = vip_plot + theme_classic() + facet_wrap(~ plsda_plot$labels$title) + ggtitle(title_PLSDA_VIP)

# The files are exported

ggsave(plot = fig_PLSDA, filename = filename_PLSDA , width = 10, height = 10)
ggsave(plot = fig_PLSDA_VIP, filename = filename_PLSDA_VIP , width = 10, height = 10)

}

#################################################################################################
#################################################################################################
#################################################################################################
##### HClustering
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
##### PCoA  

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



data_RF = DE# DE_filtered
sample_name = data_RF$sample_meta$sample_id #### check
data_subset_norm_rf = data_RF$data
data_subset_norm_rf[sapply(data_subset_norm_rf, is.infinite)] = NA
data_subset_norm_rf[is.na(data_subset_norm_rf)] = 0


dist_metabo = vegdist(data_subset_norm_rf, method = "bray") # method="man" # is a bit better
# Why dont we use it then ???
D3_data_dist = cmdscale(dist_metabo, k = 3)
D3_data_dist = data.frame(D3_data_dist)
D3_data_dist$sample_name = sample_name
D3_data_dist = D3_data_dist[order(D3_data_dist$sample_name), ]
metadata_merge = data_RF$sample_meta[order(sample_name), ]

data_PCOA_merge = data.frame(cbind(D3_data_dist, metadata_merge))

cols = data_PCOA_merge[params$target$sample_metadata_header]
cols = cols[, 1]


fig_PCoA = ggplot(data_PCOA_merge, aes(x = X1, y = X2, color = cols)) +
  geom_point() +
  ggtitle(title_PCoA) +
  theme_classic()


fig_PCoA3D = plot_ly(
  x = data_PCOA_merge$X1, y = data_PCOA_merge$X2, z = data_PCOA_merge$X3,
  type = "scatter3d", mode = "markers", color = cols,
  hoverinfo = "text",
  text = ~ paste(
    "</br> name: ", data_PCOA_merge$sample_name,
    "</br> num: ", data_PCOA_merge$sample_id
  )
) %>% layout(title = title_PCoA3D,
legend = list(title=list(text=params$target$sample_metadata_header)))


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
    htmlwidgets::saveWidget(file = filename_PCoA3D, selfcontained = TRUE,libdir = "lib")
unlink("lib", recursive = FALSE)

}



#################################################################################################
#################################################################################################
#################################################################################################
##### group detection

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
##### Volcano plot and Heatmap filtered by Random Forest

message("Launching Fold Changes and Tukey’s Honest Significant Difference calculations ...")


# # prepare model sequence

# MS_heatmap = filter_smeta(mode = "include", levels = params$filters$to_include, factor_name = "sample_type") +
#  #log_transform(base = 10) +
#   filter_by_name(mode = "include", dimension = "variable", names = names_var)


# # apply model sequence

# DE_MS_heat = model_apply(MS_heatmap, DE)

# DE_MS_heat = DE_MS_heat[length(DE_MS_heat)]

######################################################
######################################################
################# heat filter

# data_RF = DE
# sample_name = paste(data_RF$sample_meta$sample_id, data_RF$sample_meta[[params$target$sample_metadata_header]], sep = "_")

# data_subset_norm_rf = data_RF$data
# data_subset_norm_rf[sapply(data_subset_norm_rf, is.infinite)] = NA
# data_subset_norm_rf[is.na(data_subset_norm_rf)] = 0
# # data_subset_norm_rf = normalize(data_subset_norm_rf)
# # data_subset_norm_rf[is.na(data_subset_norm_rf)] = 0

#############################################################################
#############################################################################
############## Volcano Plot #################################################
#############################################################################
#############################################################################


# matt_trait = data_RF$sample_meta[params$target$sample_metadata_header]
# vec_trait = matt_trait[, 1]
# data_subset_norm_rf$treatment = as.factor(vec_trait) ### select the variable
# data_subset_norm_rf = data.frame(data_subset_norm_rf)

# matt_volcano_posthoc = NULL
# matt_volcano_lm_sum = NULL
# for (i in c(1:(ncol(data_subset_norm_rf) - 1))) {
#   peak = data_subset_norm_rf[, i]
#   if (sum(peak, na.rm = T) == 0) {
#     next
#   }
#   treatment = data_subset_norm_rf$treatment
#   dat = data.frame(peak, treatment)
#   model = lm(peak ~ treatment, data = dat)
#   em = emmeans(model, list(pairwise ~ treatment), adjust = "tukey")

#   # Summary of the analysis posthoc
#   xx = data.frame(em$`pairwise differences of treatment`)
#   xx$mol = rep(colnames(data_subset_norm_rf)[i], nrow(xx))
#   matt_volcano_posthoc = rbind(matt_volcano_posthoc, xx)
# }

# matt_volcano_posthoc$log10P = 1 - (log10(matt_volcano_posthoc$p.value))

# matt_volcano_tot = matt_volcano_posthoc

# head(DE$sample_meta)



### Here we wil work on outputting pvalues and fc for time series.

# We build a for loop to iterate over the different time points
# This loop generate a set of DE results for each time point

# params = yaml.load_file('/Users/pma/Dropbox/git_repos/mapp-metabolomics-unit/biostat_toolbox/params/params.yaml')


if (params$actions$calculate_multi_series_fc == 'TRUE') {


l <- list()

for (i in params$multi_series$points) {

  print(i)
  filter_smeta_model <- filter_smeta(mode = 'include',
                            factor_name = params$multi_series$colname,
                            levels = i)

  # apply model sequence
  filter_smeta_result = model_apply(filter_smeta_model, DE)

  DE_tp = filter_smeta_result@filtered
  # assign(paste("DE_filtered", i, sep = "_"), filter_smeta_result@filtered)

  # The formula is defined externally
  formula = as.formula(paste0('y', '~', params$target$sample_metadata_header, '+' ,
  'Error(sample_id/',
  params$target$sample_metadata_header,
  ')'
  )
  )

  # DE$sample_meta

  HSDEM_model = HSDEM(alpha = params$posthoc$p_value,
  formula = formula, mtc = 'none')

  HSDEM_result = model_apply(HSDEM_model,DE_tp)

  HSDEM_result_p_value = HSDEM_result$p_value

  # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.

  colnames(HSDEM_result_p_value) = plotrix::pasteCols(sapply(strsplit(colnames(HSDEM_result_p_value), " - "), sort), sep = "_")

  # We now add a specific suffix (`_p_value`) to each of the colnames

  colnames(HSDEM_result_p_value) = paste0("tp_", i, "_", colnames(HSDEM_result_p_value), "_p_value")

  p_value_column = colnames(HSDEM_result_p_value)

  # We set the row names as columns row_id to be able to merge the two dataframes

  HSDEM_result_p_value$row_id = rownames(HSDEM_result_p_value)


  # We build a fold change model

  fold_change_model = fold_change(
    factor_name = params$target$sample_metadata_header,
    paired = FALSE,
    sample_name = character(0),
    threshold = 0.5,
    control_group = character(0),
    method = "mean",
    conf_level = 0.95
    )


  fold_change_result = model_apply(fold_change_model, DE_tp)

  #view(DE$data)
  #DE$data[,2]  <- c(-500,-500,-500,-500,500,500,500,500)

  # We suffix the column name of the dataframe with `_fold_change`, using dplyr rename function

  fold_change_result_fold_change = fold_change_result$fold_change

  # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.
  #n !!!! We need to make sure that the header of metadata variable is not in the colnames of the fold change result


  colnames(fold_change_result_fold_change) = plotrix::pasteCols(sapply(strsplit(colnames(fold_change_result_fold_change), "/"), sort), sep = "_")


  # We now add a specific suffix (`_p_value`) to each of the colnames

  colnames(fold_change_result_fold_change) = paste0("tp_", i, "_", colnames(fold_change_result_fold_change), "_fold_change")


  fc_column = colnames(fold_change_result_fold_change)


  # We set the row names as columns row_id to be able to merge the two dataframes

  fold_change_result_fold_change$row_id = rownames(fold_change_result_fold_change)

  # # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

  # fold_change_result_fold_change = pivot_longer(fold_change_result_fold_change, cols = -row_id, names_to = "pairs", values_to = "fold_change")


  # We merge the two dataframes according to both the row_id and the pairs columns. 

  DE_foldchange_pvalues = merge(HSDEM_result_p_value, fold_change_result_fold_change,  by = "row_id")


  # We add columns corresponding to the Log2 of the fold change column (suffix by fold_change). For this we use mutate_at function from dplyr package. We save the results in new columns with a novel suffix `_log2_FC`.

  message("Calculating logs ...")

  DE_foldchange_pvalues = DE_foldchange_pvalues %>%
      mutate( across(contains('_fold_change'), 
                      .fns = list(log2 = ~log2(.)),
                      .names = "{col}_{fn}" ) )  %>% 
      mutate( across(contains('_p_value'), 
                      .fns = list(minus_log10 = ~-log10(.)),
                      .names = "{col}_{fn}" ) )


  l[[i]] <- DE_foldchange_pvalues

}

# We now merge the different dataframes in the list l

DE_foldchange_pvalues = Reduce(function(x, y) merge(x, y, by = "row_id"), l)

} else {


# The formula is defined externally
formula = as.formula(paste0('y', '~', params$target$sample_metadata_header, '+' ,
'Error(sample_id/',
 params$target$sample_metadata_header,
 ')'
)
)


# DE$sample_meta

HSDEM_model = HSDEM(alpha = params$posthoc$p_value,
formula = formula, mtc = 'none')

HSDEM_result = model_apply(HSDEM_model,DE)

HSDEM_result_p_value = HSDEM_result$p_value


# We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.

colnames(HSDEM_result_p_value) = plotrix::pasteCols(sapply(strsplit(colnames(HSDEM_result_p_value), " - "), sort), sep = "_vs_")

# We now add a specific suffix (`_p_value`) to each of the colnames

colnames(HSDEM_result_p_value) = paste0(colnames(HSDEM_result_p_value), "_p_value")

p_value_column = colnames(HSDEM_result_p_value)


# We set the row names as columns row_id to be able to merge the two dataframes

HSDEM_result_p_value$row_id = rownames(HSDEM_result_p_value)

# # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

# HSDEM_result_p_value_long = pivot_longer(HSDEM_result_p_value, cols = -row_id, names_to = "pairs", values_to = "p_value")



fold_change_model = fold_change(
  factor_name = params$target$sample_metadata_header,
  paired = FALSE,
  sample_name = character(0),
  threshold = 0.5,
  control_group = character(0),
  method = "mean",
  conf_level = 0.95
  )

# Check if this can be important to apply.

DE_fc = DE

DE_fc$data = DE_fc$data + 1

fold_change_result = model_apply(fold_change_model, DE_fc)

#view(DE$data)
#DE$data[,2]  <- c(-500,-500,-500,-500,500,500,500,500)

# We suffix the column name of the dataframe with `_fold_change`, using dplyr rename function

fold_change_result_fold_change = fold_change_result$fold_change

# We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.
#n !!!! We need to make sure that the header of metadata variable is not in the colnames of the fold change result


colnames(fold_change_result_fold_change) = plotrix::pasteCols(sapply(strsplit(colnames(fold_change_result_fold_change), "/"), sort), sep = "_vs_")


# We now add a specific suffix (`_p_value`) to each of the colnames

colnames(fold_change_result_fold_change) = paste0(colnames(fold_change_result_fold_change), "_fold_change")


fc_column = colnames(fold_change_result_fold_change)


# We set the row names as columns row_id to be able to merge the two dataframes

fold_change_result_fold_change$row_id = rownames(fold_change_result_fold_change)

# # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

# fold_change_result_fold_change = pivot_longer(fold_change_result_fold_change, cols = -row_id, names_to = "pairs", values_to = "fold_change")


# We merge the two dataframes according to both the row_id and the pairs columns. 

DE_foldchange_pvalues = merge(HSDEM_result_p_value, fold_change_result_fold_change,  by = "row_id")


# We add columns corresponding to the Log2 of the fold change column (suffix by fold_change). For this we use mutate_at function from dplyr package. We save the results in new columns with a novel suffix `_log2_FC`.

message("Calculating logs ...")

DE_foldchange_pvalues = DE_foldchange_pvalues %>%
     mutate( across(contains('_fold_change'), 
                    .fns = list(log2 = ~log2(.)),
                    .names = "{col}_{fn}" ) )  %>% 
     mutate( across(contains('_p_value'), 
                    .fns = list(minus_log10 = ~-log10(.)),
                    .names = "{col}_{fn}" ) )
}

# We now merge the DE_foldchange_pvalues with the variable metadata using the row_ID column and the rownames of the variable metadata

DE_foldchange_pvalues = merge(DE_foldchange_pvalues, DE$variable_meta, by.x = "row_id", by.y = "row.names")


# The file is exported

write.table(DE_foldchange_pvalues, file = filename_foldchange_pvalues, sep = ",")


# glimpse(DE_foldchange_pvalues)

# # Using this datatable we prepare a Volcano plot using the volcano_plot function from the plotly package. We use the p_value_log10 and the fold_change_log2_FC columns as x and y axis respectively. We use the row_ID column as the label of the points.

# vp <- EnhancedVolcano(DE_foldchange_pvalues,
#     lab = row_ID,
#     x = 'old_young_fold_change_log2',
#     y = 'old_young_p_value',
#     title = 'Old versus Young',
#     pCutoff = 10e-2,
#     FCcutoff = 0.5,
#     pointSize = 3.0,
#     labSize = 6.0)
# DE_foldchange_pvalues

# x <- filter( DE_foldchange_pvalues, !is.na(old_young_p_value))
# p <- ggplot(data=x, aes(x=old_young_fold_change_log2, y= -log10(old_young_p_value), text=row_ID )) +
#      geom_point(alpha=0.3, size=1.5, color="blue") +
#      xlab("Log2 fold change") + ylab("-Log10 p-value") +xlim(-6,6)

# y <- filter(x, old_young_p_value < 1e-5)
# p + geom_text( data=y, aes(x=old_young_fold_change_log2, y= -log10(old_young_p_value), label=row_ID),
#        hjust="left", nudge_x=.1)

# ggplotly(p)



# # keep only the fields needed for the plot
# diff_df <- DE_foldchange_pvalues[c("row_ID", "old_young_fold_change_log2", "old_young_p_value")]

# # add a grouping column; default value is "not significant"
# diff_df["group"] <- "NotSignificant"

# # for our plot, we want to highlight 
# # p value < 0.05
# # Fold Change > 1

# # change the grouping for the entries with significance but not a large enough Fold change
# diff_df[which(diff_df['old_young_p_value'] < 0.05 & abs(diff_df['old_young_fold_change_log2']) < 1 ),"group"] <- "Significant"

# # change the grouping for the entries a large enough Fold change but not a low enough p value
# diff_df[which(diff_df['old_young_p_value'] > 0.05 & abs(diff_df['old_young_fold_change_log2']) > 1 ),"group"] <- "FoldChange"

# # change the grouping for the entries with both significance and large enough fold change
# diff_df[which(diff_df['old_young_p_value'] < 0.05 & abs(diff_df['old_young_fold_change_log2']) > 1 ),"group"] <- "Significant&FoldChange"

# # make the Plot.ly plot
# p <- plot_ly(data = diff_df, x = ~old_young_fold_change_log2, y = ~-log10(old_young_p_value), text = row_ID, mode = "markers", color = ~group) %>% 
#   layout(title ="DiffBind Volcano Plot")
# p


# C = fold_change_plot(number_features = 20, orientation = 'landscape')
# chart_plot(C,fold_change_result)



# cols = c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey")
# sizes = c("up" = 2, "down" = 2, "ns" = 1)
# alphas = c("up" = 1, "down" = 1, "ns" = 0.5)

# volc_annot = data_RF$variable_meta[c("row_ID", "NPC.superclass_canopus", "NPC.pathway_canopus")]


# matt_volcano_tot$mol = gsub("X", "", matt_volcano_tot$mol)
# matt_volcano_plot = merge(matt_volcano_tot, volc_annot, by.x = "mol", by.y = "row_ID")


# matt_volcano_plot$lab_plotly = matt_volcano_plot$NPC.superclass_canopus
# matt_volcano_plot$lab_plotly[matt_volcano_plot$p.value > 0.05] = NA
# matt_volcano_plot$col_points = rep("darkred", nrow(matt_volcano_plot))
# matt_volcano_plot$col_points[matt_volcano_plot$p.value < 0.05] = "darkgreen"



# fig_volcano = ggplot(data = matt_volcano_plot, aes(x = estimate, y = -log10(p.value), color = X1, label = NPC.superclass_canopus)) +
#   geom_point() + # color = matt_volcano_plot$col_points
#   theme_minimal() +
#   geom_label_repel(max.overlaps = 10) + ### control the number of include annoation
#   # geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#   geom_hline(yintercept = -log10(0.05), col = "red") +
#   # scale_color_manual(values= wes_palette("Darjeeling1")) +
#   ggtitle(title_volcano)






# colnames(DE_foldchange_pvalues)



# fig_volcano = ggplot(data = DE_foldchange_pvalues, aes(x = lo, y = -log10(p.value), color = X1, label = NPC.superclass_canopus)) +
#   geom_point() + # color = matt_volcano_plot$col_points
#   theme_minimal() +
#   geom_label_repel(max.overlaps = 10) + ### control the number of include annoation
#   # geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#   geom_hline(yintercept = -log10(0.05), col = "red") +
#   # scale_color_manual(values= wes_palette("Darjeeling1")) +
#   ggtitle(title_volcano)



# # We need to have the feature id displayed on the Volcano plots !
# # The label = paste(lab_plotly, mol, sep ="_")) works but is not ideal

# fig_volcano_interactive = ggplot(data = matt_volcano_plot, aes(x = estimate, y = -log10(p.value), color = X1, label = paste(lab_plotly, mol, sep ="_"))) +
#   geom_point(alpha = 0.2) + # color = matt_volcano_plot$col_points
#   theme_minimal() +
#   scale_color_discrete(name = "Compared groups") +
#   # geom_text(position=position_jitter(width=0.2),size=1,col="black") +### control the number of include annoation
#   # geom_vline(xintercept=c(-0.6, 0.6), col="red") +
#   geom_hline(yintercept = -log10(0.05), col = "red") +
#   # scale_color_manual(values= wes_palette("Darjeeling1"))+
#   ggtitle(title_volcano)

# fig_volcano_interactive = ggplotly(fig_volcano_interactive)


# # The files are exported

# ggsave(plot = fig_volcano, filename = filename_volcano , width = 10, height = 10)

# fig_volcano_interactive %>%
#     htmlwidgets::saveWidget(file = filename_volcano_interactive , selfcontained = TRUE)
##############################################################################
##############################################################################
############ treemap fold 

################################### function 
################################# npc_classifier_parser



treat_npclassifier_json <- function(taxonomy) {
     taxonomy_classes <- taxonomy$Class %>%
     rbind()
     rownames(taxonomy_classes) <- "id_class"
     taxonomy_classes <- taxonomy_classes %>%
     t() %>%
     data.frame() %>%
     mutate(
     class = rownames(.),
     id_class = as.numeric(id_class))

  taxonomy_superclasses <- taxonomy$Superclass %>%
      rbind()
      rownames(taxonomy_superclasses) <- "id_superclass"
      taxonomy_superclasses <- taxonomy_superclasses %>%
      t() %>%
      data.frame() %>%
      mutate(
      superclass = rownames(.),
      id_superclass = as.numeric(id_superclass))

  taxonomy_pathways <- taxonomy$Pathway %>%
      rbind()
      rownames(taxonomy_pathways) <- "id_pathway"
      taxonomy_pathways <- taxonomy_pathways %>%
      t() %>%
      data.frame() %>%
      mutate(
      pathway = rownames(.),
      id_pathway = as.numeric(id_pathway))

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

dt_for_treemap = function(datatable, parent_value, value, count) {
  parent_value = enquo(parent_value)
  value = enquo(value)
  count = enquo(count)

  datatable = data.frame(datatable %>%
    group_by(!!parent_value, !!value, ) %>%
    summarise(count = sum(as.numeric(!!count),na.rm=T)))

  datatable = datatable %>%
    select(!!parent_value, !!value, count) %>% # create id labels for each row # Notre the !! to pass aruguments to a dplyr function
    rename(
      parent.value = !!parent_value,
      value = !!value
    ) %>%
    mutate(ids = ifelse(parent.value == "", value,
      paste0(value, "-", parent.value) # Notre that here we are passing argument to a non dplyr function call
    )) %>%
    select(ids, everything())

  par_info = datatable %>% dplyr::group_by(parent.value) %>% # group by parent
    dplyr::summarise(count = sum(as.numeric(count),na.rm=T)) %>% # parent total
    rename(value = parent.value) %>% # parent labels for the item field
    mutate(parent.value = "", ids = value) %>% # add missing fields for my_data
    select(names(datatable)) # put cols in same order as my_data

  data_for_plot = rbind(datatable, par_info)

  return(data_for_plot)
}
# ###################################################################################
# ###################################################################################

dt_for_treemap_mean = function(datatable, parent_value, value, count) {
  parent_value = enquo(parent_value)
  value = enquo(value)
  count = enquo(count)

  datatable = data.frame(datatable %>%
    group_by(!!parent_value, !!value, ) %>%
    summarise(count = mean(as.numeric(!!count),na.rm=T)))

  datatable = datatable %>%
    select(!!parent_value, !!value, count) %>% # create id labels for each row # Notre the !! to pass aruguments to a dplyr function
    rename(
      parent.value = !!parent_value,
      value = !!value
    ) %>%
    mutate(ids = ifelse(parent.value == "", value,
      paste0(value, "-", parent.value) # Notre that here we are passing argument to a non dplyr function call
    )) %>%
    select(ids, everything())

  par_info = datatable %>% dplyr::group_by(parent.value) %>% # group by parent
    dplyr::summarise(count = mean(as.numeric(count),na.rm=T)) %>% # parent total
    rename(value = parent.value) %>% # parent labels for the item field
    mutate(parent.value = "", ids = value) %>% # add missing fields for my_data
    select(names(datatable)) # put cols in same order as my_data

  data_for_plot = rbind(datatable, par_info)

  return(data_for_plot)
}

if (params$actions$run_fc_treemaps == 'TRUE') {

  cat("Great ! You decided to launch the fc treemaps calculations :) :\n")
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
  colnames(npclassifier_newpath) <-  c("NPC.class_canopus","NPC.superclass_canopus","NPC.pathway_canopus")
  npclassifier_newpath$NPC.superclass_canopus[grep(" x ",npclassifier_newpath$NPC.pathway_canopus)] <- paste(npclassifier_newpath$NPC.superclass_canopus[grep(" x ",npclassifier_newpath$NPC.pathway_canopus)],"x")

  # Here we list the distinct values in the npclassifier_newpath$NPC.pathway_canopus and order them alphabetically

  # npclassifier_newpath  %>%  
  # distinct(NPC.pathway_canopus)  %>% 
  # arrange(NPC.pathway_canopus)


  index <- sort(unique(paste(npclassifier_newpath$NPC.superclass_canopus,npclassifier_newpath$NPC.pathway_canopus)))


  # DE_foldchange_pvalues_signi <- DE_foldchange_pvalues[DE_foldchange_pvalues$C_WT_p_value < 0.05,]


  # Extract prefixes of columns with "_p_value" suffix
  conditions <- sub("_p_value$", "", grep("_p_value$", names(DE_foldchange_pvalues), value = TRUE))

  # Print message before iterating over conditions
  cat("Iterating over the following conditions:\n")


  # Iterate over the prefixes
  for (condition in conditions) {

    # Perform filtering using the prefix as a condition
    DE_foldchange_pvalues_signi <- DE_foldchange_pvalues %>%
      filter(!!sym(paste0(condition, "_p_value")) < 0.05)
    
    # Print the filtered data
    cat("Filtered data for condition", condition, ":\n")
    # print(head(DE_foldchange_pvalues_signi))
    cat("\n")

    condition_parts <- strsplit(condition, "_vs_")[[1]]
    first_part <- condition_parts[1]
    second_part <- condition_parts[2]


    mydata_meta <- select(DE_foldchange_pvalues_signi, "InChIkey2D_sirius", "row_id","name_sirius","smiles_sirius",
    "cluster.index_gnps","feature_rt","feature_mz", "adduct_sirius", "chebiasciiname_sirius", "chebiid_sirius", "molecularFormula_sirius", "componentindex_gnps", "GNPSLinkout_Cluster_gnps", "LibraryID_gnps")
    mydata_meta$name_comp <- "unknown"
    mydata_meta$name_comp[!is.na(mydata_meta$InChIkey2D_sirius)] <- mydata_meta$name_sirius[!is.na(mydata_meta$InChIkey2D_sirius)]



    mydata1 <- select(DE_foldchange_pvalues_signi,
    !!sym(paste0(condition, "_fold_change_log2")), "name_sirius", "row_id",
    "NPC.class_canopus")  %>% 
    # this line remove rows with NA in the NPC.class_canopus column using the filter function
    filter(!is.na(NPC.class_canopus))

    mydata1 <- merge(mydata1,npclassifier_newpath,by="NPC.class_canopus")


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

    mydata1 <- mydata1[!is.na(mydata1$NPC.pathway_canopus), ]
    mydata1$counter <- 1

    # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
    mydata1_neg <- mydata1_neg[!is.na(mydata1_neg$NPC.pathway_canopus), ]
    mydata1_neg$counter <- 1
    mydata1_neg$fold_dir <- paste("neg", mydata1_neg$NPC.superclass_canopus, sep = "_")
    # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
    mydata1_pos <- mydata1_pos[!is.na(mydata1_pos$NPC.superclass_canopus), ]
    mydata1_pos$counter <- 1
    mydata1_pos$fold_dir <- paste("pos", mydata1_pos$NPC.superclass_canopus, sep = "_")

    #####################################################################
    #####################################################################


    dt_se_prop_prep_count_tot = dt_for_treemap(
      datatable = mydata1,
      parent_value = NPC.pathway_canopus,
      value = NPC.superclass_canopus,
      count = counter
    )


    dt_se_prop_prep_fold_tot = dt_for_treemap_mean(
      datatable = mydata1,
      parent_value = NPC.pathway_canopus,
      value = NPC.superclass_canopus,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_tot <- dt_se_prop_prep_fold_tot %>% 
    select(-c("value","parent.value"))
    matt_class_fig_tot <- merge(dt_se_prop_prep_count_tot,dt_se_prop_prep_fold_tot,by="ids")

    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_pos = dt_for_treemap(
      datatable = mydata1_pos,
      parent_value = NPC.superclass_canopus,
      value = fold_dir,
      count = counter
    )

    dt_se_prop_prep_fold_pos = dt_for_treemap_mean(
      datatable = mydata1_pos,
      parent_value = NPC.superclass_canopus,
      value = fold_dir,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_pos <- dt_se_prop_prep_fold_pos %>% 
    select(-c("value","parent.value"))
    matt_class_fig_pos_dir <- merge(dt_se_prop_prep_count_pos,dt_se_prop_prep_fold_pos,by="ids")


    matt_class_fig_pos_dir <- matt_class_fig_pos_dir[!(matt_class_fig_pos_dir$parent.value == ""),]
    matt_class_fig_pos_dir <- na.omit(matt_class_fig_pos_dir)

    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_neg = dt_for_treemap(
      datatable = mydata1_neg,
      parent_value = NPC.superclass_canopus,
      value = fold_dir,
      count = counter
    )

    dt_se_prop_prep_fold_neg = dt_for_treemap_mean(
      datatable = mydata1_neg,
      parent_value = NPC.superclass_canopus,
      value = fold_dir,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_neg <- dt_se_prop_prep_fold_neg %>% 
    select(-c("value","parent.value"))
    matt_class_fig_neg_dir <- merge(dt_se_prop_prep_count_neg,dt_se_prop_prep_fold_neg,by="ids")

    matt_class_fig_neg_dir <- matt_class_fig_neg_dir[!(matt_class_fig_neg_dir$parent.value == ""),]
    matt_class_fig_neg_dir <- na.omit(matt_class_fig_neg_dir)

    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_pos_sirius = dt_for_treemap(
      datatable = mydata1_pos,
      parent_value = fold_dir,
      value = row_id,
      count = counter
    )

    dt_se_prop_prep_fold_pos_sirius = dt_for_treemap_mean(
      datatable = mydata1_pos,
      parent_value = fold_dir,
      value = row_id,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_pos_sirius <- dt_se_prop_prep_fold_pos_sirius %>% 
    select(-c("value","parent.value"))
    matt_class_fig_pos_dir_sirius <- merge(dt_se_prop_prep_count_pos_sirius,dt_se_prop_prep_fold_pos_sirius,by="ids")

    matt_class_fig_pos_dir_sirius <- matt_class_fig_pos_dir_sirius[!(matt_class_fig_pos_dir_sirius$parent.value == ""),]
    matt_class_fig_pos_dir_sirius <- na.omit(matt_class_fig_pos_dir_sirius)


    #####################################################################
    #####################################################################

    dt_se_prop_prep_count_neg_sirius = dt_for_treemap(
      datatable = mydata1_neg,
      parent_value = fold_dir,
      value = row_id,
      count = counter
    )

    dt_se_prop_prep_fold_neg_sirius = dt_for_treemap_mean(
      datatable = mydata1_neg,
      parent_value = fold_dir,
      value = row_id,
      count = !!sym(paste0(condition, "_fold_change_log2"))
    )

    dt_se_prop_prep_fold_neg_sirius <- dt_se_prop_prep_fold_neg_sirius %>% 
    select(-c("value","parent.value"))
    matt_class_fig_neg_dir_sirius <- merge(dt_se_prop_prep_count_neg_sirius,dt_se_prop_prep_fold_neg_sirius,by="ids")

    matt_class_fig_neg_dir_sirius <- matt_class_fig_neg_dir_sirius[!(matt_class_fig_neg_dir_sirius$parent.value == ""),]
    matt_class_fig_neg_dir_sirius <- na.omit(matt_class_fig_neg_dir_sirius)



    #####################################################################
    #####################################################################
    #####################################################################
    #####################################################################

    matttree <- rbind(matt_class_fig_tot,matt_class_fig_pos_dir,matt_class_fig_neg_dir,matt_class_fig_pos_dir_sirius,matt_class_fig_neg_dir_sirius)
    matttree$labels_adjusted <- matttree$value
    matttree$labels_adjusted[grep("pos_",matttree$labels_adjusted)] <- ""
    matttree$labels_adjusted[grep("neg_",matttree$labels_adjusted)] <- ""
    matttree$labels_adjusted <- gsub(" x"," ",matttree$labels_adjusted)

  # We rename the count.x column as count and the count.y column as foldchange_log2
    matttree <- matttree %>% 
    rename(count = count.x)  %>% 
    rename(foldchange_log2 = count.y)



    matttree <- merge(matttree,mydata_meta,by.x="labels_adjusted",by.y="row_id",all.x =T)

    matttree$labels_adjusted[!is.na(matttree$name_comp)] <- matttree$name_comp[!is.na(matttree$name_comp)]
    matttree$value[matttree$labels_adjusted== "unknown"] <- ""
    matttree$value[matttree$labels_adjusted== "unknown"] <- ""


    #####################################################################

    # The follow function creates a new hyperlink column based on the labels_adjusted columns

    # matttree$hl <- paste0("https://en.wikipedia.org/wiki/", matttree$labels_adjusted)

    # # <a href='https://example.com/box1' target='_blank'>Box 1</a>
    # matttree$full_hl <- paste0("<a href='", matttree$hl, "' target='_blank'>", matttree$labels_adjusted, "</a>")
    # matttree$full_hl <- paste0(
    #   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
    # )

    # matttree$hl <- paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", matttree$InChIkey2D_sirius, "&sort=annothitcnt")

    # # <a href='https://example.com/box1' target='_blank'>Box 1</a>
    # matttree$full_hl <- paste0(
    #   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
    # )

    # <a href='https://example.com/box1' target='_blank'>Box 1</a>
    matttree$smiles_url <- paste0(
      "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=", matttree$smiles_sirius, "&zoom=2.0&annotate=cip"
    )

      # Generate hl URL only if InChIkey2D_sirius is not NA
    matttree$hl <- ifelse(!is.na(matttree$InChIkey2D_sirius),
                          paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", matttree$InChIkey2D_sirius, "&sort=annothitcnt"),
                          NA)

    # Generate full_hl hyperlink only if hl is not NA
    matttree$full_hl <- paste0(
      "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
    )

      # Generate hl URL only if InChIkey2D_sirius is not NA
    matttree$chebi_hl <- ifelse(!is.na(matttree$chebiid_sirius),
                          paste0("https://www.ebi.ac.uk/chebi/searchId.do?chebiId=", matttree$chebiid_sirius),
                          NA)

    # Generate full_hl hyperlink only if hl is not NA
    matttree$chebi_hl_formatted <- ifelse(!is.na(matttree$chebiid_sirius),
    paste0(
      "<a href='", matttree$chebi_hl, "' target='_blank' style='color: black;'>", matttree$chebiid_sirius, "</a>"
    ), "")

    # Generate full_hl hyperlink only if hl is not NA
    matttree$gnps_hl_formatted <- ifelse(!is.na(matttree$cluster.index_gnps),
    paste0(
      "<a href='", matttree$GNPSLinkout_Cluster_gnps, "' target='_blank' style='color: black;'>", matttree$cluster.index_gnps, "</a>"
    ), "")

    # Generate smiles_url only if smiles_sirius is not NA
    matttree$smiles_url <- ifelse(!is.na(matttree$smiles_sirius),
                                  paste0("https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=", matttree$smiles_sirius, "&zoom=2.0&annotate=cip"),
                                  NA)
    # Generate clickable smiles_url only if smiles_url is not NA
    matttree$smiles_clickable_url <- ifelse(!is.na(matttree$smiles_url),
                                  paste0("<a href='", matttree$smiles_url, "' target='_blank' style='color: black;'>", matttree$smiles_sirius, "</a>"),
                                  NA)


        # "molecularFormula_sirius", "componentindex_gnps", "GNPSLinkout_Cluster_gnps", "LibraryID_GNPS"
    # mattree$smiles_clickable_url <- paste0("<a href=", matttree$smiles_url, " target='_blank' rel='noopener noreferrer'>", matttree$smiles_sirius, "</a>")


    # Here we replace all NA by empty cells in the matttree$smiles_clickable_url column

    matttree$smiles_clickable_url[is.na(matttree$smiles_clickable_url)] <- ""
    matttree$chebiid_sirius[is.na(matttree$chebiid_sirius)] <- ""
    matttree$chebiasciiname_sirius[is.na(matttree$chebiasciiname_sirius)] <- ""
    matttree$LibraryID_gnps[is.na(matttree$LibraryID_gnps)] <- ""


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
    matttree <- matttree[order(matttree$value),]


    #########################################################
    #########################################################

    txt <- as.character(paste0
    ("feature id: ", matttree$cluster.index_gnps,"<br>",
    "component id: ", matttree$componentindex_gnps,"<br>",
    "name: ", matttree$labels_adjusted,"<br>",
    "m/z: ", round(matttree$feature_mz,4),"<br>",
    "RT: ", round(matttree$feature_rt,2),"<br>",
    "MF: ", matttree$molecularFormula_sirius,"<br>",
    "adduct: ", matttree$adduct_sirius,"<br>",
    "FC (log 2): ", round(matttree$foldchange_log2,2),
    "<extra></extra>"
    ))


    matttree$txt <- txt



    fig_treemap_qual = plot_ly(
      data = matttree,
      type = "treemap",
      ids = ~value,
      labels = ~paste0("<b>", matttree$full_hl, "</b><br>", matttree$smiles_clickable_url, "<br><b>", matttree$gnps_hl_formatted, "<br>", matttree$matttree$LibraryID_gnps, "<br>", matttree$chebiasciiname_sirius, "</b><br>", matttree$chebi_hl_formatted, "<br>", "</a>"),
      parents = ~parent.value,
      values = ~count,
      branchvalues = "total",
      maxdepth=3,
      hovertemplate = ~txt,
      marker = list(
        colors = matttree$colors  # Use the colors column from the data frame
      )
      ) %>% 
      layout(
        title = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters: [", filter_sample_metadata_status, "]"),
                   margin = list(
    l = 100,  # Left margin in pixels, adjust as needed
    r = 100,  # Right margin in pixels, adjust as needed
    t = 100,  # Top margin in pixels, adjust as needed
    b = 100   # Bottom margin in pixels, adjust as needed
  )
      )

    fig_treemap_quan <- plot_ly(
      data = matttree,
      type = "treemap",
      ids = ~value,
      labels = ~paste0("<b>", matttree$full_hl, "</b><br>", matttree$smiles_clickable_url, "<br><b>", matttree$gnps_hl_formatted, "<br>", matttree$matttree$LibraryID_gnps, "<br>", matttree$chebiasciiname_sirius, "</b><br>", matttree$chebi_hl_formatted, "<br>", "</a>"),
      parents = ~parent.value,
      values = ~count,
      branchvalues = "total",
      maxdepth = 4,
      hovertemplate = ~txt,
      marker = list(
        colors = matttree$foldchange_log2,
        colorscale = list(
          c(0, 0.5, 1),
          c("#A89639", "#FFFFFF", "#337AB7")),
        cmin = max(abs(matttree$foldchange_log2)) * (-1),
        cmax = max(abs(matttree$foldchange_log2)),
        showscale = TRUE,
        colorbar = list(
        # the title html is set to add a line return
        title = '',
          tickmode = "array",
          tickvals = c((quantile(abs(matttree$foldchange_log2),probs=0.75)*(-1)), 0, (quantile(abs(matttree$foldchange_log2),probs=0.75))),
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
        reversescale = FALSE  # Set to FALSE to maintain the color gradient order
      )
    )%>% 
      layout(
        title = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters: [", filter_sample_metadata_status, "]"),
            margin = list(
    l = 100,  # Left margin in pixels, adjust as needed
    r = 100,  # Right margin in pixels, adjust as needed
    t = 100,  # Top margin in pixels, adjust as needed
    b = 100   # Bottom margin in pixels, adjust as needed
  )
      )

    # We now save the treempa as a html file locally

  if (params$operating_system$system == "unix") {
    ###linux version
    htmlwidgets::saveWidget(fig_treemap_qual, file = paste0("Treemap_", first_part, "_vs_", second_part, "_qual.html"), selfcontained = TRUE) # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_qual.html")
   
    htmlwidgets::saveWidget(fig_treemap_quan, file = paste0("Treemap_", first_part, "_vs_", second_part, "_quan.html"), selfcontained = TRUE) # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_quan.html")
    }


if (params$operating_system$system == "windows") {
    ###windows version
    Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
    htmlwidgets::saveWidget(fig_treemap_qual, file = paste0("Treemap_", first_part, "_vs_", second_part, "_qual.html"), selfcontained = TRUE,libdir = "lib") # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_qual.html")
    unlink("lib", recursive = FALSE)
    
    htmlwidgets::saveWidget(fig_treemap_quan, file = paste0("Treemap_", first_part, "_vs_", second_part, "_quan.html"), selfcontained = TRUE,libdir = "lib") # paste0(file_prefix, "_", first_part, "_vs_", second_part, "_treemap_qual.html")
    unlink("lib", recursive = FALSE)
    }
    }
  }


# mydata_meta <- select(DE_foldchange_pvalues, "InChIkey2D_sirius", "row_id","name_sirius","smiles_sirius",
# "cluster.index_gnps","feature_rt","feature_mz")  
# mydata_meta$name_comp<- "unknown"
# mydata_meta$name_comp[!is.na(mydata_meta$InChIkey2D_sirius)] <- mydata_meta$name_sirius[!is.na(mydata_meta$InChIkey2D_sirius)]




# mydata1 <- select(DE_foldchange_pvalues,
# "C_WT_fold_change_log2","name_sirius", "row_id",
# "NPC.class_canopus")  %>% 
# # this line remove rows with NA in the NPC.class_canopus column using the filter function
# filter(!is.na(NPC.class_canopus))  %>% 
# filter(!is.na(C_WT_fold_change_log2))


# mydata1 <- merge(mydata1,npclassifier_newpath,by="NPC.class_canopus")



# mydata1 <- mydata1 %>% 
# mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x)) %>% 
# mutate_if(is.numeric, function(x) ifelse(is.nan(x), 0, x))


# mydata1_neg <- mydata1[mydata1$C_WT_fold_change_log2 < 0,]
# #mydata1_neg$C_WT_fold_change_log2 <- abs(mydata1_neg$C_WT_fold_change_log2)
# mydata1_pos <- mydata1[mydata1$C_WT_fold_change_log2 >= 0,]

# # Aggregate the data
# ####
# mydata1 = mydata1[!is.na(mydata1$NPC.pathway_canopus), ]
# mydata1$counter = 1

# # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
# mydata1_neg = mydata1_neg[!is.na(mydata1_neg$NPC.pathway_canopus), ]
# mydata1_neg$counter = 1
# mydata1_neg$fold_dir <- paste("neg",mydata1_neg$NPC.superclass_canopus,sep="_")
# # matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
# mydata1_pos = mydata1_pos[!is.na(mydata1_pos$NPC.superclass_canopus), ]
# mydata1_pos$counter = 1
# mydata1_pos$fold_dir <- paste("pos",mydata1_pos$NPC.superclass_canopus,sep="_")

# #####################################################################
# #####################################################################

# dt_se_prop_prep_count_tot = dt_for_treemap(
#   datatable = mydata1,
#   parent_value = NPC.pathway_canopus,
#   value = NPC.superclass_canopus,
#   count = counter
# )


# dt_se_prop_prep_fold_tot = dt_for_treemap_mean(
#   datatable = mydata1,
#   parent_value = NPC.pathway_canopus,
#   value = NPC.superclass_canopus,
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
#   parent_value = NPC.superclass_canopus,
#   value = fold_dir,
#   count = counter
# )

# dt_se_prop_prep_fold_pos = dt_for_treemap_mean(
#   datatable = mydata1_pos,
#   parent_value = NPC.superclass_canopus,
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
#   parent_value = NPC.superclass_canopus,
#   value = fold_dir,
#   count = counter
# )

# dt_se_prop_prep_fold_neg = dt_for_treemap_mean(
#   datatable = mydata1_neg,
#   parent_value = NPC.superclass_canopus,
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

# matttree$hl <- paste0("https://pubchem.ncbi.nlm.nih.gov/#query=", matttree$InChIkey2D_sirius, "&sort=annothitcnt")

# # <a href='https://example.com/box1' target='_blank'>Box 1</a>
# matttree$full_hl <- paste0(
#   "<a href='", matttree$hl, "' target='_blank' style='color: black;'>", matttree$labels_adjusted, "</a>"
# )

# # <a href='https://example.com/box1' target='_blank'>Box 1</a>
# matttree$smiles_url <- paste0(
#   "https://www.simolecule.com/cdkdepict/depict/bow/svg?smi=", matttree$smiles_sirius, "&zoom=2.0&annotate=cip"
# )


# matttree <- matttree[order(matttree$value),]
# ####################### annoation structure 
# #########################################################

# ############# ad mol
# #mol_list_2d <- list()
# #for (i in c(1:length(matttree$smiles_sirius))) {
# #  psml <- parse.smiles(matttree$smiles_sirius[i], omit.nulls = TRUE)
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
# ("feature id: ",matttree$cluster.index_gnps,"<br>",
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

matt_donust = DE_foldchange_pvalues %>%
  filter(if_any(ends_with('_p_value'), ~ .x < params$posthoc$p_value))

# matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
matt_donust2 = matt_donust[!is.na(matt_donust$NPC.superclass_canopus), ]
matt_donust2$counter = 1

 

dt_for_treemap = function(datatable, parent_value, value, count) {
  parent_value = enquo(parent_value)
  value = enquo(value)
  count = enquo(count)

  datatable = data.frame(datatable %>%
    group_by(!!parent_value, !!value, ) %>%
    summarise(count = sum(as.numeric(!!count))))

  datatable = datatable %>%
    select(!!parent_value, !!value, count) %>% # create id labels for each row # Notre the !! to pass aruguments to a dplyr function
    rename(
      parent.value = !!parent_value,
      value = !!value
    ) %>%
    mutate(ids = ifelse(parent.value == "", value,
      paste0(value, "-", parent.value) # Notre that here we are passing argument to a non dplyr function call
    )) %>%
    select(ids, everything())

  par_info = datatable %>% dplyr::group_by(parent.value) %>% # group by parent
    dplyr::summarise(count = sum(as.numeric(count))) %>% # parent total
    rename(value = parent.value) %>% # parent labels for the item field
    mutate(parent.value = "", ids = value) %>% # add missing fields for my_data
    select(names(datatable)) # put cols in same order as my_data

  data_for_plot = rbind(datatable, par_info)

  return(data_for_plot)
}

dt_se_prop_prep_tm = dt_for_treemap(
  datatable = matt_donust2,
  parent_value = NPC.superclass_canopus,
  value = NPC.class_canopus,
  count = counter
)


fig_treemap = plot_ly(
  data = dt_se_prop_prep_tm,
  type = "treemap",
  labels = ~value,
  parents = ~parent.value,
  values = ~count,
  branchvalues = "total"
)

# Why "significant ? According to what ?

fig_treemap = fig_treemap %>% 
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
    htmlwidgets::saveWidget(file = filename_treemap, selfcontained = TRUE,libdir = "lib")
unlink("lib", recursive = FALSE)

}


#############################################################################
#############################################################################
############## Random Forest ################################################
#############################################################################
#############################################################################

message("Launching Random Forest calculations ...")

sink(filename_random_forest_model)

# Here we traduce to fit Manu's inputs ... to be updated later 

features_of_importance = DE_foldchange_pvalues %>%
  filter((!!as.symbol(p_value_column)) < params$posthoc$p_value)  %>% 
  select(feature_id) %>%
  # we output the data as a vector
  pull()


# We select all columns except the params$target$sample_metadata_header columns in 
# data_subset_norm_rf_filter and we prefix the column names with an X.
# We use the dplyr syntax to do this and the rename function to rename the columns
# We then subset the data to keep only the columns that are in the imp_filter1 variable

data_subset_for_RF = DE$data %>%
  select(all_of(as.character(features_of_importance))) %>% 
  rename_all(~ paste0("X", .))  %>% 
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names")  %>% 
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names)  %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X"))  %>% 
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header)))

# We define the formula externally to inject the external variable # params$target$sample_metadata_header

formula = as.formula(paste0(params$target$sample_metadata_header," ~ ."))

# We launch the rfPermute function

data.rp = rfPermute(formula, data = data_subset_for_RF, na.action = na.omit, ntree = 500, num.rep = 500)

imp_table_rf = data.frame(data.rp$pval)
imp_table_rf = importance(data.rp)
imp_table_rf = data.frame(imp_table_rf)

summary(data.rp)

#f = plotImportance(data.rp, plot.type = "bar", plot = FALSE)

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
  y = ~reorder(row.names(imp_table_rf), -MeanDecreaseGini,decreasing = TRUE),  # Use reorder to maintain sorting order
  type = "bar",
  orientation = "h"
) %>%
  layout(
    title = title_random_forest,
    xaxis = list(title = "Importance", tickfont = list(size = 12)),   # Adjust the label size here (e.g., size = 12)
    yaxis = list(title = "Features", tickfont = list(size = 5)),      # Adjust the label size here (e.g., size = 10)
    margin = list(l = 100, r = 20, t = 50, b = 70),
    showlegend = FALSE
  )
fig_rf
# The file is exported
# The title should be updated !!! 


if (params$operating_system$system == "unix") {
### linux version
fig_rf %>%
    htmlwidgets::saveWidget(file = filename_random_forest , selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
### windows version
Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
fig_rf %>%
    htmlwidgets::saveWidget(file = filename_random_forest, selfcontained = TRUE,libdir = "lib")
unlink("lib", recursive = FALSE)

}

#############################################################################
#############################################################################
############## Random Forest selected Box Plots #############################
#############################################################################
#############################################################################

message("Preparing RF-selected Box plots ...")


imp.scaled = rfPermute::importance(data.rp, scale = TRUE)
imp.scaled = data.frame(imp.scaled)
imp.scaled = imp.scaled[order(imp.scaled$MeanDecreaseGini, decreasing = TRUE), ]
# boxplot_top_N = row.names(imp.scaled)[1:params$boxplot$topN]
# This below using head is safer in case the number of features is smaller than the number of features to plot
boxplot_top_N = row.names(head(imp.scaled, n = params$boxplot$topN))

data_subset_boxplot = data_subset_for_RF %>%
  select(all_of(boxplot_top_N), params$target$sample_metadata_header)

# We now establish a side by side box plot for each columns of the data_subset_norm_boxplot
# We use the melt function to reshape the data to a long format
# We then use the ggplot2 syntax to plot the data and the facet_wrap function to plot the data side by side

# Gather value columns into key-value pairs
df_long <- tidyr::gather(data_subset_boxplot, key = "variable", value = "value", -params$target$sample_metadata_header)


# Create boxplots faceted by variable and colored by age
# Note how we use the get function to access the variable name
# See here for details  https://stackoverflow.com/a/22309328/4908629
p = ggplot(df_long, aes(x = get(params$target$sample_metadata_header), y = value, fill = get(params$target$sample_metadata_header))) +
  geom_boxplot() +
  facet_wrap(~ variable, ncol = 4) +
  theme_minimal()+
  ggtitle(title_box_plots) 


fig_boxplot = p + facet_wrap(~variable, scales = "free", dir = "v") + theme(
  legend.position = "top",
  legend.title = element_blank()
)

# fig_boxplotly = data_subset_boxplot %>%
#   split(data_subset_boxplot$variable) %>%
#   map(~ {
#     plot_ly(data = .x, x = .x$treatment, y = .x$value, color = .x$treatment, type = "box") %>%
#       layout(showlegend = F, xaxis = list(title = .x$variable[1])) %>%
#       layout(xaxis = list(titlefont = list(size = 10), tickfont = list(size = 10))) %>%
#       layout(yaxis = list(titlefont = list(size = 12), tickfont = list(size = 12)))
#   }) %>%
#   subplot(margin = 0.02, nrows = 4, titleX = TRUE)  %>% 
#   layout(title = title_box_plots,
#   legend = list(title = list(text = paste("<b>class </b>")))) # Find a way to define the legend title

# The files are exported

ggsave(plot = fig_boxplot, filename = filename_box_plots, width = 10, height = 10)
# fig_boxplotly %>%
#     htmlwidgets::saveWidget(file = filename_box_plots_interactive, selfcontained = TRUE)



#############################################################################
#############################################################################
############## p-Value selected Box Plots #############################
#############################################################################
#############################################################################

message("Preparing p-value selected Box plots ...")

features_of_importance_boxplots = DE_foldchange_pvalues %>%
  # we keep only the features that have a p-value lower than the threshold
  # filter((!!as.symbol(p_value_column)) < params$posthoc$p_value)  %>% 
  # We keep only the lowest top n = params$boxplot$topN in the p_value_column
  top_n(-params$boxplot$topN, !!as.symbol(p_value_column))  %>%
  # we order the features by increasing p-value
  arrange(!!as.symbol(p_value_column))  %>%
  select(feature_id) %>%
  # we output the data as a vector
  pull()


data_subset_for_boxplots = DE$data %>%
  select(all_of(as.character(features_of_importance_boxplots))) %>% 
  rename_all(~ paste0("X", .))  %>% 
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names")  %>% 
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names)  %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X"))  %>% 
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header)))  %>% 
  # Finally we remove the X from the columns names
  rename_all(~ gsub("X", "", .))

# We now establish a side by side box plot for each columns of the data_subset_norm_boxplot
# We use the melt function to reshape the data to a long format
# We then use the ggplot2 syntax to plot the data and the facet_wrap function to plot the data side by side

# Gather value columns into key-value pairs
df_long <- tidyr::gather(data_subset_for_boxplots, key = "variable", value = "value", -params$target$sample_metadata_header)

# Here we merge the df_long with the DE$variable_meta data frame to get the variable type

df_long_informed <- merge(df_long, DE_foldchange_pvalues, by.x = "variable", by.y = "feature_id")


p = ggplot(df_long_informed, aes(x = !!sym(params$target$sample_metadata_header), y = value, fill = !!sym(params$target$sample_metadata_header))) +
  geom_boxplot() +
  facet_wrap(~feature_id_full_annotated , ncol = 4) +
  # theme_minimal()+
  ggtitle(title_box_plots) 

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

fig_boxplot = p + facet_wrap2(~ chebiasciiname_sirius + feature_id_full, labeller = label_value, strip = ridiculous_strips) + theme(
  legend.position = "top",
  legend.title = element_blank()
)

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
    geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +  # Add data points with jitter for better visibility
    # ggtitle(paste("Box Plot for", "\n", 
    # "Compound name: ", data_for_plot$chebiasciiname_sirius[1], "\n",
    # "Feature details: ", data_for_plot$feature_id_full[1]))
labs(x=params$target$sample_metadata_header,
       y="Normalized Intensity",
       title = paste("Compared intensities for feature:", var),
       subtitle = paste("\n",
       "Compound name: ", data_for_plot$chebiasciiname_sirius[1], "\n",
       "Feature details: ", data_for_plot$feature_id_full[1]),
       caption  = paste("Calculated p-value is ~ ", rounded_p_value)) +
  theme(plot.caption = element_text(hjust = 0, face= "italic"), #Default is hjust=1
        plot.title.position = "plot", #NEW parameter. Apply for subtitle too.
        plot.caption.position =  "plot") #NEW parameter

  # Save the plot to a file with a unique filename for each variable
  filename <- paste(output_directory_bp, "boxplot_", gsub(" ", "_", var), ".png", sep = "")
  ggsave(plot = p, filename = filename, width = 8, height = 8)
}

#############################################################################
############## Pvalue filtered Heat Map  #############################
#############################################################################
#############################################################################

message("Preparing p-value filtered Heatmap ...")

features_of_importance = DE_foldchange_pvalues %>%
  filter((!!as.symbol(p_value_column)) < params$posthoc$p_value)  %>% 
  select(feature_id) %>%
  # we output the data as a vector
  pull()


# data_subset_for_Pval = DE$data %>%
#   select(all_of(as.character(features_of_importance))) %>% 
#   rename_all(~ paste0("X", .))  %>% 
#   # here we join the data with the associated sample metadata using the row.names as index
#   merge(DE$sample_meta, ., by = "row.names")  %>% 
#   # We keep the row.names columnn as row.names
#   transform(row.names = Row.names)  %>%
#   # We keep the params$target$sample_metadata_header column and the columns that start with X
#   select(params$target$sample_metadata_header, starts_with("X"))  %>% 
#   # We set the params$target$sample_metadata_header column as a factor
#   mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header)))


data_subset_for_pval_hm = DE$data %>%
  select(all_of(as.character(features_of_importance))) %>% 
  rename_all(~ paste0("X", .))  %>% 
  # here we join the data with the associated sample metadata using the row.names as index
  merge(DE$sample_meta, ., by = "row.names")  %>% 
  # We keep the row.names columnn as row.names
  transform(row.names = Row.names)  %>%
  # We keep the params$target$sample_metadata_header column and the columns that start with X
  select(params$target$sample_metadata_header, starts_with("X"))  %>% 
  # We set the params$target$sample_metadata_header column as a factor
  mutate(!!as.symbol(params$target$sample_metadata_header) := factor(!!as.symbol(params$target$sample_metadata_header)))  %>% 
  # Finally we remove the X from the columns names
  rename_all(~ gsub("X", "", .))

data_subset_for_pval_hm_sel = data_subset_for_pval_hm %>%
  select(params$target$sample_metadata_header)

# imp_filter2X = row.names(imp_table_rf)
# imp_filter2 = gsub("X", "", imp_filter2X)


# features_of_importance_heatmap = DE_foldchange_pvalues %>%
#   # we keep only the features that have a p-value lower than the threshold
#   # filter((!!as.symbol(p_value_column)) < params$posthoc$p_value)  %>% 
#   # We keep only the lowest top n = params$boxplot$topN in the p_value_column
#   top_n(-params$heatmap$topN, !!as.symbol(p_value_column))  %>%
#   # we order the features by increasing p-value
#   arrange(!!as.symbol(p_value_column))  %>%
#   select(feature_id) %>%
#   # we output the data as a vector
#   pull()



data_subset_for_pval_hm = data_subset_for_pval_hm[, colnames(data_subset_for_pval_hm) %in% features_of_importance]
# my_sample_col = DE$sample_meta$sample_id

# data_subset_for_Pval = data_subset_for_Pval[, colnames(data_subset_for_Pval) %in% imp_filter2X]
# # my_sample_col = DE$sample_meta$sample_id

my_sample_col = paste(DE$sample_meta$sample_id, DE$sample_meta[[params$target$sample_metadata_header]], sep = "_")

# annot_col = data.frame(paste(DE$variable_meta$NPC.pathway_canopus, DE$variable_meta$NPC.superclass_canopus, sep = "_"), DE$variable_meta$NPC.pathway_canopus)

# colnames(annot_col) = c("Superclass", "Pathway")

# rownames(annot_col) = DE$variable_meta$feature_id
# annot_col_filter = annot_col[rownames(annot_col) %in% imp_filter2, ]

# We filter the annotation table (DE$variable_meta) to keep only the features of interest identified in the (features_of_importance). We use dplyr

selected_variable_meta = DE$variable_meta %>%
  filter(feature_id %in% features_of_importance) 
  # %>%
  # select(feature_id, NPC.pathway_canopus, NPC.superclass_canopus) %>%
  # mutate(NPC.pathway_canopus = paste(NPC.pathway_canopus, NPC.superclass_canopus, sep = "_")) %>%
  # select(feature_id, NPC.pathway_canopus) %>%
  # column_to_rownames("feature_id")

selected_variable_meta_NPC = DE$variable_meta %>%
  filter(feature_id %in% features_of_importance)  %>% 
  select(feature_id, NPC.superclass_canopus, NPC.pathway_canopus, NPC.class_canopus) %>%
  mutate(NPC.superclass_merged_canopus = paste(NPC.pathway_canopus, NPC.superclass_canopus, sep = "_")) %>%
  mutate(NPC.class_merged_canopus = paste(NPC.superclass_merged_canopus, NPC.class_canopus, sep = "_")) %>%
  select(NPC.class_merged_canopus, NPC.superclass_merged_canopus, NPC.pathway_canopus)

ByPal = colorRampPalette(c(wes_palette("Zissou1")))

# data_subset_for_Pval = apply(data_subset_for_Pval, 2, as.numeric)
# # heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")


data_subset_for_pval_hm = apply(data_subset_for_pval_hm, 2, as.numeric)
# heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")


heatmap_filtered_pval = heatmaply(
  percentize(data_subset_for_pval_hm),
  seriate = "none", # none , GW , mean, OLO
  col_side_colors = data.frame(selected_variable_meta_NPC, check.names = FALSE),
  col_side_palette = ByPal,
  row_side_colors = data_subset_for_pval_hm_sel,
  # row_side_palette = ByPal,
  labRow = as.vector(as.character(my_sample_col)), # [vec_plot]
  labCol = selected_variable_meta$feature_id_full_annotated,
  subplot_margin = 0.01,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "lightsteelblue2",
    # mid = "goldenrod1",
    high = "firebrick3",
    midpoint = 0.5,
    limits = c(0, 1),
    position = "left"
  ),
  hide_colorbar = TRUE,
  branches_lwd = 0.3,
  k_row = 4,
  distfun_row = "pearson",
  distfun_col = "pearson",
  fontsize_row = 9,
  fontsize_col = 9,
  Rowv = FALSE,
  side_color_colorbar_len = 1,
  # # # # Colv = NULL,
  plot_method = "plotly"
)  %>% layout(
  title = list(text = title_heatmap_pval, font = list(size = 14), x = 0.1),
  margin = list(t = 150, b = 20) # Adjust the top margin value (e.g., 80) to move the title to the top
)

heatmap_filtered_pval

https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html


# x  <- as.matrix(datasets::mtcars)
# rc <- colorspace::rainbow_hcl(nrow(x))


# heatmaply(
#   x[, -c(8, 9)],
#   seriate = "mean",
#   col_side_colors = c(rep(0, 5), rep(1, 4)),
#   row_side_colors = x[, 8:9]
# )

# c(rep(0, 8), rep(1, 8))

# as.character(data_subset_for_pval_hm_sel[[1]])

# The file is exported


if (params$operating_system$system == "unix") {
### linux version

heatmap_filtered_pval %>%
    htmlwidgets::saveWidget(file = filename_heatmap_pval, selfcontained = TRUE) 
}

if (params$operating_system$system == "windows") {
### windows version
Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
heatmap_filtered_pval %>%
    htmlwidgets::saveWidget(file = filename_heatmap_pval, selfcontained = TRUE,libdir = "lib")
unlink("lib", recursive = FALSE)

}

library(iheatmapr)
data(measles, package = "iheatmapr")

main_heatmap(measles, name = "Measles<br>Cases", x_categorical = FALSE,
             layout = list(font = list(size = 8))) %>%
  add_col_groups(ifelse(1930:2001 < 1961,"No","Yes"),
                  side = "bottom", name = "Vaccine<br>Introduced?",
                  title = "Vaccine?",
                  colors = c("lightgray","blue")) %>%
  add_col_labels(ticktext = seq(1930,2000,10),font = list(size = 8)) %>%
  add_row_labels(size = 0.3,font = list(size = 6)) %>% 
  add_col_summary(layout = list(title = "Average<br>across<br>states"),
                  yname = "summary")  %>%                 
  add_col_title("Measles Cases from 1930 to 2001", side= "top") %>%
  add_row_summary(groups = TRUE, 
                  type = "bar",
                  layout = list(title = "Average<br>per<br>year",
                                font = list(size = 8)))
                                

#############################################################################
#############################################################################
############## Random Forest filtered Heat Map  #############################
#############################################################################
#############################################################################

message("Preparing Random Forest filtered Heatmap ...")

imp_table_rf_order = imp_table_rf[order(imp_table_rf$MeanDecreaseGini, decreasing = TRUE), ] # imp_table_rf[order(imp_table_rf$MeanDecreaseGini,decreasing=TRUE),]

imp_filter2X = row.names(imp_table_rf_order)[1:params$heatmap$topN]
imp_filter2 = gsub("X", "", imp_filter2X)

data_subset_for_RF = data_subset_for_RF[, colnames(data_subset_for_RF) %in% imp_filter2X]
# my_sample_col = DE$sample_meta$sample_id

my_sample_col = paste(DE$sample_meta$sample_id, DE$sample_meta[[params$target$sample_metadata_header]], sep = "_")

annot_col = data.frame(paste(DE$variable_meta$NPC.pathway_canopus, DE$variable_meta$NPC.superclass_canopus, sep = "_"), DE$variable_meta$NPC.pathway_canopus)

colnames(annot_col) = c("Superclass", "Pathway")

rownames(annot_col) = DE$variable_meta$feature_id
annot_col_filter = annot_col[rownames(annot_col) %in% imp_filter2, ]


ByPal = colorRampPalette(c(wes_palette("Zissou1")))

data_subset_for_RF = apply(data_subset_for_RF, 2, as.numeric)
# heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")


heatmap_filtered_rf = heatmaply(
  percentize(data_subset_for_RF),
  seriate = "mean", # none , GW , mean, OLO
  col_side_colors = data.frame(annot_col_filter, check.names = FALSE),
  col_side_palette = ByPal,
  labRow = as.vector(as.character(my_sample_col)), # [vec_plot]
  subplot_margin = 0.01,
  scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
    low = "lightsteelblue2",
    # mid = "goldenrod1",
    high = "firebrick3",
    midpoint = 0.5,
    limits = c(0, 1)
  ),
  fontsize_col = 5,
  branches_lwd = 0.3,
  k_row = 4,
  distfun_row = "pearson",
  distfun_col = "pearson"
)

heatmap_filtered_rf = heatmap_filtered_rf %>% layout(title = list(text = title_heatmap_rf, y = 0.05))

# The file is exported


if (params$operating_system$system == "unix") {
### linux version

heatmap_filtered_rf %>%
    htmlwidgets::saveWidget(file = filename_heatmap_rf, selfcontained = TRUE)
}

if (params$operating_system$system == "windows") {
### windows version
Sys.setenv(RSTUDIO_PANDOC = params$operating_system$pandoc)
heatmap_filtered_rf %>%
    htmlwidgets::saveWidget(file = filename_heatmap_rf, selfcontained = TRUE,libdir = "lib")
unlink("lib", recursive = FALSE)

}

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
# name_sirius = DE$variable_meta$name_sirius
# smiles_sirius = DE$variable_meta$smiles_sirius
# InChI_sirius = DE$variable_meta$InChI_sirius
# NPC.pathway_canopus = DE$variable_meta$NPC.pathway_canopus
# NPC.superclass_canopus = DE$variable_meta$NPC.superclass_canopus
# Annotation_merge = data.frame(feature_id, sample_raw_id, name_sirius, smiles_sirius, InChI_sirius, NPC.pathway_canopus, NPC.superclass_canopus)

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

summary_stat_output_full = DE_foldchange_pvalues

# We filter the DE_foldchange_pvalues table to only keep the top N features (any column ending with _p_value string should have a value < 0.05)
# We use the dplyr synthax to filter the table
# We need to make sure to remove the rwonames() before exporting


summary_stat_output_selected = DE_foldchange_pvalues %>% 
  filter(if_any(ends_with("_p_value"), ~ . < params$posthoc$p_value))  %>% 
  arrange(across(ends_with("_p_value"))) %>%
  select(
  feature_id_full,
  feature_id,
  feature_mz,
  feature_rt,
  contains("p_value"), 
  contains("fold"), 
  NPC.pathway_canopus,
  NPC.superclass_canopus,
  NPC.class_canopus,
  name_sirius,
  LibraryID_gnps, 
  contains("smiles", ignore.case = TRUE), 
  contains("inchi_", ignore.case = TRUE),
  contains("inchikey", ignore.case = TRUE)
  )


# glimpse(summary_stat_output_selected)


# We also prepare Metaboverse outputs from the fc and pvalues tables


metaboverse_table = DE_foldchange_pvalues

# We then keep the keep the first occurence of the chebiasciiname_sirius

metaboverse_table = metaboverse_table %>%
  distinct(chebiasciiname_sirius, .keep_all = TRUE)

# We now format the table for Metaboverse
# For this we apply the foillowing steps:
# 1. We select the columns we want to keep (chebiasciiname_sirius, Co_KO_p_value, Co_KO_fold_change_log2)
# 2. We rename the columns to the names Metaboverse expects. Using the rename_with and gsub we replace the the _p_value and _fold_change_log2 suffixes to _stat and _fc
# 3. We reorganize the columns to the order Metaboverse expects (chebiasciiname_sirius, _stat, _fc) 
# 4. We remove any rows containing NA values in the dataframe
# 5. We replace the name of the `chebiasciiname_sirius` column by an empty string


metaboverse_table = metaboverse_table %>%
  select(chebiasciiname_sirius, ends_with('_fold_change_log2'), ends_with('_p_value')) %>%
  # rename_with(~gsub("_p_value", "_stat", .)) %>%
  # rename_with(~gsub("_fold_change_log2", "_fc", .)) %>%
  # select(chebiasciiname_sirius, Co_KO_fc, Co_KO_stat) %>%
  # We remove row containing the `Inf` value  
  # filter(!grepl('Inf', ends_with('_fold_change_log2')))  %>% 
  # We remove any row containing the `Inf` value across all columns of the dataframe
  filter(if_any(everything(), ~!str_detect(., "Inf")))  %>% 
  # filter(!grepl('Inf', ends_with('_fold_change_log2')))  %>% 
  na.omit()

colnames(metaboverse_table)[1] <-""

# We now sort columns alphabetically

metaboverse_table = metaboverse_table[,order(colnames(metaboverse_table))]




# The file is exported

write.table(summary_stat_output_full, file = filename_summary_stats_table_full, sep = ",", row.names = FALSE)
write.table(summary_stat_output_selected, file = filename_summary_stats_table_selected, sep = ",", row.names = FALSE)
write.table(metaboverse_table, file = filename_metaboverse_table, sep = "\t", row.names = FALSE, quote = FALSE)


#############################################################################
#############################################################################
######################################## summmary table with structure


summary_stat_output_selected_simple = DE_foldchange_pvalues %>% 
  filter(if_any(ends_with("_p_value"), ~ . < params$posthoc$p_value))  %>% 
  arrange(across(ends_with("_p_value"))) %>%
  select(
  feature_id,
  feature_id_full,
  contains("p_value"), 
  NPC.pathway_canopus,
  NPC.superclass_canopus,
  NPC.class_canopus,
  name_sirius,
  smiles_sirius
  )


# mol_list_2d <- list()

# for (i in c(1:nrow(summary_stat_output_selected_simple))) {
#   psml <- parse.smiles(summary_stat_output_selected_simple$smiles_sirius[i], omit.nulls = TRUE)
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

graphml_file = file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "gnps_molecular_network_graphml", params$filenames$gnps_graphml)


g = read.graph(file = graphml_file, format = "graphml")
# net_gnps = igraph::simplify(g, remove.multiple = FALSE, edge.attr.comb = "ignore")

df_from_graph_edges_original = igraph::as_data_frame(g, what = c("edges"))
df_from_graph_vertices_original = igraph::as_data_frame(g, what = c("vertices"))


# We define drop the from and to columns from the edges dataframe
# And then rename the node 1 and node 2 columns to from and to, respectively
# These columns are placed at the beginning of the dataframe
# and converted to numerics

df_from_graph_edges = df_from_graph_edges_original  %>%
  select(-from, -to) %>%
  rename(from = node1, to = node2) %>%
  select(from, to, everything()) %>%
  mutate_at(vars(from, to), as.numeric)

# the id column of the vertices dataframe is converted to numerics

df_from_graph_vertices = df_from_graph_vertices_original %>%
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

summary_stat_output_red = summary_stat_output_full %>%
  select(-contains("_sirius")) %>%
  select(-contains("_canopus")) %>%
  select(-contains("_metannot")) %>%
  select(-contains("_gnps"))  %>% 
  #select(-ends_with("_id"))  %>% 
  select(-ends_with("_mz"))  %>%
  select(-ends_with("_rt"))


df_from_graph_vertices_plus = merge(DE_original$variable_meta, summary_stat_output_red, by.x = "feature_id", by.y = "feature_id", all.x = T)


# We merge the data from the DE$data dataframe with the DE$sample_meta dataframe using rownames as the key

merged_D_SM = merge(DE_original$sample_meta, DE_original$data, by = 'row.names', all = TRUE)

# We replace NA values with 0 in the merged dataframe
# Check why we dont do this before ??
merged_D_SM[is.na(merged_D_SM)] <- 0


# The function below allows to group data by multiple factors and return a dataframe with the mean of each group


dfList <- list()

for (i in params$colnames$to_output) {
  dfList[[i]] <- merged_D_SM %>%
    group_by(!!as.symbol(i)) %>%
    summarise(across(colnames(DE_original$data), mean),
      .groups = "drop"
    ) %>%
    select(!!all_of(i), colnames(DE_original$data)) %>%
    pivot_longer(-!!i) %>%
    pivot_wider(names_from = all_of(i), values_from = value)  %>% 
    # We prefix all columns with the factor name
    rename_with(.cols=-name, ~paste0("mean_int", "_", i, "_", .x))
}


flat_dfList = reduce(dfList, full_join, by = "name")

# We now add the raw feature list to the dataframe

full_flat_dfList = merge(flat_dfList, t(DE_original$data), by.x = 'name', by.y = 'row.names', all = TRUE)


# We add the raw feature list

df_from_graph_vertices_plus_plus = merge(df_from_graph_vertices_plus, full_flat_dfList, by.x = "feature_id", by.y = "name", all.x = T)


# We set back the id column as the first column of the dataframe

df_from_graph_vertices_plus_plus = df_from_graph_vertices_plus_plus %>%
  select(feature_id, everything())


# We then add the attributes to the edges dataframe and generate the igraph object

# In the case when we have been filtering the X data we will add the filtered X data to the vertices dataframe prior to merging. 




generated_g = graph_from_data_frame(df_from_graph_edges, directed = FALSE, vertices = df_from_graph_vertices_plus_plus)


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

file.copy(script_name, file.path(output_directory,filename_R_script), overwrite = TRUE)


message("Done !")
