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
usePackage("dplyr")
usePackage("emmeans")
usePackage("EnhancedVolcano")
usePackage("fpc")
usePackage("ggdendro")
usePackage("ggplot2")
usePackage("ggraph")
usePackage("ggrepel")
usePackage("ggtree")
usePackage("graphlayouts")
usePackage("gridExtra")
usePackage("here")
usePackage("heatmaply")
usePackage("igraph")
usePackage("manhattanly")
usePackage("openxlsx")
usePackage("plotly")
usePackage("pmp")
usePackage("pls")
usePackage("randomcoloR")
usePackage("randomForest")
usePackage("readr")
usePackage("rfPermute")
usePackage("rgl")
usePackage("ropls")
# usePackage("structToolbox")
usePackage("this.path")
usePackage("tidyverse")
usePackage("vegan")
usePackage("viridis")
usePackage("wesanderson")
usePackage("yaml")

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

# We call the external params

path_to_params = "./params/params.yaml"

params = yaml.load_file(path_to_params)



# We set the working directory

working_directory = file.path(params$path$docs, params$mapp_project, params$mapp_batch, params$polarity)

# We set the output directory

if (params$actions$scale_data == "TRUE") {
scaling_status = "scaled"
} else { scaling_status = "raw" }

if (params$actions$filter_sample_metadata_one == "TRUE" & params$actions$filter_sample_metadata_two == "TRUE") {
filter_sample_metadata_status = paste(params$filter_sample_metadata_one$mode,
params$filter_sample_metadata_one$factor_name,
params$filter_sample_metadata_one$levels,
params$filter_sample_metadata_two$mode,
params$filter_sample_metadata_two$factor_name,
params$filter_sample_metadata_two$levels,
sep = "_")
} else if (params$actions$filter_sample_metadata_one == "TRUE") {
filter_sample_metadata_status = paste(params$filter_sample_metadata_one$mode,
params$filter_sample_metadata_one$factor_name,
params$filter_sample_metadata_one$levels,
sep = "_") 
} else { filter_sample_metadata_status = "no_sm_filter" }


if (params$actions$filter_variable_metadata_one == "TRUE" & params$actions$filter_variable_metadata_two == "TRUE") {
filter_variable_metadata_status = paste(params$filter_variable_metadata_one$mode,
params$filter_variable_metadata_one$factor_name,
params$filter_variable_metadata_one$levels,
params$filter_variable_metadata_two$mode,
params$filter_variable_metadata_two$factor_name,
params$filter_variable_metadata_two$levels,
sep = "_")
} else if (params$actions$filter_variable_metadata_one == "TRUE") {
filter_variable_metadata_status = paste(params$filter_variable_metadata_one$mode,
params$filter_variable_metadata_one$factor_name,
params$filter_variable_metadata_one$levels,
sep = "_") 
} else { filter_variable_metadata_status = "no_vm_filter" }


output_directory = file.path(working_directory, "results", "stats", paste(params$mapp_batch, params$filters$metadata_variable, filter_variable_metadata_status, filter_sample_metadata_status, scaling_status, sep = '_'), sep = "")


dir.create(output_directory)


#################################################################################################
#################################################################################################
################### Filename and paths establishment ##########################################
#################################################################################################


# The Figures titles are conditionally defined according to the user's choices and option in the parameters file


title_PCA = paste("PCA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$filters$metadata_variable, sep = " ") 
title_PCA3D = paste("PCA3D", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$filters$metadata_variable, sep = " ")
title_PCoA = paste("PCoA", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$filters$metadata_variable, sep = " ") 
title_PCoA3D = paste("PCoA3D", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.","Colored according to", params$filters$metadata_variable, sep = " ")
title_volcano = paste("Volcano plot", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
title_treemap = paste("Treemap", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
title_random_forest = paste("Random Forest results", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
title_box_plots = paste("Top", params$boxplot$topN, "boxplots", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")
title_heatmap = paste("Heatmap of","top", params$heatmap$topN,"Random Forest filtered features", "for dataset", filter_variable_metadata_status, "and", filter_sample_metadata_status, "level.", sep = " ")


# The Figures filename is conditionally defined according to the user's choice of filtering the dataset according to CANOPUS NPClassifier classifications or not.


file_prefix = paste(params$mapp_batch, 
                    params$filters$metadata_variable, 
                    filter_variable_metadata_status, 
                    filter_sample_metadata_status, 
                    params$polarity, 
                    scaling_status, 
                    sep = "_")


filename_PCA <- paste(file_prefix, "_PCA.html", sep = "")
filename_PCA3D <- paste(file_prefix, "_PCA3D.html", sep = "")
filename_PCoA <- paste(file_prefix, "_PCoA.pdf", sep = "")
filename_PCoA3D <- paste(file_prefix, "_PCoA3D.html", sep = "")
filename_volcano <- paste(file_prefix, "_Volcano.pdf", sep = "")
filename_volcano_interactive <- paste(file_prefix, "_Volcano_interactive.html", sep = "")
filename_treemap <- paste(file_prefix, "_Treemap_interactive.html", sep = "")
filename_random_forest <- paste(file_prefix, "_RF_importance.html", sep = "")
filename_random_forest_model <- paste(file_prefix, "_RF_model.txt", sep = "")
filename_box_plots <- paste(file_prefix, "_Boxplots.pdf", sep = "")
filename_box_plots_interactive <- paste(file_prefix, "_Boxplots_interactive.html", sep = "")
filename_heatmap <- paste(file_prefix, "_Heatmap.html", sep = "")
filename_summary_stats_table <- paste(file_prefix, "_summary_stats_table.csv", sep = "")
filename_graphml <- paste(file_prefix, "_graphml.graphml", sep = "")
filename_params <- paste(file_prefix, "_params.yaml", sep = "")
filename_session_info <- paste(file_prefix, "_session_info.txt", sep = "")
filename_R_script <- paste(file_prefix, "_R_script_backup.R", sep = "")
filename_DE_model <- paste(file_prefix, "_DE_description.txt", sep = "")
filename_formatted_peak_table <- paste(file_prefix, "_formatted_peak_table.txt", sep = "")
filename_formatted_annotation_table <- paste(file_prefix, "_formatted_annotation_table.txt", sep = "")
filename_formatted_metadata_table <- paste(file_prefix, "_formatted_metadata_table.txt", sep = "")
filename_foldchange_pvalues <- paste(file_prefix, "_foldchange_pvalues.csv", sep = "")


## We save the used params.yaml

message("Writing params.yaml ...")


file.copy(path_to_params, file.path(output_directory,filename_params), overwrite = TRUE)


# We move to the output directory

setwd(output_directory)



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
  select(feature_id_full, contains(" Peak area")) %>%
  rename_with(~gsub(" Peak area", "", .x)) %>%
  column_to_rownames(var = "feature_id_full") %>%
  as.data.frame() %>%
  t()

# We keep the feature_table_intensities dataframe in a separate variable

X = feature_table_intensities

# We keep the feature metadata in a separate dataframe

feature_metadata = feature_table %>%
  select(feature_id_full, feature_id, feature_mz, feature_rt)

############################### load annotation tables #####################################
############################################################################################


# The Sirius data is loaded

data_sirius = read_delim(file.path(working_directory, "results", "sirius", params$filenames$sirius_annotations),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

# The column names are modified to include the source of the data

colnames(data_sirius) = paste(colnames(data_sirius), "sirius", sep = "_")

# We now build a unique feature_id for each feature in the Sirius data

data_sirius$feature_id = sub("^.*_([[:alnum:]]+)$", "\\1", data_sirius$id_sirius)
data_sirius$feature_id = as.numeric(data_sirius$feature_id)

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

# We now convert the VM tibble into a dataframe and set the `feature_id_full` column as the rownames

VM = as.data.frame(VM)
row.names(VM) = VM$feature_id_full


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
  stop("Some rownames in X are not present in the rownames of SMDF. Please check the rownames in X and the rownames of SMDF.")
}


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

if (params$actions$filter_sample_type == "TRUE") {

filter_smeta_model <- filter_smeta(mode = params$filter_sample_type$mode,
                          factor_name = params$filter_sample_type$factor_name,
                          levels = params$filter_sample_type$levels)

# apply model sequence
filter_smeta_result = model_apply(filter_smeta_model, DE_original)

DE_original = DE_original_filtered@filtered

}

if (params$actions$filter_sample_metadata_one == "TRUE") {

MS_filter <- filter_smeta(mode = params$filter_sample_metadata_one$mode,
                          factor_name = params$filter_sample_metadata_one$factor_name,
                          levels = params$filter_sample_metadata_one$levels)

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

DE_original = DE_original_filtered@filtered

}

if (params$actions$filter_sample_metadata_two == "TRUE") {

MS_filter <- filter_smeta(mode = params$filter_sample_metadata_two$mode,
                          factor_name = params$filter_sample_metadata_two$factor_name,
                          levels = params$filter_sample_metadata_two$levels)

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

DE_original = DE_original_filtered@filtered

}


if (params$actions$filter_variable_metadata_one == "TRUE") {

MS_filter <- filter_vmeta(mode = params$filter_variable_metadata_one$mode,
                          factor_name = params$filter_variable_metadata_one$factor_name,
                          levels = params$filter_variable_metadata_one$levels)

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

DE_original = DE_original_filtered@filtered

}

if (params$actions$filter_variable_metadata_two == "TRUE") {

MS_filter <- filter_vmeta(mode = params$filter_variable_metadata_two$mode,
                          factor_name = params$filter_variable_metadata_two$factor_name,
                          levels = params$filter_variable_metadata_two$levels)

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

DE_original = DE_original_filtered@filtered

}


if (params$actions$scale_data == "FALSE") {

DE = DE_original

# We display the properties of the DatasetExperiment object to the user.
message("DatasetExperiment object properties: ")

sink(filename_DE_model)

print(DE)

sink() } else if (params$actions$scale_data == "TRUE") {

# Overall Pareto scaling (test)

M = pareto_scale()
M = model_train(M,DE_original)
M = model_predict(M,DE_original)
DE = M$scaled

# We display the properties of the DatasetExperiment object to the user.
message("DatasetExperiment object properties: ")

sink(filename_DE_model)

print(DE)

sink() 
} else {
  stop("Please check the value of the 'scale_data' parameter in the params file.")
}


################################################################################################
################################################################################################
######################## structool box formatted data export 

formatted_peak_table <- DE$data

formatted_annotation_table <- DE$variable_meta ### need to be filter with only usefull output
col_filter <- c("feature_id", "row_ID" ,"row_mz_full" ,"row_rt_full","molecularFormula_sirius","InChIkey2D_sirius","InChI_sirius",
"name_sirius","smiles_sirius", "pubchemids_sirius", "molecularFormula_canopus", "NPC.pathway_canopus","NPC.pathway.Probability_canopus",
"NPC.superclass_canopus", "NPC.class_canopus","ClassyFire.most.specific.class_canopus","ClassyFire.most.specific.class.Probability_canopus",
"ClassyFire.level.5_canopus","ClassyFire.subclass_canopus","ClassyFire.class_canopus","ClassyFire.superclass_canopus","ClassyFire.all.classifications_canopus",
"...1_metannot","structure_wikidata_metannot","structure_inchikey_metannot","structure_inchi_metannot","structure_smiles_metannot","structure_molecular_formula_metannot",
"short_inchikey_metannot","structure_taxonomy_npclassifier_01pathway_metannot","structure_taxonomy_npclassifier_02superclass_metannot",
"structure_taxonomy_npclassifier_03class_metannot","organism_wikidata_metannot","organism_name_metannot","organism_taxonomy_ottid_metannot","organism_taxonomy_01domain_metannot",
"organism_taxonomy_02kingdom_metannot","organism_taxonomy_03phylum_metannot","organism_taxonomy_04class_metannot","organism_taxonomy_05order_metannot",
"organism_taxonomy_06family_metannot","organism_taxonomy_07tribe_metannot","organism_taxonomy_08genus_metannot","organism_taxonomy_09species_metannot","organism_taxonomy_10varietas_metannot",
"matched_domain_metannot","matched_kingdom_metannot","matched_phylum_metannot","matched_class_metannot","matched_order_metannot","matched_family_metannot","matched_tribe_metannot",
"matched_genus_metannot","matched_species_metannot","score_taxo_metannot","score_max_consistency_metannot","final_score_metannot","rank_final_metannot","component_id_metannot",
"structure_taxonomy_npclassifier_01pathway_consensus_metannot","freq_structure_taxonomy_npclassifier_01pathway_metannot","structure_taxonomy_npclassifier_02superclass_consensus_metannot",
"freq_structure_taxonomy_npclassifier_02superclass_metannot","structure_taxonomy_npclassifier_03class_consensus_metannot","freq_structure_taxonomy_npclassifier_03class_metannot")
formatted_annotation_table_filtered <- formatted_annotation_table[col_filter]

formatted_metadata_table <- DE$sample_meta

write.table(formatted_peak_table, file = filename_formatted_peak_table, sep = ",")
write.table(formatted_annotation_table_filtered, file = filename_formatted_annotation_table, sep = ",")
write.table(formatted_annotation_table_filtered, file = filename_formatted_metadata_table, sep = ",")


#################################################################################################
#################################################################################################
#################################################################################################
##### PCA filtered data #######################################################################

message("Launching PCA calculations ...")


MS_PCA <- 
  filter_na_count(threshold = 1, factor_name = "sample_type") +
  knn_impute(neighbours = 5) +
  vec_norm() +
 #log_transform(base = 10) +
  mean_centre() +
  PCA(number_components = 3)


# apply model sequence
DE_PCA = model_apply(MS_PCA, DE)

# Fetching the PCA data object
DATA_PCA = DE_PCA[length(DE_PCA)]

# PCA scores plot

C = pca_scores_plot(
  factor_name = params$filters$metadata_variable,
  label_factor = "sample_id",
  ellipse_type = "t",
  ellipse_confidence = 0.9,
  points_to_label = "all"
)

# plot
PCA = chart_plot(C, DE_PCA[length(DE_PCA)])


fig_PCA = ggplotly(PCA + theme_classic() + facet_wrap(~ PCA$labels$title) + ggtitle(title_PCA))

# We merge PCA scores and metadata info in a single df

PCA_meta = merge(x = DATA_PCA$scores$sample_meta, y = DATA_PCA$scores$data, by = 0, all = TRUE)


fig_PCA3D = plot_ly(PCA_meta, x = ~PC1, y = ~PC2, z = ~PC3, color = PCA_meta[,params$filters$metadata_variable])
fig_PCA3D = fig_PCA3D %>% add_markers()
fig_PCA3D = fig_PCA3D %>% layout(scene = list(
  xaxis = list(title = "PC1"),
  yaxis = list(title = "PC2"),
  zaxis = list(title = "PC3")
),
legend = list(title=list(text=params$filters$metadata_variable)),
title = title_PCA3D
)


# The files are exported

fig_PCA %>%
    htmlwidgets::saveWidget(file = filename_PCA, selfcontained = TRUE)
fig_PCA3D %>%
    htmlwidgets::saveWidget(file = filename_PCA3D, selfcontained = TRUE)



# #################################################################################################
# #################################################################################################
# #################################################################################################
# ##### PLSDA filtered data #######################################################################

# # prepare model sequence
# M = autoscale() + PLSDA(factor_name='genotype')
# M = model_apply(M,DE)


# C = pls_scores_plot(factor_name = 'genotype')
# chart_plot(C,M[2])



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

data_RF = DE_original
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

cols = data_PCOA_merge[params$filters$metadata_variable]
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
legend = list(title=list(text=params$filters$metadata_variable)))


# The files are exported

ggsave(plot = fig_PCoA, filename = filename_PCoA, width = 10, height = 10)

fig_PCoA3D %>%
    htmlwidgets::saveWidget(file = filename_PCoA3D, selfcontained = TRUE)


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

message("Launching Volcano Plots calculations ...")


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
# sample_name = paste(data_RF$sample_meta$sample_id, data_RF$sample_meta[[params$filters$metadata_variable]], sep = "_")

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





# matt_trait = data_RF$sample_meta[params$filters$metadata_variable]
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


# The formula is defined externally
formula = as.formula(paste0('y', '~', params$filters$metadata_variable, '+' ,
'Error(sample_id/',
 params$filters$metadata_variable,
 ')'
)
)

DE = DE_original


HSDEM_model = HSDEM(formula = formula, mtc = 'none')

HSDEM_result = model_apply(HSDEM_model,DE)

HSDEM_result_p_value = HSDEM_result$p_value


glimpse(HSDEM_result_p_value)

# We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.

colnames(HSDEM_result_p_value) = plotrix::pasteCols(sapply(strsplit(colnames(HSDEM_result_p_value), " - "), sort), sep = "_")

# We now add a specific suffix (`_p_value`) to each of the colnames

colnames(HSDEM_result_p_value) = paste0(colnames(HSDEM_result_p_value), "_p_value")

str(HSDEM_result_p_value)

# We set the row names as columns row_id to be able to merge the two dataframes

HSDEM_result_p_value$row_id = rownames(HSDEM_result_p_value)

# # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

# HSDEM_result_p_value_long = pivot_longer(HSDEM_result_p_value, cols = -row_id, names_to = "pairs", values_to = "p_value")

# 

str(HSDEM_result_p_value_long)



fold_change_model = fold_change(
  factor_name= params$filters$metadata_variable,
  paired = FALSE,
  sample_name = character(0),
  threshold = 0.5,
  control_group = character(0),
  method = "geometric",
  conf_level = 0.95
  )

fold_change_result = model_apply(fold_change_model, DE)

# We suffix the column name of the dataframe with `_fold_change`, using dplyr rename function

fold_change_result_fold_change = fold_change_result$fold_change

glimpse(fold_change_result_fold_change)
# We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.
#n !!!! We need to make sure that the header of metadata variable is not in the colnames of the fold change result



colnames(fold_change_result_fold_change) = plotrix::pasteCols(sapply(strsplit(colnames(fold_change_result_fold_change), "/"), sort), sep = "_")


# We now add a specific suffix (`_p_value`) to each of the colnames

colnames(fold_change_result_fold_change) = paste0(colnames(fold_change_result_fold_change), "_fold_change")

# We set the row names as columns row_id to be able to merge the two dataframes

fold_change_result_fold_change$row_id = rownames(fold_change_result_fold_change)

# # We pivot the data from wide to long using the row_id as identifier and the colnames as variable

# fold_change_result_fold_change = pivot_longer(fold_change_result_fold_change, cols = -row_id, names_to = "pairs", values_to = "fold_change")



# We merge the two dataframes according to both the row_id and the pairs columns. 

DE_foldchange_pvalues = merge(HSDEM_result_p_value, fold_change_result_fold_change,  by = "row_id")

glimpse(DE_foldchange_pvalues)

summary(DE_foldchange_pvalues)


# We add columns corresponding to the Log2 of the fold change column (suffix by fold_change). For this we use mutate_at function from dplyr package. We save the results in new columns with a novel suffix `_log2_FC`.


DE_foldchange_pvalues = DE_foldchange_pvalues %>%
     mutate( across(contains('_fold_change'), 
                    .fns = list(log2 = ~log2(.)),
                    .names = "{col}_{fn}" ) )  %>% 
     mutate( across(contains('_p_value'), 
                    .fns = list(minus_log10 = ~-log10(.)),
                    .names = "{col}_{fn}" ) )



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



#############################################################################
#############################################################################
############## Tree Map #####################################################
#############################################################################
#############################################################################

message("Preparing Tree Map ...")


matt_donust = matt_volcano_plot[matt_volcano_plot$p.value < params$posthoc$p_value, ]
matt_donust2 = matt_donust[!is.na(matt_donust$NPC.superclass_canopus), ]
matt_donust2$counter = 1


DE_foldchange_pvalues




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
  parent_value = NPC.pathway_canopus,
  value = NPC.superclass_canopus,
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

fig_treemap %>%
    htmlwidgets::saveWidget(file = filename_treemap, selfcontained = TRUE)


#############################################################################
#############################################################################
############## Random Forest ################################################
#############################################################################
#############################################################################

message("Launching Random Forest calculations ...")

sink(filename_random_forest_model)

imp_filter1 = paste("X", matt_volcano_tot$mol[matt_volcano_tot$p.value < params$posthoc$p_value ], sep = "")

data_subset_norm_rf_filter = data_subset_norm_rf[, colnames(data_subset_norm_rf) %in% c("treatment", imp_filter1)]

ozone.rp = rfPermute(treatment ~ ., data = data_subset_norm_rf_filter, na.action = na.omit, ntree = 500, num.rep = 500)
imp_table_rf = data.frame(ozone.rp$pval)
imp_table_rf = importance(ozone.rp)
imp_table_rf = data.frame(imp_table_rf)
summary(ozone.rp)

f = plotImportance(ozone.rp, plot.type = "bar", plot = FALSE)

sink()

########### plot importance

fig_rf = ggplotly(f[[(length(f) - 1)]] + theme_classic() + facet_wrap(~ f[[(length(f) - 1)]]$labels$title))
fig_rf = subplot(fig_rf) %>%
  layout(title = title_random_forest)


# The file is exported
# The title should be updated !!! 


fig_rf %>%
    htmlwidgets::saveWidget(file = filename_random_forest , selfcontained = TRUE)


#############################################################################
#############################################################################
############## Random Forest selected Box Plots #############################
#############################################################################
#############################################################################

message("Preparing Box plots ...")


imp.scaled = rfPermute::importance(ozone.rp, scale = TRUE)
imp.scaled = data.frame(imp.scaled)
imp.scaled = imp.scaled[order(imp.scaled$MeanDecreaseGini, decreasing = TRUE), ]
boxplot_top_N = row.names(imp.scaled)[1:params$boxplot$topN]
data_subset_norm_boxplot = data_subset_norm_rf[, colnames(data_subset_norm_rf) %in% boxplot_top_N]
data_subset_norm_boxplot$treatment = as.factor(vec_trait)
data_subset_norm_boxplot = reshape2::melt(data_subset_norm_boxplot, id = c("treatment"))

p = ggplot(data = data_subset_norm_boxplot, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = treatment)) +
  theme_bw() +
  ggtitle(title_box_plots)

fig_boxplot = p + facet_wrap(~variable, scales = "free", dir = "v") + theme(
  legend.position = "top",
  legend.title = element_blank()
)


fig_boxplotly = data_subset_norm_boxplot %>%
  split(data_subset_norm_boxplot$variable) %>%
  map(~ {
    plot_ly(data = .x, x = .x$treatment, y = .x$value, color = .x$treatment, type = "box") %>%
      layout(showlegend = F, xaxis = list(title = .x$variable[1])) %>%
      layout(xaxis = list(titlefont = list(size = 10), tickfont = list(size = 10))) %>%
      layout(yaxis = list(titlefont = list(size = 12), tickfont = list(size = 12)))
  }) %>%
  subplot(margin = 0.02, nrows = 4, titleX = TRUE)  %>% 
  layout(title = title_box_plots,
  legend = list(title = list(text = paste("<b>class </b>")))) # Find a way to define the legend title


# The files are exported

ggsave(plot = fig_boxplot, filename = filename_box_plots, width = 10, height = 10)
fig_boxplotly %>%
    htmlwidgets::saveWidget(file = filename_box_plots_interactive, selfcontained = TRUE)


#############################################################################
#############################################################################
############## Random Forest filtered Heat Map  #############################
#############################################################################
#############################################################################

message("Preparing Random Forest filtered Heatmap ...")

imp_table_rf_order = imp.scaled[order(imp.scaled$MeanDecreaseGini, decreasing = TRUE), ] # imp_table_rf[order(imp_table_rf$MeanDecreaseGini,decreasing=TRUE),]


imp_filter2X = row.names(imp_table_rf_order)[1:params$heatmap$topN]
imp_filter2 = gsub("X", "", imp_filter2X)

data_subset_norm_rf_filtered = data_subset_norm_rf[, colnames(data_subset_norm_rf) %in% imp_filter2X]
my_sample_col = sample_name
annot_col = data.frame(paste(data_RF$variable_meta$NPC.pathway_canopus, data_RF$variable_meta$NPC.superclass_canopus, sep = "_"), data_RF$variable_meta$NPC.pathway_canopus)
colnames(annot_col) = c("Superclass", "Pathway")

rownames(annot_col) = data_RF$variable_meta$row_ID
annot_col_filter = annot_col[rownames(annot_col) %in% imp_filter2, ]


ByPal = colorRampPalette(c(wes_palette("Zissou1")))

data_subset_norm_rf_filtered = apply(data_subset_norm_rf_filtered, 2, as.numeric)
# heatmap(as.matrix(data_subset_norm_rf_filtered), scale="column")


heatmap_filtered = heatmaply(
  percentize(data_subset_norm_rf_filtered),
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

heatmap_filtered = heatmap_filtered %>% layout(title = list(text = title_heatmap, y = 0.05))

# The file is exported

heatmap_filtered %>%
    htmlwidgets::saveWidget(file = filename_heatmap, selfcontained = TRUE)


#############################################################################
#############################################################################
############## Summary Table ################################################
#############################################################################
#############################################################################

message("Outputing Summary Table ...")


# Output is not clean. Feature id are repeated x times. 
# To tidy ---

feature_id = data_RF$variable_meta$feature_id
sample_raw_id = paste("X", data_RF$variable_meta$row_ID, sep = "")
name_sirius = data_RF$variable_meta$name_sirius
smiles_sirius = data_RF$variable_meta$smiles_sirius
InChI_sirius = data_RF$variable_meta$InChI_sirius
NPC.pathway_canopus = data_RF$variable_meta$NPC.pathway_canopus
NPC.superclass_canopus = data_RF$variable_meta$NPC.superclass_canopus
Annotation_merge = data.frame(feature_id, sample_raw_id, name_sirius, smiles_sirius, InChI_sirius, NPC.pathway_canopus, NPC.superclass_canopus)


RF_importance = imp_table_rf$MeanDecreaseGini
sample_raw_id = row.names(imp_table_rf)
RF_importance_merge = data.frame(sample_raw_id, RF_importance)

summary_matt1 = merge(Annotation_merge, RF_importance_merge, by = "sample_raw_id", all = TRUE)

sample_raw_id = paste("X", matt_volcano_tot$mol, sep = "")
group = matt_volcano_tot$X1
posthoc_Pvalue = matt_volcano_tot$p.value
log10P = matt_volcano_tot$log10P
posthoc_result_merge = data.frame(sample_raw_id, group, posthoc_Pvalue, log10P)
split_group = split(posthoc_result_merge, group)
matt_split = do.call("cbind", split_group)
matt_split$sample_raw_id = matt_split[, 1]

summary_stat_output = merge(summary_matt1, matt_split, by = "sample_raw_id", all = TRUE)


# The file is exported

write.table(summary_stat_output, file = filename_summary_stats_table, sep = ",")


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


# glimpse(df_from_graph_vertices)

# We then add the attributes to the vertices dataframe
# For this we merge the vertices dataframe with the VM output using the id column and the feature_id column, respectively

df_from_graph_vertices_plus = merge(df_from_graph_vertices, VM, by.x = "id", by.y = "feature_id", all.x = T)

# Now we will add the results of the statistical outputs to the vertices dataframe
# For this we merge the vertices dataframe with the summary_stat_output using the id column and the feature_id column, respectively

# First we clean the summary_stat_output dataframe
# For this we remove columns that are not needed. The one containing the sirius and canopus pattern in the column names. Indeed they arr already present in the VM dataframe

summary_stat_output = summary_stat_output %>%
  select(-contains("sirius")) %>%
  select(-contains("canopus")) %>% 
  select(-contains("_raw_id"))


df_from_graph_vertices_plus = merge(df_from_graph_vertices_plus, summary_stat_output, by.x = "id", by.y = "feature_id", all.x = T)


# We merge the data from the DE$data dataframe with the DE$sample_meta dataframe using rownames as the key

merged_D_SM = merge(DE_original$sample_meta, DE_original$data, by = 'row.names', all = TRUE)

# We replace NA values with 0 in the merged dataframe
# Check why we dont do this before ??
merged_D_SM[is.na(merged_D_SM)] <- 0


# The function below weill allow to group data by multiple factors and return a dataframe with the mean of each group

dfList <- list()

for (i in params$colnames$to_output) {
  dfList[[i]] <- merged_D_SM %>%
    group_by(!!as.symbol(i)) %>%
    summarise(across(contains("_peak"), mean),
      .groups = "drop"
    ) %>%
    select(!!all_of(i), contains("_peak")) %>%
    pivot_longer(-!!i) %>%
    pivot_wider(names_from = all_of(i), values_from = value)  %>% 
    # We prefix all columns with the factor name
    rename_with(.cols=-name, ~paste0("mean_int", "_", i, "_", .x))
}


flat_dfList = reduce(dfList, full_join, by = "name")

# We now add the raw feature list to the dataframe

full_flat_dfList = merge(flat_dfList, t(DE_original$data), by.x = 'name', by.y = 'row.names', all = TRUE)


# We add the raw feature list

df_from_graph_vertices_plus_plus = merge(df_from_graph_vertices_plus, full_flat_dfList, by.x = "row_ID", by.y = "name", all.x = T)

# We set back the id column as the first column of the dataframe

df_from_graph_vertices_plus_plus = df_from_graph_vertices_plus_plus %>%
  select(id, everything())


# We then add the attributes to the edges dataframe and generate the igraph object

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
