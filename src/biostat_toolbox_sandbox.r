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


library(MAPPstructToolbox)


path_to_params = "./params/params.yaml"

params = yaml.load_file(path_to_params)



working_directory = file.path(params$path$docs, params$mapp_project, params$mapp_batch, params$polarity)



data = read_delim(file.path(working_directory,  "results", "mzmine", paste0(params$mapp_batch, "_quant.csv")),
  delim = ",", escape_double = FALSE,
  trim_ws = TRUE
)

data$"row m/z" = round(data$"row m/z", digits = 2)
data$"row retention time" = round(data$"row retention time", digits = 1)

feature_id = data$"row ID"
row_mz_full = data$"row m/z"
row_rt_full = data$"row retention time"

data$"row ID" = paste(data$"row ID",
  data$"row m/z",
  data$"row retention time",
  "peak",
  sep = "_"
)

row_ID = data$"row ID"

data = data %>%
  select(-"row m/z", -"row retention time")

data = data[1:(length(data) - 1)]

names(data) = gsub(" Peak area", "", names(data))

data_prep = data %>%
  remove_rownames() %>%
  column_to_rownames(var = "row ID")

data_t = as.data.frame(t(data_prep))

feature_list = data.frame(feature_id, row_ID, row_mz_full, row_rt_full)



data_sirius = read_delim(file.path(working_directory, "results", "sirius", params$filenames$sirius_annotations),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)


colnames(data_sirius) = paste(colnames(data_sirius), "sirius", sep = "_")

data_canopus = read_delim(file.path(working_directory, "results", "sirius", params$filenames$canopus_annotations),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

colnames(data_canopus) = paste(colnames(data_canopus), "canopus", sep = "_")


data_metannot = read_delim(file.path(working_directory, "results", "met_annot_enhancer", params$met_annot_enhancer_folder, paste0(params$met_annot_enhancer_folder, "_spectral_match_results_repond.tsv")),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)


colnames(data_metannot) = paste(colnames(data_metannot), "metannot", sep = "_")

colnames(data_metannot)[colnames(data_metannot) == "feature_id_metannot"] = "feature_id"



VM_sir_can = merge(x = data_sirius, y = data_canopus, by.x = "id_sirius", by.y = "id_canopus", all = TRUE)

VM_sir_can$feature_id = sub("^.*_([[:alnum:]]+)$", "\\1", VM_sir_can$id_sirius)

VM_sir_can_metannot = merge(x = VM_sir_can, y = data_metannot, by = "feature_id", all = TRUE)

VM_sir_can_metannot_full = merge(x = feature_list, y = VM_sir_can_metannot, by = "feature_id", all = TRUE)

VM = VM_sir_can_metannot_full

row.names(VM) = VM$row_ID

X = data_t[, 1:ncol(data_t)]


metadata = read_delim(file.path(working_directory, "metadata", "treated", paste(params$mapp_batch,  "metadata.txt", sep = "_")),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

SM = data.frame(metadata)

df = SM %>% 
  filter(sample_type == "sample")


cols = c(params$colnames$to_combine)


for (n in 1:length(cols)) {
  combos = combn(cols, n, simplify = FALSE)
  for (combo in combos) {
    new_col_name = paste(combo, collapse = "_")
    df[new_col_name] = apply(df[combo], 1, paste, collapse = "_")
  }
}


SM = merge(x = SM, y = df,  all.x = TRUE)


SM[is.na(SM)] = "ND"



SMDF = as.data.frame(SM)

SMDF = SMDF %>%
  remove_rownames() %>%
  column_to_rownames(var = "filename")

XSM = merge(x = X, y = SMDF, by = "row.names", all = TRUE)



X_pond = XSM %>%
  select(Row.names, grep("peak", colnames(XSM))) %>%
  column_to_rownames(var = "Row.names")

X_pond = X_pond[order(row.names(X_pond)), ]
X_pond = X_pond[, order(colnames(X_pond))]
SMDF = SMDF[order(row.names(SMDF)), ]
VM = VM[order(row.names(VM)), ]


if (any(colnames(X_pond) != row.names(VM))) {
  stop("The order of the columns in X_pond is not the same as the order of the rows in VM. Please check the order of the columns in X_pond and the order of the rows in VM.")
}

# We repeat for row.names(SMDF) == row.names(X_pond)

if (any(row.names(SMDF) != row.names(X_pond))) {
  stop("The order of the rows in SMDF is not the same as the order of the rows in X_pond. Please check the order of the rows in SMDF and the order of the rows in X_pond.")
}




# The DatasetExperiment object is created using the X_pond, SMDF and VM objects.

DE_original = DatasetExperiment(
  data = X_pond,
  sample_meta = SMDF,
  variable_meta = VM,
  name = params$dataset_experiment$name,
  description = params$dataset_experiment$description
)


## Filtering steps

if (params$actions$filter_sample_type == "TRUE") {

MS_filter <- filter_smeta(mode = params$filter_sample_type$mode,
                          factor_name = params$filter_sample_type$factor_name,
                          levels = params$filter_sample_type$levels)

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

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



D = MTBLS79_DatasetExperiment()
M = fold_change(factor_name='Class')
M = model_apply(M,D)

DE_original

M = fold_change(factor_name='genotype')
M = model_apply(M,DE_original)


D = iris_DatasetExperiment()
D$sample_meta$class = D$sample_meta$Species

M = fold_change(factor_name='Species')
M = model_apply(M,D)


M$fold_change

M

C = fold_change_plot(number_features = 2, orientation = 'landscape')

chart_plot(C,M)

C = pca_scores_plot(factor_name='class') # colour by class
chart_plot(C,M[2])



------


M = fold_change(factor_name='genotype')
M = model_apply(M,D)



MS_filter <- filter_smeta(mode = 'include',
                          factor_name = 'sample_type',
                          levels = 'sample')

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

DE_original = DE_original_filtered@filtered


head(DE_original$sample_meta)


M = fold_change(factor_name='genotype')
M = model_apply(M,DE_original)


C = fold_change_plot(number_features = 300, orientation = 'landscape')

chart_plot(C,M)



M = fold_change_int(factor_name=c('genotype','brainregion'))
M = model_apply(M,DE_original)

chart_plot(C,M)


glimpse(M$fold_change)
glimpse(M$significant)
