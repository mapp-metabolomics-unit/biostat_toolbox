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



suppressPackageStartupMessages({
    # Bioconductor packages
    library(MAPPstructToolbox)
    library(pmp)
    library(ropls)
    library(BiocFileCache)
  
    # CRAN libraries
    library(ggplot2)
    library(gridExtra)
    library(cowplot)
    library(openxlsx)
    library(plotly)

})


path_to_params = "./params/params.yaml"

params = yaml.load_file(path_to_params)



working_directory = file.path(params$path$docs, params$mapp_project, params$mapp_batch, params$polarity)



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


# We here load the sample metadata

sample_metadata = read_delim(file.path(working_directory, "metadata", "treated", paste(params$mapp_batch,  "metadata.txt", sep = "_")),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

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

F = fold_change(factor_name='Species')
F = model_apply(F,D)



M$threshold

M$


F$fold_change


T = ttest(factor_name='Species')
T = model_apply(T,D)



D = iris_DatasetExperiment()
A = ANOVA(formula=y~Species)
A = model_apply(A,D)


D = iris_DatasetExperiment()
D$sample_meta$id=rownames(D) # dummy id column
M = HSDEM(formula = y~Species+ Error(id/Species))
M = model_apply(M,D)

M$p_value

head(M$significant)



M$description


DE_original$sample_meta$id=rownames(DE_original) # dummy id column


MS_filter <- filter_smeta(mode = 'include',
                          factor_name = 'sample_type',
                          levels = 'sample')

# apply model sequence
DE_original_filtered = model_apply(MS_filter, DE_original)

DE_original = DE_original_filtered@filtered




M = HSDEM(formula = y~genotype + Error(id/genotype))

M = model_apply(M,DE_original)

M$p_value


A$p_value




M$significant


M$conf_level

M$upper_ci



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




####### Testing PCA and PLS-DA

DE = iris_DatasetExperiment()



M = autoscale() + PCA(number_components = 3)
# apply model sequence to dataset
M = model_apply(M,DE)

# pca scores plots
g=list()
for (k in colnames(DE$sample_meta)) {
    C = pca_scores_plot(factor_name = k)
    g[[k]] = chart_plot(C,M[2])
}
# plot using cowplot
plot_grid(plotlist=g, nrow=1, align='vh', labels=c('A','B','C'))



# prepare model sequence
M = autoscale() + PLSDA(factor_name='Species')
M = model_apply(M,DE)


C = pls_scores_plot(factor_name = 'Species')
chart_plot(C,M[2])



str(DE_original$sample_meta)

DE_original$sample_meta$age = as.factor(DE_original$sample_meta$age)
DE_original$sample_meta$genotype = as.factor(DE_original$sample_meta$genotype)

# prepare model sequence
M = autoscale() + PLSDA(factor_name='age')
M = model_apply(M,DE_original)


C = pls_scores_plot(factor_name = 'age')
chart_plot(C,M[2])


C = pls_vip_plot(ycol='young')
vip_plot <- chart_plot(C,M[2])

ggplotly(vip_plot)


C = plsda_feature_importance_plot(n_features=30,metric='sr')
sr_plot <- chart_plot(C,M[2])

ggplotly(sr_plot)

C = plsda_roc_plot(factor_name='age')
chart_plot(C,M[2])


D = iris_DatasetExperiment()
DE_original$sample_meta$run_order=1:nrow(DE_original)
C = tic_chart(factor_name='age',run_order='run_order')
chart_plot(C,DE_original)

D = iris_DatasetExperiment()
M = HCA(factor_name='Species')
M = model_apply(M,D)

C = hca_dendrogram()
chart_plot(C,M)


D = DE_original
M = HCA(factor_name='genotype')
M = model_apply(M,D)

C = hca_dendrogram()
chart_plot(C,M)


D = DE
C = DatasetExperiment_boxplot(factor_name='age',number=50,per_class=TRUE)
chart_plot(C,D)

D = DE
C = DatasetExperiment_dist(factor_name='age')
chart_plot(C,D)


D = iris_DatasetExperiment()
C = DatasetExperiment_factor_boxplot(factor_names='Species',feature_to_plot='Petal.Width')
chart_plot(C,D)

D = DE
C = DatasetExperiment_factor_boxplot(factor_names='age',feature_to_plot='102_537.17_3.8_peak')
c <- chart_plot(C,D)
ggplotly(c)

D = DE
C = DatasetExperiment_heatmap()
chart_plot(C,D)


