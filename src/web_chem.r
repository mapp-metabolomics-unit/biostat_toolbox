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
usePackage("funModeling")
usePackage("rcdk")
usePackage("gt")
usePackage("purrr")


usePackage("webchem")



# We use the MAPPstructToolbox package 
# Uncomment the lines below to download the MAPPstructToolbox package from github

#library(devtools)
#install_github("mapp-metabolomics-unit/MAPPstructToolbox", force = TRUE)
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

path_to_params = "./params/params_webchem.yaml"

params = yaml.load_file(path_to_params)



# We set the working directory

working_directory = file.path(params$path$docs, params$mapp_project, params$mapp_batch, params$polarity)



# We fetch the sirius data

data_sirius = read_delim(file.path(working_directory, "results", "sirius", params$filenames$sirius_annotations),
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

# The column names are modified to include the source of the data

colnames(data_sirius) = paste(colnames(data_sirius), "sirius", sep = "_")

# We now build a unique feature_id for each feature in the Sirius data

data_sirius$feature_id = sub("^.*_([[:alnum:]]+)$", "\\1", data_sirius$id_sirius)
data_sirius$feature_id = as.numeric(data_sirius$feature_id)


data_sirius_sub = head(data_sirius, 100)



cts_convert("XEFQLINVKFYRCS-UHFFFAOYSA-N", "inchikey", "Chemical Name")
### multiple inputs


keys <- c("XEFQLINVKFYRCS-UHFFFAOYSA-N", "VLKZOEOYAKHREP-UHFFFAOYSA-N")
cts_convert(keys, "inchikey", "cas")


get_chebiid('BPGDAMSIGCZZLK')

comp <- c('Aspirin')

comp <- data_sirius$smiles_sirius


chebi_ids = get_chebiid(comp, from = 'smiles', to = 'chebiid', match = "best")


View(chebi_ids)

write.table(chebi_ids, file = 'chebi_ids.csv', sep = ",", row.names = FALSE)


# We first remove duplicates form the Sirius smiles columns

for_chembiid_smiles = unique(data_sirius$smiles_sirius)

chebi_ids = get_chebiid(for_chembiid_smiles, from = 'smiles', to = 'chebiid', match = "best")


# We first remove duplicates form the Sirius smiles columns

for_chembiid_inchikey2d = unique(data_sirius$InChIkey2D_sirius)

chebi_ids_from_ik = get_chebiid(for_chembiid_inchikey2d, from = 'inchikey', to = 'chebiid', match = "best")

write.table(chebi_ids_from_ik, file = 'chebi_ids_from_ik.csv', sep = ",", row.names = FALSE)

length(for_chembiid_inchikey2d)


chebiids = head(chebi_ids_from_ik$chebiid, 10)


chebiids_infos = chebi_comp_entity(chebiids)

chebi_comp_entity

glimpse(chebiids_infos)

# We merge the chebi ids with the sirius data according to two columns

data_sirius_chebi = merge(data_sirius, chebi_ids_from_ik, by.x = "InChIkey2D_sirius", by.y = "query")

data_sirius_chebi = merge(data_sirius_chebi, chebi_ids, by.x = "smiles_sirius", by.y = "query")

glimpse(data_sirius_chebi)


write.table(data_sirius_chebi, file = 'data_sirius_chebi.csv', sep = ",", row.names = FALSE)