############################################################################################
############################################################################################
#####################################     PACKAGES     #####################################
############################################################################################
############################################################################################

suppressPackageStartupMessages({
  library("crosstalk")
  library("digest")
  library("dplyr")
  library("DT")
  library("emmeans")
  library("forcats")
  library("ggh4x")
  library("ggplot2")
  library("ggrepel")
  library("htmltools")
  library("iheatmapr")
  library("janitor")
  library("jsonlite")
  library("microshades")
  library("nlme")
  library("plotly")
  library("pls")
  library("pmp")
  library("purrr")
  library("readr")
  library("rfPermute")
  library("rlang")
  library("rockchalk")
  library("stringr")
  library("svglite")
  library("tibble")
  library("tidyr")
  library("vegan")
  library("httr")
  library("webchem")
  library("wesanderson")
  library("WikidataQueryServiceR")
  library("optparse")
  library("yaml")
  library(MAPPstructToolbox)
})

############################################################################################
############################################################################################
################################ LOAD REQUIRED FUNCTIONS  ##################################
############################################################################################
############################################################################################

# We load the required functions from the MAPPstructToolbox package
# these are in the helpers.r file

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
if (length(script_path)) {
  script_path <- normalizePath(script_path[1])
} else {
  script_path <- normalizePath(file.path(getwd(), "src", "biostat_toolbox.r"), mustWork = FALSE)
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

source(file.path(script_dir, "helpers.r"))

option_list <- list(
  make_option(c("-p", "--params"), default = file.path(repo_root, "params", "params.yaml"), help = "Path to params.yaml [default repo params/params.yaml]"),
  make_option(c("-u", "--params-user"), default = file.path(repo_root, "params", "params_user.yaml"), help = "Path to params_user.yaml [default repo params/params_user.yaml]"),
  make_option(c("--params-overlay"), default = NULL, help = "Optional YAML recursively merged over --params")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

normalize_option_name <- function(opt, underscore_name, hyphen_name) {
  if (is.null(opt[[underscore_name]]) && !is.null(opt[[hyphen_name]])) {
    opt[[underscore_name]] <- opt[[hyphen_name]]
  }
  opt
}

opt <- normalize_option_name(opt, "params_user", "params-user")
opt <- normalize_option_name(opt, "params_overlay", "params-overlay")

resolve_relative_path <- function(path_value, fallback_dir) {
  if (is.null(path_value) || !length(path_value) || !nzchar(path_value[1])) {
    return(path_value)
  }
  path_value <- path_value[1]
  if (grepl("^/", path_value)) {
    return(path_value)
  }
  if (file.exists(path_value)) {
    return(normalizePath(path_value))
  }
  candidate <- file.path(fallback_dir, path_value)
  if (file.exists(candidate)) {
    return(normalizePath(candidate))
  }
  normalizePath(candidate, mustWork = FALSE)
}

path_to_params <- resolve_relative_path(opt$params, repo_root)
path_to_params_user <- resolve_relative_path(opt$params_user, repo_root)
path_to_params_overlay <- resolve_relative_path(opt$params_overlay, repo_root)

if (!file.exists(path_to_params)) {
  stop(sprintf("params.yaml not found: %s", path_to_params))
}
if (!file.exists(path_to_params_user)) {
  stop(sprintf("params_user.yaml not found: %s", path_to_params_user))
}
if (!is.null(path_to_params_overlay) && !file.exists(path_to_params_overlay)) {
  stop(sprintf("params overlay YAML not found: %s", path_to_params_overlay))
}


############################################################################################
############################################################################################
################################ LOAD & FORMAT  DATA  ######################################
############################################################################################
############################################################################################


print(script_path)

if (!exists("params") || !exists("my_path_params")) {
  my_path_params <- getwd()
} ### conserve the path after multiple run
if (exists("params") && exists("my_path_params")) {
  setwd(my_path_params)
} ### conserve the path after multiple run

# Load the params.yaml file

params <- yaml.load_file(path_to_params)
if (!is.null(path_to_params_overlay)) {
  params <- utils::modifyList(params, yaml.load_file(path_to_params_overlay), keep.null = TRUE)
}
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

if (params$actions$scale_method == "none") {
  scaling_status <- "raw"
} else {
  scaling_status <- "scaled"
}

possible_modes <- c("exclude", "include", "above", "below", "activated", "deactivated")



filter_sample_type_status = formatted_filter_status(params$filter_sample_type)
filter_sample_metadata_one_status = formatted_filter_status(params$filter_sample_metadata_one)
filter_sample_metadata_two_status = formatted_filter_status(params$filter_sample_metadata_two)
filter_sample_metadata_three_status = formatted_filter_status(params$filter_sample_metadata_three)

filter_variable_metadata_one_status = formatted_filter_status(params$filter_variable_metadata_one)
filter_variable_metadata_two_status = formatted_filter_status(params$filter_variable_metadata_two)
filter_variable_metadata_annotated_status = formatted_filter_status(params$filter_variable_metadata_annotated)
filter_variable_metadata_num_status = formatted_filter_status(params$filter_variable_metadata_num)

filter_sample_metadata_status = paste(filter_sample_type_status, filter_sample_metadata_one_status, filter_sample_metadata_two_status, filter_sample_metadata_three_status, sep = "_")

filter_variable_metadata_status = paste(filter_variable_metadata_one_status, filter_variable_metadata_two_status, filter_variable_metadata_annotated_status, filter_variable_metadata_num_status, sep = "_")

# We make sure that no multiple _ exists using the sanitize_string function

filter_sample_metadata_status = sanitize_string(filter_sample_metadata_status)
filter_variable_metadata_status = sanitize_string(filter_variable_metadata_status)



#################################################################################################
#################################################################################################
################### Filename and paths establishment ##########################################
#################################################################################################


file_prefix <- paste("")


filename_box_plots <- paste(file_prefix, "Boxplots.pdf", sep = "")
filename_box_plots_svg <- paste(file_prefix, "Boxplots.svg", sep = "")
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
dir_npc_summed_intensity <- paste(file_prefix, "NPC_summed_intensity", sep = "")
dir_npc_summed_intensity_filtered <- file.path(dir_npc_summed_intensity, "filtered")
dir_npc_summed_intensity_raw <- file.path(dir_npc_summed_intensity, "raw")
filename_npc_summed_intensity_pdf <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity.pdf")
filename_npc_summed_intensity_png <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity.png")
filename_npc_summed_intensity_table <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity.tsv")
filename_npc_summed_intensity_stats_table <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity_stats.tsv")
filename_npc_summed_intensity_driver_table <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity_feature_drivers.tsv")
filename_npc_summed_intensity_ratio_pdf <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity_ratio.pdf")
filename_npc_summed_intensity_ratio_png <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity_ratio.png")
filename_npc_summed_intensity_ratio_table <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity_ratio.tsv")
filename_npc_summed_intensity_ratio_stats_table <- file.path(dir_npc_summed_intensity_filtered, "NPC_summed_intensity_ratio_stats.tsv")
filename_npc_summed_intensity_raw_pdf <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_raw.pdf")
filename_npc_summed_intensity_raw_png <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_raw.png")
filename_npc_summed_intensity_raw_table <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_raw.tsv")
filename_npc_summed_intensity_raw_stats_table <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_raw_stats.tsv")
filename_npc_summed_intensity_raw_driver_table <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_feature_drivers_raw.tsv")
filename_npc_summed_intensity_ratio_raw_pdf <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_ratio_raw.pdf")
filename_npc_summed_intensity_ratio_raw_png <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_ratio_raw.png")
filename_npc_summed_intensity_ratio_raw_table <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_ratio_raw.tsv")
filename_npc_summed_intensity_ratio_raw_stats_table <- file.path(dir_npc_summed_intensity_raw, "NPC_summed_intensity_ratio_raw_stats.tsv")
filename_params <- paste(file_prefix, "params.yaml", sep = "")
filename_params_user <- paste(file_prefix, "params_user.yaml", sep = "")
filename_PCA <- paste(file_prefix, "PCA.pdf", sep = "")
filename_PCA_svg <- paste(file_prefix, "PCA.svg", sep = "")
filename_PCA_scores <- paste(file_prefix, "PCA_scores.tsv", sep = "")
filename_PCA_loadings <- paste(file_prefix, "PCA_loadings.tsv", sep = "")
filename_PCA3D <- paste(file_prefix, "PCA3D.html", sep = "")
filename_PCoA <- paste(file_prefix, "PCoA.pdf", sep = "")
filename_PCoA_svg <- paste(file_prefix, "PCoA.svg", sep = "")
filename_PCoA3D <- paste(file_prefix, "PCoA3D.html", sep = "")
filename_PLSDA <- paste(file_prefix, "PLSDA.pdf", sep = "")
filename_PLSDA_svg <- paste(file_prefix, "PLSDA.svg", sep = "")
filename_PLSDA_loadings <- paste(file_prefix, "PLSDA_loadings.tsv", sep = "")
filename_PLSDA_scores <- paste(file_prefix, "PLSDA_scores.tsv", sep = "")
filename_PLSDA_VIP_plot <- paste(file_prefix, "PLSDA_VIP.pdf", sep = "")
filename_PLSDA_VIP_table <- paste(file_prefix, "PLSDA_VIP.tsv", sep = "")
filename_DFA <- paste(file_prefix, "DFA.pdf", sep = "")
filename_DFA_loadings <- paste(file_prefix, "DFA_loadings.tsv", sep = "")
filename_DFA_scores <- paste(file_prefix, "DFA_scores.tsv", sep = "")
filename_DFA_eigenvalues <- paste(file_prefix, "DFA_eigenvalues.tsv", sep = "")
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



###################################################################################################
######################### rename main folder - short version




# common_df_path <- file.path(params$paths$output, "params_log.rds")
common_tsv_path <- file.path(params$paths$output, "params_log.tsv")


  # Convert the YAML content to a dataframe row
new_row_df <- convert_yaml_to_single_row_df_with_hash(params)

# Append this row to the common dataframe and save
append_to_common_df_and_save(new_row_df, common_tsv_path)



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

# feature_table_intensities <- feature_table %>%
#   select(feature_id, contains(" Peak height")) %>%
#   rename_with(~ gsub(" Peak height", "", .x)) %>%
#   column_to_rownames(var = "feature_id") %>%
#   as.data.frame() %>%
#   t()

# We check if both ` Peak area` and ` Peak height` patterns are present in the dataframe

if (any(grepl(" Peak area", colnames(feature_table))) && any(grepl(" Peak height", colnames(feature_table)))) {
  warning("Both ` Peak area` and ` Peak height` patterns are present in the dataframe. Keeping only the ` Peak area` pattern.")
}


# We make the same operation but we make it work both for ` Peak area` and ` Peak height` patterns. If both are present in the dataframe we
# raise a warning and keep only the ` Peak area` pattern

feature_table_intensities <- feature_table %>%
  select(feature_id, contains(" Peak area"), contains(" Peak height")) %>%
  rename_with(~ gsub(" Peak area", "", .x)) %>%
  rename_with(~ gsub(" Peak height", "", .x)) %>%
  column_to_rownames(var = "feature_id") %>%
  as.data.frame() %>%
  t()





# We keep the feature_table_intensities dataframe in a separate variable

X <- feature_table_intensities


# We order the X by rownames and by column names

X <- X[order(row.names(X)), ]
X <- X[, order(colnames(X))]

X <- as.data.frame(X)



# min(X)

# Uncomment for testing purposes
# X <- X[,1:100]


# We keep the feature metadata in a separate dataframe

feature_metadata <- feature_table %>%
  select(feature_id_full, feature_id, feature_mz, feature_rt)

############################### load annotation tables #####################################
############################################################################################

# Sirius data is treated
# Determin sirius version. If structure_identifications.tsv exists in dir then version 6, else version 5

sirius_version <- if (file.exists(file.path(working_directory, "results", "sirius", "structure_identifications.tsv"))) {
  "6"
} else {
  "5"
}

# Sirius filenames

if (sirius_version == "6") {
  sirius_annotations_filename = "structure_identifications.tsv"
  canopus_annotations_filename = "canopus_structure_summary.tsv"
} else {
  sirius_annotations_filename = "compound_identifications.tsv"
  canopus_annotations_filename = "canopus_compound_summary.tsv"
}

# Check if a chebied version exists, if not we create it

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
  # We first remove duplicates from the Sirius smiles columns

  for_chembiid_smiles <- unique(data_sirius$smiles)

  print("Getting ChEBI IDs from SMILES via ChEBI REST API ...")

  # Query the new ChEBI backend REST API (replaces the defunct SOAP/webchem route)
  # Endpoint: https://www.ebi.ac.uk/chebi/backend/api/public/structure_search/
  get_chebi_from_smiles <- function(smiles_vec) {
    base_url <- "https://www.ebi.ac.uk/chebi/backend/api/public/structure_search/"
    n        <- length(smiles_vec)
    pb       <- txtProgressBar(min = 0, max = n, style = 3, width = 60)

    results <- lapply(seq_along(smiles_vec), function(i) {
      smi <- smiles_vec[[i]]
      setTxtProgressBar(pb, i)

      if (is.na(smi) || nchar(trimws(smi)) == 0) {
        return(data.frame(query = smi, chebiid = NA_character_, chebiasciiname = NA_character_,
                          stringsAsFactors = FALSE))
      }
      res <- tryCatch(
        httr::GET(base_url, query = list(
          smiles          = smi,
          search_type     = "connectivity",
          three_star_only = "false",
          size            = 1
        ), httr::timeout(30)),
        error = function(e) NULL
      )
      if (is.null(res) || httr::status_code(res) != 200) {
        return(data.frame(query = smi, chebiid = NA_character_, chebiasciiname = NA_character_,
                          stringsAsFactors = FALSE))
      }
      parsed <- tryCatch(jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"),
                                            simplifyVector = FALSE),
                         error = function(e) NULL)
      hits <- parsed$results
      if (is.null(hits) || length(hits) == 0) {
        return(data.frame(query = smi, chebiid = NA_character_, chebiasciiname = NA_character_,
                          stringsAsFactors = FALSE))
      }
      best <- hits[[1]]$`_source`
      data.frame(
        query          = smi,
        chebiid        = best$chebi_accession,
        chebiasciiname = best$ascii_name,
        stringsAsFactors = FALSE
      )
    })
    close(pb)
    do.call(rbind, results)
  }

  chebi_ids <- get_chebi_from_smiles(for_chembiid_smiles)
  message(sprintf("ChEBI lookup complete: %d/%d SMILES matched.", sum(!is.na(chebi_ids$chebiid)), nrow(chebi_ids)))


  # And we merge the data_sirius dataframe with the chebi_ids dataframe
  data_sirius <- merge(data_sirius, chebi_ids, by.x = "smiles", by.y = "query")

  # The column names are modified to include the source of the data

  colnames(data_sirius) <- paste("sirius", colnames(data_sirius),  sep = "_")


  # We now build a unique feature_id for each feature in the Sirius data

  # data_sirius$feature_id <- sub("^.*_([[:alnum:]]+)$", "\\1", data_sirius$sirius_id)
  # Previous line is now deprecated with the new Sirius outputs

  if (sirius_version == "6") {
    data_sirius$feature_id <- as.numeric(data_sirius$sirius_mappingFeatureId)
  } else {
    data_sirius$feature_id <- as.numeric(data_sirius$sirius_featureId)
  }

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

#data_canopus$feature_id <- sub("^.*_([[:alnum:]]+)$", "\\1", data_canopus$canopus_id)
# Previous line is now deprecated with the new Sirius outputs

if (sirius_version == "6") {
  data_canopus$feature_id <- as.numeric(data_canopus$canopus_mappingFeatureId)
} else {
  data_canopus$feature_id <- as.numeric(data_canopus$canopus_featureId)
}


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
  # data_gnps <- merge(data_gnps_mn, data_gnps_lib, by = "#Scan#")
  list_df <- list(data_gnps_mn, data_gnps_lib)
  data_gnps <- list_df %>% reduce(full_join, by = "#Scan#")


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

# We sanitize VM column names using the same snake_case convention used elsewhere

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


sample_metadata <- read_delim(file.path(working_directory, "metadata", "treated", paste(params$mapp_batch, "metadata.tsv", sep = "_")),
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

if (!"sample_id" %in% colnames(sample_metadata) && "samp.mzMLle_id" %in% colnames(sample_metadata)) {
  warning("Renaming legacy metadata column 'samp.mzMLle_id' to 'sample_id'.")
  sample_metadata <- sample_metadata %>% rename(sample_id = `samp.mzMLle_id`)
}

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
  col_name <- make_clean_names(col_name, case = "snake")

  # Initialize aggregated groups with original condition variable
  SM[paste(col_name, "simplified", sep = "_")] <- as.factor(SM[[col_name]])

  # Iterate over each group in params$tocomb
  for (group in names(col_info$groups)) {
    group_info <- col_info$groups[[group]]
    levels <- group_info$levels

    # now make sure to sort the levels
    levels <- sort(levels, decreasing = FALSE)
    # Use an explicit group label when provided; otherwise retain the historical joined-level label.
    new_label <- if (is.null(group_info$label)) paste(levels, collapse = "_") else group_info$label


    # Combine levels for the current group
    SM[paste(col_name, "simplified", sep = "_")] <- combineLevels(SM[[paste(col_name, "simplified", sep = "_")]], levs = levels, newLabel = c(new_label))
  }
}

# The function below is used to create metadata combinations
# Then we work horizontally (across SM columns)


if (!is.null(params$to_combine_horizontally$factor_name)) {
  df <- SM %>%
    filter(sample_type == "sample")

  # This line allows us to make sure that the columns will be combined in alphabetical order
  cols <- sort(c(params$to_combine_horizontally$factor_name), decreasing = FALSE)

  # Here we normalize the column names present in cols to snake_case

  cols <- make_clean_names(cols, case = "snake")


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


# Ponderation stage.

# First we check from the params that the apply_ponderation is set to TRUE and that the factor_name is not NULL

if (params$actions$ponderate_data$run == "TRUE" && !is.null(params$actions$ponderate_data$factor_name)) {
  # Prepare the data: convert params$actions$ponderate_data$factor_name to numeric, handling non-numeric "ND" values
  # We print a message to the console to inform the user that the ponderation is being applied and the factor_name used
  print(paste("Ponderation is being applied using the factor:", params$actions$ponderate_data$factor_name))
  
  factor_name <- make_clean_names(params$actions$ponderate_data$factor_name, case = "snake")
   
  X <- X %>%
    as.data.frame() %>%
    rownames_to_column("SampleID") %>%  # Temporarily create a SampleID column from rownames
    left_join(SM %>% select(filename, all_of(factor_name)) %>% 
                rename(SampleID = filename) %>% 
                mutate(across(all_of(factor_name), as.numeric)), by = "SampleID") %>%  # Join on SampleID and convert the factor_name to numeric
    mutate(across(-c(SampleID, all_of(factor_name)), ~ ifelse(is.na(.data[[factor_name]]), .x, .x / .data[[factor_name]]))) %>%  # Perform row-wise division, keep original if the factor is NA
    select(-all_of(factor_name)) %>%  # Drop factor_name column after ponderation
    column_to_rownames("SampleID")  # Convert SampleID back to rownames
} else {
  # If ponderation is not required, we keep the original data
  X <- X
}

# Pruning stage.

# Check if the pruning threshold is set in the params
if (!is.null(params$actions$prune_data$threshold)) {
  threshold <- params$actions$prune_data$threshold  # Define the threshold value from params

  threshold <- as.numeric(threshold)
  
  # Print a message to the console to inform the user that pruning is being applied
  print(paste("Pruning is being applied using the threshold:", threshold))
  
  # Prune the X dataframe by dropping columns where the maximum value doesn't reach the threshold
  pruned_columns <- X %>%
    as.data.frame() %>%
    select(where(~ max(.x, na.rm = TRUE) >= threshold)) %>%
    colnames()
  
  # Update X to keep only the pruned columns
  X <- X[, pruned_columns, drop = FALSE]
  
  # Prune the VM dataframe by keeping only the rows corresponding to the pruned columns of X
  VM <- VM[rownames(VM) %in% pruned_columns, , drop = FALSE]
  
} else {
  # If pruning is not required, we keep the original data
  X <- X
  VM <- VM
}

# Min value imputation (to be checked !!!)

half_min <- min(X[X > 0], na.rm = TRUE) / 2
min <- min(X[X > 0], na.rm = TRUE)


X[X == 0] <- min


missing_variable_metadata <- setdiff(colnames(X), row.names(VM))
if (length(missing_variable_metadata)) {
  stop(sprintf(
    "Some columns in X are not present in the rownames of VM: %s",
    paste(missing_variable_metadata, collapse = ", ")
  ))
}
VM <- VM[colnames(X), , drop = FALSE]

# We repeat for row.names(SMDF) == row.names(X_pond)

missing_sample_metadata <- setdiff(row.names(X), row.names(SM))
if (length(missing_sample_metadata)) {
  stop(sprintf(
    "Some rownames in X are not present in the rownames of SM: %s",
    paste(missing_sample_metadata, collapse = ", ")
  ))
}

metadata_only_samples <- setdiff(row.names(SM), row.names(X))
if (length(metadata_only_samples)) {
  message(sprintf(
    "Dropping %d metadata-only sample(s) absent from the quantification matrix: %s",
    length(metadata_only_samples),
    paste(metadata_only_samples, collapse = ", ")
  ))
}
SM <- SM[row.names(X), , drop = FALSE]

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
    factor_name = make_clean_names(params$filter_blank$factor_name, case = "snake"),
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
    factor_name = make_clean_names(params$filter_sample_type$factor_name, case = "snake"),
    levels = params$filter_sample_type$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
} 

if (params$filter_sample_metadata_one$mode %in% possible_modes){
  filter_smeta_model <- filter_smeta(
    mode = params$filter_sample_metadata_one$mode,
    factor_name = make_clean_names(params$filter_sample_metadata_one$factor_name, case = "snake"),
    levels = params$filter_sample_metadata_one$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
}

if (params$filter_sample_metadata_two$mode %in% possible_modes) {
  filter_smeta_model <- filter_smeta(
    mode = params$filter_sample_metadata_two$mode,
    factor_name = make_clean_names(params$filter_sample_metadata_two$factor_name, case = "snake"),
    levels = params$filter_sample_metadata_two$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
}

if (params$filter_sample_metadata_three$mode %in% possible_modes) {
  filter_smeta_model <- filter_smeta(
    mode = params$filter_sample_metadata_three$mode,
    factor_name = make_clean_names(params$filter_sample_metadata_three$factor_name, case = "snake"),
    levels = params$filter_sample_metadata_three$levels
  )

  # apply model sequence
  filter_smeta_result <- model_apply(filter_smeta_model, DE_filtered)

  DE_filtered <- filter_smeta_result@filtered
}


if (params$filter_variable_metadata_one$mode %in% possible_modes) {
  filter_vmeta_model <- filter_vmeta(
    mode = params$filter_variable_metadata_one$mode,
    factor_name = make_clean_names(params$filter_variable_metadata_one$factor_name, case = "snake"),
    levels = params$filter_variable_metadata_one$levels
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}

if (params$filter_variable_metadata_two$mode %in% possible_modes) {
  filter_vmeta_model <- filter_vmeta(
    mode = params$filter_variable_metadata_two$mode,
    factor_name = make_clean_names(params$filter_variable_metadata_two$factor_name, case = "snake"),
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
    factor_name = make_clean_names(params$filter_variable_metadata_annotated$factor_name, case = "snake"),
    levels = as.character(params$filter_variable_metadata_annotated$levels)
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}


if (params$filter_variable_metadata_num$mode %in% possible_modes) {
  filter_vmeta_model <- filter_vmeta_num(
    mode = params$filter_variable_metadata_num$mode,
    factor_name = make_clean_names(params$filter_variable_metadata_num$factor_name, case = "snake"),
    level = params$filter_variable_metadata_num$level
  )

  # apply model sequence
  filter_vmeta_result <- model_apply(filter_vmeta_model, DE_filtered)

  DE_filtered <- filter_vmeta_result@filtered
}


if (params$actions$scale_method == "none") {
  DE <- DE_filtered

} else if (params$actions$scale_method == "pareto") {
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
  DE$data <- apply(DE$data, 2, range01)


  # Min value imputation also after the scaling stage (to be checked !!!)

  # half_min_sec = min(DE$data[DE$data > 0], na.rm = TRUE) / 2

  # min_sec = min(DE$data[DE$data > 0], na.rm = TRUE)

  # DE$data[DE$data == 0] = min_sec


} else if (params$actions$scale_method == "autoscale") 
{ 
  # Overall Pareto scaling (test)

  M <- autoscale()
  M <- model_train(M, DE_filtered)
  M <- model_predict(M, DE_filtered)
  DE <- M$autoscaled

  # We use the filter_na_count function to filter out features with a number of NAs greater than the threshold

  M <- filter_na_count(threshold = 1, factor_name = "sample_type")
  M <- model_apply(M, DE)

  DE <- M$filtered

  ##### we range all feature from 0 to 1

  # @Manu !!! Why do we have this ?!
  DE$data <- apply(DE$data, 2, range01)

} else {
stop("Please check the value of the 'scale_method' parameter in the params file.")
}


# We make sure that params$target$sample_metadata_header matches cleaned snake_case metadata headers

params$target$sample_metadata_header <- make_clean_names(params$target$sample_metadata_header, case = "snake")

# Here we check if the params$paths$out value exist and use it else we use the default output_directory

target_name = paste(as.vector(sort(as.character(unique(DE$sample_meta[[params$target$sample_metadata_header]])), decreasing = FALSE)), collapse = "_vs_")


################################################################################################
################################################################################################
######################## structool box formatted data export

message("Outputting X, VM and SM ...")

formatted_peak_table <- DE$data

formatted_variable_metadata <- DE$variable_meta ### need to be filter with only usefull output


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


title_PCA <- paste(
  paste("PCA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

title_PLSDA <- paste(
  paste("PLSDA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

title_PLSDA_VIP <- paste("PLSDA selected Features of Importance", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")

title_DFA <- paste(
  paste("DFA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)


title_PCA3D <- paste(
  paste("PCA3D", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)


title_PCoA <- paste(
  paste("PCoA", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)


title_PCoA3D <- paste(
  paste("PCoA3D", "for dataset", params$mapp_batch),
  paste("Comparison across:", params$target$sample_metadata_header, target_name),
  paste("Filter Sample Metadata Status:", filter_sample_metadata_status),
  paste("Filter Variable Metadata Status:", filter_variable_metadata_status),
  paste("Scaling Status:", scaling_status),
  sep = "\n"
)

title_treemap <- paste("Treemap", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")
title_random_forest <- paste("Random Forest results", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")
title_box_plots <- paste("Top", params$boxplot$topN, "boxplots", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")
title_heatmap_rf <- paste("Heatmap of", "top", params$heatmap$topN, "Random Forest filtered features", "for dataset", params$target$sample_metadata_header, target_name, filter_variable_metadata_status, scaling_status, sep = " ")


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

# Sample and sort unique color values from wes_palettes
wes_palettes_vec <- sample(sort(unique(unlist(wes_palettes[names(wes_palettes)]))))

# Extract unique factor names from metadata
factor_name_meta <- unlist(unique(DE$sample_meta[params$target$sample_metadata_header]))

# Check that all members of params$colors$all$key are present in factor_name_meta.
# If not, return the values of params$colors$all$key that are not present in factor_name_meta.
if (!all(params$colors$all$key %in% factor_name_meta)) {
  missing_colors <- params$colors$all$key[!params$colors$all$key %in% factor_name_meta]
  factor_name_meta_str <- paste(unique(factor_name_meta), collapse=", ")
  stop(paste("The following values in params$colors$all$key are not present in the sample metadata:", 
             paste(missing_colors, collapse=", "), 
             "Check the spelling of values in params$colors$all$key, they should match the following available values:", 
             factor_name_meta_str))
}

# We establish a named vector for the whole dataset
if (params$colors$continuous) {
  # Apply numerical sorting to the keys if continuous color scale is requested
  sorted_keys <- as.character(sort(as.numeric(params$colors$all$key)))
  
  # Apply the Viridis color scale
  viridis_colors <- viridis(length(sorted_keys))
  
  # Assign the Viridis colors to the sorted keys
  custom_colors <- setNames(viridis_colors, sorted_keys)
} else {
  if (length(params$colors$all$key) > 0) {
    custom_colors <- setNames(c(params$colors$all$value), c(params$colors$all$key))
  } else {
    custom_colors <- wes_palettes_vec[sample(c(1:length(wes_palettes_vec)), length(factor_name_meta))]
    names(custom_colors) <- factor_name_meta
  }
}




#################################################################################################
#################################################################################################
#################################################################################################
##### PCA filtered data #######################################################################

message("Launching PCA calculations ...")

get_param_scalar <- function(value, default) {
  if (is.null(value) || !length(value) || is.na(value[1])) {
    return(default)
  }
  value[1]
}

normalize_param_vector <- function(value) {
  if (is.null(value) || !length(value)) {
    return(character(0))
  }
  value <- as.character(unlist(value, use.names = FALSE))
  value <- trimws(value)
  value[!is.na(value) & nzchar(value)]
}

ordination_point_size <- as.numeric(get_param_scalar(params$ordination$point_size, get_param_scalar(params$pca$point_size, 3)))
ordination_axis_text_size <- as.numeric(get_param_scalar(params$ordination$axis_text_size, get_param_scalar(params$pca$axis_text_size, 13)))
ordination_axis_title_size <- as.numeric(get_param_scalar(params$ordination$axis_title_size, get_param_scalar(params$pca$axis_title_size, 15)))
ordination_legend_text_size <- as.numeric(get_param_scalar(params$ordination$legend_text_size, get_param_scalar(params$pca$legend_text_size, 13)))
ordination_legend_title_size <- as.numeric(get_param_scalar(params$ordination$legend_title_size, get_param_scalar(params$pca$legend_title_size, 14)))
ordination_title_size <- as.numeric(get_param_scalar(params$ordination$title_size, get_param_scalar(params$pca$title_size, 14)))
ordination_export_width <- as.numeric(get_param_scalar(params$ordination$export_width, get_param_scalar(params$pca$export_width, 11)))
ordination_export_height <- as.numeric(get_param_scalar(params$ordination$export_height, get_param_scalar(params$pca$export_height, 9)))
ordination_point_alpha <- as.numeric(get_param_scalar(params$ordination$point_alpha, 0.9))
ordination_points_to_label <- as.character(get_param_scalar(params$ordination$points_to_label, get_param_scalar(params$pca$points_to_label, "none")))
if (!ordination_points_to_label %in% c("none", "all", "outliers")) {
  stop("params$ordination$points_to_label must be one of: none, all, outliers")
}
ordination_label_size <- as.numeric(get_param_scalar(params$ordination$label_size, 3.88))

publication_ordination_theme <- function() {
  theme_classic() +
    theme(
      plot.title = element_text(size = ordination_title_size, face = "bold", hjust = 0.5, lineheight = 1.05),
      axis.title = element_text(size = ordination_axis_title_size, face = "bold"),
      axis.text = element_text(size = ordination_axis_text_size, colour = "black"),
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks.length = grid::unit(0.22, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = ordination_axis_title_size, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = ordination_legend_title_size, face = "bold"),
      legend.text = element_text(size = ordination_legend_text_size),
      legend.key.size = grid::unit(0.7, "cm"),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 16, 12, 12)
    )
}

apply_ordination_point_style <- function(plot_obj) {
  if (length(plot_obj$layers) > 0) {
    for (layer_index in seq_along(plot_obj$layers)) {
      if (inherits(plot_obj$layers[[layer_index]]$geom, "GeomPoint")) {
        plot_obj$layers[[layer_index]]$aes_params$size <- ordination_point_size
        plot_obj$layers[[layer_index]]$aes_params$alpha <- ordination_point_alpha
        plot_obj$layers[[layer_index]]$aes_params$shape <- 16
        plot_obj$layers[[layer_index]]$aes_params$stroke <- 0
      }
    }
  }
  plot_obj
}

ordination_colour_scale <- function() {
  list(
    scale_colour_manual(name = "Groups", values = custom_colors),
    guides(colour = guide_legend(override.aes = list(size = ordination_point_size + 1, alpha = 1, shape = 16, stroke = 0)))
  )
}

build_manual_plsda_scores_plot <- function(plsda_object) {
  plsda_scores_df <- data.frame(plsda_object$scores$data, check.names = FALSE)
  plsda_scores_df$sample_row <- rownames(plsda_scores_df)
  plsda_meta_df <- data.frame(plsda_object$scores$sample_meta, check.names = FALSE)
  plsda_meta_df$sample_row <- rownames(plsda_meta_df)
  plsda_plot_df <- merge(plsda_scores_df, plsda_meta_df, by = "sample_row", all.x = TRUE)
  plsda_numeric_columns <- colnames(plsda_scores_df)[vapply(plsda_scores_df, is.numeric, logical(1))]
  plsda_axis_columns <- intersect(c("LV1", "LV2"), plsda_numeric_columns)
  if (length(plsda_axis_columns) < 2) {
    plsda_axis_columns <- head(plsda_numeric_columns, 2)
  }
  if (length(plsda_axis_columns) < 2) {
    stop("PLSDA scores plot cannot be built because fewer than two numeric score columns were found.")
  }
  plsda_group_column <- params$target$sample_metadata_header
  if (!plsda_group_column %in% colnames(plsda_plot_df)) {
    stop(sprintf("PLSDA scores plot cannot be built because '%s' is missing from PLSDA sample metadata.", plsda_group_column))
  }
  plsda_label_column <- if ("sample_id" %in% colnames(plsda_plot_df)) "sample_id" else "sample_row"
  plsda_plot_df[[plsda_group_column]] <- factor(plsda_plot_df[[plsda_group_column]])

  plot_obj <- ggplot(
    plsda_plot_df,
    aes(x = .data[[plsda_axis_columns[1]]], y = .data[[plsda_axis_columns[2]]], colour = .data[[plsda_group_column]])
  ) +
    geom_point() +
    labs(
      x = plsda_axis_columns[1],
      y = plsda_axis_columns[2],
      colour = "Groups"
    )

  if (ordination_points_to_label != "none") {
    plsda_label_data <- plsda_plot_df
    if (ordination_points_to_label == "outliers") {
      plsda_distance_to_center <- sqrt(
        (plsda_plot_df[[plsda_axis_columns[1]]] - mean(plsda_plot_df[[plsda_axis_columns[1]]], na.rm = TRUE))^2 +
          (plsda_plot_df[[plsda_axis_columns[2]]] - mean(plsda_plot_df[[plsda_axis_columns[2]]], na.rm = TRUE))^2
      )
      plsda_outlier_threshold <- stats::quantile(plsda_distance_to_center, probs = 0.9, na.rm = TRUE)
      plsda_label_data <- plsda_plot_df[plsda_distance_to_center >= plsda_outlier_threshold, , drop = FALSE]
    }
    plot_obj <- plot_obj +
      ggrepel::geom_text_repel(
        data = plsda_label_data,
        aes(label = .data[[plsda_label_column]]),
        size = ordination_label_size,
        show.legend = FALSE,
        max.overlaps = Inf
      )
  }

  plot_obj
}

#################################################################################################
#################################################################################################
##### NPC summed-intensity plots ################################################################

npc_summed_intensity_params <- params$npc_summed_intensity
if (is.null(npc_summed_intensity_params)) {
  npc_summed_intensity_params <- list()
}
npc_plot_pathway <- normalize_param_vector(npc_summed_intensity_params$pathway)
npc_plot_superclass <- normalize_param_vector(npc_summed_intensity_params$superclass)
npc_plot_class <- normalize_param_vector(npc_summed_intensity_params$class)
npc_expand_all <- as.logical(get_param_scalar(npc_summed_intensity_params$expand_all, FALSE))
npc_expand_pathway <- normalize_param_vector(npc_summed_intensity_params$expand_pathway)
npc_expand_levels <- normalize_param_vector(npc_summed_intensity_params$expand_levels)
if (!length(npc_expand_levels)) {
  npc_expand_levels <- "class"
}
if ("all" %in% npc_expand_levels) {
  npc_expand_levels <- c("superclass", "class")
}
npc_expand_levels <- intersect(npc_expand_levels, c("superclass", "class"))
npc_plot_terms <- data.frame(
  npc_level = c(
    rep("pathway", length(npc_plot_pathway)),
    rep("superclass", length(npc_plot_superclass)),
    rep("class", length(npc_plot_class))
  ),
  npc_term = c(npc_plot_pathway, npc_plot_superclass, npc_plot_class),
  stringsAsFactors = FALSE
)

if (nrow(npc_plot_terms) > 0 || length(npc_expand_pathway) > 0 || isTRUE(npc_expand_all)) {
  message("Preparing NPC summed-intensity plots ...")

  npc_min_probability <- as.numeric(get_param_scalar(npc_summed_intensity_params$min_probability, 0))
  npc_transform <- as.character(get_param_scalar(npc_summed_intensity_params$transform, "log10"))
  npc_static_export_top_n <- as.integer(get_param_scalar(npc_summed_intensity_params$static_export_top_n, 40))
  npc_ratio_params <- npc_summed_intensity_params$ratios
  if (is.null(npc_ratio_params)) {
    npc_ratio_params <- list()
  }
  npc_ratio_enabled <- as.logical(get_param_scalar(npc_ratio_params$enabled, TRUE))
  npc_ratio_denominator_level <- as.character(get_param_scalar(npc_ratio_params$denominator_level, "pathway"))
  npc_ratio_pseudocount <- as.numeric(get_param_scalar(npc_ratio_params$pseudocount, 0))
  npc_export_raw <- as.logical(get_param_scalar(npc_summed_intensity_params$raw_export, TRUE))
  if (!npc_transform %in% c("log10", "none")) {
    stop("params$npc_summed_intensity$transform must be one of: log10, none")
  }
  if (npc_ratio_denominator_level != "pathway") {
    stop("params$npc_summed_intensity$ratios$denominator_level currently supports only: pathway")
  }

  npc_level_columns <- c(
    pathway = "canopus_npc_pathway",
    superclass = "canopus_npc_superclass",
    class = "canopus_npc_class"
  )
  npc_probability_columns <- c(
    pathway = "canopus_npc_pathway_probability",
    superclass = "canopus_npc_superclass_probability",
    class = "canopus_npc_class_probability"
  )

  npc_treat_npclassifier_json <- function(taxonomy) {
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

    taxonomy_hierarchy_by_class <- cbind(id_pathway, id_superclass, id_class) %>%
      data.frame() %>%
      mutate(id_class = as.numeric(id_class)) %>%
      unnest(id_superclass) %>%
      unnest(id_pathway)

    taxonomy_hierarchy_superclass <- taxonomy$Super_hierarchy
    id_pathway_2 <- list()
    id_superclass <- list()

    for (i in seq_len(length(taxonomy_hierarchy_superclass))) {
      id_pathway_2[[i]] <- taxonomy_hierarchy_superclass[[i]]$Pathway
      id_superclass[[i]] <- names(taxonomy_hierarchy_superclass[i])
    }

    taxonomy_hierarchy_by_superclass <- cbind(id_pathway_2, id_superclass) %>%
      data.frame() %>%
      mutate(id_superclass = as.numeric(id_superclass)) %>%
      unnest(id_pathway_2)

    full_join(taxonomy_hierarchy_by_class, taxonomy_classes) %>%
      full_join(., taxonomy_superclasses) %>%
      full_join(., taxonomy_pathways) %>%
      distinct(class, superclass, pathway)
  }

  npc_taxonomy <- tryCatch(
    {
      npc_taxonomy_url <- "https://raw.githubusercontent.com/mwang87/NP-Classifier/master/Classifier/dict/index_v1.json"
      npc_treat_npclassifier_json(jsonlite::fromJSON(npc_taxonomy_url))
    },
    error = function(e) {
      warning(sprintf("Could not load NP-Classifier taxonomy; NPC ratio plots will be skipped. Error: %s", e$message))
      NULL
    }
  )

  npc_safe_file_part <- function(value) {
    value <- tolower(as.character(value))
    value <- gsub("[^a-z0-9]+", "_", value)
    gsub("^_|_$", "", value)
  }

  npc_terms_to_pathway <- function(npc_level, npc_term) {
    if (is.null(npc_taxonomy)) {
      return(character(0))
    }
    if (npc_level == "class") {
      canonical_paths <- npc_taxonomy$pathway[npc_taxonomy$class == npc_term]
    } else if (npc_level == "superclass") {
      canonical_paths <- npc_taxonomy$pathway[npc_taxonomy$superclass == npc_term]
    } else {
      canonical_paths <- npc_term
    }
    canonical_paths <- unique(as.character(canonical_paths[!is.na(canonical_paths)]))
    canonical_paths[nzchar(canonical_paths)]
  }

  npc_normalize_requested_pathway <- function(pathway) {
    if (pathway %in% c("Lipids", "Lipid", "Lipids and lipid-like molecules")) {
      return("Fatty acids")
    }
    pathway
  }

  npc_add_taxonomy_pathway <- function(variable_meta) {
    variable_meta$npc_taxonomy_pathway <- if ("canopus_npc_pathway" %in% colnames(variable_meta)) {
      as.character(variable_meta$canopus_npc_pathway)
    } else {
      rep(NA_character_, nrow(variable_meta))
    }
    if (is.null(npc_taxonomy)) {
      return(variable_meta)
    }
    if ("canopus_npc_class" %in% colnames(variable_meta)) {
      class_lookup <- npc_taxonomy %>%
        filter(!is.na(class), !is.na(pathway)) %>%
        distinct(class, pathway) %>%
        group_by(class) %>%
        summarise(pathway = paste(unique(pathway), collapse = " x "), .groups = "drop")
      class_match <- match(variable_meta$canopus_npc_class, class_lookup$class)
      class_has_match <- !is.na(class_match)
      variable_meta$npc_taxonomy_pathway[class_has_match] <- class_lookup$pathway[class_match[class_has_match]]
    }
    if ("canopus_npc_superclass" %in% colnames(variable_meta)) {
      superclass_lookup <- npc_taxonomy %>%
        filter(!is.na(superclass), !is.na(pathway)) %>%
        distinct(superclass, pathway) %>%
        group_by(superclass) %>%
        summarise(pathway = paste(unique(pathway), collapse = " x "), .groups = "drop")
      missing_taxonomy_path <- is.na(variable_meta$npc_taxonomy_pathway) | !nzchar(variable_meta$npc_taxonomy_pathway)
      superclass_match <- match(variable_meta$canopus_npc_superclass, superclass_lookup$superclass)
      superclass_has_match <- missing_taxonomy_path & !is.na(superclass_match)
      variable_meta$npc_taxonomy_pathway[superclass_has_match] <- superclass_lookup$pathway[superclass_match[superclass_has_match]]
    }
    variable_meta
  }

  npc_pathway_contains <- function(values, pathway) {
    pathway_pattern <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", pathway)
    grepl(paste0("(^| x )", pathway_pattern, "($| x )"), values)
  }

  npc_write_table <- function(table_df, filename) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    write.table(table_df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  npc_feature_driver_table <- function(npc_data_matrix, npc_sample_meta, npc_variable_meta, feature_ids, npc_level, npc_term, npc_label, npc_source_label) {
    feature_ids <- intersect(as.character(feature_ids), colnames(npc_data_matrix))
    if (!length(feature_ids)) {
      return(data.frame())
    }
    group_values <- names(custom_colors)[names(custom_colors) %in% as.character(unique(npc_sample_meta[[params$target$sample_metadata_header]]))]
    if (!length(group_values)) {
      group_values <- sort(unique(as.character(npc_sample_meta[[params$target$sample_metadata_header]])))
    }
    group_vector <- factor(as.character(npc_sample_meta[[params$target$sample_metadata_header]]), levels = group_values)
    feature_rows <- lapply(feature_ids, function(feature_id) {
      values <- suppressWarnings(as.numeric(npc_data_matrix[, feature_id]))
      means <- tapply(values, group_vector, mean, na.rm = TRUE)
      medians <- tapply(values, group_vector, stats::median, na.rm = TRUE)
      means <- means[group_values]
      medians <- medians[group_values]
      valid_means <- means[!is.na(means)]
      if (length(valid_means)) {
        higher_group <- names(valid_means)[which.max(valid_means)]
        lower_group <- names(valid_means)[which.min(valid_means)]
        mean_difference <- unname(max(valid_means) - min(valid_means))
      } else {
        higher_group <- NA_character_
        lower_group <- NA_character_
        mean_difference <- NA_real_
      }
      row <- data.frame(
        data_source = npc_source_label,
        npc_level = npc_level,
        npc_term = npc_term,
        npc_label = npc_label,
        feature_id = feature_id,
        higher_mean_group = higher_group,
        lower_mean_group = lower_group,
        mean_difference_max_min = mean_difference,
        abs_mean_difference_max_min = abs(mean_difference),
        average_intensity = mean(values, na.rm = TRUE),
        median_intensity = stats::median(values, na.rm = TRUE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      for (group_name in group_values) {
        safe_group <- npc_safe_file_part(group_name)
        row[[paste0("mean_", safe_group)]] <- unname(means[group_name])
        row[[paste0("median_", safe_group)]] <- unname(medians[group_name])
      }
      row
    })
    driver_df <- bind_rows(feature_rows)
    driver_df <- driver_df[order(-driver_df$abs_mean_difference_max_min, -driver_df$average_intensity, driver_df$feature_id, na.last = TRUE), , drop = FALSE]
    total_abs_difference <- sum(driver_df$abs_mean_difference_max_min, na.rm = TRUE)
    total_average_intensity <- sum(driver_df$average_intensity, na.rm = TRUE)
    driver_df$contribution_fraction_of_term_difference <- if (total_abs_difference > 0) {
      driver_df$abs_mean_difference_max_min / total_abs_difference
    } else {
      NA_real_
    }
    driver_df$cumulative_contribution_fraction <- cumsum(ifelse(is.na(driver_df$contribution_fraction_of_term_difference), 0, driver_df$contribution_fraction_of_term_difference))
    driver_df$average_intensity_fraction_of_term <- if (total_average_intensity > 0) {
      driver_df$average_intensity / total_average_intensity
    } else {
      NA_real_
    }
    annotation_columns <- intersect(
      c(
        "feature_id",
        "feature_id_full_annotated",
        "sirius_chebiasciiname",
        "sirius_name",
        "sirius_adduct",
        "gnps_component",
        "gnps_componentindex",
        "canopus_npc_pathway",
        "canopus_npc_superclass",
        "canopus_npc_class",
        "canopus_npc_pathway_probability",
        "canopus_npc_superclass_probability",
        "canopus_npc_class_probability"
      ),
      colnames(npc_variable_meta)
    )
    if (length(annotation_columns)) {
      annotation_df <- npc_variable_meta[match(driver_df$feature_id, as.character(npc_variable_meta$feature_id)), annotation_columns, drop = FALSE]
      annotation_df$feature_id <- as.character(annotation_df$feature_id)
      driver_df <- left_join(driver_df, annotation_df, by = "feature_id")
    }
    driver_df
  }

  npc_format_p_value <- function(p_value) {
    if (is.na(p_value)) {
      return("NA")
    }
    if (p_value < 0.001) {
      return(formatC(p_value, format = "e", digits = 1))
    }
    formatC(p_value, format = "f", digits = 3)
  }

  npc_significance_label <- function(q_value) {
    if (is.na(q_value)) {
      return("ns")
    }
    if (q_value < 0.001) {
      return("***")
    }
    if (q_value < 0.01) {
      return("**")
    }
    if (q_value < 0.05) {
      return("*")
    }
    "ns"
  }

  npc_compute_stats <- function(plot_df, label_column, full_label_column, y_column, data_source, value_name) {
    labels <- unique(as.character(plot_df[[label_column]]))
    stats_rows <- list()
    for (label in labels) {
      label_df <- plot_df[as.character(plot_df[[label_column]]) == label, , drop = FALSE]
      label_df <- label_df[!is.na(label_df[[y_column]]) & !is.na(label_df$group), , drop = FALSE]
      groups <- levels(label_df$group)
      groups <- groups[groups %in% as.character(unique(label_df$group))]
      if (!length(groups)) {
        groups <- sort(unique(as.character(label_df$group)))
      }
      group_count <- length(groups)
      if (group_count < 2) {
        next
      }
      full_label <- as.character(label_df[[full_label_column]][1])
      group_sizes <- table(as.character(label_df$group))
      group_medians <- tapply(label_df[[y_column]], as.character(label_df$group), stats::median, na.rm = TRUE)
      group_means <- tapply(label_df[[y_column]], as.character(label_df$group), mean, na.rm = TRUE)

      overall_p <- NA_real_
      overall_test <- NA_character_
      if (group_count == 2) {
        overall_test <- "wilcoxon_rank_sum"
        overall_p <- tryCatch(
          stats::wilcox.test(label_df[[y_column]] ~ label_df$group, exact = FALSE)$p.value,
          error = function(err) NA_real_
        )
      } else {
        overall_test <- "kruskal_wallis"
        overall_p <- tryCatch(
          stats::kruskal.test(label_df[[y_column]] ~ label_df$group)$p.value,
          error = function(err) NA_real_
        )
      }
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        data_source = data_source,
        value = value_name,
        npc_label = full_label,
        plot_label = label,
        test = overall_test,
        contrast = "overall",
        group_1 = NA_character_,
        group_2 = NA_character_,
        n_group_1 = NA_integer_,
        n_group_2 = NA_integer_,
        mean_group_1 = NA_real_,
        mean_group_2 = NA_real_,
        mean_difference = NA_real_,
        abs_mean_difference = if (length(group_means) > 1) max(group_means, na.rm = TRUE) - min(group_means, na.rm = TRUE) else NA_real_,
        median_group_1 = NA_real_,
        median_group_2 = NA_real_,
        median_difference = NA_real_,
        abs_median_difference = if (length(group_medians) > 1) max(group_medians, na.rm = TRUE) - min(group_medians, na.rm = TRUE) else NA_real_,
        p_value = overall_p,
        stringsAsFactors = FALSE
      )

      for (pair in combn(groups, 2, simplify = FALSE)) {
        pair_df <- label_df[as.character(label_df$group) %in% pair, , drop = FALSE]
        pair_df$pair_group <- factor(as.character(pair_df$group), levels = pair)
        pair_p <- tryCatch(
          stats::wilcox.test(pair_df[[y_column]] ~ pair_df$pair_group, exact = FALSE)$p.value,
          error = function(err) NA_real_
        )
        median_1 <- unname(group_medians[pair[1]])
        median_2 <- unname(group_medians[pair[2]])
        mean_1 <- unname(group_means[pair[1]])
        mean_2 <- unname(group_means[pair[2]])
        stats_rows[[length(stats_rows) + 1]] <- data.frame(
          data_source = data_source,
          value = value_name,
          npc_label = full_label,
          plot_label = label,
          test = "pairwise_wilcoxon_rank_sum",
          contrast = paste(pair, collapse = "_vs_"),
          group_1 = pair[1],
          group_2 = pair[2],
          n_group_1 = unname(group_sizes[pair[1]]),
          n_group_2 = unname(group_sizes[pair[2]]),
          mean_group_1 = mean_1,
          mean_group_2 = mean_2,
          mean_difference = mean_2 - mean_1,
          abs_mean_difference = abs(mean_2 - mean_1),
          median_group_1 = median_1,
          median_group_2 = median_2,
          median_difference = median_2 - median_1,
          abs_median_difference = abs(median_2 - median_1),
          p_value = pair_p,
          stringsAsFactors = FALSE
        )
      }
    }
    if (!length(stats_rows)) {
      return(data.frame())
    }
    stats_df <- bind_rows(stats_rows)
    stats_df$q_value <- ave(
      stats_df$p_value,
      stats_df$test,
      stats_df$contrast,
      FUN = function(values) stats::p.adjust(values, method = "BH")
    )
    stats_df$significance <- vapply(stats_df$q_value, npc_significance_label, character(1))
    stats_df
  }

  npc_attach_overall_stats <- function(plot_df, stats_df, label_column) {
    if (!nrow(stats_df)) {
      plot_df$npc_overall_test <- NA_character_
      plot_df$npc_overall_p_value <- NA_real_
      plot_df$npc_overall_q_value <- NA_real_
      plot_df$npc_overall_significance <- NA_character_
      plot_df$npc_max_mean_difference <- NA_real_
      plot_df$npc_max_median_difference <- NA_real_
      return(plot_df)
    }
    overall_stats <- stats_df[stats_df$contrast == "overall", c("plot_label", "test", "p_value", "q_value", "significance", "abs_mean_difference", "abs_median_difference"), drop = FALSE]
    colnames(overall_stats) <- c(label_column, "npc_overall_test", "npc_overall_p_value", "npc_overall_q_value", "npc_overall_significance", "npc_max_mean_difference", "npc_max_median_difference")
    left_join(plot_df, overall_stats, by = label_column)
  }


  npc_selected_static_labels <- function(stats_df) {
    if (!nrow(stats_df)) {
      return(character())
    }
    overall_stats <- stats_df[stats_df$contrast == "overall", , drop = FALSE]
    if (!nrow(overall_stats)) {
      return(character())
    }
    ordered_stats <- overall_stats[order(overall_stats$q_value, overall_stats$p_value, -overall_stats$abs_mean_difference, overall_stats$plot_label, na.last = TRUE), , drop = FALSE]
    head(ordered_stats$plot_label, npc_static_export_top_n)
  }

  npc_save_plot <- function(plot_obj, filename, width = ordination_export_width, height = ordination_export_height) {
    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    if (tolower(tools::file_ext(filename)) == "png" && requireNamespace("ragg", quietly = TRUE)) {
      ggsave(plot = plot_obj, filename = filename, width = width, height = height, units = "in", dpi = 300, device = ragg::agg_png)
    } else {
      ggsave(plot = plot_obj, filename = filename, width = width, height = height, units = "in", dpi = 300)
    }
  }

  npc_apply_plot_intensity <- function(plot_df) {
    plot_df$plot_intensity <- if (npc_transform == "log10") {
      log10(plot_df$summed_intensity + 1)
    } else {
      plot_df$summed_intensity
    }
    plot_df
  }

  npc_facet_ncol <- function(label_count) {
    if (label_count >= 16) {
      return(4)
    }
    if (label_count >= 8) {
      return(3)
    }
    if (label_count >= 4) {
      return(2)
    }
    label_count
  }

  npc_combined_plot_height <- function(label_count) {
    max(ordination_export_height, ceiling(label_count / npc_facet_ncol(label_count)) * 2.65)
  }

  npc_static_stats_text <- function(plot_df) {
    if (!all(c("npc_overall_p_value", "npc_overall_q_value") %in% colnames(plot_df))) {
      return(NULL)
    }
    p_value <- plot_df$npc_overall_p_value[1]
    q_value <- plot_df$npc_overall_q_value[1]
    effect_text <- if ("npc_max_mean_difference" %in% colnames(plot_df) && !is.na(plot_df$npc_max_mean_difference[1])) {
      paste0(" | max mean delta=", signif(plot_df$npc_max_mean_difference[1], 3))
    } else {
      ""
    }
    paste0("Overall p=", npc_format_p_value(p_value), " / FDR q=", npc_format_p_value(q_value), " (", npc_significance_label(q_value), ")", effect_text)
  }

  npc_add_static_label <- function(plot_df, label_column, static_label_column) {
    if (all(c("npc_overall_p_value", "npc_overall_q_value") %in% colnames(plot_df))) {
      plot_df[[static_label_column]] <- paste0(
        plot_df[[label_column]],
        "\np=", vapply(plot_df$npc_overall_p_value, npc_format_p_value, character(1)),
        " / q=", vapply(plot_df$npc_overall_q_value, npc_format_p_value, character(1))
      )
    } else {
      plot_df[[static_label_column]] <- plot_df[[label_column]]
    }
    plot_df
  }

  npc_intensity_plot <- function(plot_df, plot_title, plot_subtitle, facet = TRUE) {
    plot_df$npc_plot_label <- if ("npc_plot_label" %in% colnames(plot_df)) {
      plot_df$npc_plot_label
    } else {
      plot_df$npc_label
    }
    plot_df <- npc_add_static_label(plot_df, "npc_plot_label", "npc_static_label")
    plot_df$plot_tooltip <- paste0(
      "NPC ", plot_df$npc_level, ": ", plot_df$npc_term,
      "<br>Group: ", plot_df$group,
      "<br>Sample: ", plot_df$sample_id,
      "<br>Summed intensity: ", signif(plot_df$summed_intensity, 4),
      "<br>Features: ", plot_df$n_features,
      if ("npc_overall_p_value" %in% colnames(plot_df)) paste0("<br>Overall p: ", vapply(plot_df$npc_overall_p_value, npc_format_p_value, character(1))) else "",
      if ("npc_overall_q_value" %in% colnames(plot_df)) paste0("<br>Overall q: ", vapply(plot_df$npc_overall_q_value, npc_format_p_value, character(1))) else ""
    )
    plot_obj <- ggplot(plot_df, aes(x = group, y = plot_intensity, color = group, fill = group)) +
      geom_boxplot(aes(group = group), width = 0.58, alpha = 0.18, outlier.shape = NA, linewidth = 0.65) +
      geom_jitter(aes(text = plot_tooltip), width = 0.1, size = ordination_point_size * 0.58, alpha = ordination_point_alpha, shape = 16, stroke = 0) +
      stat_summary(aes(group = group, fill = group), fun = mean, geom = "point", shape = 23, size = 2.8, colour = "black", stroke = 0.45, show.legend = FALSE) +
      scale_colour_manual(name = "Groups", values = custom_colors) +
      scale_fill_manual(name = "Groups", values = custom_colors) +
      labs(
        title = plot_title,
        subtitle = plot_subtitle,
        x = params$target$sample_metadata_header,
        y = if (npc_transform == "log10") "log10 summed intensity + 1" else "Summed intensity",
        caption = if (!facet) npc_static_stats_text(plot_df) else NULL
      ) +
      publication_ordination_theme() +
      theme(
        axis.text.x = element_text(size = ordination_axis_text_size, colour = "black", angle = 30, hjust = 1),
        plot.caption = element_text(size = ordination_axis_text_size, colour = "black", hjust = 0),
        legend.position = "none"
      )
    if (facet) {
      plot_obj <- plot_obj + facet_wrap(~ npc_static_label, scales = "free_y", ncol = npc_facet_ncol(length(unique(plot_df$npc_plot_label))), labeller = label_wrap_gen(width = 24))
    }
    plot_obj
  }

  npc_ratio_plot <- function(plot_df, plot_title, plot_subtitle, facet = TRUE) {
    plot_df$ratio_plot_label <- if ("ratio_plot_label" %in% colnames(plot_df)) {
      plot_df$ratio_plot_label
    } else {
      plot_df$ratio_label
    }
    plot_df <- npc_add_static_label(plot_df, "ratio_plot_label", "ratio_static_label")
    plot_df$plot_tooltip <- paste0(
      "Class: ", plot_df$numerator_term,
      "<br>Pathway: ", plot_df$denominator_term,
      "<br>Group: ", plot_df$group,
      "<br>Sample: ", plot_df$sample_id,
      "<br>Ratio: ", signif(plot_df$ratio, 4),
      "<br>Numerator features: ", plot_df$numerator_n_features,
      "<br>Denominator features: ", plot_df$denominator_n_features,
      if ("npc_overall_p_value" %in% colnames(plot_df)) paste0("<br>Overall p: ", vapply(plot_df$npc_overall_p_value, npc_format_p_value, character(1))) else "",
      if ("npc_overall_q_value" %in% colnames(plot_df)) paste0("<br>Overall q: ", vapply(plot_df$npc_overall_q_value, npc_format_p_value, character(1))) else ""
    )
    plot_obj <- ggplot(plot_df, aes(x = group, y = ratio, color = group, fill = group)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.35) +
      geom_boxplot(aes(group = group), width = 0.58, alpha = 0.18, outlier.shape = NA, linewidth = 0.65) +
      geom_jitter(aes(text = plot_tooltip), width = 0.1, size = ordination_point_size * 0.58, alpha = ordination_point_alpha, shape = 16, stroke = 0) +
      stat_summary(aes(group = group, fill = group), fun = mean, geom = "point", shape = 23, size = 2.8, colour = "black", stroke = 0.45, show.legend = FALSE) +
      scale_colour_manual(name = "Groups", values = custom_colors) +
      scale_fill_manual(name = "Groups", values = custom_colors) +
      labs(
        title = plot_title,
        subtitle = plot_subtitle,
        x = params$target$sample_metadata_header,
        y = "Fraction of corresponding NPC pathway summed intensity",
        caption = if (!facet) npc_static_stats_text(plot_df) else NULL
      ) +
      publication_ordination_theme() +
      theme(
        axis.text.x = element_text(size = ordination_axis_text_size, colour = "black", angle = 30, hjust = 1),
        plot.caption = element_text(size = ordination_axis_text_size, colour = "black", hjust = 0),
        legend.position = "none"
      )
    if (facet) {
      plot_obj <- plot_obj + facet_wrap(~ ratio_static_label, scales = "free_y", ncol = npc_facet_ncol(length(unique(plot_df$ratio_plot_label))), labeller = label_wrap_gen(width = 24))
    }
    plot_obj
  }

  npc_process_source <- function(npc_data_matrix, npc_sample_meta, npc_variable_meta, npc_source_label, filenames, file_suffix) {
    npc_variable_meta <- npc_add_taxonomy_pathway(npc_variable_meta)
    npc_summed_intensity_tables <- list()
    npc_feature_sets <- list()
    npc_source_plot_terms <- npc_plot_terms

    if (isTRUE(npc_expand_all)) {
      for (npc_expand_level in npc_expand_levels) {
        npc_expand_column <- npc_level_columns[[npc_expand_level]]
        if (!npc_expand_column %in% colnames(npc_variable_meta)) {
          next
        }
        npc_expand_probability_column <- npc_probability_columns[[npc_expand_level]]
        npc_expand_probability_values <- if (npc_expand_probability_column %in% colnames(npc_variable_meta)) {
          suppressWarnings(as.numeric(npc_variable_meta[[npc_expand_probability_column]]))
        } else {
          rep(1, nrow(npc_variable_meta))
        }
        npc_expand_keep <- !is.na(npc_expand_probability_values) & npc_expand_probability_values >= npc_min_probability
        npc_expanded_terms <- unique(as.character(npc_variable_meta[[npc_expand_column]][npc_expand_keep]))
        npc_expanded_terms <- npc_expanded_terms[!is.na(npc_expanded_terms) & nzchar(npc_expanded_terms)]
        if (length(npc_expanded_terms)) {
          npc_source_plot_terms <- bind_rows(
            npc_source_plot_terms,
            data.frame(npc_level = npc_expand_level, npc_term = sort(npc_expanded_terms), stringsAsFactors = FALSE)
          ) %>%
            distinct(npc_level, npc_term, .keep_all = TRUE)
        }
      }
    }

    if (length(npc_expand_pathway)) {
      for (npc_requested_pathway in npc_expand_pathway) {
        npc_pathway <- npc_normalize_requested_pathway(npc_requested_pathway)
        npc_pathway_keep <- !is.na(npc_variable_meta$npc_taxonomy_pathway) &
          npc_pathway_contains(npc_variable_meta$npc_taxonomy_pathway, npc_pathway)
        for (npc_expand_level in npc_expand_levels) {
          npc_expand_column <- npc_level_columns[[npc_expand_level]]
          if (!npc_expand_column %in% colnames(npc_variable_meta)) {
            next
          }
          npc_expand_probability_column <- npc_probability_columns[[npc_expand_level]]
          npc_expand_probability_values <- if (npc_expand_probability_column %in% colnames(npc_variable_meta)) {
            suppressWarnings(as.numeric(npc_variable_meta[[npc_expand_probability_column]]))
          } else {
            rep(1, nrow(npc_variable_meta))
          }
          npc_expand_keep <- npc_pathway_keep &
            !is.na(npc_expand_probability_values) &
            npc_expand_probability_values >= npc_min_probability
          npc_expanded_terms <- unique(as.character(npc_variable_meta[[npc_expand_column]][npc_expand_keep]))
          npc_expanded_terms <- npc_expanded_terms[!is.na(npc_expanded_terms) & nzchar(npc_expanded_terms)]
          if (length(npc_expanded_terms)) {
            npc_source_plot_terms <- bind_rows(
              npc_source_plot_terms,
              data.frame(npc_level = npc_expand_level, npc_term = sort(npc_expanded_terms), stringsAsFactors = FALSE)
            ) %>%
              distinct(npc_level, npc_term, .keep_all = TRUE)
          }
        }
      }
    }

    if (!nrow(npc_source_plot_terms)) {
      warning(sprintf("NPC summed-intensity plotting was requested for %s data, but no selected or expanded NPC terms were available.", npc_source_label))
      return(invisible(NULL))
    }

    for (npc_term_index in seq_len(nrow(npc_source_plot_terms))) {
      npc_level <- npc_source_plot_terms$npc_level[npc_term_index]
      npc_term <- npc_source_plot_terms$npc_term[npc_term_index]
      npc_column <- npc_level_columns[[npc_level]]
      npc_probability_column <- npc_probability_columns[[npc_level]]

      if (!npc_column %in% colnames(npc_variable_meta)) {
        warning(sprintf("Skipping %s NPC %s '%s': missing column %s.", npc_source_label, npc_level, npc_term, npc_column))
        next
      }
      npc_probability_values <- if (npc_probability_column %in% colnames(npc_variable_meta)) {
        suppressWarnings(as.numeric(npc_variable_meta[[npc_probability_column]]))
      } else {
        rep(1, nrow(npc_variable_meta))
      }
      npc_feature_keep <- !is.na(npc_variable_meta[[npc_column]]) &
        npc_variable_meta[[npc_column]] == npc_term &
        !is.na(npc_probability_values) &
        npc_probability_values >= npc_min_probability
      npc_feature_ids <- as.character(npc_variable_meta$feature_id[npc_feature_keep])
      npc_feature_ids <- intersect(npc_feature_ids, colnames(npc_data_matrix))

      if (!length(npc_feature_ids)) {
        warning(sprintf("No retained features matched %s NPC %s '%s'.", npc_source_label, npc_level, npc_term))
        next
      }

      npc_summed_values <- rowSums(npc_data_matrix[, npc_feature_ids, drop = FALSE], na.rm = TRUE)
      npc_label <- paste0("NPC ", npc_level, ": ", npc_term)
      npc_term_df <- data.frame(
        sample_id = npc_sample_meta$sample_id,
        group = as.character(npc_sample_meta[[params$target$sample_metadata_header]]),
        data_source = npc_source_label,
        npc_level = npc_level,
        npc_term = npc_term,
        npc_label = npc_label,
        npc_plot_label = npc_term,
        n_features = length(npc_feature_ids),
        summed_intensity = npc_summed_values,
        stringsAsFactors = FALSE
      )
      npc_summed_intensity_tables[[length(npc_summed_intensity_tables) + 1]] <- npc_term_df
      npc_feature_sets[[length(npc_feature_sets) + 1]] <- list(
        npc_level = npc_level,
        npc_term = npc_term,
        npc_label = npc_label,
        feature_ids = npc_feature_ids
      )
    }

    if (!length(npc_summed_intensity_tables)) {
      warning(sprintf("NPC summed-intensity plotting was requested for %s data, but no selected terms matched retained features.", npc_source_label))
      return(invisible(NULL))
    }

    npc_summed_intensity_df <- bind_rows(npc_summed_intensity_tables)
    npc_group_levels <- names(custom_colors)[names(custom_colors) %in% unique(npc_summed_intensity_df$group)]
    if (!length(npc_group_levels)) {
      npc_group_levels <- sort(unique(npc_summed_intensity_df$group))
    }
    npc_summed_intensity_df$group <- factor(npc_summed_intensity_df$group, levels = npc_group_levels)
    npc_summed_intensity_df <- npc_apply_plot_intensity(npc_summed_intensity_df)
    npc_intensity_stats_df <- npc_compute_stats(
      npc_summed_intensity_df,
      "npc_plot_label",
      "npc_label",
      "plot_intensity",
      npc_source_label,
      if (npc_transform == "log10") "log10_summed_intensity_plus_1" else "summed_intensity"
    )
    npc_summed_intensity_df <- npc_attach_overall_stats(npc_summed_intensity_df, npc_intensity_stats_df, "npc_plot_label")
    npc_intensity_static_labels <- npc_selected_static_labels(npc_intensity_stats_df)
    if (!length(npc_intensity_static_labels)) {
      npc_intensity_static_labels <- unique(npc_summed_intensity_df$npc_plot_label)
    }
    npc_summed_intensity_static_df <- npc_summed_intensity_df[npc_summed_intensity_df$npc_plot_label %in% npc_intensity_static_labels, , drop = FALSE]
    npc_intensity_label_count <- length(unique(npc_summed_intensity_static_df$npc_plot_label))
    npc_intensity_plot_height <- npc_combined_plot_height(npc_intensity_label_count)

    npc_summed_intensity_plot <- npc_intensity_plot(
      npc_summed_intensity_static_df,
      paste("Top", npc_intensity_label_count, "summed NPC-classified feature intensities for", params$mapp_batch, paste0("(", npc_source_label, " data)")),
      paste("Comparison across:", params$target$sample_metadata_header),
      facet = TRUE
    )

    npc_driver_dir <- file.path(dirname(filenames$intensity_table), "feature_drivers")
    npc_driver_tables <- list()
    for (npc_feature_set in npc_feature_sets) {
      npc_driver_df <- npc_feature_driver_table(
        npc_data_matrix,
        npc_sample_meta,
        npc_variable_meta,
        npc_feature_set$feature_ids,
        npc_feature_set$npc_level,
        npc_feature_set$npc_term,
        npc_feature_set$npc_label,
        npc_source_label
      )
      if (!nrow(npc_driver_df)) {
        next
      }
      npc_driver_suffix <- paste(npc_safe_file_part(npc_feature_set$npc_level), npc_safe_file_part(npc_feature_set$npc_term), sep = "_")
      npc_driver_file <- file.path(npc_driver_dir, paste0("NPC_summed_intensity_feature_drivers", file_suffix, "_", npc_driver_suffix, ".tsv"))
      npc_write_table(npc_driver_df, npc_driver_file)
      npc_driver_tables[[length(npc_driver_tables) + 1]] <- npc_driver_df
    }

    npc_write_table(npc_summed_intensity_df, filenames$intensity_table)
    npc_write_table(npc_intensity_stats_df, filenames$intensity_stats_table)
    if (length(npc_driver_tables)) {
      npc_write_table(bind_rows(npc_driver_tables), filenames$intensity_driver_table)
    }
    npc_save_plot(npc_summed_intensity_plot, filenames$intensity_pdf, height = npc_intensity_plot_height)
    npc_save_plot(npc_summed_intensity_plot, filenames$intensity_png, height = npc_intensity_plot_height)

    if (isTRUE(npc_ratio_enabled)) {
      npc_ratio_tables <- list()
      for (npc_feature_set in npc_feature_sets) {
        if (npc_feature_set$npc_level == "pathway") {
          next
        }
        npc_denominator_terms <- npc_terms_to_pathway(npc_feature_set$npc_level, npc_feature_set$npc_term)
        if (!length(npc_denominator_terms)) {
          warning(sprintf("Skipping %s NPC ratio for %s: no NP-Classifier pathway found.", npc_source_label, npc_feature_set$npc_label))
          next
        }
        npc_denominator_term <- npc_denominator_terms[1]
        if (length(npc_denominator_terms) > 1) {
          warning(sprintf(
            "NPC taxonomy maps %s to multiple pathways; using '%s' as denominator.",
            npc_feature_set$npc_label,
            npc_denominator_term
          ))
        }
        npc_pathway_probability_values <- if ("canopus_npc_pathway_probability" %in% colnames(npc_variable_meta)) {
          suppressWarnings(as.numeric(npc_variable_meta$canopus_npc_pathway_probability))
        } else {
          rep(1, nrow(npc_variable_meta))
        }
        npc_denominator_keep <- !is.na(npc_variable_meta$npc_taxonomy_pathway) &
          npc_pathway_contains(npc_variable_meta$npc_taxonomy_pathway, npc_denominator_term) &
          !is.na(npc_pathway_probability_values) &
          npc_pathway_probability_values >= npc_min_probability
        npc_denominator_feature_ids <- as.character(npc_variable_meta$feature_id[npc_denominator_keep])
        npc_denominator_feature_ids <- intersect(npc_denominator_feature_ids, colnames(npc_data_matrix))
        if (!length(npc_denominator_feature_ids)) {
          warning(sprintf("Skipping %s NPC ratio for %s: denominator pathway '%s' has no retained features.", npc_source_label, npc_feature_set$npc_label, npc_denominator_term))
          next
        }

        npc_ratio_numerator_feature_ids <- as.character(npc_feature_set$feature_ids)
        npc_ratio_numerator_feature_ids <- intersect(npc_ratio_numerator_feature_ids, colnames(npc_data_matrix))
        if (!length(npc_ratio_numerator_feature_ids)) {
          warning(sprintf("Skipping %s NPC ratio for %s: no numerator features remained inside denominator pathway '%s'.", npc_source_label, npc_feature_set$npc_label, npc_denominator_term))
          next
        }

        npc_numerator_values <- rowSums(npc_data_matrix[, npc_ratio_numerator_feature_ids, drop = FALSE], na.rm = TRUE)
        npc_denominator_values <- rowSums(npc_data_matrix[, npc_denominator_feature_ids, drop = FALSE], na.rm = TRUE)
        npc_ratio_denominator <- npc_denominator_values + npc_ratio_pseudocount
        npc_ratio_values <- ifelse(npc_ratio_denominator > 0, (npc_numerator_values + npc_ratio_pseudocount) / npc_ratio_denominator, NA_real_)
        npc_ratio_tables[[length(npc_ratio_tables) + 1]] <- data.frame(
          sample_id = npc_sample_meta$sample_id,
          group = as.character(npc_sample_meta[[params$target$sample_metadata_header]]),
          data_source = npc_source_label,
          numerator_level = npc_feature_set$npc_level,
          numerator_term = npc_feature_set$npc_term,
          denominator_level = "pathway",
          denominator_term = npc_denominator_term,
          npc_pathway = npc_denominator_term,
          ratio_label = paste0(npc_feature_set$npc_label, " / NPC pathway: ", npc_denominator_term),
          ratio_plot_label = npc_feature_set$npc_term,
          numerator_n_features = length(npc_ratio_numerator_feature_ids),
          denominator_n_features = length(npc_denominator_feature_ids),
          numerator_summed_intensity = npc_numerator_values,
          denominator_summed_intensity = npc_denominator_values,
          ratio = npc_ratio_values,
          stringsAsFactors = FALSE
        )
      }

      if (length(npc_ratio_tables)) {
        npc_ratio_df <- bind_rows(npc_ratio_tables)
        npc_ratio_df$group <- factor(npc_ratio_df$group, levels = npc_group_levels)
        npc_ratio_stats_df <- npc_compute_stats(
          npc_ratio_df,
          "ratio_plot_label",
          "ratio_label",
          "ratio",
          npc_source_label,
          "pathway_normalized_ratio"
        )
        npc_ratio_df <- npc_attach_overall_stats(npc_ratio_df, npc_ratio_stats_df, "ratio_plot_label")
        npc_ratio_static_labels <- npc_selected_static_labels(npc_ratio_stats_df)
        if (!length(npc_ratio_static_labels)) {
          npc_ratio_static_labels <- unique(npc_ratio_df$ratio_plot_label)
        }
        npc_ratio_static_df <- npc_ratio_df[npc_ratio_df$ratio_plot_label %in% npc_ratio_static_labels, , drop = FALSE]
        npc_ratio_label_count <- length(unique(npc_ratio_static_df$ratio_plot_label))
        npc_ratio_plot_height <- npc_combined_plot_height(npc_ratio_label_count)
        npc_ratio_combined_plot <- npc_ratio_plot(
          npc_ratio_static_df,
          paste("Top", npc_ratio_label_count, "NPC class/superclass pathway-normalized intensities for", params$mapp_batch, paste0("(", npc_source_label, " data)")),
          paste("Comparison across:", params$target$sample_metadata_header),
          facet = TRUE
        )
        npc_write_table(npc_ratio_df, filenames$ratio_table)
        npc_write_table(npc_ratio_stats_df, filenames$ratio_stats_table)
        npc_save_plot(npc_ratio_combined_plot, filenames$ratio_pdf, height = npc_ratio_plot_height)
        npc_save_plot(npc_ratio_combined_plot, filenames$ratio_png, height = npc_ratio_plot_height)
      }
    }
  }

  npc_filtered_data_matrix <- DE$data
  npc_filtered_sample_meta <- DE$sample_meta[rownames(npc_filtered_data_matrix), , drop = FALSE]
  npc_filtered_variable_meta <- DE$variable_meta[as.character(colnames(npc_filtered_data_matrix)), , drop = FALSE]
  npc_process_source(
    npc_filtered_data_matrix,
    npc_filtered_sample_meta,
    npc_filtered_variable_meta,
    "filtered",
    list(
      intensity_table = filename_npc_summed_intensity_table,
      intensity_stats_table = filename_npc_summed_intensity_stats_table,
      intensity_driver_table = filename_npc_summed_intensity_driver_table,
      intensity_pdf = filename_npc_summed_intensity_pdf,
      intensity_png = filename_npc_summed_intensity_png,
      ratio_table = filename_npc_summed_intensity_ratio_table,
      ratio_stats_table = filename_npc_summed_intensity_ratio_stats_table,
      ratio_pdf = filename_npc_summed_intensity_ratio_pdf,
      ratio_png = filename_npc_summed_intensity_ratio_png
    ),
    ""
  )

  if (isTRUE(npc_export_raw)) {
    npc_raw_sample_ids <- intersect(rownames(DE$data), rownames(DE_original$data))
    npc_raw_data_matrix <- DE_original$data[npc_raw_sample_ids, , drop = FALSE]
    npc_process_source(
      npc_raw_data_matrix,
      DE$sample_meta[npc_raw_sample_ids, , drop = FALSE],
      DE_original$variable_meta[as.character(colnames(DE_original$data)), , drop = FALSE],
      "raw",
      list(
        intensity_table = filename_npc_summed_intensity_raw_table,
        intensity_stats_table = filename_npc_summed_intensity_raw_stats_table,
        intensity_driver_table = filename_npc_summed_intensity_raw_driver_table,
        intensity_pdf = filename_npc_summed_intensity_raw_pdf,
        intensity_png = filename_npc_summed_intensity_raw_png,
        ratio_table = filename_npc_summed_intensity_ratio_raw_table,
        ratio_stats_table = filename_npc_summed_intensity_ratio_raw_stats_table,
        ratio_pdf = filename_npc_summed_intensity_ratio_raw_pdf,
        ratio_png = filename_npc_summed_intensity_ratio_raw_png
      ),
      "_raw"
    )
  }
}

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
  label_size = ordination_label_size,
  ellipse_type = "t",
  ellipse_confidence = 0.9,
  points_to_label = ordination_points_to_label
)

# We keep the PCA scores

pca_scores = pca_object$scores$data

# We keep the PCA loadings

pca_loadings = pca_object$loadings

# plot
pca_plot <- chart_plot(pca_scores_plot, pca_object)



fig_PCA <- pca_plot +
  facet_wrap(~ pca_plot$labels$title) +
  ggtitle(title_PCA) +
  publication_ordination_theme()


fig_PCA <- fig_PCA + ordination_colour_scale()
fig_PCA <- apply_ordination_point_style(fig_PCA)



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

ggsave(plot = fig_PCA, filename = filename_PCA, width = ordination_export_width, height = ordination_export_height, units = "in")
ggsave(plot = fig_PCA, filename = filename_PCA_svg, width = ordination_export_width, height = ordination_export_height, units = "in")


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

# We export the PCA scores

pca_scores <- pca_scores %>%
  rownames_to_column(var = "samples")

write.table(pca_scores, file = filename_PCA_scores, sep = "\t", row.names = FALSE)

# We export the PCA loadings

# First we add the missing column name


pca_loadings <- pca_loadings %>%
  rownames_to_column(var = "features")  %>% 
  # Kept as numeric
  mutate(features = as.numeric(features))


write.table(pca_loadings, file = filename_PCA_loadings, sep = "\t", row.names = FALSE)


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

  # We keep the PLSDA scores

  plsda_scores = plsda_object$scores$data

  # We keep the PLSDA loadings

  plsda_loadings = plsda_object$loadings


  # We merge the plsda_object$vip object with the DE$variable_meta object. We use dplyr syntax and keep a new object called variable_meta_plsda. We keep the rownames of the plsda_object$vip.

  plsda_object_vip <- plsda_object$vip %>%
    rownames_to_column(var = "feature_id") %>%
    mutate(feature_id = as.numeric(feature_id))

  vip_variable_meta <- plsda_object_vip %>%
    left_join(DE$variable_meta, by = "feature_id") %>%
    select(feature_id, feature_id_full_annotated)

  rownames(plsda_object$vip) <- vip_variable_meta$feature_id_full_annotated


  C <- pls_scores_plot(
    factor_name = params$target$sample_metadata_header,
    label_factor = "sample_id",
    label_size = ordination_label_size,
    points_to_label = ordination_points_to_label
  )

  plsda_plot <- tryCatch(
    chart_plot(C, plsda_object),
    error = function(err) {
      warning(sprintf("Could not build PLSDA chart with chart_plot(): %s", conditionMessage(err)))
      NULL
    }
  )

  if (inherits(plsda_plot, "ggplot")) {
    fig_PLSDA <- plsda_plot +
      ggtitle(title_PLSDA) +
      publication_ordination_theme()
  } else {
    warning("Could not build PLSDA chart with chart_plot(); using a direct scores plot fallback.")
    fig_PLSDA <- build_manual_plsda_scores_plot(plsda_object) +
      ggtitle(title_PLSDA) +
      publication_ordination_theme()
  }


  fig_PLSDA <- fig_PLSDA + ordination_colour_scale()
  fig_PLSDA <- apply_ordination_point_style(fig_PLSDA)


  # We output the feature importance

  C <- plsda_feature_importance_plot(n_features = 30, metric = "vip")

  vip_plot <- tryCatch(
    chart_plot(C, plsda_object),
    error = function(err) {
      warning(sprintf("Could not build PLSDA VIP chart with chart_plot(): %s", conditionMessage(err)))
      NULL
    }
  )



  fig_PLSDA_VIP <- NULL
  if (inherits(vip_plot, "ggplot")) {
    fig_PLSDA_VIP <- vip_plot + theme_classic() + ggtitle(title_PLSDA_VIP)
  } else {
    warning("Skipping PLSDA VIP PDF because the VIP chart object was not available; the VIP TSV will still be exported.")
  }

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

  ggsave(plot = fig_PLSDA, filename = filename_PLSDA, width = ordination_export_width, height = ordination_export_height, units = "in")
  if (inherits(fig_PLSDA_VIP, "ggplot")) {
    ggsave(plot = fig_PLSDA_VIP, filename = filename_PLSDA_VIP_plot, width = 20, height = 10)
  }

  # We export the loadings

  # write.table(loadings, file = filename_PLSDA_loadings, sep = "\t", row.names = FALSE)

  # We export the vip

  write.table(vip, file = filename_PLSDA_VIP_table, sep = "\t", row.names = FALSE)


  # We export the PLSDA scores

  plsda_scores <- plsda_scores %>%
    rownames_to_column(var = "samples")

  write.table(plsda_scores, file = filename_PLSDA_scores, sep = "\t", row.names = FALSE)

  # We export the PLSDA loadings

  # First we add the missing column name


  plsda_loadings <- plsda_loadings %>%
    rownames_to_column(var = "features")  %>% 
    # Kept as numeric
    mutate(features = as.numeric(features))


  write.table(plsda_loadings, file = filename_PLSDA_loadings, sep = "\t", row.names = FALSE)

}

#################################################################################################
#################################################################################################
#################################################################################################
# ##### DFA filtered data #########################################################################


if (!is.null(params$actions$run_DFA) && params$actions$run_DFA == "TRUE") {
  message("Launching DFA calculations ...")

  DE$sample_meta[, params$target$sample_metadata_header] <- as.factor(DE$sample_meta[, params$target$sample_metadata_header])

  dfa_seq_model <- # autoscale() +
    filter_na_count(threshold = 1, factor_name = params$target$sample_metadata_header) +
    DFA(factor_name = params$target$sample_metadata_header, number_components = 2)

  dfa_seq_result <- tryCatch(
    model_apply(dfa_seq_model, DE),
    error = function(e) {
      warning(paste("DFA failed on the full feature space:", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(dfa_seq_result)) {
    dfa_object <- dfa_seq_result[length(dfa_seq_result)]

    dfa_scores <- dfa_object$scores$data
    dfa_loadings <- dfa_object$loadings
    dfa_eigenvalues <- dfa_object$eigenvalues

    C <- dfa_scores_plot(factor_name = params$target$sample_metadata_header)

    dfa_plot <- chart_plot(C, dfa_object)

    fig_DFA <- dfa_plot + theme_classic() + facet_wrap(~ dfa_plot$labels$title) + ggtitle(title_DFA)

    fig_DFA <- fig_DFA +
      scale_colour_manual(name = "Groups", values = custom_colors)

    ggsave(plot = fig_DFA, filename = filename_DFA, width = 10, height = 10)

    dfa_scores <- dfa_scores %>%
      rownames_to_column(var = "samples")

    write.table(dfa_scores, file = filename_DFA_scores, sep = "\t", row.names = FALSE)

    dfa_loadings <- dfa_loadings %>%
      rownames_to_column(var = "features")

    write.table(dfa_loadings, file = filename_DFA_loadings, sep = "\t", row.names = FALSE)
    write.table(dfa_eigenvalues, file = filename_DFA_eigenvalues, sep = "\t", row.names = TRUE)
  }
}

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
  publication_ordination_theme()

if (ordination_points_to_label != "none") {
  pcoa_label_data <- data_PCOA_merge
  if (ordination_points_to_label == "outliers") {
    pcoa_distance_to_center <- sqrt((data_PCOA_merge$X1 - mean(data_PCOA_merge$X1, na.rm = TRUE))^2 + (data_PCOA_merge$X2 - mean(data_PCOA_merge$X2, na.rm = TRUE))^2)
    pcoa_outlier_threshold <- stats::quantile(pcoa_distance_to_center, probs = 0.9, na.rm = TRUE)
    pcoa_label_data <- data_PCOA_merge[pcoa_distance_to_center >= pcoa_outlier_threshold, , drop = FALSE]
  }
  fig_PCoA <- fig_PCoA +
    ggrepel::geom_text_repel(
      data = pcoa_label_data,
      aes(label = sample_name),
      size = ordination_label_size,
      show.legend = FALSE,
      max.overlaps = Inf
    )
}



fig_PCoA <- fig_PCoA + ordination_colour_scale()
fig_PCoA <- apply_ordination_point_style(fig_PCoA)



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

ggsave(plot = fig_PCoA, filename = filename_PCoA, width = ordination_export_width, height = ordination_export_height, units = "in")


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

  if (params$actions$scale_method == "pareto" | params$actions$scale_method == "none") {
    # The formula is defined externally
    formula <- as.formula(paste0(
      "y", "~", params$target$sample_metadata_header, "+",
      "Error(sample_id/",
      params$target$sample_metadata_header,
      ")"
    ))


    model <- HSDEM(
      alpha = params$posthoc$p_value,
      formula = formula, mtc = "none"
    )
  } else if (params$actions$scale_method == "autoscale") {
    # The formula is defined externally
    formula <- as.formula(paste0(
      "y", "~", params$target$sample_metadata_header
    ))

    model <- HSD(
      alpha = params$posthoc$p_value,
      formula = formula, mtc = "none", unbalanced = FALSE
    )
  }



  HSDEM_result <- model_apply(model, DE)

  HSDEM_result_p_value <- HSDEM_result$p_value


  # We split each colnames according to the `-` character. We then rebuild the colnames, alphabetically ordered.

  colnames(HSDEM_result_p_value) <- plotrix::pasteCols(sapply(strsplit(colnames(HSDEM_result_p_value), " - "), sort), sep = "_vs_")

  # Additionally we make sure to remove the headers name from the colnames (this one can be added when the data are numerics.)

  colnames(HSDEM_result_p_value) <- gsub(params$target$sample_metadata_header, "", colnames(HSDEM_result_p_value))

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

write.table(DE_foldchange_pvalues, file = filename_foldchange_pvalues, sep = ",", row.names = FALSE)



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
      contains("sirius_confidencescore"),
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
      contains("sirius_confidencescore"),
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
    )
    # We set the type of the sirius_confidencescoreapproximate column to numeric
    if ("sirius_confidencescoreapproximate" %in% colnames(DE_foldchange_pvalues)) {
      de4dt <- de4dt %>%
        mutate(sirius_confidencescoreapproximate = as.numeric(sirius_confidencescoreapproximate))
    } else if ("sirius_confidencescore" %in% colnames(DE_foldchange_pvalues)) {
      de4dt <- de4dt %>%
        mutate(sirius_confidencescore = as.numeric(sirius_confidencescore))
    }

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
      contains("sirius_confidencescore"),
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
      contains("sirius_confidencescore"),
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
    )
    # We set the type of the sirius_confidencescoreapproximate column to numeric
    if ("sirius_confidencescoreapproximate" %in% colnames(DE_foldchange_pvalues)) {
      de4dt <- de4dt %>%
        mutate(sirius_confidencescoreapproximate = as.numeric(sirius_confidencescoreapproximate))
    } else if ("sirius_confidencescore" %in% colnames(DE_foldchange_pvalues)) {
      de4dt <- de4dt %>%
        mutate(sirius_confidencescore = as.numeric(sirius_confidencescore))
    }
}
# We output a generic DT for data exploration of the whole set


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
    "canopus_npc_class_probability"
    # "sirius_confidencescoreapproximate"
  ), digits = 2)

    # Add conditional formatting for sirius_confidencescore columns
  if ("sirius_confidencescoreapproximate" %in% colnames(DE_foldchange_pvalues)) {
    DT_volcano <- DT_volcano %>%
      formatRound(c("sirius_confidencescoreapproximate"), digits = 2)
  } else if ("sirius_confidencescore" %in% colnames(DE_foldchange_pvalues)) {
    DT_volcano <- DT_volcano %>%
      formatRound(c("sirius_confidencescore"), digits = 2)
  }



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
      "canopus_npc_class_probability"
    # "sirius_confidencescoreapproximate"
  ), digits = 2)

    # Add conditional formatting for sirius_confidencescore columns
  if ("sirius_confidencescoreapproximate" %in% colnames(DE_foldchange_pvalues)) {
    DT_volcano <- DT_volcano %>%
      formatRound(c("sirius_confidencescoreapproximate"), digits = 2)
  } else if ("sirius_confidencescore" %in% colnames(DE_foldchange_pvalues)) {
    DT_volcano <- DT_volcano %>%
      formatRound(c("sirius_confidencescore"), digits = 2)
  }

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
    # We save the plot in a svg file
    tryCatch(
      {
        ggsave(plot = gg_volcano, filename = paste0("Volcano_", first_part, "_vs_", second_part, "_", label_column, ".svg"), width = 10, height = 10)
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

  # condition = "blastogenesis_vs_healing"
  # Iterate over the prefixes
  for (condition in conditions) {
    message("Generating treemaps for condition: ", condition)
    # glimpse(DE_foldchange_pvalues)
    # Perform filtering using the prefix as a condition
    DE_foldchange_pvalues_signi <- DE_foldchange_pvalues %>%
      filter(!!sym(paste0(condition, "_p_value")) < params$posthoc$p_value)
    # Check that the filtered data frame has at least 10 rows.
    # Else exit the loop and print a message.
    if (nrow(DE_foldchange_pvalues_signi) < 10) {
      message("Not enough significant features for condition: ", condition, ". Skipping treemap generation.\n")
      next
    }

    # Print the filtered data
    message("Filtered data for condition: ", condition, ":\n")
    # print(head(DE_foldchange_pvalues_signi))
    message("\n")

    condition_parts <- strsplit(condition, "_vs_")[[1]]
    first_part <- condition_parts[1]
    second_part <- condition_parts[2]
    # glimpse(DE_foldchange_pvalues_signi)
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
        title = list(
          text = paste0("<b>Metabolic variations across ", first_part, " vs ", second_part, "</b>", "<br>", "Sample metadata filters:", "<br>", "[", filter_sample_metadata_status, "]"),
          font = list(size = 14), # Adjust the font size
          x = 0.15, # Align title to the left
          xanchor = "left" # Ensure the title starts from the left
        ),
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
if (!nrow(matt_donust2)) {
  message("No significant CANOPUS-classified features; skipping the global treemap.")
} else {
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
  pull() %>%
  as.character() %>%
  intersect(colnames(DE$data))

imp_table_rf <- data.frame(
  MeanDecreaseGini = rep(NA_real_, length(features_of_importance)),
  row.names = paste0("X", features_of_importance)
)

if (length(features_of_importance) < 2) {
  message(sprintf(
    "Random Forest requires at least two significant predictors; found %d. Skipping Random Forest.",
    length(features_of_importance)
  ))
} else {


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

}



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
  pull() %>%
  as.character() %>%
  intersect(colnames(DE$data))

if (length(features_of_importance) < 2) {
  message(sprintf(
    "The p-value filtered heatmap requires at least two significant features; found %d. Skipping the heatmap.",
    length(features_of_importance)
  ))
} else {

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

#### Iheatmapr


target_metadata <- as.factor(DE$sample_meta[[params$target$sample_metadata_header]])


##########################
# We make sure to order the colors.
custom_colors_heatmap <- custom_colors[order(names(custom_colors))]


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
  add_row_annotation(data.frame("Classification" = mdf_selected_variable_meta_NPC_simple_resolved_colored$group),
    side = "right",
    buffer = 0.05,
    colors = list("Classification" = col_np_pathway)
  ) %>%
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


# The file is exported

iheatmap %>% save_iheatmap(file = filename_heatmap_pval) # Save interactive HTML

}




#############################################################################
#############################################################################
############## Summary Table ################################################
#############################################################################
#############################################################################

message("Outputing Summary Table ...")

# Output is not clean. Feature id are repeated x times.
# To tidy ---



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
    contains("sirius_confidencescore"),
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



#############################################################################
#############################################################################
############## GraphML output ################################################
#############################################################################
#############################################################################

# message("Generating GraphML output ...")


# # We first load the GNPS graphml file

# if (gnps2_job) {
#   graphml_dir <- file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "nf_output", "networking")
#   graphml_file <- list.files(path = graphml_dir, pattern = "network.graphml$")

# } else {
#   graphml_dir <- file.path(working_directory, "results", "met_annot_enhancer", params$gnps_job_id, "gnps_molecular_network_graphml")
#   graphml_file <- list.files(path = graphml_dir, pattern = "\\.graphml$")

# }


# graphml_file_path <- file.path(graphml_dir, graphml_file)



# g <- read.graph(file = graphml_file_path, format = "graphml")
# # net_gnps = igraph::simplify(g, remove.multiple = FALSE, edge.attr.comb = "ignore")

# df_from_graph_edges_original <- igraph::as_data_frame(g, what = c("edges"))
# df_from_graph_vertices_original <- igraph::as_data_frame(g, what = c("vertices"))


# # We define drop the from and to columns from the edges dataframe
# # And then rename the node 1 and node 2 columns to from and to, respectively
# # These columns are placed at the beginning of the dataframe
# # and converted to numerics
# if (gnps2_job) {
#   df_from_graph_edges <- df_from_graph_edges_original %>%
#     select(-from, -to) %>%
#     rename(from = scan1, to = scan2) %>%
#     select(from, to, everything()) %>%
#     mutate_at(vars(from, to), as.numeric)
# } else {
#   df_from_graph_edges <- df_from_graph_edges_original %>%
#     select(-from, -to) %>%
#     rename(from = node1, to = node2) %>%
#     select(from, to, everything()) %>%
#     mutate_at(vars(from, to), as.numeric)
# }

# # the id column of the vertices dataframe is converted to numerics

# df_from_graph_vertices <- df_from_graph_vertices_original %>%
#   mutate_at(vars(id), as.numeric)


# # We then add the attributes to the vertices dataframe
# # For this we merge the vertices dataframe with the VM output using the id column and the feature_id column, respectively

# # vm_minus_gnps = DE_original$variable_meta  %>%
# # select(-contains("_gnps"))


# # df_from_graph_vertices_plus = merge(df_from_graph_vertices, vm_minus_gnps, by.x = "id", by.y = "feature_id", all.x = T)

# # glimpse(df_from_graph_vertices_plus)


# # Now we will add the results of the statistical outputs to the vertices dataframe
# # For this we merge the vertices dataframe with the summary_stat_output using the id column and the feature_id column, respectively

# # First we clean the summary_stat_output dataframe
# # For this we remove columns that are not needed. The one containing the sirius and canopus pattern in the column names. Indeed they arr already present in the VM dataframe

# # summary_stat_output_red = summary_stat_output_full %>%
# #   select(-contains("_sirius")) %>%
# #   select(-contains("_canopus")) %>%
# #   select(-contains("_met_annot")) %>%
# #   select(-contains("_gnps"))  %>%
# #   #select(-ends_with("_id"))  %>%
# #   select(-ends_with("_mz"))  %>%
# #   select(-ends_with("_rt"))

# DE_original_features <- DE_original$variable_meta 
# # %>%
# #   select(feature_id)

# # We merge DE_original_features and summary_stat_output_selected_cytoscape but make sure to drop duplicated columns from the summary_stat_output_selected_cytoscape dataframe

# # Identify the names of columns that are duplicated between the two data frames, excluding the merging key column.
# common_cols <- setdiff(intersect(names(DE_original_features), names(summary_stat_output_selected_cytoscape)), "feature_id")

# # Remove Duplicated Columns from One DataFrame
# summary_stat_output_selected_cytoscape <- summary_stat_output_selected_cytoscape %>%
#   select(-all_of(common_cols))


# df_from_graph_vertices_plus <- DE_original_features %>%
#   left_join(summary_stat_output_selected_cytoscape, by = "feature_id")



# # We merge the data from the DE$data dataframe with the DE$sample_meta dataframe using rownames as the key

# merged_D_SM <- merge(DE_original$sample_meta, DE_original$data, by = "row.names", all = TRUE)

# # We replace NA values with 0 in the merged dataframe
# # Check why we dont do this before ??
# merged_D_SM[is.na(merged_D_SM)] <- 0


# # The function below allows to group data by multiple factors and return a dataframe with the mean of each group


# dfList <- list()

# for (i in params$to_mean$factor_name) {
#   dfList[[i]] <- merged_D_SM %>%
#     group_by(!!as.symbol(i)) %>%
#     summarise(across(colnames(DE_original$data), mean),
#       .groups = "drop"
#     ) %>%
#     select(!!all_of(i), colnames(DE_original$data)) %>%
#     pivot_longer(-!!i) %>%
#     pivot_wider(names_from = all_of(i), values_from = value) %>%
#     # We prefix all columns with the factor name
#     rename_with(.cols = -name, ~ paste0("mean_int", "_", i, "_", .x))
# }


# flat_dfList <- reduce(dfList, full_join, by = "name")

# # We now add the raw feature list to the dataframe

# full_flat_dfList <- merge(flat_dfList, t(DE_original$data), by.x = "name", by.y = "row.names", all = TRUE)


# #  We add the raw feature list

# df_from_graph_vertices_plus_plus <- merge(df_from_graph_vertices_plus, full_flat_dfList, by.x = "feature_id", by.y = "name", all.x = T)


# # We set back the id column as the first column of the dataframe

# df_from_graph_vertices_plus_plus <- df_from_graph_vertices_plus_plus %>%
#   select(feature_id, everything())

# node_size <- df_from_graph_vertices_plus_plus$MeanDecreaseGini
# node_size[is.na(node_size)] <- min(node_size, na.rm = T)

# df_from_graph_vertices_plus_plus$node_size <- node_size
# # We then add the attributes to the edges dataframe and generate the igraph object

# # In the case when we have been filtering the X data we will add the filtered X data to the vertices dataframe prior to merging.


# # We then make sure to have the df_from_graph_vertices_plus_plus dataframe ordered by decreasing value of the MeanDecreaseGini column

# df_from_graph_vertices_plus_plus <- df_from_graph_vertices_plus_plus %>%
#   arrange(desc(MeanDecreaseGini))


# generated_g <- graph_from_data_frame(df_from_graph_edges, directed = FALSE, vertices = df_from_graph_vertices_plus_plus)


################################################################################
################################################################################
##### add annotations to igraph


# The file is exported

# Not outputted by default
# write_graph(generated_g, file = filename_graphml, format = "graphml")


message("... the R session info file ...")

sink(filename_session_info)
sessionInfo()
sink()

message("... and the R script file !")

# script_path is the absolute path resolved at startup — safe to use after any setwd()
file.copy(script_path, file.path(output_directory, filename_R_script), overwrite = TRUE)


message("Done !")
