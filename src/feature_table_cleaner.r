############################################################################################
############################################################################################
##################################### PACKAGES #####################################
############################################################################################
############################################################################################


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

usePackage("base")
usePackage("data.table")
usePackage("dplyr")
usePackage("tidyverse")
usePackage("yaml")




############################################################################################
############################################################################################
################################LOAD & FORMATDATA#################################
############################################################################################
############################################################################################


# We call the external params


params = yaml.load_file("./params/params_ft_cleaner.yaml")


# Make this more generic


# We load the data from the params file (file path is in params$paths$file_to_clean)


data = fread(params$paths$file_to_clean)

# We apply the filtering according to the params file (filtering conditions are in params$filters)

data_sel = data %>% 
filter(`old_Wt - young_Wt.posthoc_Pvalue` <= 0.400 &
`old_Arg_II - old_Wt.posthoc_Pvalue` <= 0.600 &
componentindex == 24)

data_sel = data %>% 
filter(`young_Arg_II - young_Wt.posthoc_Pvalue` <= 0.03)

data_sel = data %>% 
filter(`old_Wt - young_Wt.posthoc_Pvalue` <= 0.125)




# We define the columns to keep

columns_to_keep_full = c(
"adduct_sirius",
"#adducts_sirius",
"charge",
"cluster index",
"componentindex",
"Compound_Name",
"GNPSLibraryURL",
"id_sirius",
"INCHI", 
"InChI_sirius",
"InChIkey2D_sirius",
"libname_metannot",
"links_sirius",
"matched_class_metannot",
"matched_domain_metannot",
"matched_family_metannot",
"matched_genus_metannot",
"matched_kingdom_metannot",
"matched_order_metannot",
"matched_phylum_metannot",
"matched_species_metannot",
"molecularFormula_canopus",
"molecularFormula_sirius",
"msms_score_metannot",
"MZErrorPPM",
"name",
"name_sirius",
"NPC#class_canopus",
"NPC#pathway_canopus",
"NPC#superclass_canopus",
"organism_name_metannot",
"organism_taxonomy_01domain_metannot",
"organism_taxonomy_02kingdom_metannot",
"organism_taxonomy_03phylum_metannot",
"organism_taxonomy_04class_metannot",
"organism_taxonomy_05order_metannot",
"organism_taxonomy_06family_metannot",
"organism_taxonomy_07tribe_metannot",
"organism_taxonomy_08genus_metannot",
"organism_taxonomy_09species_metannot",
"organism_taxonomy_ottid_metannot",
"organism_wikidata_metannot",
"parent mass",
"PI",
"precursor mass",
"pubchemids_sirius",
"rank_final_metannot",
"rank_sirius",
"rank_spec_metannot",
"retentionTimeInSeconds_sirius",
"row_ID",
"row_mz_full",
"row_rt_full",
"shared name",
"short_inchikey_metannot",
"Smiles",
"smiles_sirius",
"structure_exact_mass_metannot",
"structure_inchi_metannot",
"structure_inchikey_metannot",
"structure_molecular_formula_metannot",
"structure_nameTraditional_metannot",
"structure_smiles_metannot",
"structure_taxonomy_npclassifier_01pathway_consensus_metannot",
"structure_taxonomy_npclassifier_01pathway_metannot",
"structure_taxonomy_npclassifier_02superclass_consensus_metannot",
"structure_taxonomy_npclassifier_02superclass_metannot",
"structure_taxonomy_npclassifier_03class_consensus_metannot",
"structure_taxonomy_npclassifier_03class_metannot",
"structure_wikidata_metannot")


columns_to_keep_min = c(
"cluster index",
"componentindex",
"row_ID",
"row_mz_full",
"row_rt_full",
"shared name",
"adduct_sirius",
"Compound_Name",
"INCHI", 
"InChI_sirius",
"InChIkey2D_sirius",
"name",
"name_sirius",
"NPC#class_canopus",
"NPC#pathway_canopus",
"NPC#superclass_canopus",
"organism_wikidata_metannot",
"parent mass",
"precursor mass",
"pubchemids_sirius",
"short_inchikey_metannot",
"Smiles",
"smiles_sirius",
"structure_exact_mass_metannot",
"structure_inchi_metannot",
"structure_inchikey_metannot",
"structure_molecular_formula_metannot",
"structure_nameTraditional_metannot",
"structure_smiles_metannot",
"structure_taxonomy_npclassifier_01pathway_metannot",
"structure_taxonomy_npclassifier_02superclass_metannot",
"structure_taxonomy_npclassifier_03class_metannot",
"structure_wikidata_metannot")



columns_to_keep_ultra_min = c(
"cluster index",
"componentindex",
"row_ID",
"row_mz_full",
"row_rt_full",
"Compound_Name",
"name_sirius",
"structure_nameTraditional_metannot",
"Smiles",
"smiles_sirius",
"structure_smiles_metannot",
"INCHI", 
"InChI_sirius",
"InChIkey2D_sirius",
"structure_inchi_metannot",
"structure_inchikey_metannot",
"NPC#class_canopus",
"NPC#pathway_canopus",
"NPC#superclass_canopus",
"organism_wikidata_metannot",
"pubchemids_sirius",
"short_inchikey_metannot",
"structure_exact_mass_metannot",
"structure_molecular_formula_metannot",
"structure_taxonomy_npclassifier_01pathway_metannot",
"structure_taxonomy_npclassifier_02superclass_metannot",
"structure_taxonomy_npclassifier_03class_metannot",
"structure_wikidata_metannot")


# We keep only the columns defined above, using dplyr syntax

data_sel = data_sel %>% select(columns_to_keep_ultra_min)


data_sel


# We then save the data in the params file (file path is in params$paths$cleaned_file)

fwrite(data_sel, params$paths$cleaned_file)


