# Candidate rows are not joined across annotation engines: their ranks and classifications
# describe different evidence and a join would multiply candidates for the same feature.
normalize_annotations <- function(dataset) {
  feature_ids <- as.character(colnames(dataset$raw))
  if (is.null(feature_ids)) feature_ids <- as.character(dataset$variable_metadata$feature_id)
  if (anyNA(feature_ids) || anyDuplicated(feature_ids)) stop("Dataset feature IDs must be unique and non-missing.", call. = FALSE)

  columns <- c("feature_id", "source", "candidate_id", "label", "confidence", "confidence_metric",
               "npc_pathway", "npc_superclass", "npc_class", "component_id", "feature_in_dataset",
               "rank", "molecular_formula", "adduct", "inchikey", "smiles", "library", "msms_score",
               "final_score", "class_probability", "superclass_probability", "pathway_probability",
               "organism", "identification_links")
  empty <- as.data.frame(setNames(rep(list(character()), length(columns)), columns), stringsAsFactors = FALSE)
  empty$confidence <- empty$msms_score <- empty$final_score <- empty$class_probability <-
    empty$superclass_probability <- empty$pathway_probability <- numeric()
  empty$feature_in_dataset <- logical()

  pick <- function(table, keys, numeric = FALSE) {
    result <- if (numeric) rep(NA_real_, nrow(table)) else rep(NA_character_, nrow(table))
    for (column in intersect(keys, names(table))) {
      value <- table[[column]]
      if (numeric) {
        value <- suppressWarnings(as.numeric(value))
      } else {
        value <- trimws(as.character(value))
        value[is.na(value) | !nzchar(value) | value %in% c("NA", "None", "null")] <- NA_character_
      }
      replace <- is.na(result) & !is.na(value)
      result[replace] <- value[replace]
    }
    result
  }
  append_source <- function(table, source, id_columns, label_columns, confidence_columns, metric,
                            pathway_columns, superclass_columns, class_columns, component_columns = character(),
                            rank_columns = character(), formula_columns = character(), adduct_columns = character(),
                            inchikey_columns = character(), smiles_columns = character(), library_columns = character(),
                            msms_columns = character(), final_columns = character(), organism_columns = character(),
                            links_columns = character(), class_probability_columns = character(),
                            superclass_probability_columns = character(), pathway_probability_columns = character()) {
    if (is.null(table)) return(empty)
    if (!is.data.frame(table)) stop("Annotation source must be a table: ", source, call. = FALSE)
    if (!nrow(table)) return(empty)
    if (!length(intersect(id_columns, names(table)))) stop("Annotation source lacks feature IDs: ", source, call. = FALSE)
    ids <- pick(table, id_columns)
    result <- data.frame(
      feature_id = ids, source = source, candidate_id = paste0(source, ":", seq_len(nrow(table))),
      label = pick(table, label_columns), confidence = pick(table, confidence_columns, TRUE),
      confidence_metric = metric, npc_pathway = pick(table, pathway_columns),
      npc_superclass = pick(table, superclass_columns), npc_class = pick(table, class_columns),
      component_id = pick(table, component_columns), feature_in_dataset = ids %in% feature_ids,
      rank = pick(table, rank_columns), molecular_formula = pick(table, formula_columns),
      adduct = pick(table, adduct_columns), inchikey = pick(table, inchikey_columns),
      smiles = pick(table, smiles_columns), library = pick(table, library_columns),
      msms_score = pick(table, msms_columns, TRUE), final_score = pick(table, final_columns, TRUE),
      class_probability = pick(table, class_probability_columns, TRUE),
      superclass_probability = pick(table, superclass_probability_columns, TRUE),
      pathway_probability = pick(table, pathway_probability_columns, TRUE),
      organism = pick(table, organism_columns), identification_links = pick(table, links_columns),
      stringsAsFactors = FALSE
    )
    # Retain unjoinable candidates for provenance instead of dropping them silently.
    result
  }
  tables <- dataset$annotations
  if (is.null(tables)) tables <- list()
  structures <- tables$structures
  if (is.null(structures)) structures <- tables$sirius
  canopus <- tables$canopus
  enhancer <- tables$met_annot_enhancer
  rows <- list(
    append_source(structures, "sirius", c("mappingFeatureId", "featureId", "feature_id", "sirius_mappingFeatureId"),
                  c("name", "chebiasciiname", "sirius_name", "molecularFormula"),
                  c("ConfidenceScoreExact", "ConfidenceScoreApproximate", "sirius_ConfidenceScoreExact"),
                  "SIRIUS confidence score", character(), character(), character(),
                  rank_columns = c("structurePerIdRank", "sirius_structurePerIdRank"),
                  formula_columns = c("molecularFormula", "sirius_molecularFormula"),
                  adduct_columns = c("adduct", "sirius_adduct"),
                  inchikey_columns = c("InChIkey2D", "sirius_InChIkey2D"),
                  smiles_columns = c("smiles", "sirius_smiles"), links_columns = c("links", "pubchemids")),
    append_source(canopus, "canopus", c("mappingFeatureId", "featureId", "feature_id", "canopus_mappingFeatureId"),
                  c("NPC#class", "NPC#superclass", "NPC#pathway", "molecularFormula"),
                  c("NPC#class Probability", "npc_class_probability"),
                  "CANOPUS NPC class probability", c("NPC#pathway", "npc_pathway"),
                  c("NPC#superclass", "npc_superclass"), c("NPC#class", "npc_class"),
                  rank_columns = c("formulaRank"), formula_columns = c("molecularFormula"),
                  adduct_columns = c("adduct"),
                  class_probability_columns = c("NPC#class Probability", "npc_class_probability"),
                  superclass_probability_columns = c("NPC#superclass Probability", "npc_superclass_probability"),
                  pathway_probability_columns = c("NPC#pathway Probability", "npc_pathway_probability")),
    append_source(enhancer, "met_annot_enhancer", c("feature_id", "mappingFeatureId"),
                  c("structure_nameTraditional", "structure_name", "structure_inchikey", "short_inchikey"),
                  c("final_score"), "met-annot-enhancer final score (not probability)",
                  c("structure_taxonomy_npclassifier_01pathway", "structure_taxonomy_npclassifier_01pathway_consensus"),
                  c("structure_taxonomy_npclassifier_02superclass", "structure_taxonomy_npclassifier_02superclass_consensus"),
                  c("structure_taxonomy_npclassifier_03class", "structure_taxonomy_npclassifier_03class_consensus"),
                  component_columns = c("component_id"), rank_columns = c("rank_spec"),
                  formula_columns = c("structure_molecular_formula"), adduct_columns = c("adduct"),
                  inchikey_columns = c("structure_inchikey", "short_inchikey"),
                  smiles_columns = c("structure_smiles"), library_columns = c("libname"),
                  msms_columns = c("msms_score"), final_columns = c("final_score"),
                  organism_columns = c("organism_name"), links_columns = c("structure_wikidata"))
  )
  candidates <- do.call(rbind, rows)
  if (is.null(candidates)) candidates <- empty
  unannotated <- setdiff(feature_ids, candidates$feature_id[candidates$feature_in_dataset])
  if (length(unannotated)) {
    missing <- empty[rep(NA_integer_, length(unannotated)), , drop = FALSE]
    missing$feature_id <- unannotated
    missing$feature_in_dataset <- TRUE
    candidates <- rbind(candidates, missing)
  }
  rownames(candidates) <- NULL
  candidates
}
