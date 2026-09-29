if (!requireNamespace("shiny", quietly = TRUE)) stop("The V2 explorer requires shiny.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The V2 explorer requires ggplot2.")
if (!requireNamespace("yaml", quietly = TRUE)) stop("The V2 explorer requires yaml.")
if (!requireNamespace("digest", quietly = TRUE)) stop("The V2 explorer requires digest.")

library(shiny)

# The administrator supplies a root; visitors never supply filesystem paths.
configured_root <- Sys.getenv("MAPP_STATS_ROOT", unset = "")
root <- if (nzchar(configured_root) && dir.exists(configured_root)) normalizePath(configured_root, mustWork = TRUE) else NULL

within_root <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  identical(resolved, root) || startsWith(resolved, paste0(root, .Platform$file.sep))
}

safe_path <- function(run_dir, relative) {
  if (is.null(root) || !grepl("^(objects|tables)/[A-Za-z0-9_.-]+[.](rds|tsv)$", relative))
    stop("Unapproved run file.", call. = FALSE)
  path <- file.path(run_dir, relative)
  if (!file.exists(path) || !within_root(path)) stop("Run file is missing or outside the approved root.", call. = FALSE)
  parts <- strsplit(substring(normalizePath(run_dir, mustWork = FALSE), nchar(root) + 2L), .Platform$file.sep, fixed = TRUE)[[1]]
  chain <- root
  for (part in c(parts, strsplit(relative, "/", fixed = TRUE)[[1]])) {
    chain <- file.path(chain, part)
    if (nzchar(Sys.readlink(chain))) stop("Linked paths are not allowed in completed runs.", call. = FALSE)
  }
  path
}

verified_file <- function(run, relative) {
  expected <- run$manifest$outputs[[relative]]
  if (!is.character(expected) || length(expected) != 1L || !grepl("^[[:xdigit:]]{64}$", expected))
    stop(paste("Run has no valid checksum for", relative), call. = FALSE)
  path <- safe_path(run$path, relative)
  actual <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(actual), tolower(expected))) stop(paste("Checksum mismatch for", relative), call. = FALSE)
  path
}

# Walk MAPP project/batch directories only, never other results trees or symlinks.
catalogue_runs <- function() {
  if (is.null(root)) return(data.frame(label = character(), path = character()))
  batches <- character()
  walk <- function(directory) {
    if (!within_root(directory) || nzchar(Sys.readlink(directory))) return()
    if (dir.exists(file.path(directory, "results"))) {
      batches <<- c(batches, directory)
      return()
    }
    children <- list.files(directory, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE)
    for (child in children[dir.exists(children) & grepl("^mapp_(project|batch)_", basename(children))]) walk(child)
  }
  walk(root)
  paths <- character()
  labels <- character()
  for (batch in batches) {
    stats_dir <- file.path(batch, "results", "stats_v2")
    if (!dir.exists(stats_dir) || !within_root(stats_dir) || nzchar(Sys.readlink(stats_dir))) next
    entries <- list.files(stats_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE)
    for (entry in entries[dir.exists(entries) & grepl("^[[:xdigit:]]{64}$", basename(entries))]) {
      if (!within_root(entry) || nzchar(Sys.readlink(entry))) next
      if (!file.exists(file.path(entry, "COMPLETE")) || !file.exists(file.path(entry, "manifest.yaml"))) next
      # Neither a marker nor a manifest alone authorizes an object read.
      if (!file.exists(file.path(entry, "objects", "processed.rds")) ||
          !file.exists(file.path(entry, "objects", "analyses.rds"))) next
      paths <- c(paths, entry)
      labels <- c(labels, paste(basename(batch), substr(basename(entry), 1, 12), sep = " / "))
    }
  }
  data.frame(label = labels, path = paths, stringsAsFactors = FALSE)
}

open_run <- function(path) {
  if (is.null(root) || !within_root(path) || !grepl("^[[:xdigit:]]{64}$", basename(path))) stop("Run is outside the catalogue.", call. = FALSE)
  marker <- file.path(path, "COMPLETE")
  manifest_file <- file.path(path, "manifest.yaml")
  if (!file.exists(marker) || !file.exists(manifest_file) || nzchar(Sys.readlink(marker)) || nzchar(Sys.readlink(manifest_file)) ||
      !within_root(marker) || !within_root(manifest_file)) stop("Run is not safely completed.", call. = FALSE)
  completion <- yaml::read_yaml(marker)
  if (!is.list(completion) || !identical(completion$run_hash, basename(path)) ||
      !is.character(completion$content_sha256) || length(completion$content_sha256) != 1L ||
      !grepl("^[[:xdigit:]]{64}$", completion$content_sha256))
    stop("Run completion marker is incomplete or does not match the selected run.", call. = FALSE)
  manifest <- yaml::read_yaml(manifest_file)
  if (!identical(manifest$schema_version, "1.0.0"))
    stop("Unsupported completed-run schema (expected 1.0.0).", call. = FALSE)
  if (!identical(manifest$status, "complete") || !identical(manifest$hashes$run_hash, basename(path)) || !is.list(manifest$outputs))
    stop("Run manifest is incomplete or does not match the selected run.", call. = FALSE)
  run <- list(path = path, manifest = manifest)
  # Check the sealed objects even though this viewer never deserializes RDS.
  verified_file(run, "objects/processed.rds")
  verified_file(run, "objects/analyses.rds")
  run
}

saved_table <- function(run, name, optional = FALSE) {
  relative <- paste0("tables/", name, ".tsv")
  if (optional && is.null(run$manifest$outputs[[relative]])) return(NULL)
  path <- verified_file(run, relative)
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"), quote = "", comment.char = "")
}

available_tables <- function(run) {
  folder <- file.path(run$path, "tables")
  if (!dir.exists(folder) || !within_root(folder) || nzchar(Sys.readlink(folder))) return(character())
  names <- list.files(folder, pattern = "^[A-Za-z0-9_.-]+[.]tsv$", full.names = FALSE)
  names[vapply(names, function(name) !is.null(run$manifest$outputs[[paste0("tables/", name)]]), logical(1))]
}

ui <- fluidPage(
  titlePanel("MAPP statistics | Completed run explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("run", "Completed batch / run", choices = character()),
      actionButton("refresh", "Refresh catalogue"),
      hr(), uiOutput("controls"), width = 3
    ),
    mainPanel(
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Overview", tableOutput("summary"), h4("Samples by selected group"), tableOutput("group_counts")),
        tabPanel("PCA", plotOutput("pca", height = 550), textOutput("pca_note")),
        tabPanel("PCoA", plotOutput("pcoa", height = 550), textOutput("pcoa_note")),
        tabPanel("Saved contrasts", plotOutput("volcano", height = 550, click = "volcano_click"),
                 tableOutput("volcano_table")),
        tabPanel("Validated PLS-DA", uiOutput("plsda_interpretation"), tableOutput("plsda_validation"),
                 h4("Held-out fold predictions"), tableOutput("plsda_predictions"),
                 h4("Saved label permutations"), plotOutput("plsda_null", height = 320), tableOutput("plsda_permutations"),
                 h4("Full-fit scores (illustration, not validation)"), plotOutput("plsda_scores", height = 500)),
        tabPanel("Feature", plotOutput("feature", height = 550), tableOutput("feature_metadata")),
        tabPanel("Annotations", h4("Candidate records and distinct features"), tableOutput("annotation_counts"),
                 textOutput("annotation_limit"), tableOutput("annotation_table"),
                 h4("Selected candidate evidence"), tableOutput("evidence")),
        tabPanel("Provenance & downloads", verbatimTextOutput("manifest"),
                 selectInput("download_table", "Existing saved table", choices = character()),
                 downloadButton("download", "Download verified TSV"))
      ), width = 9
    )
  )
)

server <- function(input, output, session) {
  catalogue <- reactiveVal(catalogue_runs())
  observeEvent(input$refresh, catalogue(catalogue_runs()))
  observe({
    listing <- catalogue()
    previous <- isolate(input$run)
    updateSelectInput(session, "run", choices = stats::setNames(listing$path, listing$label),
                      selected = if (previous %in% listing$path) previous else if (nrow(listing)) listing$path[1] else character())
  })
  run <- reactive({
    listing <- catalogue()
    req(input$run, input$run %in% listing$path)
    tryCatch(open_run(input$run), error = function(error) validate(need(FALSE, conditionMessage(error))))
  })
  table_data <- function(name, optional = FALSE) {
    tryCatch(saved_table(run(), name, optional), error = function(error) validate(need(FALSE, conditionMessage(error))))
  }
  metadata <- reactive(table_data("sample_metadata"))
  variables <- reactive(table_data("variable_metadata"))
  annotations <- reactive({
    result <- table_data("annotation_candidates", optional = TRUE)
    if (is.null(result)) return(data.frame(feature_id = character(), source = character()))
    result
  })

  output$status <- renderUI({
    if (is.null(root)) return(tags$p(class = "text-danger", "Administrator: set MAPP_STATS_ROOT to an existing batch directory or a directory containing MAPP project/batch trees."))
    if (!nrow(catalogue())) return(tags$p("No completed V2 runs under the approved root. Use Refresh catalogue after running the CLI."))
    value <- run()
    tags$p(class = "text-success", paste("Completed", value$manifest$dataset_id, "run", basename(value$path), "(saved results only)"))
  })
  output$controls <- renderUI({
    value <- run(); samples <- metadata(); features <- variables()
    group_default <- value$manifest$effective_recipe$design$group
    categorical <- names(samples)[vapply(samples, function(column) {
      count <- length(unique(stats::na.omit(column)))
      count > 0L && count <= 30L && !is.list(column)
    }, logical(1))]
    categorical <- setdiff(categorical, c("filename", "sample_id"))
    role_col <- value$manifest$effective_recipe$roles$column
    if (is.null(role_col) || !role_col %in% names(samples)) role_col <- "sample_type"
    groups <- unique(c(intersect(group_default, categorical), categorical))
    stages <- sub("^matrix_(.*)[.]tsv$", "\\1", grep("^matrix_.*[.]tsv$", available_tables(value), value = TRUE))
    contrasts <- table_data("differential", optional = TRUE)
    feature_ids <- if ("feature_id" %in% names(features)) as.character(features$feature_id) else character()
    tagList(
      selectInput("group", "Descriptive grouping / point colour", choices = groups),
      if (role_col %in% names(samples)) selectInput("role", "Sample role filter", choices = c("All roles", sort(unique(as.character(stats::na.omit(samples[[role_col]]))))), selected = "All roles"),
      uiOutput("group_levels"),
      hr(), selectInput("stage", "Saved intensity stage", choices = stages, selected = if ("transformed" %in% stages) "transformed" else stages[1]),
      selectizeInput("feature_id", "Feature ID", choices = feature_ids, options = list(placeholder = "Search feature ID")),
      if (!is.null(contrasts) && "contrast" %in% names(contrasts)) selectInput("contrast", "Saved planned contrast", choices = unique(contrasts$contrast)),
      hr(), textInput("annotation_query", "Find annotation (feature, label, InChIKey)", ""),
      selectInput("annotation_source", "Candidate source", choices = c("All sources", sort(unique(stats::na.omit(annotations()$source))))),
      uiOutput("classification_filter"),
      uiOutput("candidate_selection")
    )
  })
  output$group_levels <- renderUI({
    samples <- metadata(); req(input$group, input$group %in% names(samples))
    levels <- sort(unique(as.character(stats::na.omit(samples[[input$group]]))))
    selectInput("levels", "Groups to show", choices = levels, selected = levels, multiple = TRUE)
  })
  output$classification_filter <- renderUI({
    rows <- annotations()
    if (!all(c("npc_pathway", "npc_superclass", "npc_class") %in% names(rows))) return(NULL)
    tagList(lapply(c(pathway = "npc_pathway", superclass = "npc_superclass", class = "npc_class"), function(column) {
      selectInput(paste0("filter_", column), paste("NPC", sub("npc_", "", column)),
                  choices = c("All", sort(unique(stats::na.omit(rows[[column]])))))
    }))
  })
  sample_subset <- reactive({
    samples <- metadata()
    req(input$group, input$group %in% names(samples))
    role_col <- run()$manifest$effective_recipe$roles$column
    if (is.null(role_col) || !role_col %in% names(samples)) role_col <- "sample_type"
    if (role_col %in% names(samples) && !is.null(input$role) && input$role != "All roles")
      samples <- samples[!is.na(samples[[role_col]]) & samples[[role_col]] == input$role, , drop = FALSE]
    samples <- samples[!is.na(samples[[input$group]]) & as.character(samples[[input$group]]) %in% input$levels, , drop = FALSE]
    samples
  })
  sample_key <- function(samples) if ("filename" %in% names(samples)) as.character(samples$filename) else as.character(samples$sample)
  output$summary <- renderTable({
    value <- run(); annotated <- annotations(); validation <- value$manifest$validation
    data.frame(metric = c("Samples in dataset", "Features in dataset", "Selected samples", "Annotation candidate records", "Unique annotated dataset features", "Unmatched candidate records"),
               value = c(validation$dimensions$samples, validation$dimensions$features, nrow(sample_subset()),
                         sum(!is.na(annotated$source)),
                         length(unique(annotated$feature_id[!is.na(annotated$source) & annotated$feature_id %in% variables()$feature_id])),
                         if ("feature_in_dataset" %in% names(annotated)) sum(!is.na(annotated$source) & !annotated$feature_in_dataset) else 0))
  })
  output$group_counts <- renderTable({
    rows <- sample_subset(); req(nrow(rows))
    data.frame(group = names(table(rows[[input$group]], useNA = "ifany")), samples = as.integer(table(rows[[input$group]], useNA = "ifany")))
  })
  ordination <- function(stem, x, y) {
    rows <- table_data(paste0(stem, "_scores"), optional = TRUE)
    validate(need(!is.null(rows) && all(c("sample", x) %in% names(rows)), paste("No saved", toupper(stem), "scores for this run.")))
    samples <- sample_subset()
    index <- match(rows$sample, sample_key(samples))
    rows <- rows[!is.na(index), , drop = FALSE]
    index <- index[!is.na(index)]
    validate(need(nrow(rows), "No saved scores match the selected sample filters."))
    rows$selected_group <- as.character(samples[[input$group]][index])
    variance <- table_data(paste0(stem, "_variance"), optional = TRUE)
    axis_label <- function(axis) {
      percent <- if (!is.null(variance) && all(c("component", "variance_percent") %in% names(variance)))
        variance$variance_percent[match(axis, variance$component)] else NA_real_
      if (length(percent) && is.finite(percent)) sprintf("%s (%.1f%%)", axis, percent) else axis
    }
    if (!y %in% names(rows)) {
      rows$one_dimension <- 0
      return(ggplot2::ggplot(rows, ggplot2::aes(x = .data[[x]], y = one_dimension, colour = selected_group)) +
               ggplot2::geom_jitter(height = 0.08, width = 0, size = 3) + ggplot2::theme_classic() +
               ggplot2::labs(x = axis_label(x), y = "One-dimensional ordination", colour = input$group,
                             subtitle = "Only one component was saved; sample filters do not refit the analysis"))
    }
    ggplot2::ggplot(rows, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = selected_group)) +
      ggplot2::geom_point(size = 3) + ggplot2::theme_classic() +
      ggplot2::labs(x = axis_label(x), y = axis_label(y), colour = input$group,
                    subtitle = "Saved ordination scores; sample filters only")
  }
  output$pca <- renderPlot(ordination("pca", "PC1", "PC2"))
  output$pcoa <- renderPlot(ordination("pcoa", "PCoA1", "PCoA2"))
  output$pca_note <- renderText("PCA is a saved analysis; changing grouping or sample filters does not refit it.")
  output$pcoa_note <- renderText("PCoA is a saved analysis; changing grouping or sample filters does not recompute distances.")
  plsda_status <- reactive({
    status <- table_data("plsda_status", optional = TRUE)
    if (is.null(status) || !nrow(status)) return(list(status = "disabled", reason = "PLS-DA was not enabled for this run."))
    validate(need(all(c("status", "reason") %in% names(status)), "Saved PLS-DA status is incomplete."))
    list(status = as.character(status$status[1]), reason = as.character(status$reason[1]))
  })
  validated_plsda <- reactive({
    state <- plsda_status()
    validate(need(identical(state$status, "validated"), paste("No validated PLS-DA:", state$reason)))
    summary <- table_data("plsda_summary", optional = TRUE)
    validate(need(!is.null(summary) && nrow(summary) && all(c("balanced_accuracy", "null_mean_balanced_accuracy", "permutation_p_value") %in% names(summary)),
                  "Saved PLS-DA validation summary is incomplete."))
    summary
  })
  output$plsda_interpretation <- renderUI({
    state <- plsda_status()
    if (!identical(state$status, "validated"))
      return(tags$p(class = "text-warning", paste("No validated supervised model:", state$reason)))
    summary <- validated_plsda()
    score <- as.numeric(summary$balanced_accuracy[1])
    null <- as.numeric(summary$null_mean_balanced_accuracy[1])
    p <- as.numeric(summary$permutation_p_value[1])
    supported <- is.finite(score) && is.finite(null) && is.finite(p) && p < 0.05 && score > null + 0.05
    tags$p(class = if (supported) "text-success" else "text-warning",
           if (supported) "Held-out prediction exceeds the permutation baseline; interpret alongside the saved design and sample size."
           else "Exploratory / unsupported discrimination: held-out accuracy is within 0.05 of the permuted mean or the permutation result is not significant. Full-fit score separation alone is not evidence.")
  })
  output$plsda_validation <- renderTable({
    summary <- validated_plsda()
    keys <- intersect(c("model", "validation", "selection", "preprocessing", "samples", "groups", "folds", "repeats",
                        "permutations", "full_fit_components", "accuracy", "balanced_accuracy", "null_mean_balanced_accuracy",
                        "null_sd_balanced_accuracy", "permutation_p_value"), names(summary))
    data.frame(metric = keys, saved_value = vapply(summary[keys], function(column) as.character(column[1]), character(1)))
  })
  output$plsda_predictions <- renderTable({
    validated_plsda()
    rows <- table_data("plsda_fold_predictions", optional = TRUE)
    if (!is.null(rows)) rows[seq_len(min(100L, nrow(rows))), , drop = FALSE]
  })
  output$plsda_permutations <- renderTable({
    validated_plsda()
    rows <- table_data("plsda_permutations", optional = TRUE)
    if (!is.null(rows)) rows[seq_len(min(100L, nrow(rows))), , drop = FALSE]
  })
  output$plsda_null <- renderPlot({
    summary <- validated_plsda()
    rows <- table_data("plsda_permutations", optional = TRUE)
    validate(need(!is.null(rows) && "balanced_accuracy" %in% names(rows) && nrow(rows), "No saved label permutations."))
    ggplot2::ggplot(rows, ggplot2::aes(x = balanced_accuracy)) +
      ggplot2::geom_histogram(bins = min(20L, max(5L, nrow(rows) %/% 2L)), fill = "grey65", colour = "white") +
      ggplot2::geom_vline(xintercept = summary$balanced_accuracy[1], colour = "#C23B22", linewidth = 1) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Balanced accuracy from saved label permutations", y = "Permutation count",
                    subtitle = "Red line: observed held-out balanced accuracy; see saved permutation p-value above")
  })
  output$plsda_scores <- renderPlot({
    validated_plsda()
    rows <- table_data("plsda_scores", optional = TRUE)
    validate(need(!is.null(rows) && all(c("sample", "LV1") %in% names(rows)), "Full-fit scores are unavailable."))
    samples <- sample_subset()
    index <- match(rows$sample, sample_key(samples))
    rows <- rows[!is.na(index), , drop = FALSE]
    index <- index[!is.na(index)]
    validate(need(nrow(rows), "No full-fit scores match the selected sample filters."))
    rows$selected_group <- as.character(samples[[input$group]][index])
    one_component <- !"LV2" %in% names(rows)
    if (one_component) rows$LV2 <- 0
    ggplot2::ggplot(rows, ggplot2::aes(x = LV1, y = LV2, colour = selected_group)) +
      (if (one_component) ggplot2::geom_jitter(height = 0.08, width = 0, size = 3) else ggplot2::geom_point(size = 3)) +
      ggplot2::theme_classic() +
      ggplot2::labs(y = if (one_component) "One component saved" else "LV2", colour = input$group,
                    subtitle = "Full-fit scores illustrate the model; held-out and permutation results determine evidential support")
  })
  selected_volcano <- reactive({
    rows <- table_data("differential", optional = TRUE)
    validate(need(!is.null(rows) && "contrast" %in% names(rows) && !is.null(input$contrast), "No saved planned contrasts for this run."))
    rows <- rows[!is.na(rows$contrast) & rows$contrast == input$contrast, , drop = FALSE]
    validate(need(nrow(rows), "No saved results for this planned contrast."))
    rows
  })
  output$volcano <- renderPlot({
    rows <- selected_volcano()
    rows$minus_log10_p <- -log10(pmax(rows$p_value, .Machine$double.xmin))
    rows$significant <- !is.na(rows$q_value) & rows$q_value < 0.05
    label <- unique(rows$effect_scale)
    ggplot2::ggplot(rows, ggplot2::aes(x = effect, y = minus_log10_p, colour = significant)) +
      ggplot2::geom_point(alpha = 0.65) + ggplot2::scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#C23B22")) +
      ggplot2::theme_classic() + ggplot2::labs(title = input$contrast, x = paste(label, collapse = ", "), y = "-log10(saved p-value)", colour = "Saved BH q < 0.05")
  })
  output$volcano_table <- renderTable({
    rows <- selected_volcano()
    rows[order(rows$q_value, na.last = TRUE), intersect(c("feature_id", "effect", "effect_scale", "p_value", "q_value"), names(rows)), drop = FALSE][seq_len(min(20L, nrow(rows))), , drop = FALSE]
  })
  observeEvent(input$volcano_click, {
    rows <- selected_volcano()
    rows$minus_log10_p <- -log10(pmax(rows$p_value, .Machine$double.xmin))
    closest <- shiny::nearPoints(rows, input$volcano_click, xvar = "effect", yvar = "minus_log10_p", maxpoints = 1)
    if (nrow(closest)) updateSelectizeInput(session, "feature_id", selected = as.character(closest$feature_id[1]))
  })
  output$feature <- renderPlot({
    req(input$stage, input$feature_id, input$group)
    matrix <- table_data(paste0("matrix_", input$stage))
    validate(need(input$feature_id %in% names(matrix) && "sample" %in% names(matrix), "Feature absent from this saved stage."))
    samples <- sample_subset()
    index <- match(matrix$sample, sample_key(samples))
    selected <- !is.na(index)
    validate(need(any(selected), "No samples match the selected filters."))
    rows <- data.frame(group = as.character(samples[[input$group]][index[selected]]), intensity = as.numeric(matrix[[input$feature_id]][selected]))
    ggplot2::ggplot(rows, ggplot2::aes(x = group, y = intensity, colour = group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) + ggplot2::geom_jitter(width = 0.12, na.rm = TRUE) +
      ggplot2::theme_classic() + ggplot2::labs(title = paste("Feature", input$feature_id),
                     subtitle = paste("Saved", input$stage, "intensities; descriptive plot only"), x = input$group)
  })
  output$feature_metadata <- renderTable({
    req(input$feature_id)
    rows <- variables()
    if (!"feature_id" %in% names(rows)) return(NULL)
    rows[rows$feature_id == input$feature_id, , drop = FALSE]
  })
  filtered_annotations <- reactive({
    rows <- annotations()
    if (!nrow(rows)) return(rows)
    if (!is.null(input$annotation_source) && input$annotation_source != "All sources") rows <- rows[!is.na(rows$source) & rows$source == input$annotation_source, , drop = FALSE]
    for (column in c("npc_pathway", "npc_superclass", "npc_class")) {
      selected <- input[[paste0("filter_", column)]]
      if (column %in% names(rows) && !is.null(selected) && selected != "All")
        rows <- rows[!is.na(rows[[column]]) & rows[[column]] == selected, , drop = FALSE]
    }
    query <- tolower(trimws(if (is.null(input$annotation_query)) "" else input$annotation_query))
    if (nzchar(query)) {
      fields <- intersect(c("feature_id", "label", "inchikey", "npc_pathway", "npc_superclass", "npc_class", "component_id"), names(rows))
      matches <- Reduce(`|`, lapply(rows[fields], function(column) grepl(query, tolower(as.character(column)), fixed = TRUE)), init = rep(FALSE, nrow(rows)))
      rows <- rows[matches, , drop = FALSE]
    }
    rows
  })
  output$annotation_counts <- renderTable({
    rows <- filtered_annotations()
    in_dataset <- if ("feature_in_dataset" %in% names(rows)) !is.na(rows$feature_in_dataset) & rows$feature_in_dataset
                  else rows$feature_id %in% variables()$feature_id
    candidate <- !is.na(rows$source)
    data.frame(metric = c("Candidate records (matching filters)", "Distinct dataset features with candidates",
                          "Distinct unmatched feature IDs", "Unmatched candidate records"),
               value = c(sum(candidate),
                         length(unique(rows$feature_id[candidate & in_dataset & !is.na(rows$feature_id)])),
                         length(unique(rows$feature_id[candidate & !in_dataset & !is.na(rows$feature_id)])),
                         sum(candidate & !in_dataset)))
  })
  output$candidate_selection <- renderUI({
    rows <- filtered_annotations()
    if (!"candidate_id" %in% names(rows) || !any(!is.na(rows$candidate_id))) return(NULL)
    if (!nrow(rows)) return(NULL)
    labels <- paste(rows$source, rows$feature_id, ifelse(is.na(rows$label), "unlabelled", rows$label), rows$candidate_id, sep = " | ")
    selectizeInput("candidate", "Candidate evidence (search any matching candidate)", choices = stats::setNames(rows$candidate_id, labels),
                   options = list(maxOptions = 100))
  })
  observeEvent(input$candidate, {
    rows <- annotations()
    if (!"candidate_id" %in% names(rows)) return()
    feature <- rows$feature_id[!is.na(rows$candidate_id) & rows$candidate_id == input$candidate]
    if (length(feature) && feature[1] %in% variables()$feature_id)
      updateSelectizeInput(session, "feature_id", selected = as.character(feature[1]))
  })
  output$annotation_limit <- renderText({
    count <- nrow(filtered_annotations())
    if (count > 100L) paste("Showing the first 100 of", count, "rows. Narrow the filters or download the complete saved table.")
  })
  output$annotation_table <- renderTable({
    rows <- filtered_annotations()
    if (!nrow(rows)) return(data.frame(message = "No annotation candidates match the filters (annotations may be optional)."))
    fields <- intersect(c("feature_id", "source", "label", "confidence", "confidence_metric", "npc_pathway", "npc_superclass", "npc_class", "component_id", "feature_in_dataset"), names(rows))
    rows[seq_len(min(100L, nrow(rows))), fields, drop = FALSE]
  }, striped = TRUE, hover = TRUE)
  output$evidence <- renderTable({
    req(input$candidate)
    rows <- annotations()
    if (!"candidate_id" %in% names(rows)) return(NULL)
    record <- rows[!is.na(rows$candidate_id) & rows$candidate_id == input$candidate, , drop = FALSE]
    if (!nrow(record)) return(NULL)
    data.frame(field = names(record), value = vapply(record, function(column) as.character(column[1]), character(1)), row.names = NULL)
  })
  output$manifest <- renderPrint({
    manifest <- run()$manifest
    inputs <- manifest$inputs
    if (length(inputs)) {
      labels <- names(inputs)
      bad <- is.na(labels) | !nzchar(labels) | grepl("/", labels, fixed = TRUE) | grepl("\\", labels, fixed = TRUE)
      labels[bad] <- paste0("input_", which(bad))
      names(inputs) <- labels
    }
    str(list(dataset_id = manifest$dataset_id, created_at = manifest$created_at,
             hashes = manifest$hashes, input_checksums = inputs,
             recipe = manifest$effective_recipe, validation = manifest$validation,
             software = manifest$environment[c("schema_version", "git", "code_sha256", "renv_lock_sha256", "packages", "r_version")],
             output_checksums = manifest$outputs), max.level = 5)
  })
  observe({
    value <- run()
    files <- available_tables(value)
    updateSelectInput(session, "download_table", choices = files, selected = if (length(files)) files[1] else character())
  })
  output$download <- downloadHandler(
    filename = function() { req(input$download_table); basename(input$download_table) },
    content = function(file) {
      selected <- input$download_table
      req(selected, selected %in% available_tables(run()))
      source <- tryCatch(verified_file(run(), paste0("tables/", selected)), error = function(error) validate(need(FALSE, conditionMessage(error))))
      if (!file.copy(source, file, overwrite = TRUE)) stop("Unable to copy verified table.", call. = FALSE)
    }
  )
}

shinyApp(ui, server)
