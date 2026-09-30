if (!requireNamespace("shiny", quietly = TRUE)) stop("The V2 explorer requires shiny.")
if (!requireNamespace("htmltools", quietly = TRUE)) stop("The V2 explorer requires htmltools.")
if (!requireNamespace("yaml", quietly = TRUE)) stop("The V2 explorer requires yaml.")
if (!requireNamespace("digest", quietly = TRUE)) stop("The V2 explorer requires digest.")
if (!requireNamespace("plotly", quietly = TRUE)) stop("The V2 explorer requires plotly.")

library(shiny)

# The administrator supplies a root; visitors never supply filesystem paths.
configured_root <- Sys.getenv("MAPP_STATS_ROOT", unset = "")
root <- if (nzchar(configured_root) && dir.exists(configured_root)) normalizePath(configured_root, mustWork = TRUE) else NULL

within_root <- function(path) {
  resolved <- normalizePath(path, mustWork = TRUE)
  identical(resolved, root) || startsWith(resolved, paste0(root, .Platform$file.sep))
}

safe_path <- function(run_dir, relative) {
  if (is.null(root) || !grepl("^(objects|tables|plots)/[A-Za-z0-9_.-]+[.](rds|tsv|png|pdf)$", relative))
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
  if (length(paths)) {
    newest_first <- order(file.info(file.path(paths, "COMPLETE"))$mtime, decreasing = TRUE, na.last = TRUE)
    paths <- paths[newest_first]
    labels <- labels[newest_first]
    labels[1] <- paste0(labels[1], " (latest)")
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
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"), quote = "\"", comment.char = "")
}

available_tables <- function(run) {
  folder <- file.path(run$path, "tables")
  if (!dir.exists(folder) || !within_root(folder) || nzchar(Sys.readlink(folder))) return(character())
  names <- list.files(folder, pattern = "^[A-Za-z0-9_.-]+[.]tsv$", full.names = FALSE)
  names[vapply(names, function(name) !is.null(run$manifest$outputs[[paste0("tables/", name)]]), logical(1))]
}

# Browser point labels never expose source filenames or other path-like metadata.
hover_value <- function(value) {
  value <- as.character(value)
  if (length(value) != 1L || is.na(value) || !nzchar(value) || grepl("[/\\\\]", value)) return(NULL)
  as.character(htmltools::htmlEscape(value))
}

sample_hover <- function(samples, index, extra) {
  vapply(seq_along(index), function(i) {
    row <- index[i]
    identifier <- if ("sample_id" %in% names(samples)) samples$sample_id[row] else
      if ("sample" %in% names(samples)) samples$sample[row] else samples$filename[row]
    id <- hover_value(basename(as.character(identifier)))
    if (is.null(id)) id <- "(unavailable)"
    fields <- c(paste0("Sample: ", id), extra[i])
    columns <- setdiff(names(samples), c("sample_id", "sample", "filename"))
    columns <- columns[!grepl("path|file|dir|uri|url|[/\\\\]", columns, ignore.case = TRUE)]
    for (column in columns) {
      value <- hover_value(samples[[column]][row])
      if (!is.null(value))
        fields <- c(fields, paste0(htmltools::htmlEscape(column), ": ", value))
    }
    paste(fields, collapse = "<br>")
  }, character(1))
}

plot_controls <- function(plot, filename) {
  plotly::config(plot, displayModeBar = TRUE, displaylogo = FALSE,
                 toImageButtonOptions = list(format = "png", filename = filename, width = 1100, height = 750, scale = 2))
}

saved_axes <- function(rows, prefix) {
  if (is.null(rows)) return(character())
  axes <- grep(paste0("^", prefix, "[1-9][0-9]*$"), names(rows), value = TRUE)
  axes <- axes[vapply(rows[axes], is.numeric, logical(1))]
  axes[order(as.integer(sub(prefix, "", axes, fixed = TRUE)))]
}

score_controls <- function(stem, axes) {
  if (!length(axes)) return(NULL)
  tagList(
    if (length(axes) >= 3L) radioButtons(paste0(stem, "_mode"), "Plot view", choices = c("2D", "3D"), inline = TRUE),
    fluidRow(
      column(4, selectInput(paste0(stem, "_x"), "X axis", choices = axes, selected = axes[1])),
      if (length(axes) >= 2L) column(4, selectInput(paste0(stem, "_y"), "Y axis", choices = axes, selected = axes[2])),
      if (length(axes) >= 3L) column(4, conditionalPanel(
        paste0("input.", stem, "_mode === '3D'"),
        selectInput(paste0(stem, "_z"), "Z axis", choices = axes, selected = axes[3])
      ))
    )
  )
}

selected_axes <- function(input, stem, axes) {
  shiny::validate(shiny::need(length(axes), "No saved components are available."))
  is_3d <- length(axes) >= 3L && identical(input[[paste0(stem, "_mode")]], "3D")
  selected <- c(input[[paste0(stem, "_x")]],
                if (length(axes) >= 2L) input[[paste0(stem, "_y")]],
                if (is_3d) input[[paste0(stem, "_z")]])
  shiny::req(length(selected) == if (is_3d) 3L else if (length(axes) >= 2L) 2L else 1L)
  shiny::validate(shiny::need(all(selected %in% axes) && length(unique(selected)) == length(selected),
                             "Choose different saved components for each axis."))
  list(x = selected[1], y = if (length(axes) >= 2L) selected[2] else NULL,
       z = if (is_3d) selected[3] else NULL)
}

group_scatter <- function(rows, x, y, x_label, y_label, subtitle, filename, z = NULL, z_label = NULL) {
  groups <- sort(unique(rows$selected_group))
  colors <- stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
  plot <- plotly::plot_ly()
  for (group in groups) {
    part <- rows[rows$selected_group == group, , drop = FALSE]
    name <- as.character(htmltools::htmlEscape(group))
    if (is.null(z)) {
      plot <- plotly::add_markers(plot, x = part[[x]], y = part[[y]], name = name,
                                  text = part$hover, hoverinfo = "text",
                                  marker = list(size = 9, color = colors[[group]], opacity = 0.8))
    } else {
      plot <- plotly::add_trace(plot, type = "scatter3d", mode = "markers",
                                x = part[[x]], y = part[[y]], z = part[[z]], name = name,
                                text = part$hover, hoverinfo = "text",
                                marker = list(size = 5, color = colors[[group]], opacity = 0.8))
    }
  }
  annotation <- list(list(text = subtitle, xref = "paper", yref = "paper",
                          x = 0, y = 1.12, showarrow = FALSE, xanchor = "left"))
  legend <- list(itemclick = "toggle", itemdoubleclick = "toggleothers")
  if (!is.null(z)) {
    return(plot_controls(plotly::layout(plot,
                                       scene = list(xaxis = list(title = x_label),
                                                    yaxis = list(title = y_label),
                                                    zaxis = list(title = z_label)),
                                       annotations = annotation, legend = legend), filename))
  }
  plot_controls(plotly::layout(plot, xaxis = list(title = x_label), yaxis = list(title = y_label),
                               annotations = annotation, legend = legend, hovermode = "closest"), filename)
}

wikidata_id <- function(value) {
  if (length(value) != 1L || is.na(value)) return(NULL)
  value <- trimws(value)
  if (grepl("^Q[1-9][0-9]*$", value)) return(value)
  if (!grepl("^https?://(www[.])?wikidata[.]org/(entity|wiki)/Q[1-9][0-9]*$", value))
    return(NULL)
  sub("^.*/(Q[1-9][0-9]*)$", "\\1", value)
}

wikidata_anchor <- function(value) {
  id <- wikidata_id(value)
  if (is.null(id)) return(NULL)
  tags$a(href = paste0("https://www.wikidata.org/wiki/", id),
         target = "_blank", rel = "noopener noreferrer", paste("Wikidata", id))
}

candidate_field <- function(rows, field) {
  if (field %in% names(rows)) as.character(rows[[field]]) else rep(NA_character_, nrow(rows))
}

feature_table <- function(rows) {
  if (is.null(rows) || !nrow(rows)) return(NULL)
  tags$div(class = "table-responsive",
    tags$table(class = "table table-striped table-hover",
      tags$thead(tags$tr(lapply(names(rows), tags$th))),
      tags$tbody(lapply(seq_len(nrow(rows)), function(i) {
        tags$tr(lapply(seq_along(rows), function(j) {
          value <- as.character(rows[[j]][i])
          if (is.na(value)) value <- ""
          linked <- identical(names(rows)[j], "feature_id") ||
            (identical(names(rows)[j], "value") && "field" %in% names(rows) &&
             identical(as.character(rows$field[i]), "feature_id"))
          wikidata_field <- names(rows)[j] %in% c("identification_links", "organism_wikidata") ||
            (identical(names(rows)[j], "value") && "field" %in% names(rows) &&
             as.character(rows$field[i]) %in% c("identification_links", "organism_wikidata"))
          wikidata <- if (wikidata_field) wikidata_anchor(value) else NULL
          tags$td(if (linked && nzchar(value))
            tags$a(href = "#feature_drawer", class = "feature-link", `data-feature-id` = value, value)
            else if (!is.null(wikidata)) wikidata else value)
        }))
      }))
    )
  )
}

structure_cards <- function(rows) {
  rows <- rows[!is.na(rows$smiles) & nzchar(trimws(rows$smiles)) &
                 !tolower(trimws(rows$smiles)) %in% c("na", "none", "null"), , drop = FALSE]
  if (!nrow(rows)) return(tags$p("No saved SMILES for this feature."))
  rows$smiles <- trimws(rows$smiles)
  # Collapse identical saved strings, not chemically equivalent encodings or stereoisomers.
  tags$div(class = "structure-grid", lapply(unique(rows$smiles), function(smiles) {
    group <- rows[rows$smiles == smiles, , drop = FALSE]
    labels <- unique(group$label)
    taxa <- unique(group[!is.na(group$organism) & nzchar(group$organism),
                         c("organism", "organism_wikidata"), drop = FALSE])
    links <- Filter(Negate(is.null), lapply(unique(group$identification_links), wikidata_anchor))
    tags$figure(class = "structure-card",
      tags$figcaption(paste(labels, collapse = " · ")),
      if (nchar(smiles, type = "bytes") <= 5000L)
        tags$img(class = "structure-depiction", loading = "lazy", referrerpolicy = "no-referrer",
                 src = paste0("https://api.naturalproducts.net/latest/depict/2D?smiles=",
                              utils::URLencode(smiles, reserved = TRUE),
                              "&toolkit=rdkit&width=440&height=300"),
                 alt = paste("2D depiction of", paste(labels, collapse = ", ")))
      else tags$p("SMILES exceeds the depiction service's 5000-character limit."),
      if (nrow(taxa)) tags$div(class = "structure-taxa", tags$strong("Reported taxon (met-annot-enhancer):"),
        tags$ul(lapply(seq_len(nrow(taxa)), function(i) {
          link <- wikidata_anchor(taxa$organism_wikidata[i])
          tags$li(taxa$organism[i], if (!is.null(link)) tagList(" — ", link))
        }))),
      if (length(links)) tags$div(class = "structure-links", tags$strong("Structure Wikidata:"),
                                 tags$ul(lapply(links, tags$li))),
      tags$p(class = "structure-fallback", "Depiction unavailable; the saved SMILES is below."),
      tags$details(tags$summary("Saved SMILES"), tags$code(smiles))
    )
  }))
}

ui <- fluidPage(
  tags$head(
    tags$style(htmltools::HTML("
      :root {--feature-drawer-width: clamp(320px, 36vw, 620px);}
      #explorer_workspace {transition: margin-right .2s ease;}
      body.feature-drawer-open #explorer_workspace {margin-right: var(--feature-drawer-width);}
      #feature_drawer {position: fixed; top: 0; right: 0; bottom: 0; width: var(--feature-drawer-width);
        z-index: 1100; background: #fff; border-left: 1px solid #ccc; box-shadow: -4px 0 16px #0002;
        overflow-y: auto; padding: 18px; transform: translateX(105%); transition: transform .2s ease;}
      #feature_drawer.is-open {transform: translateX(0);}
      #feature_drawer_toggle {position: fixed; top: 40%; right: 0; z-index: 1101; writing-mode: vertical-rl;
        padding: 12px 8px; border-radius: 5px 0 0 5px; background: #22668a; color: white;}
      #feature_drawer.is-open ~ #feature_drawer_toggle {display: none;}
      .structure-grid {display: grid; grid-template-columns: repeat(auto-fit, minmax(min(250px, 100%), 1fr)); gap: 12px;}
      .structure-card {margin: 0; padding: 10px; border: 1px solid #ddd; min-width: 0;}
      .structure-card figcaption {font-weight: bold; overflow-wrap: anywhere;}
      .structure-depiction {display: block; width: 100%; height: 210px; object-fit: contain;}
      .structure-card code {white-space: normal; overflow-wrap: anywhere;}
      .structure-fallback {display: none;}
      .structure-card.is-unavailable .structure-fallback {display: block;}
      .structure-card.is-unavailable .structure-depiction {display: none;}
      @media (prefers-reduced-motion: reduce) {
        #explorer_workspace, #feature_drawer {transition: none;}
      }
    ")),
    tags$script(htmltools::HTML("
      $(document).on('shiny:connected', function() {
        Shiny.addCustomMessageHandler('featureDrawer', function(open) {
          document.getElementById('feature_drawer').classList.toggle('is-open', open);
          document.body.classList.toggle('feature-drawer-open', open);
          document.getElementById('feature_drawer_toggle').setAttribute('aria-expanded', String(open));
          setTimeout(function() { window.dispatchEvent(new Event('resize')); }, 220);
        });
      });
      document.addEventListener('click', function(event) {
        var link = event.target.closest('a.feature-link');
        if (!link) return;
        event.preventDefault();
        Shiny.setInputValue('open_feature', link.dataset.featureId, {priority: 'event'});
      });
      document.addEventListener('error', function(event) {
        if (event.target.matches && event.target.matches('img.structure-depiction'))
          event.target.closest('.structure-card').classList.add('is-unavailable');
      }, true);
    "))
  ),
  titlePanel("MAPP statistics | Completed run explorer"),
  tags$div(id = "explorer_workspace", sidebarLayout(
    sidebarPanel(
      selectInput("run", "Completed batch / run", choices = character()),
      actionButton("refresh", "Refresh catalogue"),
      hr(), uiOutput("controls"), width = 3
    ),
    mainPanel(
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Overview", tableOutput("summary"), h4("Samples by selected group"), tableOutput("group_counts")),
        tabPanel("PCA", uiOutput("pca_axes"), plotly::plotlyOutput("pca", height = 550), uiOutput("pca_downloads"), textOutput("pca_note"),
                 h4("Feature loadings (select a point for details)"),
                 tags$p("These points are features on the saved loading scale, not the samples above."),
                 plotly::plotlyOutput("pca_loadings", height = 450)),
        tabPanel("PCoA", uiOutput("pcoa_axes"), plotly::plotlyOutput("pcoa", height = 550), uiOutput("pcoa_downloads"), textOutput("pcoa_note")),
        tabPanel("Saved contrasts", plotly::plotlyOutput("volcano", height = 550), uiOutput("volcano_downloads"),
                 uiOutput("volcano_table")),
        tabPanel("Validated PLS-DA", uiOutput("plsda_interpretation"), tableOutput("plsda_validation"),
                 h4("Held-out fold predictions"), tableOutput("plsda_predictions"),
                 h4("Saved label permutations"), plotly::plotlyOutput("plsda_null", height = 320), tableOutput("plsda_permutations"),
                 h4("Full-fit scores (illustration, not validation)"), uiOutput("plsda_axes"), plotly::plotlyOutput("plsda_scores", height = 500)),
        tabPanel("Annotations", uiOutput("consensus_section"),
                 h4("Candidate records and distinct features"), tableOutput("annotation_counts"),
                 tags$p("A feature ID identifies one measured feature. Source record keys (for example, sirius:402) identify annotation input rows, not additional features."),
                 textOutput("annotation_limit"), uiOutput("annotation_table"),
                 h4("Selected candidate evidence"), uiOutput("evidence"), uiOutput("evidence_structure")),
        tabPanel("Provenance & downloads", verbatimTextOutput("manifest"),
                 selectInput("download_table", "Existing saved table", choices = character()),
                 downloadButton("download", "Download verified TSV"))
      ), width = 9
    ),
  )),
  tags$aside(id = "feature_drawer", class = "feature-drawer", role = "complementary", `aria-label` = "Feature description",
             tags$div(class = "clearfix", actionButton("close_feature_drawer", "Close", class = "pull-right")),
             h3(textOutput("feature_heading", inline = TRUE)),
             tags$p(textOutput("feature_consensus", inline = TRUE)),
             tags$p("Saved descriptive intensities: boxes show the median, quartiles and whiskers; dots are individual samples. No statistical test is recomputed."),
             uiOutput("feature_plot_region"),
             h4("Measured feature metadata"), uiOutput("feature_metadata"),
             h4("Candidate annotations (not confirmed identities)"), textOutput("feature_candidate_count"),
             tags$p("All candidates below belong to this feature. Their source record keys identify annotation input rows, not other feature IDs."),
             uiOutput("feature_candidates"),
             h4("Source-specific structures (not confirmed identities)"),
             uiOutput("feature_structures")),
  actionButton("feature_drawer_toggle", "Feature details", class = "feature-drawer-tab",
               `aria-controls` = "feature_drawer", `aria-expanded` = "false")
)

server <- function(input, output, session) {
  catalogue <- reactiveVal(catalogue_runs())
  refresh_count <- reactiveVal(0L)
  observeEvent(input$refresh, {
    catalogue(catalogue_runs())
    refresh_count(refresh_count() + 1L)
  })
  observe({
    listing <- catalogue()
    refreshed <- refresh_count()
    previous <- isolate(input$run)
    selected <- if (!nrow(listing)) character() else if (refreshed == 0L && length(previous) == 1L && previous %in% listing$path) previous else listing$path[1]
    updateSelectInput(session, "run", choices = stats::setNames(listing$path, listing$label), selected = selected)
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
  input_variables <- reactive(table_data("variable_metadata_input", optional = TRUE))
  annotations <- reactive({
    result <- table_data("annotation_candidates", optional = TRUE)
    if (is.null(result)) return(data.frame(feature_id = character(), source = character()))
    result
  })
  horizontal_annotations <- reactive(table_data("annotation_horizontal", optional = TRUE))
  raw_matrix <- reactive(table_data("matrix_raw"))
  raw_feature_ids <- reactive(setdiff(names(raw_matrix()), "sample"))
  intensity_measure <- reactive({
    measure <- run()$manifest$effective_dataset$intensity_measure
    if (is.character(measure) && length(measure) == 1L && measure %in% c("height", "area"))
      paste("peak", measure) else "intensity"
  })
  selected_feature <- reactiveVal(NULL)
  drawer_open <- reactiveVal(FALSE)
  known_feature_ids <- reactive(unique(c(raw_feature_ids(), as.character(annotations()$feature_id),
                                          as.character(horizontal_annotations()$feature_id))))
  open_feature <- function(id) {
    if (length(id) != 1L || is.na(id) || !nzchar(id) || !id %in% known_feature_ids()) return()
    selected_feature(id)
    if (id %in% raw_feature_ids()) updateSelectizeInput(session, "feature_id", selected = id)
    drawer_open(TRUE)
  }
  sidebar_initialized <- reactiveVal(FALSE)
  observeEvent(input$feature_id, {
    if (!sidebar_initialized()) { sidebar_initialized(TRUE); return() }
    if (length(input$feature_id) == 1L && input$feature_id %in% raw_feature_ids()) {
      selected_feature(input$feature_id)
      drawer_open(TRUE)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$open_feature, open_feature(input$open_feature))
  observeEvent(input$candidate, {
    rows <- annotations()
    if (!"candidate_id" %in% names(rows)) return()
    feature <- rows$feature_id[!is.na(rows$candidate_id) & rows$candidate_id == input$candidate]
    if (length(feature)) open_feature(as.character(feature[1]))
  })
  observeEvent(input$feature_drawer_toggle, {
    if (is.null(selected_feature()) && length(input$feature_id) == 1L)
      selected_feature(input$feature_id)
    drawer_open(!isTRUE(drawer_open()))
  })
  observeEvent(input$close_feature_drawer, drawer_open(FALSE))
  observeEvent(input$run, {
    selected_feature(NULL)
    drawer_open(FALSE)
    sidebar_initialized(FALSE)
  }, ignoreInit = TRUE)
  observe(session$sendCustomMessage("featureDrawer", isTRUE(drawer_open())))
  output$feature_heading <- renderText(if (is.null(selected_feature())) "Choose a feature" else paste("Feature", selected_feature()))
  output$feature_consensus <- renderText({
    req(selected_feature())
    rows <- horizontal_annotations()
    if (is.null(rows)) return("No saved InChIKey2D consensus table in this run.")
    record <- rows[!is.na(rows$feature_id) & as.character(rows$feature_id) == selected_feature(), , drop = FALSE]
    if (!nrow(record)) return("No saved InChIKey2D consensus record for this feature.")
    counts <- unique(ifelse(is.na(record$sources_number_IK2D), "not recorded",
                            as.character(record$sources_number_IK2D)))
    paste("Saved InChIKey2D agreeing sources (per horizontal row):", paste(counts, collapse = ", "))
  })

  output$status <- renderUI({
    if (is.null(root)) return(tags$p(class = "text-danger", "Administrator: set MAPP_STATS_ROOT to an existing batch directory or a directory containing MAPP project/batch trees."))
    if (!nrow(catalogue())) return(tags$p("No completed V2 runs under the approved root. Use Refresh catalogue after running the CLI."))
    value <- run()
    tags$p(class = "text-success", paste("Completed", value$manifest$dataset_id, "run", basename(value$path), "(saved results only)"))
  })
  output$controls <- renderUI({
    value <- run(); samples <- metadata()
    group_default <- value$manifest$effective_recipe$design$group
    categorical <- names(samples)[vapply(samples, function(column) {
      count <- length(unique(stats::na.omit(column)))
      count > 0L && count <= 30L && !is.list(column) &&
        !any(grepl("[/\\\\]", as.character(stats::na.omit(column))))
    }, logical(1))]
    categorical <- setdiff(categorical[!grepl("path|file|dir|uri|url|^sample$", categorical, ignore.case = TRUE)], "sample_id")
    role_col <- value$manifest$effective_recipe$roles$column
    if (is.null(role_col) || !role_col %in% names(samples)) role_col <- "sample_type"
    groups <- unique(c(intersect(group_default, categorical), categorical))
    stages <- sub("^matrix_(.*)[.]tsv$", "\\1", grep("^matrix_.*[.]tsv$", available_tables(value), value = TRUE))
    contrasts <- table_data("differential", optional = TRUE)
    feature_ids <- raw_feature_ids()
    tagList(
      selectInput("group", "Descriptive grouping / point colour", choices = groups),
      if (role_col %in% names(samples)) selectInput("role", "Sample role filter", choices = c("All roles", sort(unique(as.character(stats::na.omit(samples[[role_col]]))))), selected = "All roles"),
      uiOutput("group_levels"),
      hr(), selectInput("stage", "Saved intensity stage", choices = stages, selected = if ("transformed" %in% stages) "transformed" else stages[1]),
      selectizeInput("feature_id", "Feature ID", choices = feature_ids, options = list(placeholder = "Search feature ID")),
      if (!is.null(contrasts) && "contrast" %in% names(contrasts)) selectInput("contrast", "Saved planned contrast", choices = unique(contrasts$contrast)),
      hr(), textInput("annotation_query", "Find annotation (feature, label, InChIKey)", ""),
      selectInput("annotation_source", "Candidate source", choices = c("All sources", sort(unique(stats::na.omit(annotations()$source))))),
      if (!is.null(horizontal_annotations()))
        tagList(selectInput("annotation_consensus", "Saved InChIKey2D consensus (sources agreeing)",
                            choices = c("All", stats::setNames(as.character(0:3), paste0(0:3, " sources")))),
                tags$small("Agreement count is not the number of saved SMILES or candidate records."))
      else tags$p("Consensus selection is unavailable: this completed run has no saved horizontal annotation table."),
      uiOutput("classification_filter"),
      uiOutput("candidate_selection")
    )
  })
  output$group_levels <- renderUI({
    samples <- metadata(); req(input$group, input$group %in% names(samples))
    levels <- sort(unique(as.character(stats::na.omit(samples[[input$group]]))))
    tagList(checkboxInput("all_groups", "All groups", value = TRUE),
            conditionalPanel("input.all_groups === false",
                             selectInput("levels", "Groups to show", choices = levels, selected = levels, multiple = TRUE)))
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
    if (!isTRUE(input$all_groups))
      samples <- samples[!is.na(samples[[input$group]]) & as.character(samples[[input$group]]) %in% input$levels, , drop = FALSE]
    else
      samples <- samples[!is.na(samples[[input$group]]), , drop = FALSE]
    samples
  })
  sample_key <- function(samples) if ("filename" %in% names(samples)) as.character(samples$filename) else as.character(samples$sample)
  output$summary <- renderTable({
    value <- run(); annotated <- annotations(); validation <- value$manifest$validation
    data.frame(metric = c("Samples in dataset", "Features in dataset", "Selected samples", "Annotation candidate records", "Unique annotated retained features", "Unmatched candidate records"),
               value = c(validation$dimensions$samples, validation$dimensions$features, nrow(sample_subset()),
                         sum(!is.na(annotated$source)),
                         length(unique(annotated$feature_id[!is.na(annotated$source) & annotated$feature_id %in% variables()$feature_id])),
                         if ("feature_in_dataset" %in% names(annotated)) sum(!is.na(annotated$source) & !annotated$feature_in_dataset) else 0))
  })
  output$group_counts <- renderTable({
    rows <- sample_subset(); req(nrow(rows))
    data.frame(group = names(table(rows[[input$group]], useNA = "ifany")), samples = as.integer(table(rows[[input$group]], useNA = "ifany")))
  })
  pca_scores <- reactive(table_data("pca_scores", optional = TRUE))
  pcoa_scores <- reactive(table_data("pcoa_scores", optional = TRUE))
  output$pca_axes <- renderUI(score_controls("pca", saved_axes(pca_scores(), "PC")))
  output$pcoa_axes <- renderUI(score_controls("pcoa", saved_axes(pcoa_scores(), "PCoA")))
  ordination <- function(stem, prefix, rows) {
    validate(need(!is.null(rows) && "sample" %in% names(rows) && length(saved_axes(rows, prefix)),
                  paste("No saved", toupper(stem), "scores for this run.")))
    axes <- selected_axes(input, stem, saved_axes(rows, prefix))
    x <- axes$x
    one_dimension <- is.null(axes$y)
    y <- if (one_dimension) ".ordination_jitter" else axes$y
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
    if (one_dimension) rows[[y]] <- ((seq_len(nrow(rows)) - 1L) %% 9L - 4L) * 0.02
    extra <- paste0(x, ": ", signif(rows[[x]], 4))
    if (!one_dimension) extra <- paste0(extra, "<br>", y, ": ", signif(rows[[y]], 4))
    if (!is.null(axes$z)) extra <- paste0(extra, "<br>", axes$z, ": ", signif(rows[[axes$z]], 4))
    rows$hover <- sample_hover(samples, index, extra)
    group_scatter(rows, x, y, axis_label(x),
                  if (one_dimension) "One-dimensional ordination" else axis_label(y),
                  if (one_dimension) "Only one component was saved; sample filters do not refit the analysis"
                  else "Saved ordination scores; sample filters only", stem,
                  z = axes$z, z_label = if (!is.null(axes$z)) axis_label(axes$z))
  }
  output$pca <- plotly::renderPlotly(ordination("pca", "PC", pca_scores()))
  output$pcoa <- plotly::renderPlotly(ordination("pcoa", "PCoA", pcoa_scores()))
  output$pca_loadings <- plotly::renderPlotly({
    rows <- table_data("pca_loadings", optional = TRUE)
    validate(need(!is.null(rows) && "feature_id" %in% names(rows) && nrow(rows),
                  "No saved PCA feature loadings for this run."))
    axes <- selected_axes(input, "pca", saved_axes(pca_scores(), "PC"))
    required <- c(axes$x, axes$y, axes$z)
    validate(need(all(required %in% names(rows)), "Selected PCA components have no saved feature loadings."))
    rows <- rows[!is.na(rows$feature_id) & as.character(rows$feature_id) %in% variables()$feature_id, , drop = FALSE]
    rows <- rows[stats::complete.cases(rows[required]), , drop = FALSE]
    validate(need(nrow(rows), "No retained features have finite saved PCA loadings."))
    y <- if (is.null(axes$y)) ".loading_jitter" else axes$y
    if (is.null(axes$y)) rows[[y]] <- ((seq_len(nrow(rows)) - 1L) %% 9L - 4L) * 0.02
    hover <- vapply(seq_len(nrow(rows)), function(i) {
      id <- hover_value(rows$feature_id[i])
      if (is.null(id)) id <- "(unavailable)"
      paste0("Feature: ", id, "<br>", axes$x, " loading: ", signif(rows[[axes$x]][i], 4),
             if (!is.null(axes$y)) paste0("<br>", axes$y, " loading: ", signif(rows[[axes$y]][i], 4)),
             if (!is.null(axes$z)) paste0("<br>", axes$z, " loading: ", signif(rows[[axes$z]][i], 4)))
    }, character(1))
    plot <- plotly::plot_ly(source = "pca_loadings")
    if (is.null(axes$z)) {
      plot <- plotly::add_markers(plot, x = rows[[axes$x]], y = rows[[y]], key = as.character(rows$feature_id),
                                  text = hover, hoverinfo = "text",
                                  marker = list(size = 7, color = "#2876a3", opacity = 0.7))
      plot <- plotly::layout(plot, xaxis = list(title = paste(axes$x, "loading")),
                             yaxis = list(title = if (is.null(axes$y)) "One-dimensional loadings" else paste(axes$y, "loading")),
                             showlegend = FALSE)
    } else {
      plot <- plotly::add_trace(plot, type = "scatter3d", mode = "markers",
                                x = rows[[axes$x]], y = rows[[y]], z = rows[[axes$z]],
                                key = as.character(rows$feature_id), text = hover, hoverinfo = "text",
                                marker = list(size = 4, color = "#2876a3", opacity = 0.7))
      plot <- plotly::layout(plot, scene = list(xaxis = list(title = paste(axes$x, "loading")),
                                                yaxis = list(title = paste(axes$y, "loading")),
                                                zaxis = list(title = paste(axes$z, "loading"))),
                             showlegend = FALSE)
    }
    plot_controls(plotly::event_register(plot, "plotly_click"), "pca_loadings")
  })
  observeEvent(plotly::event_data("plotly_click", source = "pca_loadings"), {
    clicked <- plotly::event_data("plotly_click", source = "pca_loadings")
    feature <- as.character(clicked$key)
    if (length(feature) == 1L) open_feature(feature)
  })
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
  output$plsda_null <- plotly::renderPlotly({
    summary <- validated_plsda()
    rows <- table_data("plsda_permutations", optional = TRUE)
    validate(need(!is.null(rows) && "balanced_accuracy" %in% names(rows) && nrow(rows), "No saved label permutations."))
    values <- as.numeric(rows$balanced_accuracy)
    values <- values[is.finite(values)]
    validate(need(length(values), "No finite saved label permutation accuracies."))
    observed <- as.numeric(summary$balanced_accuracy[1])
    plot <- plotly::plot_ly()
    plot <- plotly::add_histogram(plot, x = values, nbinsx = min(20L, max(5L, length(values) %/% 2L)),
                                  name = "Saved label permutations",
                                  marker = list(color = "grey65", line = list(color = "white", width = 1)),
                                  hovertemplate = "Balanced accuracy: %{x:.3f}<br>Permutation count: %{y}<extra></extra>")
    shapes <- if (is.finite(observed)) list(list(type = "line", x0 = observed, x1 = observed,
                                                y0 = 0, y1 = 1, yref = "paper",
                                                line = list(color = "#C23B22", width = 2))) else list()
    plot_controls(plotly::layout(plot, xaxis = list(title = "Balanced accuracy from saved label permutations"),
                                 yaxis = list(title = "Permutation count"), shapes = shapes,
                                 showlegend = FALSE,
                                 annotations = list(list(text = "Red line: observed held-out balanced accuracy; see saved permutation p-value above",
                                                         xref = "paper", yref = "paper", x = 0, y = 1.15,
                                                         showarrow = FALSE, xanchor = "left"))), "plsda_null")
  })
  plsda_score_rows <- reactive({
    validated_plsda()
    table_data("plsda_scores", optional = TRUE)
  })
  output$plsda_axes <- renderUI({
    if (!identical(plsda_status()$status, "validated")) return(NULL)
    score_controls("plsda", saved_axes(plsda_score_rows(), "LV"))
  })
  output$plsda_scores <- plotly::renderPlotly({
    rows <- plsda_score_rows()
    validate(need(!is.null(rows) && "sample" %in% names(rows) && length(saved_axes(rows, "LV")),
                  "Full-fit scores are unavailable."))
    axes <- selected_axes(input, "plsda", saved_axes(rows, "LV"))
    x <- axes$x
    one_component <- is.null(axes$y)
    y <- if (one_component) ".score_jitter" else axes$y
    samples <- sample_subset()
    index <- match(rows$sample, sample_key(samples))
    rows <- rows[!is.na(index), , drop = FALSE]
    index <- index[!is.na(index)]
    validate(need(nrow(rows), "No full-fit scores match the selected sample filters."))
    rows$selected_group <- as.character(samples[[input$group]][index])
    if (one_component) rows[[y]] <- ((seq_len(nrow(rows)) - 1L) %% 9L - 4L) * 0.02
    extra <- paste0(x, ": ", signif(rows[[x]], 4))
    if (!one_component) extra <- paste0(extra, "<br>", y, ": ", signif(rows[[y]], 4))
    if (!is.null(axes$z)) extra <- paste0(extra, "<br>", axes$z, ": ", signif(rows[[axes$z]], 4))
    rows$hover <- sample_hover(samples, index, extra)
    group_scatter(rows, x, y, x, if (one_component) "One component saved" else y,
                  "Full-fit scores illustrate the model; held-out and permutation results determine evidential support",
                  "plsda_scores", z = axes$z, z_label = axes$z)
  })
  selected_volcano <- reactive({
    rows <- table_data("differential", optional = TRUE)
    validate(need(!is.null(rows) && "contrast" %in% names(rows) && !is.null(input$contrast), "No saved planned contrasts for this run."))
    rows <- rows[!is.na(rows$contrast) & rows$contrast == input$contrast, , drop = FALSE]
    validate(need(nrow(rows), "No saved results for this planned contrast."))
    rows
  })
  output$volcano <- plotly::renderPlotly({
    rows <- selected_volcano()
    rows$minus_log10_p <- -log10(pmax(rows$p_value, .Machine$double.xmin))
    rows$significant <- !is.na(rows$q_value) & rows$q_value < 0.05
    rows <- rows[is.finite(rows$effect) & is.finite(rows$minus_log10_p), , drop = FALSE]
    validate(need(nrow(rows), "No finite saved results to plot for this contrast."))
    format_value <- function(value) if (is.na(value)) "NA" else format(signif(value, 4), trim = TRUE)
    rows$hover <- vapply(seq_len(nrow(rows)), function(i) {
      feature <- hover_value(rows$feature_id[i])
      if (is.null(feature)) feature <- "(unavailable)"
      paste0("Feature: ", feature, "<br>Effect: ", format_value(rows$effect[i]),
             "<br>Saved p: ", format_value(rows$p_value[i]),
             "<br>Saved BH q: ", format_value(rows$q_value[i]))
    }, character(1))
    plot <- plotly::plot_ly(source = "volcano")
    for (significant in c(FALSE, TRUE)) {
      part <- rows[rows$significant == significant, , drop = FALSE]
      if (!nrow(part)) next
      plot <- plotly::add_markers(plot, x = part$effect, y = part$minus_log10_p,
                                  key = as.character(part$feature_id), text = part$hover, hoverinfo = "text",
                                  name = if (significant) "Saved BH q < 0.05" else "Other features",
                                  marker = list(size = 7, opacity = 0.7,
                                                color = if (significant) "#C23B22" else "grey70"))
    }
    label <- unique(rows$effect_scale)
    plot <- plotly::layout(plot, title = list(text = as.character(htmltools::htmlEscape(input$contrast))),
                           xaxis = list(title = paste(label, collapse = ", ")),
                           yaxis = list(title = "-log10(saved p-value)"),
                           legend = list(itemclick = "toggle", itemdoubleclick = "toggleothers"),
                           hovermode = "closest")
    plot_controls(plotly::event_register(plot, "plotly_click"), "volcano")
  })
  output$volcano_table <- renderUI({
    rows <- selected_volcano()
    fields <- intersect(c("feature_id", "effect", "effect_scale", "p_value", "q_value"), names(rows))
    feature_table(rows[order(rows$q_value, na.last = TRUE), fields, drop = FALSE][seq_len(min(20L, nrow(rows))), , drop = FALSE])
  })
  observeEvent(plotly::event_data("plotly_click", source = "volcano"), {
    clicked <- plotly::event_data("plotly_click", source = "volcano")
    feature <- as.character(clicked$key)
    if (length(feature) == 1L && feature %in% selected_volcano()$feature_id)
      open_feature(feature)
  })
  feature_matrix <- reactive({
    req(selected_feature(), input$stage)
    matrix <- table_data(paste0("matrix_", input$stage))
    if (selected_feature() %in% names(matrix) || !selected_feature() %in% raw_feature_ids())
      list(matrix = matrix, stage = input$stage)
    else list(matrix = raw_matrix(), stage = "raw")
  })
  output$feature_plot_region <- renderUI({
    req(selected_feature(), input$stage)
    data <- feature_matrix()
    if (!selected_feature() %in% names(data$matrix))
      return(tags$p("Annotation-only ID: no measured peak height or area in this dataset."))
    note <- NULL
    if (data$stage != input$stage) {
      blank <- table_data("blank_filter", optional = TRUE)
      qc <- table_data("qc_rsd", optional = TRUE)
      removed_blank <- !is.null(blank) && all(c("feature_id", "removed_as_blank") %in% names(blank)) &&
        any(blank$feature_id == selected_feature() & !is.na(blank$removed_as_blank) & blank$removed_as_blank, na.rm = TRUE)
      removed_qc <- !is.null(qc) && all(c("feature_id", "removed_by_qc_rsd") %in% names(qc)) &&
        any(qc$feature_id == selected_feature() & !is.na(qc$removed_by_qc_rsd) & qc$removed_by_qc_rsd, na.rm = TRUE)
      reason <- if (removed_blank) "The saved blank filter removed this feature. " else
        if (removed_qc) "The saved QC-RSD filter removed this feature. " else ""
      note <- tags$p(paste0(reason, "Selected ", input$stage, " intensities are unavailable; showing saved raw ",
                            intensity_measure(), " instead."))
    }
    tagList(note, plotly::plotlyOutput("feature", height = 400))
  })
  output$feature <- plotly::renderPlotly({
    req(input$stage, input$group, selected_feature())
    feature <- selected_feature()
    data <- feature_matrix()
    matrix <- data$matrix
    validate(need(feature %in% names(matrix) && "sample" %in% names(matrix),
                  "No measured peak intensity for this annotation ID."))
    samples <- sample_subset()
    index <- match(matrix$sample, sample_key(samples))
    selected <- !is.na(index)
    validate(need(any(selected), "No samples match the selected filters."))
    index <- index[selected]
    rows <- data.frame(group = as.character(samples[[input$group]][index]),
                       intensity = as.numeric(matrix[[feature]][selected]))
    rows$hover <- sample_hover(samples, index,
                               paste0("Feature: ", htmltools::htmlEscape(feature), "<br>",
                                      if (data$stage == "raw") paste0("Raw ", intensity_measure()) else "Intensity",
                                      ": ", signif(rows$intensity, 5)))
    rows <- rows[is.finite(rows$intensity), , drop = FALSE]
    validate(need(nrow(rows), "No finite saved intensities match the selected filters."))
    groups <- sort(unique(rows$group))
    colors <- stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
    plot <- plotly::plot_ly()
    for (group in groups) {
      part <- rows[rows$group == group, , drop = FALSE]
      name <- as.character(htmltools::htmlEscape(group))
      plot <- plotly::add_boxplot(plot, x = rep(name, nrow(part)), y = part$intensity,
                                  name = name, legendgroup = group, boxpoints = FALSE,
                                  marker = list(color = colors[[group]], opacity = 0.65),
                                  line = list(color = colors[[group]]))
      plot <- plotly::add_markers(plot, x = rep(name, nrow(part)), y = part$intensity,
                                  name = name, legendgroup = group, showlegend = FALSE,
                                  text = part$hover, hoverinfo = "text",
                                  marker = list(color = colors[[group]], size = 6, opacity = 0.85))
    }
    y_label <- if (data$stage == "raw") paste("Saved raw", intensity_measure()) else
      paste("Saved", data$stage, "intensity")
    plot_controls(plotly::layout(plot, title = list(text = paste("Feature", htmltools::htmlEscape(feature))),
                                 xaxis = list(title = as.character(htmltools::htmlEscape(input$group))),
                                 yaxis = list(title = y_label),
                                 boxmode = "group", legend = list(itemclick = "toggle", itemdoubleclick = "toggleothers")),
                  "feature")
  })
  output$feature_metadata <- renderUI({
    req(selected_feature())
    rows <- input_variables()
    if (is.null(rows)) rows <- variables()
    if (!"feature_id" %in% names(rows)) return(NULL)
    record <- rows[!is.na(rows$feature_id) & rows$feature_id == selected_feature(), , drop = FALSE]
    if (!nrow(record)) return(tags$p(if (selected_feature() %in% raw_feature_ids())
      "This older run did not export input measurement metadata for this feature. Raw peak intensities are available; rerun the pipeline to include its metadata."
      else "Annotation-only ID: this feature was not present in the measured quantification dataset."))
    keep <- nzchar(names(record)) & vapply(record, function(column) !all(is.na(column) | !nzchar(as.character(column))), logical(1))
    feature_table(data.frame(field = names(record)[keep],
                             value = vapply(record[keep], function(column) as.character(column[1]), character(1))))
  })
  feature_candidates <- reactive({
    req(selected_feature())
    rows <- annotations()
    rows[!is.na(rows$feature_id) & rows$feature_id == selected_feature() & !is.na(rows$source), , drop = FALSE]
  })
  output$feature_candidate_count <- renderText({
    count <- nrow(feature_candidates())
    if (!count) "No candidate records for this feature." else
      paste("Showing", min(count, 20L), "of", count, "source-specific candidates.")
  })
  output$feature_candidates <- renderUI({
    rows <- feature_candidates()
    fields <- intersect(c("feature_id", "source", "candidate_id", "label", "confidence", "confidence_metric",
                          "molecular_formula", "inchikey", "organism", "organism_wikidata",
                          "identification_links", "npc_pathway", "npc_superclass", "npc_class"), names(rows))
    display <- rows[seq_len(min(20L, nrow(rows))), fields, drop = FALSE]
    names(display)[names(display) == "candidate_id"] <- "source_record_key"
    feature_table(display)
  })
  output$feature_structures <- renderUI({
    req(selected_feature())
    cards <- data.frame(label = character(), smiles = character(), organism = character(),
                        identification_links = character(), organism_wikidata = character())
    horizontal <- horizontal_annotations()
    if (!is.null(horizontal)) {
      record <- horizontal[!is.na(horizontal$feature_id) & horizontal$feature_id == selected_feature(), , drop = FALSE]
      if (nrow(record)) for (index in seq_len(nrow(record))) {
        for (column in grep("_smiles$", names(record), value = TRUE, ignore.case = TRUE)) {
          cards <- rbind(cards, data.frame(label = if (nrow(record) == 1L) column else paste(column, "row", index),
                                           smiles = as.character(record[[column]][index]), organism = NA_character_,
                                           identification_links = NA_character_, organism_wikidata = NA_character_))
        }
      }
    }
    candidates <- feature_candidates()
    if ("smiles" %in% names(candidates) && nrow(candidates)) {
      cards <- rbind(cards, data.frame(label = paste("source record", candidates$candidate_id),
                                       smiles = candidate_field(candidates, "smiles"),
                                       organism = candidate_field(candidates, "organism"),
                                       identification_links = candidate_field(candidates, "identification_links"),
                                       organism_wikidata = candidate_field(candidates, "organism_wikidata")))
    }
    structure_cards(cards)
  })
  consensus_rows <- reactive({
    rows <- horizontal_annotations()
    if (is.null(rows)) return(NULL)
    if (!is.null(input$annotation_consensus) && input$annotation_consensus != "All")
      rows <- rows[!is.na(rows$sources_number_IK2D) &
                     as.character(rows$sources_number_IK2D) == input$annotation_consensus, , drop = FALSE]
    rows
  })
  output$consensus_section <- renderUI({
    rows <- consensus_rows()
    if (is.null(rows)) return(tags$p("No saved horizontal consensus annotations in this run."))
    query <- tolower(trimws(if (is.null(input$annotation_query)) "" else input$annotation_query))
    if (nzchar(query)) {
      candidates <- annotations()
      fields <- intersect(c("feature_id", "label", "inchikey", "npc_pathway", "npc_superclass", "npc_class", "component_id"), names(candidates))
      matched <- Reduce(`|`, lapply(candidates[fields], function(column)
        grepl(query, tolower(as.character(column)), fixed = TRUE)), init = rep(FALSE, nrow(candidates)))
      ids <- as.character(rows$feature_id)
      rows <- rows[grepl(query, tolower(ids), fixed = TRUE) |
                     ids %in% as.character(candidates$feature_id[matched]), , drop = FALSE]
    }
    fields <- intersect(c("feature_id", "sources_number_IK2D", "sources_IK2D"), names(rows))
    tagList(h4("Saved InChIKey2D consensus per feature"),
            tags$p(length(unique(rows$feature_id)), "features across", nrow(rows), "saved horizontal rows match the consensus selection",
                   if (nrow(rows) > 100L) " (showing first 100 rows)." else ".",
                   " Candidate source and NPC filters apply only to candidate records below."),
            feature_table(rows[seq_len(min(100L, nrow(rows))), fields, drop = FALSE]))
  })
  filtered_annotations <- reactive({
    rows <- annotations()
    if (!nrow(rows)) return(rows)
    if (!is.null(input$annotation_consensus) && input$annotation_consensus != "All") {
      horizontal <- consensus_rows()
      ids <- if (is.null(horizontal)) character() else as.character(horizontal$feature_id)
      rows <- rows[!is.na(rows$feature_id) & as.character(rows$feature_id) %in% ids, , drop = FALSE]
    }
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
    labels <- paste(paste("Feature", rows$feature_id), rows$source, ifelse(is.na(rows$label), "unlabelled", rows$label),
                    paste("source record", rows$candidate_id), sep = " | ")
    selectizeInput("candidate", "Candidate evidence (search by feature ID or source record)",
                   choices = c("Choose a candidate" = "", stats::setNames(rows$candidate_id, labels)),
                   selected = "", options = list(maxOptions = 100))
  })
  output$annotation_limit <- renderText({
    rows <- filtered_annotations()
    count <- if ("source" %in% names(rows)) sum(!is.na(rows$source)) else 0L
    if (count > 100L) paste("Showing the first 100 of", count, "rows. Narrow the filters or download the complete saved table.")
  })
  output$annotation_table <- renderUI({
    rows <- filtered_annotations()
    rows <- rows[!is.na(rows$source), , drop = FALSE]
    if (!nrow(rows)) return(tags$p("No annotation candidates match the filters (annotations may be optional)."))
    fields <- intersect(c("feature_id", "source", "label", "confidence", "confidence_metric", "npc_pathway", "npc_superclass", "npc_class", "component_id", "feature_in_dataset"), names(rows))
    feature_table(rows[seq_len(min(100L, nrow(rows))), fields, drop = FALSE])
  })
  output$evidence <- renderUI({
    req(input$candidate)
    rows <- annotations()
    if (!"candidate_id" %in% names(rows)) return(NULL)
    record <- rows[!is.na(rows$candidate_id) & rows$candidate_id == input$candidate, , drop = FALSE]
    if (!nrow(record)) return(NULL)
    fields <- names(record)
    fields[fields == "candidate_id"] <- "source_record_key"
    feature_table(data.frame(field = fields,
                             value = vapply(record, function(column) as.character(column[1]), character(1))))
  })
  output$evidence_structure <- renderUI({
    req(input$candidate)
    rows <- annotations()
    if (!all(c("candidate_id", "smiles") %in% names(rows))) return(NULL)
    record <- rows[!is.na(rows$candidate_id) & rows$candidate_id == input$candidate, , drop = FALSE]
    if (!nrow(record)) return(NULL)
    structure_cards(data.frame(label = paste("Source record", record$candidate_id[1],
                                             "(feature", paste0(record$feature_id[1], ")")),
                               smiles = candidate_field(record, "smiles"),
                               organism = candidate_field(record, "organism"),
                               identification_links = candidate_field(record, "identification_links"),
                               organism_wikidata = candidate_field(record, "organism_wikidata")))
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
  static_relative <- function(stem, extension) {
    if (identical(stem, "volcano")) {
      if (is.null(input$contrast) || length(input$contrast) != 1L || is.na(input$contrast)) return(NULL)
      contrasts <- table_data("differential", optional = TRUE)
      if (is.null(contrasts) || !"contrast" %in% names(contrasts) ||
          !input$contrast %in% contrasts$contrast) return(NULL)
      name <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9._-]+", "_", input$contrast))
      if (!nzchar(name)) return(NULL)
      stem <- paste0("volcano_", name)
    }
    paste0("plots/", stem, ".", extension)
  }
  static_verified <- function(stem, extension) {
    relative <- static_relative(stem, extension)
    if (is.null(relative)) return(NULL)
    value <- run()
    if (is.null(value$manifest$outputs[[relative]])) return(NULL)
    tryCatch(verified_file(value, relative), error = function(error) NULL)
  }
  static_links <- function(stem) {
    links <- lapply(c("png", "pdf"), function(extension) {
      if (is.null(static_verified(stem, extension))) return(NULL)
      downloadButton(paste0("download_", stem, "_", extension),
                     paste("Download verified saved", toupper(extension)))
    })
    if (!any(vapply(links, Negate(is.null), logical(1)))) return(NULL)
    tagList(tags$p("Saved unfiltered analysis (not the currently filtered interactive view):"), links)
  }
  output$pca_downloads <- renderUI(static_links("pca"))
  output$pcoa_downloads <- renderUI(static_links("pcoa"))
  output$volcano_downloads <- renderUI(static_links("volcano"))
  for (stem in c("pca", "pcoa", "volcano")) {
    for (extension in c("png", "pdf")) {
      local({
        selected_stem <- stem
        selected_extension <- extension
        output[[paste0("download_", selected_stem, "_", selected_extension)]] <- downloadHandler(
          filename = function() basename(static_relative(selected_stem, selected_extension)),
          content = function(file) {
            source <- static_verified(selected_stem, selected_extension)
            validate(need(!is.null(source), "Saved plot is missing or failed checksum verification."))
            if (!file.copy(source, file, overwrite = TRUE)) stop("Unable to copy verified plot.", call. = FALSE)
          }
        )
      })
    }
  }
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
