if (!requireNamespace("shiny", quietly = TRUE)) stop("The V2 explorer requires shiny. Install it in the V2 environment.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The V2 explorer requires ggplot2.")

library(shiny)

initial_run <- Sys.getenv("MAPP_STATS_RUN", unset = "")

ui <- fluidPage(
  titlePanel("MAPP statistics explorer V2"),
  sidebarLayout(
    sidebarPanel(
      textInput("run_dir", "Completed run directory", initial_run),
      actionButton("load", "Load run", class = "btn-primary"),
      hr(),
      uiOutput("controls")
    ),
    mainPanel(
      uiOutput("status"),
      tabsetPanel(
        tabPanel("Data", tableOutput("summary")),
        tabPanel("PCA", plotOutput("pca", height = 650)),
        tabPanel("PCoA", plotOutput("pcoa", height = 650)),
        tabPanel("Volcano", plotOutput("volcano", height = 650), tableOutput("volcano_table")),
        tabPanel("Feature", plotOutput("feature", height = 650)),
        tabPanel("Provenance", verbatimTextOutput("manifest"))
      )
    )
  )
)

server <- function(input, output, session) {
  run <- eventReactive(input$load, {
    path <- normalizePath(input$run_dir, mustWork = FALSE)
    required <- c("COMPLETE", "manifest.yaml", "objects/processed.rds", "objects/analyses.rds")
    missing <- required[!file.exists(file.path(path, required))]
    validate(need(!length(missing), paste("Not a completed V2 run; missing:", paste(missing, collapse = ", "))))
    list(
      path = path,
      manifest = yaml::read_yaml(file.path(path, "manifest.yaml")),
      processed = readRDS(file.path(path, "objects", "processed.rds")),
      analyses = readRDS(file.path(path, "objects", "analyses.rds"))
    )
  }, ignoreInit = !nzchar(initial_run))

  output$status <- renderUI({
    req(run())
    tags$p(class = "text-success", paste("Loaded", run()$manifest$dataset_id, "run", basename(run()$path)))
  })

  output$controls <- renderUI({
    value <- run(); req(value)
    stages <- names(value$processed$matrices)
    features <- colnames(value$processed$matrices$transformed)
    contrasts <- if (!is.null(value$analyses$differential)) unique(value$analyses$differential$contrast) else character()
    tagList(
      selectInput("stage", "Feature intensity stage", stages, selected = "transformed"),
      selectInput("feature_id", "Feature", features),
      if (length(contrasts)) selectInput("contrast", "Volcano contrast", contrasts)
    )
  })

  output$summary <- renderTable({
    value <- run(); req(value)
    validation <- value$manifest$validation
    data.frame(
      metric = c("Samples", "Features", "Raw zero fraction", names(validation$role_counts)),
      value = c(validation$dimensions$samples, validation$dimensions$features, signif(validation$zero_fraction, 4), unlist(validation$role_counts)),
      row.names = NULL
    )
  })

  ordination_plot <- function(scores, x, y, variance, group) {
    ggplot2::ggplot(scores, ggplot2::aes(x = .data[[x]], y = .data[[y]], colour = .data[[group]], label = sample_id)) +
      ggplot2::geom_point(size = 3) + ggplot2::theme_classic() +
      ggplot2::labs(x = sprintf("%s (%.1f%%)", x, variance[1]), y = sprintf("%s (%.1f%%)", y, variance[2]), colour = group)
  }

  output$pca <- renderPlot({
    value <- run(); req(value$analyses$pca)
    ordination_plot(value$analyses$pca$scores, "PC1", "PC2", value$analyses$pca$variance$variance_percent, value$manifest$effective_recipe$design$group)
  })
  output$pcoa <- renderPlot({
    value <- run(); req(value$analyses$pcoa)
    ordination_plot(value$analyses$pcoa$scores, "PCoA1", "PCoA2", value$analyses$pcoa$variance$variance_percent, value$manifest$effective_recipe$design$group)
  })
  selected_volcano <- reactive({
    value <- run(); req(value$analyses$differential, input$contrast)
    value$analyses$differential[value$analyses$differential$contrast == input$contrast, , drop = FALSE]
  })
  output$volcano <- renderPlot({
    values <- selected_volcano()
    values$significant <- values$q_value < 0.05
    ggplot2::ggplot(values, ggplot2::aes(effect, -log10(pmax(p_value, .Machine$double.xmin)), colour = significant)) +
      ggplot2::geom_point(alpha = 0.65) + ggplot2::scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#C23B22")) +
      ggplot2::theme_classic() + ggplot2::labs(x = unique(values$effect_scale), y = "-log10(p-value)", colour = "BH q < 0.05")
  })
  output$volcano_table <- renderTable({
    values <- selected_volcano()
    head(values[order(values$q_value), ], 20)
  })
  output$feature <- renderPlot({
    value <- run(); req(input$stage, input$feature_id)
    matrix <- value$processed$matrices[[input$stage]]
    validate(need(input$feature_id %in% colnames(matrix), "Feature is absent from this stage."))
    metadata <- value$processed$sample_metadata
    sample_ids <- intersect(rownames(matrix), rownames(metadata))
    group <- value$manifest$effective_recipe$design$group
    values <- data.frame(sample = sample_ids, intensity = matrix[sample_ids, input$feature_id], group = metadata[sample_ids, group])
    ggplot2::ggplot(values, ggplot2::aes(group, intensity, colour = group)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.12) +
      ggplot2::theme_classic() + ggplot2::labs(title = paste("Feature", input$feature_id, "at", input$stage, "stage"), x = group)
  })
  output$manifest <- renderPrint({ str(run()$manifest, max.level = 5) })
}

shinyApp(ui, server)

