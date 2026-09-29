export_results <- function(dataset, processed, analyses, manifest, directory) {
  dir.create(file.path(directory, "objects"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(directory, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(directory, "plots"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(dataset, file.path(directory, "objects", "dataset.rds"))
  saveRDS(processed, file.path(directory, "objects", "processed.rds"))
  saveRDS(analyses, file.path(directory, "objects", "analyses.rds"))

  write_tsv(processed$sample_metadata, file.path(directory, "tables", "sample_metadata.tsv"))
  write_tsv(processed$variable_metadata, file.path(directory, "tables", "variable_metadata.tsv"))
  if (nrow(processed$diagnostics$blank)) write_tsv(processed$diagnostics$blank, file.path(directory, "tables", "blank_filter.tsv"))
  if (nrow(processed$diagnostics$qc)) write_tsv(processed$diagnostics$qc, file.path(directory, "tables", "qc_rsd.tsv"))
  for (stage in names(processed$matrices)) {
    stage_matrix <- processed$matrices[[stage]]
    write_tsv(data.frame(sample = rownames(stage_matrix), stage_matrix, check.names = FALSE), file.path(directory, "tables", paste0("matrix_", stage, ".tsv")))
  }
  if (!is.null(analyses$pca)) {
    write_tsv(analyses$pca$scores, file.path(directory, "tables", "pca_scores.tsv"))
    write_tsv(analyses$pca$loadings, file.path(directory, "tables", "pca_loadings.tsv"))
    write_tsv(analyses$pca$variance, file.path(directory, "tables", "pca_variance.tsv"))
  }
  if (!is.null(analyses$pcoa)) {
    write_tsv(analyses$pcoa$scores, file.path(directory, "tables", "pcoa_scores.tsv"))
    write_tsv(analyses$pcoa$variance, file.path(directory, "tables", "pcoa_variance.tsv"))
  }
  if (!is.null(analyses$omnibus)) write_tsv(analyses$omnibus, file.path(directory, "tables", "omnibus.tsv"))
  if (!is.null(analyses$differential)) write_tsv(analyses$differential, file.path(directory, "tables", "differential.tsv"))
  export_analysis_plots(analyses, manifest$effective_recipe, file.path(directory, "plots"))
  write_yaml_file(manifest, file.path(directory, "manifest.yaml"))
  invisible(directory)
}

export_analysis_plots <- function(analyses, recipe, plot_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is unavailable; tabular results were saved without static plots.")
    return(invisible(FALSE))
  }
  group <- recipe$design$group
  save_plot <- function(plot, stem, width = 8, height = 6) {
    ggplot2::ggsave(file.path(plot_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 180)
    ggplot2::ggsave(file.path(plot_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  }
  if (!is.null(analyses$pca) && nrow(analyses$pca$scores)) {
    variance <- analyses$pca$variance$variance_percent
    p <- ggplot2::ggplot(analyses$pca$scores, ggplot2::aes(x = PC1, y = PC2, colour = .data[[group]])) +
      ggplot2::geom_point(size = 3) + ggplot2::theme_classic() +
      ggplot2::labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]), colour = group)
    save_plot(p, "pca")
  }
  if (!is.null(analyses$pcoa) && nrow(analyses$pcoa$scores)) {
    variance <- analyses$pcoa$variance$variance_percent
    p <- ggplot2::ggplot(analyses$pcoa$scores, ggplot2::aes(x = PCoA1, y = PCoA2, colour = .data[[group]])) +
      ggplot2::geom_point(size = 3) + ggplot2::theme_classic() +
      ggplot2::labs(x = sprintf("PCoA1 (%.1f%%)", variance[1]), y = sprintf("PCoA2 (%.1f%%)", variance[2]), colour = group)
    save_plot(p, "pcoa")
  }
  if (!is.null(analyses$differential) && nrow(analyses$differential)) {
    for (contrast in unique(analyses$differential$contrast)) {
      values <- analyses$differential[analyses$differential$contrast == contrast, , drop = FALSE]
      values$significant <- !is.na(values$q_value) & values$q_value < 0.05
      p <- ggplot2::ggplot(values, ggplot2::aes(x = effect, y = -log10(pmax(p_value, .Machine$double.xmin)), colour = significant)) +
        ggplot2::geom_point(alpha = 0.65, size = 1.5) + ggplot2::scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#C23B22")) +
        ggplot2::theme_classic() + ggplot2::labs(title = contrast, x = unique(values$effect_scale), y = "-log10(p-value)", colour = "BH q < 0.05")
      save_plot(p, paste0("volcano_", safe_name(contrast)))
    }
  }
  invisible(TRUE)
}
