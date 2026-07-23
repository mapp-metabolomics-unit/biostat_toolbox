#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(digest)
  library(dplyr)
})

args_full <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_full[grep("^--file=", args_full)])
if (!length(script_path)) {
  script_path <- file.path(getwd(), "src", "annotate_hash_table.R")
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "helpers.r"), local = TRUE)

option_list <- list(
  make_option(c("-t", "--table"), default = NULL, help = "Path to params_log.tsv or reprocess_manifest.tsv"),
  make_option(c("--params-dir-column"), default = "output_dir", help = "Column containing result directories with params.yaml"),
  make_option(c("--params-file-column"), default = "source_params", help = "Fallback column containing params.yaml paths")
)

opt <- parse_args(OptionParser(option_list = option_list))

normalize_option_name <- function(opt, underscore_name, hyphen_name) {
  if (is.null(opt[[underscore_name]]) && !is.null(opt[[hyphen_name]])) {
    opt[[underscore_name]] <- opt[[hyphen_name]]
  }
  opt
}

opt <- normalize_option_name(opt, "params_dir_column", "params-dir-column")
opt <- normalize_option_name(opt, "params_file_column", "params-file-column")

if (is.null(opt$table) || !nzchar(opt$table)) {
  stop("--table is required.")
}
if (!file.exists(opt$table)) {
  stop(sprintf("Table not found: %s", opt$table))
}

hash_table <- read.delim(opt$table, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
descriptions <- character(nrow(hash_table))

for (row_index in seq_len(nrow(hash_table))) {
  params_path <- ""
  if (opt$params_dir_column %in% names(hash_table) && nzchar(hash_table[[opt$params_dir_column]][row_index])) {
    params_path <- file.path(hash_table[[opt$params_dir_column]][row_index], "params.yaml")
  }
  if (!file.exists(params_path) && opt$params_file_column %in% names(hash_table)) {
    params_path <- hash_table[[opt$params_file_column]][row_index]
  }
  descriptions[row_index] <- if (file.exists(params_path)) {
    describe_params_for_hash_table(yaml.load_file(params_path))
  } else {
    ""
  }
}

hash_table$description <- descriptions
first_columns <- intersect(c("original_hash", "new_hash", "hash", "timestamp", "description"), names(hash_table))
hash_table <- hash_table[, c(first_columns, setdiff(names(hash_table), first_columns)), drop = FALSE]
write.table(hash_table, opt$table, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Annotated %d row(s): %s", nrow(hash_table), opt$table))
