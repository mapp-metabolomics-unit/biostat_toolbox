mapp_v2_schema_version <- "1.0.0"

`%||%` <- function(x, y) {
  if (is.null(x) || !length(x)) y else x
}

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

read_yaml_file <- function(path) {
  assert_packages("yaml")
  if (!file.exists(path)) stop("YAML file not found: ", path, call. = FALSE)
  yaml::read_yaml(path)
}

write_yaml_file <- function(value, path) {
  assert_packages("yaml")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  yaml::write_yaml(value, path)
}

canonicalize <- function(x) {
  if (is.list(x)) {
    if (!is.null(names(x))) x <- x[order(names(x))]
    return(lapply(x, canonicalize))
  }
  if (is.factor(x)) return(as.character(x))
  x
}

sha256_object <- function(x) {
  assert_packages("digest")
  digest::digest(canonicalize(x), algo = "sha256", serialize = TRUE)
}

sha256_file <- function(path) {
  assert_packages("digest")
  if (!file.exists(path)) stop("Cannot checksum missing file: ", path, call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

normalize_existing_path <- function(path, base = getwd()) {
  if (!grepl("^/", path)) path <- file.path(base, path)
  normalizePath(path, mustWork = TRUE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = TRUE, qmethod = "double", row.names = FALSE, na = "")
}

write_json <- function(x, path, pretty = TRUE) {
  assert_packages("jsonlite")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = pretty, null = "null", na = "null")
}

package_versions <- function(packages) {
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  }, character(1))
  as.list(versions)
}


v2_code_checksums <- function(repo_root) {
  module_dir <- file.path(repo_root, "src", "v2")
  modules <- sort(list.files(module_dir, pattern = "[.]R$", full.names = FALSE))
  files <- c(file.path(module_dir, modules), file.path(repo_root, "src", "mapp_stats.R"))
  if (any(!file.exists(files))) stop("Cannot fingerprint missing V2 source code.", call. = FALSE)
  checksums <- vapply(files, sha256_file, character(1))
  names(checksums) <- c(file.path("src", "v2", modules), file.path("src", "mapp_stats.R"))
  as.list(checksums)
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

atomic_publish <- function(staging_dir, final_dir) {
  if (file.exists(final_dir)) stop("Refusing to overwrite existing run: ", final_dir, call. = FALSE)
  if (!file.rename(staging_dir, final_dir)) stop("Could not publish run directory: ", final_dir, call. = FALSE)
  invisible(final_dir)
}
