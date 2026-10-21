mapp_biostats_schema_version <- "1.0.0"

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


biostats_code_checksums <- function(repo_root) {
  module_dir <- file.path(repo_root, "src", "biostats")
  modules <- sort(list.files(module_dir, pattern = "[.]R$", full.names = FALSE))
  files <- c(file.path(module_dir, modules), file.path(repo_root, "src", "mapp_stats.R"))
  if (any(!file.exists(files))) stop("Cannot fingerprint missing Biostats source code.", call. = FALSE)
  checksums <- vapply(files, sha256_file, character(1))
  names(checksums) <- c(file.path("src", "biostats", modules), file.path("src", "mapp_stats.R"))
  as.list(checksums)
}

# Runs must originate from exactly the committed toolbox, including files that
# were sourced before the run began. Recheck against the saved hashes at sealing.
committed_biostats_code <- function(repo_root, code_sha256 = NULL, commit = NULL) {
  git <- Sys.which("git")
  if (!nzchar(git)) stop("A Git checkout and committed toolbox code are required to run Biostats.", call. = FALSE)
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
  git_args <- function(...) c("-C", shQuote(repo_root), ...)
  command <- function(...) {
    output <- suppressWarnings(system2(git, git_args(...), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(output, "status")) || !length(output))
      stop("Cannot verify committed toolbox code (Git checkout or HEAD unavailable).", call. = FALSE)
    output
  }
  root <- command("rev-parse", "--show-toplevel")
  if (length(root) != 1L || !identical(normalizePath(root, mustWork = TRUE), repo_root))
    stop("Biostats toolbox must be the root of its Git checkout.", call. = FALSE)
  head <- command("rev-parse", "--verify", "HEAD^{commit}")
  if (length(head) != 1L || !grepl("^[0-9a-f]{40}([0-9a-f]{24})?$", head))
    stop("Cannot verify committed toolbox code: no valid Git HEAD.", call. = FALSE)
  if (!is.null(commit) && !identical(commit, head))
    stop("Toolbox Git HEAD changed during the Biostats run.", call. = FALSE)
  status_path <- tempfile("biostats-git-status-")
  on.exit(unlink(status_path), add = TRUE)
  status <- suppressWarnings(system2(git, git_args("status", "--porcelain=v1", "-z",
    "--untracked-files=all"), stdout = status_path, stderr = FALSE))
  if (!identical(as.integer(status), 0L))
    stop("Cannot verify clean toolbox Git status.", call. = FALSE)
  if (file.info(status_path)$size > 0L)
    stop("Biostats run requires a clean committed toolbox (tracked, staged and untracked files must be committed or removed).",
         call. = FALSE)
  if (is.null(code_sha256)) code_sha256 <- biostats_code_checksums(repo_root)
  if (!is.list(code_sha256) || !length(code_sha256) || is.null(names(code_sha256)) ||
      anyDuplicated(names(code_sha256)) || !identical(names(code_sha256),
        names(biostats_code_checksums(repo_root))))
    stop("Cannot verify saved Biostats source checksums.", call. = FALSE)
  snapshot <- tempfile("biostats-head-")
  on.exit(unlink(snapshot), add = TRUE)
  for (path in names(code_sha256)) {
    result <- suppressWarnings(system2(git, git_args("show", shQuote(paste0("HEAD:", path))),
                                       stdout = snapshot, stderr = FALSE))
    if (!identical(as.integer(result), 0L) ||
        !identical(sha256_file(snapshot), code_sha256[[path]]))
      stop("Biostats source differs from committed HEAD or changed during the run: ", path,
           call. = FALSE)
  }
  head
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
