#!/usr/bin/env Rscript

# Data-free repository safety gate. This checks the repository contents without
# reading any controlled U-BIOPRED source directory. A pass is not approval for
# release of controlled study materials.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[startsWith(args, file_arg)][1])
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

if (dir.exists(file.path(repo_root, ".git")) && nzchar(Sys.which("git"))) {
  entry_relative <- system2(
    Sys.which("git"), c("-C", repo_root, "ls-files"), stdout = TRUE
  )
  entry_relative <- entry_relative[nzchar(entry_relative)]
  all_entries <- file.path(repo_root, entry_relative)
} else {
  all_entries <- list.files(
    repo_root,
    recursive = TRUE,
    all.files = TRUE,
    full.names = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  entry_relative <- substring(all_entries, nchar(repo_root) + 2L)
  keep <- !startsWith(entry_relative, "renv/library/") &
    !startsWith(entry_relative, "renv/staging/") &
    !startsWith(entry_relative, "outputs/public/") &
    !startsWith(entry_relative, "work/")
  all_entries <- all_entries[keep]
  entry_relative <- entry_relative[keep]
}

failures <- character()

link_targets <- Sys.readlink(all_entries)
link_hits <- entry_relative[nzchar(link_targets)]
if (length(link_hits)) {
  failures <- c(
    failures,
    paste0("Symbolic links are not permitted: ", paste(link_hits, collapse = ", "))
  )
}

entry_info <- file.info(all_entries)
is_file <- !is.na(entry_info$isdir) & !entry_info$isdir
all_files <- all_entries[is_file]
relative <- entry_relative[is_file]

forbidden_components <- c(
  "data/raw/", "data/controlled/", "data/private/",
  "data/participant_level/", "outputs/private/", "outputs/controlled/",
  "results/", "tmp/",
  "99_ARCHIVE", "U-Biopred Adult DataExport"
)
for (pattern in forbidden_components) {
  hits <- relative[grepl(pattern, relative, fixed = TRUE)]
  if (length(hits)) {
    failures <- c(
      failures,
      paste0("Forbidden path '", pattern, "': ", paste(hits, collapse = ", "))
    )
  }
}

hidden_work_hits <- relative[
  grepl("(^|/)\\.[^/]*_work/", relative, perl = TRUE)
]
if (length(hidden_work_hits)) {
  failures <- c(
    failures,
    paste0(
      "Hidden working directory is not permitted: ",
      paste(hidden_work_hits, collapse = ", ")
    )
  )
}

forbidden_extensions <- c(
  "rds", "rda", "rdata", "doc", "docx", "ppt", "pptx", "enl", "enlx",
  "xls", "xlsx", "tif", "tiff", "zip", "tar", "gz", "parquet",
  "feather", "fst", "h5", "hdf5", "sqlite", "db"
)
extensions <- tolower(tools::file_ext(relative))
bad_extension <- relative[extensions %in% forbidden_extensions]
if (length(bad_extension)) {
  failures <- c(
    failures,
    paste0("Forbidden file type: ", paste(bad_extension, collapse = ", "))
  )
}

sizes <- file.info(all_files)$size
too_large <- relative[is.finite(sizes) & sizes > 50 * 1024^2]
if (length(too_large)) {
  failures <- c(
    failures,
    paste0("File exceeds 50 MiB: ", paste(too_large, collapse = ", "))
  )
}

text_extensions <- c(
  "r", "md", "txt", "csv", "tsv", "yml", "yaml", "json", "cff",
  "rproj", "gitignore", "gitattributes"
)
is_text <- extensions %in% text_extensions |
  basename(relative) %in% c("LICENSE", "CITATION.cff", ".gitignore", ".gitattributes")
local_path_pattern <- paste0(
  "(", "/", "Users/[^/]+/", "|", "/", "home/[^/]+/",
  "|[A-Za-z]:[/\\\\]Users[/\\\\][^/\\\\]+[/\\\\])"
)
for (i in which(is_text)) {
  lines <- tryCatch(
    readLines(all_files[i], warn = FALSE),
    error = function(e) character()
  )
  if (any(grepl(local_path_pattern, lines, perl = TRUE))) {
    failures <- c(failures, paste0("Personal absolute path in ", relative[i]))
  }
}

manifest_path <- file.path(
  repo_root, "validation", "reference_manifests", "PUBLIC_CONTENT_MANIFEST.csv"
)
if (!file.exists(manifest_path)) {
  failures <- c(failures, "Public-content allowlist manifest is missing")
} else {
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  inventory_roots <- c("data/schemas/", "tests/synthetic/")
  in_inventory <- Reduce(
    `|`, lapply(inventory_roots, function(root) startsWith(relative, root))
  )
  listed_inventory <- Reduce(
    `|`, lapply(inventory_roots, function(root) {
      startsWith(manifest$relative_path, root)
    })
  )
  actual <- sort(relative[in_inventory])
  allowed <- sort(manifest$relative_path[listed_inventory])
  extra <- setdiff(actual, allowed)
  missing <- setdiff(allowed, actual)
  if (length(extra)) {
    failures <- c(
      failures,
      paste0("Unlisted staged file: ", paste(extra, collapse = ", "))
    )
  }
  if (length(missing)) {
    failures <- c(
      failures,
      paste0("Allowlisted staged file is missing: ", paste(missing, collapse = ", "))
    )
  }
}

public_csv <- which(
  startsWith(relative, "outputs/public/") &
    extensions %in% c("csv", "tsv")
)
sensitive_fields <- c(
  "subject", "subject_id", "participant", "participant_id",
  "sample", "sample_id", "assay", "assay_id", "analysed_gene_symbols",
  "exact_analysed_genes", "detected_genes"
)
for (i in public_csv) {
  header <- tryCatch(
    readLines(all_files[i], n = 1L, warn = FALSE),
    error = function(e) ""
  )
  separator <- if (extensions[i] == "tsv") "\t" else ","
  fields <- if (length(header)) {
    trimws(tolower(strsplit(header, separator, fixed = TRUE)[[1L]]))
  } else {
    character()
  }
  if (any(fields %in% sensitive_fields)) {
    failures <- c(
      failures,
      paste0("Potential participant-level identifier in public table: ", relative[i])
    )
  }
}

if (length(failures)) {
  cat("REPOSITORY SAFETY CHECK: FAIL\n")
  cat(paste0("- ", unique(failures), collapse = "\n"), "\n")
  quit(status = 1L)
}

cat("REPOSITORY SAFETY CHECK: PASS\n")
cat("Files checked:", length(relative), "\n")
cat("No controlled-data paths, participant-level public tables, personal absolute paths,\n")
cat("unlisted staged artifacts, forbidden binary files, symbolic links, or files\n")
cat("above 50 MiB were detected.\n")
cat("CONTROLLED MATERIALS: distribution requires institutional and U-BIOPRED approval.\n")
