# Validate either the data-free public repository or an authorised full run.
#
# Default: data-free public validation. Use --full, or set
# execution.full_validation: true in config/analysis_parameters.yml, only when
# the controlled analysis outputs are available locally.

find_entrypoint_root <- function(start = getwd()) {
  candidates <- c(start)
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument)) {
    candidates <- c(
      dirname(normalizePath(sub("^--file=", "", file_argument[1L]))),
      candidates
    )
  }
  for (candidate in unique(candidates)) {
    current <- normalizePath(candidate, mustWork = TRUE)
    repeat {
      if (file.exists(file.path(
        current, "ubio-pred-airway-pathway-analysis.Rproj"
      ))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the repository root.")
}

setwd(find_entrypoint_root())
source(file.path("R", "00_common.R"))

arguments <- commandArgs(trailingOnly = TRUE)
unknown_arguments <- setdiff(arguments, c("--public", "--full"))
if (length(unknown_arguments)) {
  stop("Unknown argument(s): ", paste(unknown_arguments, collapse = ", "))
}
if (all(c("--public", "--full") %in% arguments)) {
  stop("Choose either --public or --full, not both.")
}

full_configured <- isTRUE(config_value(
  c("execution", "full_validation"),
  FALSE
))
full_requested <- "--full" %in% arguments ||
  (full_configured && !("--public" %in% arguments))

run_r_script <- function(relative_path) {
  script_path <- file.path(PROJECT_ROOT, relative_path)
  if (!file.exists(script_path)) stop("Required script is missing: ", relative_path)
  message("\n>>> Running ", relative_path)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = shQuote(normalizePath(script_path, mustWork = TRUE)),
    stdout = "",
    stderr = ""
  )
  if (!isTRUE(status == 0L)) {
    stop(relative_path, " failed with exit status ", status, ".")
  }
  invisible(TRUE)
}

verify_manifest <- function(manifest_relative_path, path_column) {
  manifest_path <- file.path(PROJECT_ROOT, manifest_relative_path)
  if (!file.exists(manifest_path)) {
    stop("Committed manifest is missing: ", manifest_relative_path)
  }
  manifest <- utils::read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(path_column, "sha256", "size_bytes")
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns)) {
    stop(
      "Manifest ", manifest_relative_path, " is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (!nrow(manifest)) stop("Manifest has no entries: ", manifest_relative_path)

  listed_paths <- manifest[[path_column]]
  unsafe <- is.na(listed_paths) |
    !nzchar(listed_paths) |
    grepl("^(/|[A-Za-z]:)", listed_paths) |
    grepl("(^|/)\\.\\.(/|$)", listed_paths)
  if (any(unsafe)) {
    stop("Manifest contains an unsafe or blank relative path: ", manifest_relative_path)
  }

  failures <- character()
  for (i in seq_len(nrow(manifest))) {
    relative_path <- listed_paths[i]
    target_path <- file.path(PROJECT_ROOT, relative_path)
    if (!file.exists(target_path)) {
      failures <- c(failures, paste0(relative_path, " (missing)"))
      next
    }
    actual_size <- unname(file.info(target_path)$size)
    expected_size <- suppressWarnings(as.numeric(manifest$size_bytes[i]))
    if (!is.finite(expected_size) || actual_size != expected_size) {
      failures <- c(
        failures,
        paste0(relative_path, " (size ", actual_size, "; expected ", expected_size, ")")
      )
      next
    }
    actual_hash <- digest::digest(
      file = target_path,
      algo = "sha256",
      serialize = FALSE
    )
    expected_hash <- tolower(trimws(manifest$sha256[i]))
    if (!identical(tolower(actual_hash), expected_hash)) {
      failures <- c(failures, paste0(relative_path, " (SHA-256 mismatch)"))
    }
  }
  if (length(failures)) {
    stop(
      "Manifest verification failed for ", manifest_relative_path, ":\n- ",
      paste(failures, collapse = "\n- ")
    )
  }
  message("Verified ", nrow(manifest), " entries in ", manifest_relative_path)
  nrow(manifest)
}

validate_reference_inventory <- function(manifest_relative_path, path_column,
                                         expected_rows) {
  manifest_path <- file.path(PROJECT_ROOT, manifest_relative_path)
  if (!file.exists(manifest_path)) {
    stop("Reference inventory is missing: ", manifest_relative_path)
  }
  manifest <- utils::read.csv(
    manifest_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(path_column, "sha256", "size_bytes")
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns)) {
    stop(
      "Reference inventory ", manifest_relative_path,
      " is missing columns: ", paste(missing_columns, collapse = ", ")
    )
  }
  if (nrow(manifest) != expected_rows) {
    stop(
      "Reference inventory ", manifest_relative_path, " has ", nrow(manifest),
      " rows; expected ", expected_rows, "."
    )
  }
  listed_paths <- trimws(manifest[[path_column]])
  if (anyNA(listed_paths) || any(!nzchar(listed_paths)) ||
      anyDuplicated(listed_paths)) {
    stop("Reference inventory has blank or duplicate paths: ", manifest_relative_path)
  }
  unsafe <- grepl("^(/|[A-Za-z]:)", listed_paths) |
    grepl("(^|/)\\.\\.(/|$)", listed_paths)
  valid_hash <- grepl("^[0-9a-fA-F]{64}$", trimws(manifest$sha256))
  valid_size <- is.finite(suppressWarnings(as.numeric(manifest$size_bytes))) &
    suppressWarnings(as.numeric(manifest$size_bytes)) > 0
  if (any(unsafe) || any(!valid_hash) || any(!valid_size)) {
    stop("Reference inventory contains an invalid path, hash or size: ", manifest_relative_path)
  }
  message("Validated ", nrow(manifest), " metadata rows in ", manifest_relative_path)
  nrow(manifest)
}

run_public_validation <- function() {
  # Confirm that the restored lock contains every package required by the
  # scientific runners, even though the synthetic demonstration exercises
  # only a subset of those packages.
  assert_packages()
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("R.utils is required to read and verify compressed CSV outputs.")
  }

  lockfile <- renv::lockfile_read(file.path(PROJECT_ROOT, "renv.lock"))
  missing_from_lock <- setdiff(required_packages, names(lockfile$Packages))
  if (length(missing_from_lock)) {
    stop(
      "Scientific packages are installed but absent from renv.lock: ",
      paste(missing_from_lock, collapse = ", ")
    )
  }
  message(
    "Required scientific packages recorded in renv.lock: PASS (",
    length(required_packages), ")"
  )

  gzip_fixture <- data.table::data.table(
    integer_value = 1:3,
    text_value = c("alpha", "beta", "gamma"),
    numeric_value = c(-1.25, 0, 2.5)
  )
  gzip_path <- tempfile(fileext = ".csv.gz")
  data.table::fwrite(gzip_fixture, gzip_path, compress = "gzip")
  gzip_restored <- data.table::fread(gzip_path, showProgress = FALSE)
  unlink(gzip_path)
  if (!identical(gzip_fixture, gzip_restored)) {
    stop("Compressed CSV round-trip did not preserve the synthetic fixture.")
  }
  message("Compressed CSV round-trip: PASS")

  run_r_script(file.path("validation", "check_R_syntax.R"))
  run_r_script("run_demo.R")
  run_r_script(file.path("validation", "check_public_release.R"))

  tracked_manifest_specs <- list(
    list("validation/reference_manifests/SCHEMA_FILE_MANIFEST.csv", "relative_path"),
    list("validation/reference_manifests/SYNTHETIC_FIXTURE_MANIFEST.csv", "relative_path"),
    list("validation/reference_manifests/REPOSITORY_SAFE_CONTENT_MANIFEST.csv", "relative_path")
  )
  checked_files <- vapply(
    tracked_manifest_specs,
    function(spec) verify_manifest(spec[[1L]], spec[[2L]]),
    integer(1L)
  )
  reference_inventory_specs <- list(
    list("validation/reference_manifests/LOCKED_CONTROLLED_INPUT_MANIFEST.csv", "relative_path", 13L),
    list("validation/reference_manifests/PUBLIC_CONTENT_MANIFEST.csv", "relative_path", 123L),
    list("validation/reference_manifests/PUBLIC_TABLE_MANIFEST.csv", "public_relative_path", 41L),
    list("validation/reference_manifests/AUTHORITATIVE_FIGURE_PUBLIC_MANIFEST.csv", "public_relative_path", 68L)
  )
  inventory_rows <- vapply(
    reference_inventory_specs,
    function(spec) validate_reference_inventory(spec[[1L]], spec[[2L]], spec[[3L]]),
    integer(1L)
  )

  public_status_dir <- file.path(PROJECT_ROOT, "work", "validation")
  dir.create(public_status_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(
    c(
      "DATA-FREE REPOSITORY VALIDATION: PASS",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
      "R syntax check: PASS",
      paste(
        "Required scientific packages recorded in renv.lock: PASS (",
        length(required_packages), ")",
        sep = ""
      ),
      "Compressed CSV round-trip: PASS",
      "Wholly synthetic software demonstration: PASS",
      "Repository safety gate: PASS",
      paste("Tracked files verified by SHA-256:", sum(checked_files)),
      paste("Withheld-output metadata rows validated:", sum(inventory_rows)),
      "No controlled U-BIOPRED data were required or accessed.",
      "This is software/repository validation, not scientific result reproduction.",
      "This technical validation does not authorise distribution of controlled study materials."
    ),
    file.path(public_status_dir, "PUBLIC_REPOSITORY_VALIDATION_STATUS.txt")
  )
  message("\nDATA-FREE REPOSITORY VALIDATION: PASS")
  invisible(TRUE)
}

run_full_validation <- function() {
  assert_packages()

  message("\n>>> Rechecking the currently configured controlled inputs")
  sys.source(
    file.path("R", "01_validate_locked_inputs.R"),
    envir = new.env(parent = globalenv())
  )

  full_output_markers <- c(
    file.path(paths$analysis, "nested_cv_fold_level_predictions.csv.gz"),
    file.path(paths$tables, "Table_4_repeated_nested_CV_model_performance.csv"),
    file.path(paths$final_figures, "FINAL_DISSERTATION_FIGURE_EXPORT_MANIFEST.csv")
  )
  missing_markers <- full_output_markers[!file.exists(full_output_markers)]
  if (length(missing_markers)) {
    stop(
      "Full validation was requested, but controlled/full-analysis outputs are missing:\n- ",
      paste(missing_markers, collapse = "\n- "),
      "\nRun run_analysis.R and run_outputs.R in an authorised environment first, ",
      "or omit --full to run the data-free public checks."
    )
  }

  for (script_name in c(
    "10_validate_nested_cv.R",
    "26_validate_final_dissertation_figure_set_correct_numbering.R"
  )) {
    message("\n>>> Running R/", script_name)
    sys.source(file.path("R", script_name), envir = new.env(parent = globalenv()))
  }

  required_status_files <- c(
    "LOCKED_INPUT_HASH_STATUS.txt",
    "INPUT_AUDIT_STATUS.txt",
    "EXPRESSION_PREPARATION_STATUS.txt",
    "SSGSEA_STATUS.txt",
    "CLINICAL_MERGE_STATUS.txt",
    "DESCRIPTIVE_SPEARMAN_STATUS.txt",
    "CURATED_MODELS_STATUS.txt",
    "MATCHED_CONCORDANCE_STATUS.txt",
    "NESTED_CV_COMPUTATION_STATUS.txt",
    "NESTED_CV_QA_STATUS.txt",
    "PREDICTIVE_INFERENCE_STATUS.txt",
    "FIGURE_FORMAT_AND_SOURCE_QA_STATUS.txt",
    "SIGNATURE_PROVENANCE_OUTPUT_STATUS.txt"
  )
  status_paths <- file.path(paths$validation, required_status_files)
  missing_status <- status_paths[!file.exists(status_paths)]
  if (length(missing_status)) {
    stop(
      "Required validation records are missing: ",
      paste(missing_status, collapse = "; ")
    )
  }
  invalid_status <- vapply(
    status_paths,
    function(status_path) {
      status_text <- paste(readLines(status_path, warn = FALSE), collapse = "\n")
      !grepl("PASS", status_text, fixed = TRUE) ||
        grepl("FAIL", status_text, fixed = TRUE)
    },
    logical(1L)
  )
  if (any(invalid_status)) {
    stop(
      "Required validation records do not contain a clean PASS: ",
      paste(status_paths[invalid_status], collapse = "; ")
    )
  }

  final_figure_status <- file.path(
    paths$final_figures, "FINAL_DISSERTATION_FIGURE_VALIDATION.txt"
  )
  table_status <- file.path(paths$final_tables, "TABLE_1_R_GENERATION_STATUS.txt")
  if (!file.exists(final_figure_status) || !file.exists(table_status)) {
    stop("Final figure or Table 1 validation status is missing.")
  }
  final_status_paths <- c(final_figure_status, table_status)
  invalid_final_status <- vapply(
    final_status_paths,
    function(status_path) {
      status_text <- paste(readLines(status_path, warn = FALSE), collapse = "\n")
      !grepl("PASS", status_text, fixed = TRUE) ||
        grepl("FAIL", status_text, fixed = TRUE)
    },
    logical(1L)
  )
  if (any(invalid_final_status)) {
    stop(
      "Final figure or Table 1 status is not a clean PASS: ",
      paste(final_status_paths[invalid_final_status], collapse = "; ")
    )
  }
  table_stem <- file.path(
    paths$final_tables, "01_Table_1_sample_omics_and_preprocessing"
  )
  table_artifacts <- paste0(
    table_stem, c(".pdf", ".svg", "_PNG300.png", "_TIFF600.tiff")
  )
  missing_or_empty_table <- table_artifacts[
    !file.exists(table_artifacts) |
      is.na(file.info(table_artifacts)$size) |
      file.info(table_artifacts)$size <= 0
  ]
  if (length(missing_or_empty_table)) {
    stop(
      "Table 1 validation artifacts are missing or empty: ",
      paste(missing_or_empty_table, collapse = "; ")
    )
  }

  writeLines(
    c(
      "FULL PORTABLE PIPELINE VALIDATION: PASS",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
      "Nested cross-validation leakage and reconstruction checks: PASS",
      "Figures 1-9 and S1-S12 four-format computational checks: PASS",
      "Figure S12 provenance/count consistency check: PASS",
      "Table 1 R-generation check: PASS",
      "Internal validation only; this is not external validation."
    ),
    file.path(paths$validation, "FULL_PIPELINE_VALIDATION_STATUS.txt")
  )
  message("\nFULL PORTABLE PIPELINE VALIDATION: PASS")
  invisible(TRUE)
}

if (full_requested) {
  run_full_validation()
} else {
  run_public_validation()
}
