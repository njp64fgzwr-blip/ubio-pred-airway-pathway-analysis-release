# Validate the exact processed TSV exports used by the dissertation.
#
# The tracked manifest contains only logical relative paths, byte sizes and
# SHA-256 values. Controlled file contents and local absolute paths remain
# outside the repository.

if (!exists("PROJECT_ROOT") || !exists("datasets")) {
  source(file.path("R", "00_common.R"))
}
assert_packages("digest")
require_controlled_inputs()

manifest_path <- file.path(
  PROJECT_ROOT, "validation", "reference_manifests",
  "LOCKED_CONTROLLED_INPUT_MANIFEST.csv"
)
if (!file.exists(manifest_path)) {
  stop("Locked controlled-input manifest is missing: ", manifest_path)
}

manifest <- utils::read.csv(
  manifest_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_columns <- c(
  "dataset", "file_role", "relative_path", "size_bytes", "sha256"
)
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns)) {
  stop("Locked input manifest is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (nrow(manifest) != 13L || anyDuplicated(manifest$relative_path)) {
  stop("Locked input manifest must contain 13 unique processed TSV records.")
}

resolve_manifest_input <- function(dataset_name, file_role, relative_path) {
  if (identical(dataset_name, "clinical") && identical(file_role, "clinical")) {
    return(clinical_file)
  }
  if (!dataset_name %in% names(datasets)) {
    stop("Unknown dataset in locked input manifest: ", dataset_name)
  }
  field <- switch(
    file_role,
    expression = "expression_file",
    samples = "samples_file",
    platform = "platform_file",
    stop("Unknown file role in locked input manifest: ", file_role)
  )
  resolved <- datasets[[dataset_name]][[field]]
  if (!is.na(SOURCE_ROOT)) {
    expected_from_relative <- normalizePath(
      file.path(SOURCE_ROOT, relative_path), mustWork = FALSE
    )
    if (!identical(normalizePath(resolved, mustWork = FALSE), expected_from_relative)) {
      stop("Configured dataset path differs from locked logical path: ", relative_path)
    }
  }
  resolved
}

failures <- character()
for (i in seq_len(nrow(manifest))) {
  input_path <- resolve_manifest_input(
    manifest$dataset[i], manifest$file_role[i], manifest$relative_path[i]
  )
  if (is.na(input_path) || !file.exists(input_path)) {
    failures <- c(failures, paste0(manifest$relative_path[i], " (missing)"))
    next
  }
  actual_size <- unname(file.info(input_path)$size)
  expected_size <- as.numeric(manifest$size_bytes[i])
  if (!isTRUE(actual_size == expected_size)) {
    failures <- c(
      failures,
      paste0(
        manifest$relative_path[i], " (size ", actual_size,
        "; expected ", expected_size, ")"
      )
    )
    next
  }
  actual_hash <- digest::digest(
    file = input_path, algo = "sha256", serialize = FALSE
  )
  if (!identical(tolower(actual_hash), tolower(manifest$sha256[i]))) {
    failures <- c(failures, paste0(manifest$relative_path[i], " (SHA-256 mismatch)"))
  }
}

signature_hash <- digest::digest(
  file = SIGNATURE_FILE, algo = "sha256", serialize = FALSE
)
if (!identical(
  tolower(signature_hash), tolower(signature_settings$frozen_sha256)
)) {
  failures <- c(failures, "frozen 39-signature input (SHA-256 mismatch)")
}

if (length(failures)) {
  stop(
    "Controlled inputs do not match the dissertation source lock:\n- ",
    paste(failures, collapse = "\n- ")
  )
}

writeLines(
  c(
    "LOCKED CONTROLLED INPUT VALIDATION: PASS",
    "Processed TSV records verified by size and SHA-256: 13/13",
    "Frozen 39-signature file verified by SHA-256: PASS",
    "No controlled file content or local absolute path is stored in the manifest."
  ),
  file.path(paths$validation, "LOCKED_INPUT_HASH_STATUS.txt")
)
message("LOCKED CONTROLLED INPUT VALIDATION: PASS (13/13)")
