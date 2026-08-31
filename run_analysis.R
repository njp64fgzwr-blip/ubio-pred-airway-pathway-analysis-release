# Run the complete scientific analysis from authorised processed inputs.

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
assert_packages()
require_controlled_inputs()

run_expensive_nested_cv <- isTRUE(config_value(
  c("execution", "run_expensive_nested_cv"), TRUE
))
pipeline <- c(
  "01_validate_locked_inputs.R",
  "02_audit_inputs.R",
  "03_prepare_expression.R",
  "04_run_ssgsea.R",
  "05_merge_clinical_and_baseline.R",
  "06_descriptive_and_spearman.R",
  "07_curated_models.R",
  "08_matched_concordance.R",
  if (run_expensive_nested_cv) "09_nested_validation.R",
  "10_validate_nested_cv.R",
  "11_predictive_inference_and_sensitivity.R"
)

if (!run_expensive_nested_cv &&
    !file.exists(file.path(paths$analysis, "nested_cv_fold_level_predictions.csv.gz"))) {
  stop(
    "Nested cross-validation was disabled, but no existing fold-level outputs ",
    "were found under the configured output root."
  )
}

for (script_name in pipeline) {
  message("\n>>> Running R/", script_name)
  sys.source(file.path("R", script_name), envir = new.env(parent = globalenv()))
}

message(
  "\nScientific analysis complete. Generated controlled outputs are under: ",
  paths$root
)
