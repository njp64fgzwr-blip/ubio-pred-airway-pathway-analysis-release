# Generate dissertation figures and Table 1 from validated analysis outputs.

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

required_analysis_outputs <- c(
  file.path(paths$tables, "Table_2_pathway_clinical_spearman_associations.csv"),
  file.path(paths$tables, "Table_3_curated_model_coefficients_95CI.csv"),
  file.path(paths$tables, "Table_4_repeated_nested_CV_model_performance.csv"),
  file.path(paths$tables, "Table_5_matched_pathway_concordance.csv"),
  file.path(paths$tables, "Table_7_pathway_vs_biomarker_model_comparison.csv")
)
missing_outputs <- required_analysis_outputs[!file.exists(required_analysis_outputs)]
if (length(missing_outputs)) {
  stop(
    "Required analysis outputs are missing. Run run_analysis.R first: ",
    paste(missing_outputs, collapse = "; ")
  )
}

message("\n>>> Running R/25_generate_final_dissertation_figure_set_correct_numbering.R")
sys.source(
  file.path("R", "25_generate_final_dissertation_figure_set_correct_numbering.R"),
  envir = new.env(parent = globalenv())
)

run_pdf_render_qa <- isTRUE(config_value(
  c("execution", "run_pdf_render_qa"), TRUE
))
if (run_pdf_render_qa) {
  for (script_name in c(
    "15_render_all_figure_qa.R",
    "17_validate_figure_sources_and_formats.R"
  )) {
    message("\n>>> Running R/", script_name)
    sys.source(file.path("R", script_name), envir = new.env(parent = globalenv()))
  }
} else {
  warning(
    "PDF render/source QA was skipped by configuration. validate_results.R ",
    "will not report a complete PASS until this QA has been run."
  )
}

for (script_name in c(
  "28_generate_signature_provenance_outputs.R",
  "27_regenerate_table1_correct_numbering.R"
)) {
  message("\n>>> Running R/", script_name)
  sys.source(file.path("R", script_name), envir = new.env(parent = globalenv()))
}

message(
  "\nOutput generation complete. Run validate_results.R for final checks."
)
