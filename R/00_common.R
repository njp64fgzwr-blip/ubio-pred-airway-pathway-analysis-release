# Shared configuration and helper functions for the portable U-BIOPRED
# airway-pathway analysis. All scientific scripts source this file first.

options(stringsAsFactors = FALSE, warn = 1)

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(
      current, "ubio-pred-airway-pathway-analysis.Rproj"
    ))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not locate ubio-pred-airway-pathway-analysis.Rproj from: ",
        start
      )
    }
    current <- parent
  }
}

PROJECT_ROOT <- find_project_root()

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required to read config/analysis_parameters.yml.")
}

config_file <- Sys.getenv(
  "UBIOPRED_ANALYSIS_CONFIG",
  unset = file.path(PROJECT_ROOT, "config", "analysis_parameters.yml")
)
if (!file.exists(config_file)) {
  stop("Analysis configuration file is missing: ", config_file)
}
analysis_config <- yaml::read_yaml(config_file)

config_value <- function(path, default = NULL) {
  value <- analysis_config
  for (key in path) {
    if (is.null(value) || is.null(value[[key]])) return(default)
    value <- value[[key]]
  }
  value
}

resolve_config_path <- function(environment_variable, configured_value,
                                default = "") {
  value <- Sys.getenv(environment_variable, unset = "")
  if (!nzchar(value)) value <- configured_value %||% default
  if (is.null(value) || !length(value) || !nzchar(as.character(value)[1L])) {
    return(NA_character_)
  }
  value <- path.expand(as.character(value)[1L])
  is_absolute <- grepl("^/", value) || grepl("^[A-Za-z]:[/\\\\]", value)
  if (!is_absolute) value <- file.path(PROJECT_ROOT, value)
  normalizePath(value, mustWork = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

SOURCE_ROOT <- resolve_config_path(
  "UBIOPRED_SOURCE_ROOT",
  config_value(c("inputs", "ubio_pred_source_root"), "")
)
SIGNATURE_FILE <- resolve_config_path(
  "AIRWAY_SIGNATURE_FILE",
  config_value(c("inputs", "signature_file"), "")
)
SIGNATURE_PROVENANCE_FILE <- resolve_config_path(
  "AIRWAY_SIGNATURE_PROVENANCE_FILE",
  config_value(c("inputs", "signature_provenance_file"), "")
)
OUTPUT_ROOT <- resolve_config_path(
  "UBIOPRED_OUTPUT_ROOT",
  config_value(c("execution", "output_root"), "work/full_analysis"),
  default = "work/full_analysis"
)
if (is.na(OUTPUT_ROOT)) {
  stop("An output root must be supplied in the configuration or environment.")
}

# Prevent a full-data run from writing participant-level derivatives into a
# tracked repository location. Within this project, controlled outputs may be
# written only beneath the ignored work/ directory. An external authorised
# location remains permitted.
project_path <- normalizePath(PROJECT_ROOT, mustWork = TRUE)
output_path <- normalizePath(OUTPUT_ROOT, mustWork = FALSE)
project_prefix <- paste0(project_path, .Platform$file.sep)
if (startsWith(output_path, project_prefix)) {
  relative_output <- substring(output_path, nchar(project_prefix) + 1L)
  work_prefix <- paste0("work", .Platform$file.sep)
  if (!identical(relative_output, "work") && !startsWith(relative_output, work_prefix)) {
    stop(
      "Within the repository, UBIOPRED_OUTPUT_ROOT must be beneath the ignored ",
      "work/ directory. Use an authorised external directory otherwise."
    )
  }
}

paths <- list(
  project = PROJECT_ROOT,
  root = OUTPUT_ROOT,
  scripts = file.path(PROJECT_ROOT, "R"),
  config = file.path(PROJECT_ROOT, "config"),
  manifest = file.path(OUTPUT_ROOT, "data_manifest"),
  results = file.path(OUTPUT_ROOT, "results"),
  processed = file.path(OUTPUT_ROOT, "results", "processed_expression"),
  scores = file.path(OUTPUT_ROOT, "results", "ssgsea"),
  analysis = file.path(OUTPUT_ROOT, "results", "analysis_data"),
  tables = file.path(OUTPUT_ROOT, "results", "tables"),
  models = file.path(OUTPUT_ROOT, "results", "models"),
  diagnostics = file.path(OUTPUT_ROOT, "results", "diagnostics"),
  figures_main = file.path(OUTPUT_ROOT, "figures", "main"),
  figures_supp = file.path(OUTPUT_ROOT, "figures", "supplementary"),
  figures_qa = file.path(OUTPUT_ROOT, "figures", "qa_renders"),
  validation = file.path(OUTPUT_ROOT, "validation"),
  reproducibility = file.path(OUTPUT_ROOT, "reproducibility"),
  logs = file.path(OUTPUT_ROOT, "logs"),
  final_figures = file.path(
    OUTPUT_ROOT, "FINAL_DISSERTATION_FIGURES_R_ONLY_CORRECT_NUMBERING"
  ),
  final_tables = file.path(
    OUTPUT_ROOT, "FINAL_DISSERTATION_TABLES_R_ONLY_CORRECT_NUMBERING"
  )
)

output_directories <- unique(unlist(paths[c(
  "root", "manifest", "results", "processed", "scores", "analysis",
  "tables", "models", "diagnostics", "figures_main", "figures_supp",
  "figures_qa", "validation", "reproducibility", "logs", "final_figures",
  "final_tables"
)]))
for (directory in output_directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

transcriptomics_root <- file.path(
  SOURCE_ROOT, "OMICS+DATA", "Transcriptomics", "Baseline+Visit"
)

# Canonical display and processing order used throughout every figure, table
# and export that shows more than one airway dataset. Keep the two GPL570
# airway compartments together before the RNA-seq technical comparison.
canonical_dataset_order <- unname(unlist(config_value(
  c("datasets", "order"),
  c("sputum", "brushing_gpl570", "biopsy", "brushing_rnaseq")
)))
expected_sample_counts <- config_value(
  c("datasets", "expected_sample_counts"),
  list(sputum = 120L, brushing_gpl570 = 149L, biopsy = 108L,
       brushing_rnaseq = 118L)
)

datasets <- list(
  sputum = list(
    label = "Sputum GPL570",
    short_label = "Sputum",
    compartment = "Sputum",
    platform = "GPL570",
    omic_type = "Transcriptomics",
    technology = "microarray",
    expected_n = as.integer(expected_sample_counts$sputum),
    directory = file.path(transcriptomics_root, "Sputum", "GPL570_SPUTUM")
  ),
  brushing_gpl570 = list(
    label = "Bronchial brushing GPL570",
    short_label = "Brushing GPL570",
    compartment = "Bronchial brushing",
    platform = "GPL570",
    omic_type = "Transcriptomics",
    technology = "microarray",
    expected_n = as.integer(expected_sample_counts$brushing_gpl570),
    directory = file.path(
      transcriptomics_root, "Bronchial+Brushings", "GPL570_BRBRUSHING"
    )
  ),
  biopsy = list(
    label = "Bronchial biopsy GPL570",
    short_label = "Biopsy",
    compartment = "Bronchial biopsy",
    platform = "GPL570",
    omic_type = "Transcriptomics",
    technology = "microarray",
    expected_n = as.integer(expected_sample_counts$biopsy),
    directory = file.path(
      transcriptomics_root, "Bronchial+Biopsy", "GPL570_BIOPSY"
    )
  ),
  brushing_rnaseq = list(
    label = "Bronchial brushing RNA-seq",
    short_label = "Brushing RNA-seq",
    compartment = "Bronchial brushing",
    platform = "RNA-seq",
    omic_type = "Transcriptomics",
    technology = "rnaseq",
    expected_n = as.integer(expected_sample_counts$brushing_rnaseq),
    directory = file.path(
      transcriptomics_root, "Bronchial+Brushings",
      "BRONCHIAL_BRUSHING_RNASEQ"
    )
  )
)

if (!identical(names(datasets), canonical_dataset_order)) {
  stop("Dataset definitions no longer match the canonical display order.")
}

for (dataset_name in names(datasets)) {
  datasets[[dataset_name]]$expression_file <- file.path(
    datasets[[dataset_name]]$directory, "data_mrna.tsv"
  )
  datasets[[dataset_name]]$samples_file <- file.path(
    datasets[[dataset_name]]$directory, "samples.tsv"
  )
  datasets[[dataset_name]]$platform_file <- file.path(
    datasets[[dataset_name]]$directory, "platform.tsv"
  )
}

clinical_file <- resolve_config_path(
  "UBIOPRED_CLINICAL_FILE",
  config_value(c("inputs", "clinical_file"), "")
)
if (is.na(clinical_file) && !is.na(SOURCE_ROOT)) {
  clinical_file <- file.path(SOURCE_ROOT, "data_clinical.tsv")
}

required_controlled_input_files <- function() {
  dataset_files <- unlist(lapply(datasets, function(x) {
    c(x$expression_file, x$samples_file, x$platform_file)
  }), use.names = FALSE)
  c(clinical_file, SIGNATURE_FILE, dataset_files)
}

require_controlled_inputs <- function() {
  files <- required_controlled_input_files()
  missing_configuration <- is.na(files) | !nzchar(files)
  if (any(missing_configuration)) {
    stop(
      "Controlled input paths are not configured. Set UBIOPRED_SOURCE_ROOT ",
      "and AIRWAY_SIGNATURE_FILE, or edit config/analysis_parameters.yml."
    )
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) {
    stop(
      "Required controlled input files are unavailable: ",
      paste(missing_files, collapse = "; ")
    )
  }
  invisible(TRUE)
}

require_signature_provenance_input <- function() {
  if (is.na(SIGNATURE_PROVENANCE_FILE) ||
      !nzchar(SIGNATURE_PROVENANCE_FILE) ||
      !file.exists(SIGNATURE_PROVENANCE_FILE)) {
    stop(
      "The audited signature-provenance CSV is unavailable. Set ",
      "AIRWAY_SIGNATURE_PROVENANCE_FILE or inputs.signature_provenance_file."
    )
  }
  invisible(TRUE)
}

outcome_specs <- list(
  FEV1 = list(
    column = "fev1_percent_predicted",
    label = "FEV1 % predicted",
    short_label = "FEV1",
    transform = "identity",
    original_unit = "% predicted"
  ),
  FeNO = list(
    column = "feno_ppb",
    label = "FeNO",
    short_label = "FeNO",
    transform = "log1p",
    original_unit = "ppb"
  ),
  blood_eosinophils = list(
    column = "blood_eosinophils_x10_3_uL",
    label = "Blood eosinophils",
    short_label = "Blood eosinophils",
    transform = "log1p",
    original_unit = "x10^3/uL"
  ),
  sputum_eosinophils = list(
    column = "sputum_eosinophils_percent",
    label = "Sputum eosinophils",
    short_label = "Sputum eosinophils",
    transform = "log1p",
    original_unit = "%"
  ),
  exacerbation_frequency = list(
    column = "exacerbations_previous_12m",
    label = "Exacerbations in previous 12 months",
    short_label = "Prior-year exacerbations",
    transform = "log1p",
    original_unit = "events/year"
  )
)

# These pathway subsets were prespecified before examining the present rerun.
curated_pathway_models <- list(
  FEV1 = c(
    "IL13_signalling_ref26", "Eosinophils", "IL17_signalling_ref27",
    "Neutrophil", "ECM_organisation"
  ),
  FeNO = c(
    "IL13_signalling_ref26", "IL5_signalling_ref24", "Eosinophils",
    "Mast_cell", "IL33_signalling_ref28"
  ),
  blood_eosinophils = c(
    "Eosinophils", "IL5_signalling_ref24", "TAC1",
    "ILC2_up_PMID26878113", "Severe_asthma_eosinophil_signature"
  ),
  sputum_eosinophils = c(
    "Eosinophils", "TAC1", "IL5_signalling_ref24",
    "IL13_signalling_ref26", "Mast_cell"
  ),
  exacerbation_frequency = c(
    "IL13_signalling_ref26", "IL17_signalling_ref27", "Neutrophil",
    "Macrophage_GMCSF_TNFa_HS_IVS_UP", "ECM_organisation", "Inflammasome"
  )
)

biomarker_models <- list(
  FEV1 = c("feno_ppb", "blood_eosinophils_x10_3_uL", "sputum_eosinophils_percent"),
  FeNO = c("blood_eosinophils_x10_3_uL", "sputum_eosinophils_percent"),
  blood_eosinophils = c("feno_ppb", "sputum_eosinophils_percent"),
  sputum_eosinophils = c("feno_ppb", "blood_eosinophils_x10_3_uL"),
  exacerbation_frequency = c(
    "feno_ppb", "blood_eosinophils_x10_3_uL", "sputum_eosinophils_percent"
  )
)

model_types <- c(
  "Biomarker-only", "Curated pathways", "Elastic Net pathways",
  "Biomarker + curated pathways", "Biomarker + Elastic Net pathways"
)

cv_settings <- list(
  outer_repeats = as.integer(config_value(
    c("cross_validation", "outer_repeats"), 10L
  )),
  outer_folds = as.integer(config_value(
    c("cross_validation", "outer_folds"), 5L
  )),
  inner_folds = as.integer(config_value(
    c("cross_validation", "inner_folds"), 5L
  )),
  alpha_grid = as.numeric(unlist(config_value(
    c("cross_validation", "alpha_grid"),
    c(0.10, 0.25, 0.50, 0.75, 1.00)
  ))),
  lambda_fraction_grid = {
    exponents <- as.numeric(unlist(config_value(
      c("cross_validation", "lambda_fraction_exponents"), c(0, -4)
    )))
    count <- as.integer(config_value(
      c("cross_validation", "lambda_fraction_count"), 100L
    ))
    10^seq(exponents[1L], exponents[2L], length.out = count)
  },
  seed = as.integer(config_value(
    c("cross_validation", "seed"), 20260813L
  )),
  tuning_rule = paste(
    "For each alpha, choose the largest lambda fraction within one standard",
    "error of its minimum mean inner-fold MSE; then choose the alpha with",
    "the lowest pooled MSE at its one-SE lambda. Losses within a relative",
    "floating-point tolerance of 1e-12 are treated as ties, resolved toward",
    "larger alpha and then stronger penalisation."
  )
)

signature_settings <- list(
  expected_count = as.integer(config_value(
    c("signatures", "expected_count"), 39L
  )),
  frozen_sha256 = as.character(config_value(
    c("signatures", "frozen_sha256"),
    "00aaf4a958729605fd2005294b36f944699067f10cc3e160809f21ac12c8a3d8"
  )),
  provenance_sha256 = as.character(config_value(
    c("signatures", "provenance_sha256"),
    "5b572e7eecf2c410d8ad677fe6c4a0301fc2efc001cf9c32bd2389890c70cefc"
  ))
)

ssgsea_settings <- list(
  min_size = as.integer(config_value(c("ssgsea", "min_size"), 3L)),
  alpha = as.numeric(config_value(c("ssgsea", "alpha"), 0.25)),
  normalize = isTRUE(config_value(c("ssgsea", "normalize"), TRUE))
)

preprocessing_settings <- list(
  rnaseq_cpm_threshold = as.numeric(config_value(
    c("expression_preprocessing", "rnaseq_cpm_threshold"), 1
  )),
  rnaseq_minimum_sample_fraction = as.numeric(config_value(
    c("expression_preprocessing", "rnaseq_minimum_sample_fraction"), 0.10
  ))
)

predictive_inference_settings <- list(
  bootstrap_replicates = as.integer(config_value(
    c("predictive_inference", "conditional_participant_bootstrap_replicates"),
    2000L
  )),
  bootstrap_seed = as.integer(config_value(
    c("predictive_inference", "bootstrap_seed"), 20260814L
  )),
  calibration_bins = as.integer(config_value(
    c("predictive_inference", "calibration_bins"), 5L
  ))
)

pathway_labels <- c(
  IL13_stimulated_HBECs = "IL-13-stimulated HBECs",
  ILC2_up_PMID26878113 = "ILC2",
  Th2_activated = "Activated Th2",
  Eosinophils = "Eosinophils",
  TAC1 = "TAC1",
  Severe_asthma_eosinophil_signature = "Severe-asthma eosinophil",
  IL5_signalling_ref24 = "IL-5 signalling",
  IL13_signalling_ref26 = "IL-13 signalling",
  ILC1_cells = "ILC1",
  Inflammasome = "Inflammasome",
  OXPHOS = "Oxidative phosphorylation",
  Neutrophil = "Neutrophils",
  LPS_stimulated_macrophage = "LPS-stimulated macrophage",
  TAC2 = "TAC2",
  NETosis = "NETosis",
  IL6_TS = "IL-6 trans-signalling",
  IL1R_family_related = "IL-1 receptor family",
  Mast_cell = "Mast cells",
  IgE_FcepsilonR1_Activated_MC_ref_3 = "IgE/FcER1-activated mast cells",
  Repeated_FcepsilonR1_Stimulation_ref_4 = "Repeated FcER1 stimulation",
  Mast_cell_stimulated_by_IL33_ref7 = "IL-33-stimulated mast cells",
  Mast_cell_asthma_MCscbb_ref6 = "Asthma mast-cell state",
  Th17_Zhang = "Th17",
  IL17A = "IL-17A response",
  IL17_signalling_ref27 = "IL-17 signalling",
  ILC3_up_PMID26878113 = "ILC3",
  M1_macrophage_ref17 = "M1 macrophage",
  M2_macrophage_ref17 = "M2 macrophage",
  Macrophage_GMCSF_IL4 = "GM-CSF/IL-4 macrophage",
  Macrophage_GMCSF_TNFa_HS_IVS_UP = "GM-CSF/TNF-alpha macrophage",
  Treg_activated = "Activated Treg",
  CD4 = "CD4 T cells",
  TAC3 = "TAC3",
  ECM_organisation = "ECM organisation",
  Mucus = "Mucus",
  P53_pathway = "p53 pathway",
  SASP = "SASP",
  Ageing = "Ageing",
  IL33_signalling_ref28 = "IL-33 signalling"
)

required_packages <- c(
  "yaml", "data.table", "R.utils", "digest", "GSVA", "BiocParallel", "AnnotationDbi",
  "hgu133plus2.db", "org.Hs.eg.db", "glmnet", "ggplot2", "patchwork",
  "openxlsx", "MASS", "broom", "car", "lmtest", "svglite", "ragg",
  "gridExtra", "magick", "scales", "stringr"
)

assert_packages <- function(packages = required_packages) {
  missing <- packages[!vapply(
    packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
  )]
  if (length(missing)) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

message_rule <- function(title) {
  message("\n", paste(rep("=", 76), collapse = ""))
  message(title)
  message(paste(rep("=", 76), collapse = ""))
}

sha256_file <- function(file) {
  digest::digest(file = file, algo = "sha256", serialize = FALSE)
}

safe_zscore_rows <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1, stats::sd, na.rm = TRUE)
  sds[!is.finite(sds) | sds == 0] <- 1
  sweep(sweep(x, 1, means, "-"), 1, sds, "/")
}

apply_outcome_transform <- function(x, method) {
  x <- suppressWarnings(as.numeric(x))
  if (identical(method, "identity")) return(x)
  if (identical(method, "log1p")) {
    if (any(x < 0, na.rm = TRUE)) stop("Negative value supplied to log1p.")
    return(log1p(x))
  }
  stop("Unknown outcome transformation: ", method)
}

inverse_outcome_transform <- function(x, method) {
  if (identical(method, "identity")) return(x)
  if (identical(method, "log1p")) return(expm1(x))
  stop("Unknown outcome transformation: ", method)
}

matrix_to_subject_table <- function(x) {
  data.table::as.data.table(
    data.frame(Subject_ID = colnames(x), t(x), check.names = FALSE)
  )
}

find_unique_column <- function(column_names, required, excluded = character()) {
  matches <- column_names
  for (fragment in required) {
    matches <- matches[grepl(fragment, matches, fixed = TRUE)]
  }
  for (fragment in excluded) {
    matches <- matches[!grepl(fragment, matches, fixed = TRUE)]
  }
  if (length(matches) != 1L) {
    stop(
      "Expected exactly one column containing [",
      paste(required, collapse = " | "), "] and excluding [",
      paste(excluded, collapse = " | "), "]; found ", length(matches), "."
    )
  }
  matches
}

calculate_metrics <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- observed[valid]
  predicted <- predicted[valid]
  denominator <- sum((observed - mean(observed))^2)
  c(
    r_squared = if (denominator > 0) {
      1 - sum((observed - predicted)^2) / denominator
    } else NA_real_,
    rmse = sqrt(mean((observed - predicted)^2)),
    mae = mean(abs(observed - predicted))
  )
}

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  n <- length(y)
  if (n < k) stop("Fewer observations than requested folds.")
  breaks <- unique(stats::quantile(
    y, probs = seq(0, 1, length.out = min(6L, length(unique(y)) + 1L)),
    na.rm = TRUE, type = 2
  ))
  if (length(breaks) < 3L) return(sample(rep(seq_len(k), length.out = n)))
  bins <- cut(y, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  folds <- integer(n)
  for (bin in sort(unique(bins))) {
    indexes <- which(bins == bin)
    folds[indexes] <- sample(rep(seq_len(k), length.out = length(indexes)))
  }
  folds
}

dataset_metadata_table <- function() {
  data.table::rbindlist(lapply(names(datasets), function(nm) {
    x <- datasets[[nm]]
    data.table::data.table(
      dataset = nm, dataset_label = x$label, compartment = x$compartment,
      platform = x$platform, omic_type = x$omic_type,
      technology = x$technology, expected_n = x$expected_n
    )
  }))
}
