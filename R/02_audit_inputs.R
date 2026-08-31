# Audit signatures, clinical variables, sample identifiers and expression scales.
# This script produces metadata/QC only; it does not calculate scientific results.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "digest"))
require_controlled_inputs()
message_rule("02: Auditing signatures, clinical metadata and expression inputs")

# -----------------------------------------------------------------------------
# Fixed 39-signature library
# -----------------------------------------------------------------------------

signature_lines <- readLines(SIGNATURE_FILE, warn = FALSE, encoding = "UTF-8")
signature_fields <- lapply(
  signature_lines,
  function(line) trimws(strsplit(line, "\t", fixed = TRUE)[[1L]])
)
signature_fields <- lapply(signature_fields, function(x) x[nzchar(x)])

if (length(signature_fields) != signature_settings$expected_count) {
  stop(
    "Expected ", signature_settings$expected_count,
    " signature records; found ", length(signature_fields), "."
  )
}
signature_names <- vapply(signature_fields, `[`, character(1), 1L)
if (anyDuplicated(signature_names)) stop("Signature names are not unique.")
if (!setequal(signature_names, names(pathway_labels))) {
  stop("Reader-facing label dictionary does not match the fixed signatures.")
}

gene_sets_source <- stats::setNames(
  lapply(signature_fields, function(x) x[-1L]), signature_names
)
gene_sets_analysis <- lapply(gene_sets_source, unique)

signature_long <- data.table::rbindlist(lapply(seq_along(gene_sets_source), function(i) {
  genes <- gene_sets_source[[i]]
  data.table::data.table(
    pathway_order = i,
    pathway = names(gene_sets_source)[i],
    pathway_label = unname(pathway_labels[names(gene_sets_source)[i]]),
    gene_order = seq_along(genes),
    gene_symbol = genes
  )
}))

signature_audit <- data.table::rbindlist(lapply(seq_along(gene_sets_source), function(i) {
  genes <- gene_sets_source[[i]]
  data.table::data.table(
    pathway_order = i,
    pathway = names(gene_sets_source)[i],
    pathway_label = unname(pathway_labels[names(gene_sets_source)[i]]),
    source_gene_entries = length(genes),
    unique_gene_symbols = length(unique(genes)),
    duplicate_entries = sum(duplicated(genes)),
    blank_entries = sum(!nzchar(genes))
  )
}))

saveRDS(
  gene_sets_analysis,
  file.path(paths$analysis, "airway_39_signatures_analysis.rds"),
  compress = "gzip"
)
data.table::fwrite(
  signature_long,
  file.path(paths$tables, "Table_S_signature_gene_definitions.csv")
)
data.table::fwrite(
  signature_audit,
  file.path(paths$validation, "SIGNATURE_INTEGRITY_AUDIT.csv")
)

if (sha256_file(SIGNATURE_FILE) != signature_settings$frozen_sha256) {
  stop("Signature SHA-256 changed after source lock.")
}

# -----------------------------------------------------------------------------
# Clinical variables
# -----------------------------------------------------------------------------

clinical_raw <- data.table::fread(
  clinical_file, check.names = FALSE,
  na.strings = c("", "NA", "NaN", "null", "NULL"), showProgress = FALSE
)
if (!"Subject ID" %in% names(clinical_raw)) stop("Clinical Subject ID missing.")
if (anyDuplicated(clinical_raw[["Subject ID"]])) {
  stop("Clinical Subject ID is not unique.")
}

clinical_columns <- list(
  cohort = find_unique_column(names(clinical_raw), c("Study Groups", "cohort")),
  age = find_unique_column(
    names(clinical_raw), c("Demographic Data", "\\Age\\"), "Estimated Age"
  ),
  sex = find_unique_column(names(clinical_raw), c("Demographic Data", "\\Sex\\")),
  smoking_status = find_unique_column(
    names(clinical_raw), c("Smoking History", "Smoking Status")
  ),
  bmi = find_unique_column(
    names(clinical_raw), c("Subject Body Measurements", "Body Mass Index")
  ),
  ocs_frequency = find_unique_column(
    names(clinical_raw), c("Medication", "Baseline", "Oral Corticosteroids")
  ),
  ocs_dose_mg = find_unique_column(
    names(clinical_raw), c("Medication", "Baseline", "OCS Normalised Dose")
  ),
  fev1_percent_predicted = find_unique_column(
    names(clinical_raw),
    c("Clinical Data", "Lung Function", "Spirometry", "Baseline", "FEV1 % (L)")
  ),
  feno_ppb = find_unique_column(
    names(clinical_raw),
    c("Exhaled nitric oxide (NO)", "Baseline", "NO Standard Flow Rate")
  ),
  blood_eosinophils_x10_3_uL = find_unique_column(
    names(clinical_raw),
    c("Haematology and biochemistry tests", "Screening", "eosinophils (x10^3/uL)")
  ),
  blood_eosinophils_percent = find_unique_column(
    names(clinical_raw),
    c("Haematology and biochemistry tests", "Screening", "Eosinophils %")
  ),
  sputum_eosinophils_percent = find_unique_column(
    names(clinical_raw),
    c("Sputum Data", "Baseline Visit", "Baseline Visit Master list for analyses", "% Eosinophils")
  ),
  exacerbations_previous_12m = find_unique_column(
    names(clinical_raw),
    c("Recent Asthma Exacerbation History", "Baseline", "Exacerbation Per Year"),
    "Severe Exacerbation Per Year"
  )
)

clinical_clean <- data.table::data.table(
  Subject_ID = as.character(clinical_raw[["Subject ID"]]),
  cohort_source = as.character(clinical_raw[[clinical_columns$cohort]]),
  age = suppressWarnings(as.numeric(clinical_raw[[clinical_columns$age]])),
  sex = as.character(clinical_raw[[clinical_columns$sex]]),
  smoking_status = as.character(clinical_raw[[clinical_columns$smoking_status]]),
  bmi = suppressWarnings(as.numeric(clinical_raw[[clinical_columns$bmi]])),
  ocs_frequency_source = as.character(clinical_raw[[clinical_columns$ocs_frequency]]),
  ocs_dose_mg = suppressWarnings(as.numeric(clinical_raw[[clinical_columns$ocs_dose_mg]])),
  fev1_percent_predicted = suppressWarnings(as.numeric(
    clinical_raw[[clinical_columns$fev1_percent_predicted]]
  )),
  feno_ppb = suppressWarnings(as.numeric(clinical_raw[[clinical_columns$feno_ppb]])),
  blood_eosinophils_x10_3_uL = suppressWarnings(as.numeric(
    clinical_raw[[clinical_columns$blood_eosinophils_x10_3_uL]]
  )),
  blood_eosinophils_percent = suppressWarnings(as.numeric(
    clinical_raw[[clinical_columns$blood_eosinophils_percent]]
  )),
  sputum_eosinophils_percent = suppressWarnings(as.numeric(
    clinical_raw[[clinical_columns$sputum_eosinophils_percent]]
  )),
  exacerbations_previous_12m = suppressWarnings(as.numeric(
    clinical_raw[[clinical_columns$exacerbations_previous_12m]]
  ))
)

cohort_labels <- c(
  cohort_a = "Severe asthma, <5 pack-years",
  cohort_b = "Severe asthma, >=5 pack-years/current or ex-smoking",
  cohort_c = "Mild-to-moderate asthma",
  cohort_d = "Healthy comparator",
  cohort_v = "Source-coded group V (unclassified)"
)
clinical_clean[, cohort_label := unname(cohort_labels[cohort_source])]

current_ocs_levels <- c(
  "daily", "at_least_twice_a_day", "more_than_twice_week",
  "weekly", "more_than_once_a_month", "once_a_month"
)
clinical_clean[, ocs_current := data.table::fcase(
  ocs_frequency_source %in% current_ocs_levels, "Yes",
  ocs_frequency_source %in% c("never", "Previous"), "No",
  default = NA_character_
)]
clinical_clean[, ocs_dose_recorded := ifelse(is.finite(ocs_dose_mg), "Yes", "No")]

numeric_nonnegative <- c(
  "age", "bmi", "ocs_dose_mg", "fev1_percent_predicted", "feno_ppb",
  "blood_eosinophils_x10_3_uL", "blood_eosinophils_percent",
  "sputum_eosinophils_percent", "exacerbations_previous_12m"
)
for (variable in numeric_nonnegative) {
  if (any(clinical_clean[[variable]] < 0, na.rm = TRUE)) {
    stop("Negative clinical values found for ", variable, ".")
  }
}

data.table::fwrite(
  clinical_clean, file.path(paths$analysis, "clinical_analysis_data.csv")
)

clinical_dictionary <- data.table::data.table(
  analysis_variable = names(clinical_columns),
  source_column = unlist(clinical_columns, use.names = FALSE)
)
clinical_dictionary <- data.table::rbindlist(list(
  clinical_dictionary,
  data.table::data.table(
    analysis_variable = c("cohort_label", "ocs_current", "ocs_dose_recorded"),
    source_column = c(
      "Derived from Study Groups/cohort using documented U-BIOPRED group labels",
      paste0(
        "Derived from baseline Oral Corticosteroids frequency; Yes for: ",
        paste(current_ocs_levels, collapse = ", "),
        "; No for never/Previous; missing and Free_text remain missing"
      ),
      "Yes when baseline OCS Normalised Dose (mg) is present"
    )
  )
), fill = TRUE)
data.table::fwrite(
  clinical_dictionary,
  file.path(paths$tables, "Table_S_clinical_variable_dictionary.csv")
)

clinical_audit <- data.table::rbindlist(lapply(
  setdiff(names(clinical_clean), "Subject_ID"),
  function(variable) {
    x <- clinical_clean[[variable]]
    data.table::data.table(
      variable = variable,
      storage_class = class(x)[1L],
      total_n = length(x),
      available_n = sum(!is.na(x)),
      missing_n = sum(is.na(x)),
      unique_nonmissing = data.table::uniqueN(x, na.rm = TRUE),
      minimum = if (is.numeric(x)) suppressWarnings(min(x, na.rm = TRUE)) else NA_real_,
      maximum = if (is.numeric(x)) suppressWarnings(max(x, na.rm = TRUE)) else NA_real_
    )
  }
), fill = TRUE)
clinical_audit[!is.finite(minimum), minimum := NA_real_]
clinical_audit[!is.finite(maximum), maximum := NA_real_]
data.table::fwrite(
  clinical_audit, file.path(paths$validation, "CLINICAL_VARIABLE_AUDIT.csv")
)

# -----------------------------------------------------------------------------
# Sample metadata and expression scale
# -----------------------------------------------------------------------------

sample_rows <- list()
expression_rows <- list()
all_subject_sets <- list()

for (dataset_name in names(datasets)) {
  specification <- datasets[[dataset_name]]
  samples <- data.table::fread(
    specification$samples_file, check.names = FALSE, showProgress = FALSE
  )
  samples[, `Assay ID` := as.character(`Assay ID`)]
  samples[, `Subject ID` := as.character(`Subject ID`)]

  if (nrow(samples) != specification$expected_n) {
    stop(dataset_name, ": unexpected sample count.")
  }
  if (anyDuplicated(samples[["Assay ID"]]) || anyDuplicated(samples[["Subject ID"]])) {
    stop(dataset_name, ": duplicate assay or participant ID.")
  }
  unmatched <- setdiff(samples[["Subject ID"]], clinical_clean$Subject_ID)
  if (length(unmatched)) {
    stop(dataset_name, ": sample IDs do not match clinical IDs: ",
         paste(unmatched, collapse = ", "))
  }

  all_subject_sets[[dataset_name]] <- samples[["Subject ID"]]
  sample_rows[[dataset_name]] <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type,
    sample_n = nrow(samples),
    unique_assay_n = data.table::uniqueN(samples[["Assay ID"]]),
    unique_subject_n = data.table::uniqueN(samples[["Subject ID"]]),
    sample_type = paste(unique(samples[["Sample Type"]]), collapse = ";"),
    time_point = paste(unique(samples[["Time Point"]]), collapse = ";"),
    tissue_type = paste(unique(samples[["Tissue Type"]]), collapse = ";"),
    unmatched_clinical_ids = length(unmatched),
    status = "PASS"
  )

  message("Auditing expression scale for ", specification$label, "...")
  expression <- data.table::fread(
    specification$expression_file,
    select = c("Assay ID", "VALUE", "LOG2E", "PROBE", "GENE SYMBOL"),
    na.strings = c("", "NA", "NaN", "null", "NULL"),
    showProgress = FALSE
  )
  expression[, `Assay ID` := as.character(`Assay ID`)]
  relationship_error <- abs(expression$LOG2E - log2(expression$VALUE))
  expression_rows[[dataset_name]] <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    technology = specification$technology,
    input_rows = nrow(expression),
    assay_n = data.table::uniqueN(expression[["Assay ID"]]),
    feature_n = data.table::uniqueN(expression$PROBE),
    supplied_symbol_n = data.table::uniqueN(expression[["GENE SYMBOL"]], na.rm = TRUE),
    missing_symbol_rows = sum(is.na(expression[["GENE SYMBOL"]])),
    value_min = min(expression$VALUE, na.rm = TRUE),
    value_max = max(expression$VALUE, na.rm = TRUE),
    log2e_min = min(expression$LOG2E, na.rm = TRUE),
    log2e_max = max(expression$LOG2E, na.rm = TRUE),
    max_abs_log2e_minus_log2_value = max(relationship_error, na.rm = TRUE),
    assay_ids_match_samples = setequal(unique(expression[["Assay ID"]]), samples[["Assay ID"]]),
    selected_analysis_input = if (specification$technology == "microarray") {
      "Supplied LOG2E; median collapse of probes to gene symbol"
    } else {
      "Sparse raw VALUE counts; reconstruct zeros, aggregate symbols, CPM filter, log2(CPM+1)"
    },
    supplied_zscore_used = FALSE
  )
  rm(expression, samples)
  invisible(gc())
}

sample_audit <- data.table::rbindlist(sample_rows)
expression_audit <- data.table::rbindlist(expression_rows)
data.table::fwrite(
  sample_audit, file.path(paths$validation, "SAMPLE_ID_STRUCTURE_AUDIT.csv")
)
data.table::fwrite(
  expression_audit, file.path(paths$validation, "EXPRESSION_SCALE_AUDIT.csv")
)

all_subjects <- sort(unique(unlist(all_subject_sets, use.names = FALSE)))
membership <- data.table::data.table(Subject_ID = all_subjects)
for (dataset_name in names(all_subject_sets)) {
  membership[[dataset_name]] <- membership$Subject_ID %in% all_subject_sets[[dataset_name]]
}
membership[, exact_pattern := apply(
  as.matrix(.SD), 1L, function(x) paste(names(all_subject_sets)[x], collapse = "+")
), .SDcols = names(all_subject_sets)]
data.table::fwrite(
  membership, file.path(paths$analysis, "dataset_participant_membership.csv")
)

exact_overlap <- membership[, .(participants = .N), by = exact_pattern][
  order(-participants, exact_pattern)
]
data.table::fwrite(
  exact_overlap, file.path(paths$tables, "Table_S_exact_dataset_overlap.csv")
)

pair_rows <- list()
dataset_names <- names(all_subject_sets)
for (i in seq_along(dataset_names)) {
  for (j in seq_along(dataset_names)) {
    if (j <= i) next
    left <- dataset_names[i]
    right <- dataset_names[j]
    pair_rows[[length(pair_rows) + 1L]] <- data.table::data.table(
      dataset_1 = left,
      dataset_2 = right,
      matched_participants = length(intersect(all_subject_sets[[left]], all_subject_sets[[right]]))
    )
  }
}
pair_overlap <- data.table::rbindlist(pair_rows)
data.table::fwrite(
  pair_overlap, file.path(paths$tables, "Table_S_pairwise_dataset_overlap.csv")
)

raw_instrument_files <- list.files(
  SOURCE_ROOT,
  pattern = "[.](CEL|fastq|fq|bam|cram)([.]gz)?$",
  recursive = TRUE, full.names = TRUE, ignore.case = TRUE
)

input_status <- c(
  "INPUT AUDIT: PASS",
  paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste("Clinical participants:", nrow(clinical_clean)),
  paste("Unique participants with >=1 airway transcriptome:", length(all_subjects)),
  paste("Fixed signatures:", length(gene_sets_analysis)),
  paste("Raw CEL/FASTQ/BAM/CRAM files found:", length(raw_instrument_files)),
  "All airway sample IDs matched unique clinical Subject IDs.",
  "Microarrays: supplied LOG2E selected; supplied clipped ZSCORE not used.",
  "RNA-seq: sparse VALUE counts selected for R-based CPM normalization.",
  "Eicosanoid and other non-transcriptomic data are excluded from ssGSEA."
)
writeLines(input_status, file.path(paths$validation, "INPUT_AUDIT_STATUS.txt"))
message("Input audit completed successfully.")
