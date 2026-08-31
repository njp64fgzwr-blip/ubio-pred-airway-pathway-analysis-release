# Merge normalized ssGSEA scores with clinical metadata and build baseline tables.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "openxlsx"))
message_rule("05: Merging clinical data and building baseline tables")

clinical <- data.table::fread(
  file.path(paths$analysis, "clinical_analysis_data.csv"),
  check.names = FALSE
)
pathways <- names(readRDS(
  file.path(paths$analysis, "airway_39_signatures_harmonised.rds")
))

sample_rows <- list()
missingness_rows <- list()
baseline_raw_rows <- list()
baseline_formatted_rows <- list()
cohort_rows <- list()

format_mean_sd <- function(x, digits = 1L) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  sprintf(paste0("%.", digits, "f (%.", digits, "f); n=%d"),
          mean(x), stats::sd(x), length(x))
}

format_median_iqr <- function(x, digits = 1L) {
  x <- x[is.finite(x)]
  if (!length(x)) return("NA")
  q <- stats::quantile(x, c(0.25, 0.50, 0.75), names = FALSE, type = 2)
  sprintf(paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]; n=%d"),
          q[2L], q[1L], q[3L], length(x))
}

format_n_percent <- function(condition, denominator = length(condition)) {
  n <- sum(condition, na.rm = TRUE)
  available <- sum(!is.na(condition))
  if (!available) return("NA")
  sprintf("%d (%.1f%%); n=%d", n, 100 * n / available, available)
}

add_baseline_row <- function(dataset_name, variable, statistic, value, n,
                             mean = NA_real_, sd = NA_real_,
                             median = NA_real_, q1 = NA_real_, q3 = NA_real_,
                             count = NA_integer_, percent = NA_real_) {
  baseline_raw_rows[[length(baseline_raw_rows) + 1L]] <<- data.table::data.table(
    dataset = dataset_name,
    variable = variable,
    statistic = statistic,
    value = value,
    n = n,
    mean = mean,
    sd = sd,
    median = median,
    q1 = q1,
    q3 = q3,
    count = count,
    percent = percent
  )
}

for (dataset_name in names(datasets)) {
  specification <- datasets[[dataset_name]]
  raw_scores <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")
  ))
  z_scores <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_pathway_zscores.rds")
  ))
  if (!identical(colnames(raw_scores), colnames(z_scores))) {
    stop(dataset_name, ": normalized and z-score subject order differs.")
  }
  unmatched <- setdiff(colnames(raw_scores), clinical$Subject_ID)
  if (length(unmatched)) {
    stop(dataset_name, ": unmatched clinical IDs: ", paste(unmatched, collapse = ", "))
  }

  score_table <- matrix_to_subject_table(raw_scores)
  merged <- merge(score_table, clinical, by = "Subject_ID", all.x = TRUE, sort = FALSE)
  merged <- merged[match(colnames(raw_scores), Subject_ID)]
  merged[, `:=`(
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type
  )]
  data.table::setcolorder(
    merged,
    c("dataset", "dataset_label", "compartment", "platform", "omic_type",
      "Subject_ID", pathways, setdiff(names(merged), c(
        "dataset", "dataset_label", "compartment", "platform", "omic_type",
        "Subject_ID", pathways
      )))
  )
  saveRDS(
    merged,
    file.path(paths$analysis, paste0(dataset_name, "_pathway_clinical_analysis.rds")),
    compress = "gzip"
  )
  data.table::fwrite(
    merged,
    file.path(paths$analysis, paste0(dataset_name, "_pathway_clinical_analysis.csv"))
  )

  z_table <- matrix_to_subject_table(z_scores)
  data.table::fwrite(
    z_table,
    file.path(paths$analysis, paste0(dataset_name, "_pathway_zscores_for_PCA.csv"))
  )

  sample_rows[[dataset_name]] <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type,
    transcriptomic_profiles = ncol(raw_scores),
    unique_participants = data.table::uniqueN(colnames(raw_scores)),
    matched_clinical_ids = ncol(raw_scores),
    unmatched_clinical_ids = length(unmatched)
  )

  for (outcome_name in names(outcome_specs)) {
    outcome <- outcome_specs[[outcome_name]]
    x <- merged[[outcome$column]]
    missingness_rows[[length(missingness_rows) + 1L]] <- data.table::data.table(
      dataset = dataset_name,
      dataset_label = specification$label,
      compartment = specification$compartment,
      platform = specification$platform,
      omic_type = specification$omic_type,
      outcome = outcome_name,
      outcome_label = outcome$label,
      total_n = nrow(merged),
      available_n = sum(is.finite(x)),
      missing_n = sum(!is.finite(x)),
      missing_percent = 100 * mean(!is.finite(x)),
      outcome_transform = outcome$transform
    )
  }

  cohort_rows[[dataset_name]] <- merged[, .(
    participants = .N,
    percent = 100 * .N / nrow(merged)
  ), by = .(cohort_source, cohort_label)]
  cohort_rows[[dataset_name]][, dataset := dataset_name]

  baseline_formatted_rows[[dataset_name]] <- data.table::data.table(
    characteristic = c(
      "Participants",
      "Age, years",
      "Female",
      "Never smoker",
      "Former smoker",
      "Current smoker",
      "BMI, kg/m2",
      "Current/ongoing OCS exposure",
      "Baseline OCS normalised dose recorded",
      "FEV1, % predicted",
      "FeNO, ppb",
      "Blood eosinophils, x10^3/uL",
      "Sputum eosinophils, %",
      "Exacerbations in previous 12 months"
    ),
    value = c(
      as.character(nrow(merged)),
      format_mean_sd(merged$age),
      format_n_percent(merged$sex == "female"),
      format_n_percent(merged$smoking_status == "non_smoker"),
      format_n_percent(merged$smoking_status == "ex_smoker"),
      format_n_percent(merged$smoking_status == "current_smoker"),
      format_mean_sd(merged$bmi),
      format_n_percent(merged$ocs_current == "Yes"),
      format_n_percent(merged$ocs_dose_recorded == "Yes"),
      format_mean_sd(merged$fev1_percent_predicted),
      format_median_iqr(merged$feno_ppb),
      format_median_iqr(merged$blood_eosinophils_x10_3_uL, 2L),
      format_median_iqr(merged$sputum_eosinophils_percent),
      format_median_iqr(merged$exacerbations_previous_12m)
    )
  )
  baseline_formatted_rows[[dataset_name]][, dataset := dataset_name]

  continuous_specs <- list(
    age = list(label = "Age, years", x = merged$age, type = "mean_sd"),
    bmi = list(label = "BMI, kg/m2", x = merged$bmi, type = "mean_sd"),
    fev1 = list(label = "FEV1, % predicted", x = merged$fev1_percent_predicted, type = "mean_sd"),
    feno = list(label = "FeNO, ppb", x = merged$feno_ppb, type = "median_iqr"),
    blood_eos = list(label = "Blood eosinophils, x10^3/uL", x = merged$blood_eosinophils_x10_3_uL, type = "median_iqr"),
    sputum_eos = list(label = "Sputum eosinophils, %", x = merged$sputum_eosinophils_percent, type = "median_iqr"),
    exacerbations = list(label = "Exacerbations in previous 12 months", x = merged$exacerbations_previous_12m, type = "median_iqr")
  )
  for (item in continuous_specs) {
    x <- item$x[is.finite(item$x)]
    q <- if (length(x)) stats::quantile(x, c(.25, .5, .75), names = FALSE, type = 2) else rep(NA_real_, 3L)
    add_baseline_row(
      dataset_name, item$label, item$type,
      if (item$type == "mean_sd") mean(x) else q[2L],
      length(x),
      mean = if (length(x)) mean(x) else NA_real_,
      sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
      median = q[2L], q1 = q[1L], q3 = q[3L]
    )
  }
  category_specs <- list(
    female = list(label = "Female", condition = merged$sex == "female"),
    never_smoker = list(label = "Never smoker", condition = merged$smoking_status == "non_smoker"),
    former_smoker = list(label = "Former smoker", condition = merged$smoking_status == "ex_smoker"),
    current_smoker = list(label = "Current smoker", condition = merged$smoking_status == "current_smoker"),
    current_ocs = list(label = "Current/ongoing OCS exposure", condition = merged$ocs_current == "Yes"),
    ocs_dose_recorded = list(label = "Baseline OCS normalised dose recorded", condition = merged$ocs_dose_recorded == "Yes")
  )
  for (item in category_specs) {
    available <- sum(!is.na(item$condition))
    count <- sum(item$condition, na.rm = TRUE)
    add_baseline_row(
      dataset_name, item$label, "n_percent",
      if (available) 100 * count / available else NA_real_,
      available, count = count,
      percent = if (available) 100 * count / available else NA_real_
    )
  }
}

sample_table <- data.table::rbindlist(sample_rows)
missingness_table <- data.table::rbindlist(missingness_rows)
baseline_raw <- data.table::rbindlist(baseline_raw_rows, fill = TRUE)
baseline_formatted_long <- data.table::rbindlist(baseline_formatted_rows)
cohort_table <- data.table::rbindlist(cohort_rows)

characteristic_order <- unique(baseline_formatted_long$characteristic)
baseline_formatted <- data.table::dcast(
  baseline_formatted_long,
  characteristic ~ dataset,
  value.var = "value"
)
baseline_formatted[, characteristic := factor(
  characteristic, levels = characteristic_order
)]
data.table::setorder(baseline_formatted, characteristic)
baseline_formatted[, characteristic := as.character(characteristic)]

data.table::fwrite(
  sample_table, file.path(paths$tables, "Table_S_sample_availability.csv")
)
data.table::fwrite(
  missingness_table, file.path(paths$tables, "Table_S_outcome_missingness.csv")
)
data.table::fwrite(
  baseline_raw, file.path(paths$tables, "Table_1_clinical_characteristics_numeric.csv")
)
data.table::fwrite(
  baseline_formatted, file.path(paths$tables, "Table_1_clinical_characteristics.csv")
)
data.table::fwrite(
  cohort_table, file.path(paths$tables, "Table_S_cohort_composition.csv")
)

workbook <- openxlsx::createWorkbook()
openxlsx::addWorksheet(workbook, "Baseline characteristics")
openxlsx::writeData(workbook, "Baseline characteristics", baseline_formatted)
openxlsx::setColWidths(
  workbook, "Baseline characteristics", cols = seq_len(ncol(baseline_formatted)),
  widths = c(38, rep(27, ncol(baseline_formatted) - 1L))
)
openxlsx::freezePane(workbook, "Baseline characteristics", firstRow = TRUE, firstCol = TRUE)
openxlsx::addStyle(
  workbook, "Baseline characteristics",
  style = openxlsx::createStyle(
    textDecoration = "bold", fgFill = "#D9EAF2", halign = "center",
    valign = "center", wrapText = TRUE
  ), rows = 1L, cols = seq_len(ncol(baseline_formatted)), gridExpand = TRUE
)
openxlsx::addWorksheet(workbook, "Numeric source")
openxlsx::writeData(workbook, "Numeric source", baseline_raw)
openxlsx::addWorksheet(workbook, "Outcome missingness")
openxlsx::writeData(workbook, "Outcome missingness", missingness_table)
openxlsx::addWorksheet(workbook, "Cohort composition")
openxlsx::writeData(workbook, "Cohort composition", cohort_table)
openxlsx::saveWorkbook(
  workbook,
  file.path(paths$tables, "Clean_Table_1_clinical_characteristics.xlsx"),
  overwrite = TRUE
)

if (any(sample_table$unmatched_clinical_ids != 0L)) stop("Clinical ID mismatch.")
if (any(sample_table$transcriptomic_profiles != dataset_metadata_table()$expected_n)) {
  stop("Sample counts differ from expected dataset metadata.")
}

writeLines(
  c(
    "CLINICAL MERGE AND BASELINE TABLES: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "All pathway-score participant IDs matched unique clinical Subject IDs.",
    "No clinical outcome, pathway score or covariate was imputed.",
    "Exacerbations are source baseline history (events in the previous 12 months), not future events.",
    paste(
      "Current/ongoing OCS exposure is derived from the supplied baseline frequency field;",
      "never and Previous are No, ongoing frequency categories are Yes, and ambiguous text is missing."
    ),
    "A separate row reports whether a normalized OCS dose was recorded; dose presence is not called current use.",
    "Baseline characteristics describe each transcriptomic sample set and explain outcome-specific n differences."
  ),
  file.path(paths$validation, "CLINICAL_MERGE_STATUS.txt")
)
message("Clinical merge and baseline tables completed successfully.")
