# Run ssGSEA independently for all four airway transcriptomic datasets.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "GSVA", "BiocParallel"))
message_rule("04: Running ssGSEA across four airway datasets")

gene_sets <- readRDS(
  file.path(paths$analysis, "airway_39_signatures_harmonised.rds")
)
source_gene_sets <- readRDS(
  file.path(paths$analysis, "airway_39_signatures_analysis.rds")
)
if (length(gene_sets) != signature_settings$expected_count) {
  stop("Expected ", signature_settings$expected_count, " harmonised signatures.")
}

coverage_rows <- list()
score_range_rows <- list()

for (dataset_name in names(datasets)) {
  specification <- datasets[[dataset_name]]
  expression_file <- file.path(
    paths$processed, paste0(dataset_name, "_expression_gene_level.rds")
  )
  if (!file.exists(expression_file)) stop("Missing expression matrix: ", expression_file)
  expression <- readRDS(expression_file)

  coverage <- data.table::rbindlist(lapply(names(gene_sets), function(pathway) {
    detected <- intersect(gene_sets[[pathway]], rownames(expression))
    data.table::data.table(
      dataset = dataset_name,
      dataset_label = specification$label,
      compartment = specification$compartment,
      platform = specification$platform,
      omic_type = specification$omic_type,
      pathway = pathway,
      pathway_label = unname(pathway_labels[pathway]),
      source_gene_entries = length(source_gene_sets[[pathway]]),
      harmonised_unique_genes = length(gene_sets[[pathway]]),
      detected_gene_count = length(detected),
      detected_percent = 100 * length(detected) / length(gene_sets[[pathway]]),
      low_coverage_flag = length(detected) < 5L ||
        length(detected) / length(gene_sets[[pathway]]) < 0.20,
      detected_genes = paste(detected, collapse = ";")
    )
  }))
  coverage_rows[[dataset_name]] <- coverage
  if (any(coverage$detected_gene_count < ssgsea_settings$min_size)) {
    stop(
      dataset_name, ": signatures below minSize=",
      ssgsea_settings$min_size, ": ",
      paste(
        coverage[detected_gene_count < ssgsea_settings$min_size, pathway],
        collapse = ", "
      )
    )
  }

  message("ssGSEA: ", specification$label, " (", nrow(expression),
          " genes x ", ncol(expression), " participants)")
  parameter <- GSVA::ssgseaParam(
    exprData = expression,
    geneSets = gene_sets,
    minSize = ssgsea_settings$min_size,
    maxSize = Inf,
    alpha = ssgsea_settings$alpha,
    normalize = ssgsea_settings$normalize,
    checkNA = "yes",
    use = "everything",
    verbose = FALSE
  )
  scores_raw <- GSVA::gsva(
    parameter,
    verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  )
  scores_raw <- scores_raw[names(gene_sets), colnames(expression), drop = FALSE]
  scores_z <- safe_zscore_rows(scores_raw)

  if (!identical(rownames(scores_raw), names(gene_sets))) {
    stop(dataset_name, ": ssGSEA pathway order mismatch.")
  }
  if (!identical(colnames(scores_raw), colnames(expression))) {
    stop(dataset_name, ": ssGSEA participant order mismatch.")
  }
  if (anyNA(scores_raw) || any(!is.finite(scores_raw))) {
    stop(dataset_name, ": non-finite ssGSEA scores.")
  }

  saveRDS(
    scores_raw,
    file.path(paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")),
    compress = "gzip"
  )
  saveRDS(
    scores_z,
    file.path(paths$scores, paste0(dataset_name, "_ssgsea_pathway_zscores.rds")),
    compress = "gzip"
  )
  data.table::fwrite(
    matrix_to_subject_table(scores_raw),
    file.path(paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.csv"))
  )
  data.table::fwrite(
    matrix_to_subject_table(scores_z),
    file.path(paths$scores, paste0(dataset_name, "_ssgsea_pathway_zscores.csv"))
  )

  score_range_rows[[dataset_name]] <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    n_participants = ncol(scores_raw),
    n_pathways = nrow(scores_raw),
    normalized_score_minimum = min(scores_raw),
    normalized_score_maximum = max(scores_raw),
    normalized_score_mean = mean(scores_raw),
    normalized_score_sd = stats::sd(as.numeric(scores_raw)),
    zscore_minimum = min(scores_z),
    zscore_maximum = max(scores_z),
    zscore_row_mean_max_abs_error = max(abs(rowMeans(scores_z))),
    zscore_row_sd_max_abs_error = max(abs(apply(scores_z, 1L, stats::sd) - 1))
  )

  rm(expression, parameter, scores_raw, scores_z)
  invisible(gc())
}

coverage_table <- data.table::rbindlist(coverage_rows)
score_range_table <- data.table::rbindlist(score_range_rows)
data.table::fwrite(
  coverage_table,
  file.path(paths$tables, "Table_S_signature_gene_coverage_by_dataset.csv")
)
data.table::fwrite(
  score_range_table,
  file.path(paths$tables, "Table_S_ssgsea_score_scaling_audit.csv")
)

if (any(score_range_table$n_pathways != signature_settings$expected_count)) {
  stop("Not all datasets have the configured number of signature scores.")
}
if (any(score_range_table$zscore_row_mean_max_abs_error > 1e-10)) {
  stop("Pathway-wise z-score centring check failed.")
}
if (any(score_range_table$zscore_row_sd_max_abs_error > 1e-10)) {
  stop("Pathway-wise z-score scaling check failed.")
}

writeLines(
  c(
    "ssGSEA: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "Method: GSVA::ssgseaParam followed by GSVA::gsva.",
    paste0(
      "Parameters: minSize=", ssgsea_settings$min_size,
      ", maxSize=Inf, alpha=", ssgsea_settings$alpha,
      ", normalize=", toupper(as.character(ssgsea_settings$normalize)), "."
    ),
    "Normalized ssGSEA scores are retained as the primary pathway-activity values.",
    "Pathway-wise z-scores are separate derivative values for PCA and relative displays.",
    "No eicosanoid or other non-transcriptomic dataset was scored.",
    paste("Low-coverage flags:", coverage_table[low_coverage_flag == TRUE, .N])
  ),
  file.path(paths$validation, "SSGSEA_STATUS.txt")
)
message("ssGSEA completed successfully for all four datasets.")
