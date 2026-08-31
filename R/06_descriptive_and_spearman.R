# Descriptive pathway analyses and univariable clinical associations.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table"))
message_rule("06: Descriptive pathway analyses and Spearman associations")

pathways <- names(readRDS(
  file.path(paths$analysis, "airway_39_signatures_harmonised.rds")
))

variation_rows <- list()
pca_variance_rows <- list()
pca_score_rows <- list()
pca_loading_rows <- list()
correlation_rows <- list()
spearman_rows <- list()

for (dataset_name in names(datasets)) {
  specification <- datasets[[dataset_name]]
  raw_scores <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")
  ))
  z_scores <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_pathway_zscores.rds")
  ))
  merged <- readRDS(file.path(
    paths$analysis, paste0(dataset_name, "_pathway_clinical_analysis.rds")
  ))

  if (!identical(rownames(raw_scores), pathways)) stop("Pathway order mismatch.")
  if (!identical(colnames(raw_scores), merged$Subject_ID)) stop("Subject order mismatch.")

  # Participant-level variation for each pathway, retaining both score scales.
  for (score_scale in c("normalized_ssgsea", "pathway_zscore")) {
    score_matrix <- if (score_scale == "normalized_ssgsea") raw_scores else z_scores
    variation_rows[[paste(dataset_name, score_scale, sep = "__")]] <-
      data.table::rbindlist(lapply(pathways, function(pathway) {
        x <- as.numeric(score_matrix[pathway, ])
        q <- stats::quantile(x, c(.25, .50, .75), names = FALSE, type = 2)
        data.table::data.table(
          dataset = dataset_name,
          dataset_label = specification$label,
          compartment = specification$compartment,
          platform = specification$platform,
          omic_type = specification$omic_type,
          pathway = pathway,
          pathway_label = unname(pathway_labels[pathway]),
          score_scale = score_scale,
          n = length(x),
          mean = mean(x),
          sd = stats::sd(x),
          median = q[2L],
          q1 = q[1L],
          q3 = q[3L],
          iqr = q[3L] - q[1L],
          minimum = min(x),
          maximum = max(x)
        )
      }))
  }

  # PCA is descriptive and is run on pathway-wise z-standardised scores.
  pca <- stats::prcomp(t(z_scores), center = FALSE, scale. = FALSE)
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_variance_rows[[dataset_name]] <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type,
    component = paste0("PC", seq_along(variance)),
    variance_percent = variance,
    cumulative_variance_percent = cumsum(variance),
    n_participants = nrow(pca$x),
    n_pathways = ncol(pca$x)
  )
  pca_scores <- data.table::as.data.table(pca$x, keep.rownames = "Subject_ID")
  pca_scores[, `:=`(
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type
  )]
  pca_score_rows[[dataset_name]] <- pca_scores
  pca_loadings <- data.table::as.data.table(pca$rotation, keep.rownames = "pathway")
  pca_loadings[, `:=`(
    pathway_label = unname(pathway_labels[pathway]),
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type
  )]
  pca_loading_rows[[dataset_name]] <- pca_loadings

  # Spearman pathway-pathway correlation. With complete score matrices this is
  # identical on normalized and pathway-z scales; normalized scores are used.
  pathway_correlation <- stats::cor(
    t(raw_scores), method = "spearman", use = "everything"
  )
  correlation_long <- data.table::as.data.table(
    as.table(pathway_correlation), keep.rownames = FALSE
  )
  data.table::setnames(correlation_long, c("pathway_1", "pathway_2", "spearman_rho"))
  correlation_long[, `:=`(
    dataset = dataset_name,
    dataset_label = specification$label,
    compartment = specification$compartment,
    platform = specification$platform,
    omic_type = specification$omic_type,
    n = ncol(raw_scores)
  )]
  correlation_rows[[dataset_name]] <- correlation_long
  data.table::fwrite(
    data.table::as.data.table(pathway_correlation, keep.rownames = "pathway"),
    file.path(paths$tables, paste0(
      "Table_S_pathway_correlation_matrix_", dataset_name, ".csv"
    ))
  )

  # All pathway-clinical associations are univariable Spearman tests. Outcome
  # transformations are not needed because Spearman correlation is rank based.
  for (outcome_name in names(outcome_specs)) {
    outcome <- outcome_specs[[outcome_name]]
    y <- merged[[outcome$column]]
    for (pathway in pathways) {
      x <- merged[[pathway]]
      complete <- is.finite(x) & is.finite(y)
      n <- sum(complete)
      if (n < 10L) {
        rho <- p_value <- NA_real_
      } else {
        test <- suppressWarnings(stats::cor.test(
          x[complete], y[complete], method = "spearman", exact = FALSE
        ))
        rho <- unname(test$estimate)
        p_value <- test$p.value
      }
      spearman_rows[[length(spearman_rows) + 1L]] <- data.table::data.table(
        dataset = dataset_name,
        dataset_label = specification$label,
        compartment = specification$compartment,
        platform = specification$platform,
        omic_type = specification$omic_type,
        outcome = outcome_name,
        outcome_label = outcome$label,
        outcome_source_scale = outcome$original_unit,
        pathway = pathway,
        pathway_label = unname(pathway_labels[pathway]),
        n = n,
        spearman_rho = rho,
        p_value = p_value,
        test = "Spearman rank correlation; asymptotic p value",
        inference = "Univariable association; not an independent or causal effect"
      )
    }
  }
}

variation_table <- data.table::rbindlist(variation_rows)
pca_variance_table <- data.table::rbindlist(pca_variance_rows)
pca_scores_table <- data.table::rbindlist(pca_score_rows, fill = TRUE)
pca_loadings_table <- data.table::rbindlist(pca_loading_rows, fill = TRUE)
correlation_table <- data.table::rbindlist(correlation_rows)
spearman_table <- data.table::rbindlist(spearman_rows)
spearman_table[, fdr_bh := stats::p.adjust(p_value, method = "BH"),
                by = .(dataset, outcome)]
spearman_table[, significance := data.table::fcase(
  fdr_bh < 0.001, "***",
  fdr_bh < 0.01, "**",
  fdr_bh < 0.05, "*",
  default = ""
)]

data.table::fwrite(
  variation_table,
  file.path(paths$tables, "Table_S_pathway_score_variation_SD_IQR.csv")
)
data.table::fwrite(
  pca_variance_table,
  file.path(paths$tables, "Table_S_PCA_variance_explained.csv")
)
data.table::fwrite(
  pca_scores_table,
  file.path(paths$analysis, "PCA_participant_scores.csv")
)
data.table::fwrite(
  pca_loadings_table,
  file.path(paths$tables, "Table_S_PCA_pathway_loadings.csv")
)
data.table::fwrite(
  correlation_table,
  file.path(paths$tables, "Table_S_pathway_pathway_spearman_correlations_long.csv")
)
data.table::fwrite(
  spearman_table,
  file.path(paths$tables, "Table_2_pathway_clinical_spearman_associations.csv")
)

if (nrow(spearman_table) !=
    length(datasets) * length(outcome_specs) * signature_settings$expected_count) {
  stop("Expected 780 pathway-clinical associations; found ", nrow(spearman_table), ".")
}
if (any(!is.finite(spearman_table$spearman_rho))) {
  stop("Non-finite Spearman estimates found.")
}

writeLines(
  c(
    "DESCRIPTIVE AND SPEARMAN ANALYSES: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Pathway-clinical tests:", nrow(spearman_table)),
    paste("BH-supported tests at q<0.05:", sum(spearman_table$fdr_bh < 0.05)),
    "BH correction was applied separately within each dataset-outcome family of 39 pathways.",
    "PCA used pathway-wise z-standardised scores and is descriptive only.",
    "The PCA and heatmap displays do not test for, prove, or exclude discrete latent classes/endotypes.",
    "Spearman estimates are univariable associations and do not establish causality."
  ),
  file.path(paths$validation, "DESCRIPTIVE_SPEARMAN_STATUS.txt")
)
message("Descriptive and Spearman analyses completed successfully.")
