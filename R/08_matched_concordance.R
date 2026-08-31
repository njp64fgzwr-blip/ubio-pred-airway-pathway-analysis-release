# Matched-participant pathway-score concordance across compartments/platforms.

source(file.path("R", "00_common.R"))
assert_packages("data.table")
message_rule("08: Matched-participant compartment and platform concordance")

pathways <- names(readRDS(
  file.path(paths$analysis, "airway_39_signatures_harmonised.rds")
))

comparisons <- list(
  brushing_gpl570_vs_biopsy = list(
    left = "brushing_gpl570", right = "biopsy",
    comparison_type = "Same-platform cross-compartment",
    main_figure_panel = "A"
  ),
  sputum_vs_biopsy = list(
    left = "sputum", right = "biopsy",
    comparison_type = "Same-platform cross-compartment",
    main_figure_panel = "A"
  ),
  sputum_vs_brushing_gpl570 = list(
    left = "sputum", right = "brushing_gpl570",
    comparison_type = "Same-platform cross-compartment",
    main_figure_panel = "A"
  ),
  brushing_gpl570_vs_brushing_rnaseq = list(
    left = "brushing_gpl570", right = "brushing_rnaseq",
    comparison_type = "Same-compartment cross-platform",
    main_figure_panel = "B"
  ),
  biopsy_vs_brushing_rnaseq = list(
    left = "biopsy", right = "brushing_rnaseq",
    comparison_type = "Cross-compartment and cross-platform supplementary",
    main_figure_panel = NA_character_
  ),
  sputum_vs_brushing_rnaseq = list(
    left = "sputum", right = "brushing_rnaseq",
    comparison_type = "Cross-compartment and cross-platform supplementary",
    main_figure_panel = NA_character_
  )
)

scores <- lapply(names(datasets), function(dataset_name) readRDS(file.path(
  paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")
)))
names(scores) <- names(datasets)

concordance_rows <- list()
participant_rows <- list()
for (comparison_name in names(comparisons)) {
  comparison <- comparisons[[comparison_name]]
  left_name <- comparison$left
  right_name <- comparison$right
  matched <- sort(intersect(
    colnames(scores[[left_name]]), colnames(scores[[right_name]])
  ))
  if (length(matched) < 10L) stop(comparison_name, ": too few matched participants.")
  participant_rows[[comparison_name]] <- data.table::data.table(
    comparison = comparison_name,
    dataset_1 = left_name,
    dataset_2 = right_name,
    Subject_ID = matched
  )

  pathway_rows <- lapply(pathways, function(pathway) {
    x <- scores[[left_name]][pathway, matched]
    y <- scores[[right_name]][pathway, matched]
    test <- suppressWarnings(stats::cor.test(
      x, y, method = "spearman", exact = FALSE
    ))
    data.table::data.table(
      comparison = comparison_name,
      comparison_type = comparison$comparison_type,
      main_figure_panel = comparison$main_figure_panel,
      dataset_1 = left_name,
      dataset_1_label = datasets[[left_name]]$label,
      compartment_1 = datasets[[left_name]]$compartment,
      platform_1 = datasets[[left_name]]$platform,
      dataset_2 = right_name,
      dataset_2_label = datasets[[right_name]]$label,
      compartment_2 = datasets[[right_name]]$compartment,
      platform_2 = datasets[[right_name]]$platform,
      omic_type = "Transcriptomics",
      pathway = pathway,
      pathway_label = unname(pathway_labels[pathway]),
      n_matched = length(matched),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      interpretation = paste(
        "Same-participant score concordance; not external replication",
        "and not evidence of causal equivalence"
      )
    )
  })
  table <- data.table::rbindlist(pathway_rows)
  table[, fdr_bh := stats::p.adjust(p_value, method = "BH")]
  concordance_rows[[comparison_name]] <- table
}

concordance <- data.table::rbindlist(concordance_rows)
matched_participants <- data.table::rbindlist(participant_rows)
coverage <- data.table::fread(
  file.path(paths$tables, "Table_S_signature_gene_coverage_by_dataset.csv")
)
rnaseq_coverage <- coverage[dataset == "brushing_rnaseq", .(
  pathway, rnaseq_detected_gene_count = detected_gene_count,
  rnaseq_detected_percent = detected_percent,
  rnaseq_low_coverage_flag = low_coverage_flag
)]
concordance <- merge(concordance, rnaseq_coverage, by = "pathway", all.x = TRUE)

summary_table <- concordance[, .(
  n_matched = unique(n_matched),
  pathway_count = .N,
  median_spearman_rho = stats::median(spearman_rho),
  minimum_spearman_rho = min(spearman_rho),
  maximum_spearman_rho = max(spearman_rho),
  fdr_supported_pathways = sum(fdr_bh < 0.05)
), by = .(
  comparison, comparison_type, main_figure_panel,
  dataset_1, dataset_1_label, compartment_1, platform_1,
  dataset_2, dataset_2_label, compartment_2, platform_2
)]

data.table::fwrite(
  concordance,
  file.path(paths$tables, "Table_5_matched_pathway_concordance.csv")
)
data.table::fwrite(
  summary_table,
  file.path(paths$tables, "Table_S_matched_concordance_summary.csv")
)
data.table::fwrite(
  matched_participants,
  file.path(paths$analysis, "matched_concordance_participant_ids.csv")
)

expected_n <- c(
  brushing_gpl570_vs_biopsy = 99L,
  sputum_vs_biopsy = 23L,
  sputum_vs_brushing_gpl570 = 31L,
  brushing_gpl570_vs_brushing_rnaseq = 118L,
  biopsy_vs_brushing_rnaseq = 80L,
  sputum_vs_brushing_rnaseq = 25L
)
observed_n <- stats::setNames(summary_table$n_matched, summary_table$comparison)
if (!identical(as.integer(observed_n[names(expected_n)]), as.integer(expected_n))) {
  stop("Matched-participant counts differ from the source membership audit.")
}

writeLines(
  c(
    "MATCHED CONCORDANCE: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Comparison-specific matched n:", paste(
      paste0(names(expected_n), "=", expected_n), collapse = "; "
    )),
    "Main cross-compartment comparisons hold platform constant at GPL570.",
    "The brushing GPL570/RNA-seq analysis holds compartment and participants constant while changing platform.",
    "These are same-participant concordance analyses, not independent replication or external validation."
  ),
  file.path(paths$validation, "MATCHED_CONCORDANCE_STATUS.txt")
)
message("Matched concordance analyses completed successfully.")
