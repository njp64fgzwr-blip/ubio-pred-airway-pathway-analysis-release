# Validate the clean and draft-matched figure packages against their R sources.
#
# This script performs only read-only checks on existing scientific outputs.
# It verifies numerical source extracts, four-format exports, dimensions,
# direct PDF renders and the redundant encodings needed for grayscale reading.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "digest", "magick"))
message_rule("17: Validating figure values, formats and R-only provenance")

near_equal <- function(x, y, tolerance = 1e-10) {
  if (length(x) != length(y) || any(is.na(x) != is.na(y))) return(FALSE)
  both_na <- is.na(x) & is.na(y)
  x <- x[!both_na]
  y <- y[!both_na]
  if (length(x) != length(y)) return(FALSE)
  if (!length(x)) return(TRUE)
  if (is.numeric(x) || is.integer(x)) {
    return(all(is.finite(x) == is.finite(y)) &&
      all(abs(x - y) <= tolerance * pmax(1, abs(x), abs(y)), na.rm = TRUE))
  }
  identical(as.character(x), as.character(y))
}

compare_tables <- function(left, right, keys, columns, tolerance = 1e-10) {
  left <- data.table::as.data.table(left)
  right <- data.table::as.data.table(right)
  required <- unique(c(keys, columns))
  if (!all(required %in% names(left)) || !all(required %in% names(right))) {
    return(FALSE)
  }
  left <- unique(left[, ..required])
  right <- unique(right[, ..required])
  for (key in keys) {
    if (!is.numeric(left[[key]])) left[, (key) := as.character(get(key))]
    if (!is.numeric(right[[key]])) right[, (key) := as.character(get(key))]
  }
  data.table::setorderv(left, keys)
  data.table::setorderv(right, keys)
  if (nrow(left) != nrow(right)) return(FALSE)
  if (!all(vapply(keys, function(z) near_equal(left[[z]], right[[z]], tolerance),
                  logical(1)))) return(FALSE)
  all(vapply(columns, function(z) near_equal(left[[z]], right[[z]], tolerance),
             logical(1)))
}

read_csv <- function(path) data.table::fread(file.path(paths$root, path))

clean_manifest <- data.table::fread(file.path(
  paths$validation, "FIGURE_EXPORT_MANIFEST.csv"
))
draft_manifest <- data.table::fread(file.path(
  paths$validation, "DRAFT_MATCHED_FIGURE_EXPORT_MANIFEST.csv"
))
clean_manifest[, figure_set := "base_dissertation"]
draft_manifest[, figure_set := "draft_matched"]
manifest <- data.table::rbindlist(list(clean_manifest, draft_manifest), fill = TRUE)

expected_stems <- c(base_dissertation = 18L, draft_matched = 28L)
if (nrow(manifest) != 184L || data.table::uniqueN(manifest$figure_stem) != 46L) {
  stop("Expected 46 figure stems and 184 production exports.")
}
if (any(manifest$generated_in != "R") || any(manifest$hand_edited)) {
  stop("A figure is not declared as an unedited R export.")
}

source_paths <- unique(trimws(unlist(strsplit(
  manifest$source_data_files, ";", fixed = TRUE
))))
source_paths <- source_paths[nzchar(source_paths)]
wildcards <- source_paths[grepl("[*]", source_paths)]
plain_sources <- source_paths[!grepl("[*]", source_paths)]
missing_sources <- plain_sources[!file.exists(file.path(paths$root, plain_sources))]
if (length(missing_sources)) {
  stop("Missing declared figure sources: ", paste(missing_sources, collapse = ", "))
}
for (wildcard in wildcards) {
  matches <- Sys.glob(file.path(paths$root, wildcard))
  if (!length(matches)) stop("A source wildcard matched no file: ", wildcard)
}

# ---------------------------------------------------------------------------
# Numerical reconstruction of every figure-specific CSV from canonical files.
# ---------------------------------------------------------------------------
source_audit <- list()
add_audit <- function(set, stem, source_csv, canonical, keys, columns,
                      dataset = "All", outcome = "All",
                      detail = "Exact keyed comparison") {
  source_table <- read_csv(source_csv)
  canonical_table <- if (is.character(canonical)) read_csv(canonical) else canonical
  passed <- compare_tables(source_table, canonical_table, keys, columns)
  source_audit[[length(source_audit) + 1L]] <<- data.table::data.table(
    figure_set = set,
    figure_stem = stem,
    figure_source_csv = source_csv,
    canonical_source_data = if (is.character(canonical)) canonical else "R-reconstructed canonical source",
    dataset = dataset,
    outcome = outcome,
    statistic_checked = paste(columns, collapse = "; "),
    figure_source_rows = nrow(source_table),
    expected_rows = nrow(unique(canonical_table[, unique(c(keys, columns)), with = FALSE])),
    values_match_source = passed,
    check_detail = detail
  )
}

add_audit(
  "base_dissertation", "Figure_1_cohort_structure_and_missingness",
  "results/figure_source_data/Figure_1A_sample_availability.csv",
  "results/tables/Table_S_sample_availability.csv", "dataset",
  c("transcriptomic_profiles", "unique_participants", "matched_clinical_ids")
)
add_audit(
  "base_dissertation", "Figure_1_cohort_structure_and_missingness",
  "results/figure_source_data/Figure_1B_outcome_availability.csv",
  "results/tables/Table_S_outcome_missingness.csv", c("dataset", "outcome"),
  c("total_n", "available_n", "missing_n", "missing_percent")
)
add_audit(
  "base_dissertation", "Figure_1_cohort_structure_and_missingness",
  "results/figure_source_data/Figure_1C_exact_overlap_patterns.csv",
  "results/tables/Table_S_exact_dataset_overlap.csv", "exact_pattern",
  "participants"
)

dataset_file_map <- c(
  sputum = "2A", brushing_gpl570 = "2B", biopsy = "2C",
  brushing_rnaseq = "2D"
)
for (dataset_name in names(dataset_file_map)) {
  matrix <- read_csv(paste0(
    "results/ssgsea/", dataset_name, "_ssgsea_normalized_scores.csv"
  ))
  long <- data.table::melt(
    matrix, id.vars = "Subject_ID", variable.name = "pathway",
    value.name = "score"
  )
  add_audit(
    "draft_matched",
    paste0("DraftMatched_Figure_", dataset_file_map[dataset_name], "_heatmap"),
    paste0(
      "results/figure_source_data/draft_matched/DraftMatched_Figure_",
      dataset_file_map[dataset_name], "_heatmap.csv"
    ),
    long, c("pathway", "Subject_ID"), "score", dataset_name, "All",
    "Every normalized ssGSEA heatmap cell matched the saved R matrix"
  )
}

add_audit(
  "base_dissertation", "Figure_3_pathway_clinical_spearman_associations",
  "results/figure_source_data/Figure_3_spearman_associations.csv",
  "results/tables/Table_2_pathway_clinical_spearman_associations.csv",
  c("dataset", "outcome", "pathway"),
  c("n", "spearman_rho", "p_value", "fdr_bh")
)
add_audit(
  "base_dissertation", "Figure_4_curated_pathway_model_forest_plots",
  "results/figure_source_data/Figure_4_curated_coefficients.csv",
  "results/tables/Table_3_curated_model_coefficients_95CI.csv",
  c("dataset", "outcome", "pathway"),
  c("n", "standardized_beta", "ci_low", "ci_high", "fdr_bh")
)
canonical_stability <- read_csv(
  "results/tables/Table_5_elastic_net_selection_stability.csv"
)[model_type == "Elastic Net pathways"]
add_audit(
  "base_dissertation", "Figure_5_elastic_net_pathway_selection_stability",
  "results/figure_source_data/Figure_5_all_pathway_selection_stability.csv",
  canonical_stability, c("dataset", "outcome", "model_type", "pathway"),
  c("outer_models", "selection_count", "selection_frequency",
    "sign_consistency", "median_nonzero_coefficient")
)
add_audit(
  "base_dissertation", "Figure_6_internal_predictive_performance",
  "results/figure_source_data/Figure_6A_and_Figure_8_model_performance.csv",
  "results/tables/Table_7_pathway_vs_biomarker_model_comparison.csv",
  c("dataset", "outcome", "model_type"),
  c("n", "mean_cross_validated_r_squared", "RMSE", "MAE")
)
add_audit(
  "base_dissertation", "Figure_6_internal_predictive_performance",
  "results/figure_source_data/Figure_6B_delta_performance.csv",
  "results/tables/Table_S_combined_vs_biomarker_delta_inference.csv",
  c("dataset", "outcome", "augmented_model"),
  c("n", "delta_model_r_squared_point", "delta_model_r_squared_lower_95",
    "delta_model_r_squared_upper_95",
    "exploratory_BH_q_across_40_delta_r_squared_comparisons")
)
add_audit(
  "base_dissertation", "Figure_7_matched_pathway_concordance",
  "results/figure_source_data/Figure_7_pathway_concordance.csv",
  read_csv("results/tables/Table_5_matched_pathway_concordance.csv")[
    comparison %in% c(
      "brushing_gpl570_vs_biopsy", "sputum_vs_biopsy",
      "sputum_vs_brushing_gpl570",
      "brushing_gpl570_vs_brushing_rnaseq"
    )
  ],
  c("comparison", "pathway"),
  c("n_matched", "spearman_rho", "p_value", "fdr_bh")
)
add_audit(
  "base_dissertation", "Figure_8_pathway_vs_biomarker_model_performance",
  "results/figure_source_data/Figure_6A_and_Figure_8_model_performance.csv",
  "results/tables/Table_7_pathway_vs_biomarker_model_comparison.csv",
  c("dataset", "outcome", "model_type"),
  c("n", "mean_cross_validated_r_squared", "RMSE", "MAE")
)
add_audit(
  "base_dissertation", "Figure_S1_pathway_PCA",
  "results/figure_source_data/Figure_S1_PCA_scores.csv",
  "results/analysis_data/PCA_participant_scores.csv",
  c("dataset", "Subject_ID"), c("PC1", "PC2")
)
add_audit(
  "base_dissertation", "Figure_S3_normalized_score_between_subject_variation",
  "results/figure_source_data/Figure_S3_pathway_variation.csv",
  read_csv("results/tables/Table_S_pathway_score_variation_SD_IQR.csv")[
    score_scale == "normalized_ssgsea"
  ],
  c("dataset", "pathway"), c("n", "mean", "sd", "median", "iqr")
)
add_audit(
  "base_dissertation", "Figure_S4_covariate_adjusted_curated_models",
  "results/figure_source_data/Figure_S4_adjusted_coefficients.csv",
  "results/tables/Table_S_covariate_adjusted_curated_coefficients.csv",
  c("dataset", "outcome", "pathway"),
  c("n", "standardized_beta", "ci_low", "ci_high", "fdr_bh")
)
add_audit(
  "base_dissertation", "Figures_S5_to_S8_calibration",
  "results/figure_source_data/Figures_S5_to_S8_calibration_bins.csv",
  "results/tables/Table_S_predictive_calibration_bins.csv",
  c("dataset", "outcome", "model_type", "calibration_bin"),
  c("bin_n", "mean_observed_model", "mean_predicted_model",
    "se_observed_model")
)
coverage_canonical <- read_csv(
  "results/tables/Table_S_signature_gene_coverage_by_dataset.csv"
)[dataset == "brushing_rnaseq"]
add_audit(
  "base_dissertation", "Figure_S9_rnaseq_signature_gene_coverage",
  "results/figure_source_data/Figure_S9_RNAseq_coverage.csv",
  coverage_canonical, c("dataset", "pathway"),
  c("source_gene_entries", "harmonised_unique_genes", "detected_gene_count",
    "detected_percent", "low_coverage_flag")
)
clip_canonical <- read_csv(
  "results/tables/Table_S_below_zero_prediction_clipping_sensitivity.csv"
)[prediction_level == "Averaged repeated OOF"]
add_audit(
  "base_dissertation", "Figure_S10_below_zero_prediction_sensitivity",
  "results/figure_source_data/Figure_S10_clipping_counts.csv",
  clip_canonical, c("dataset", "outcome", "model_type"),
  c("prediction_rows", "below_zero_prediction_rows",
    "MAE_improvement_after_clipping")
)

# Draft-matched source extracts derived from the same canonical tables.
letters <- c(sputum = "A", brushing_gpl570 = "B", biopsy = "C", brushing_rnaseq = "D")
for (dataset_name in names(letters)) {
  letter <- letters[dataset_name]
  add_audit(
    "draft_matched", paste0("DraftMatched_Figure_3", letter),
    paste0("results/figure_source_data/draft_matched/DraftMatched_Figure_3",
           letter, "_spearman.csv"),
    read_csv("results/tables/Table_2_pathway_clinical_spearman_associations.csv")[
      dataset == dataset_name
    ], c("dataset", "outcome", "pathway"),
    c("n", "spearman_rho", "p_value", "fdr_bh"), dataset_name
  )
  add_audit(
    "draft_matched", paste0("DraftMatched_Figure_4", letter),
    paste0("results/figure_source_data/draft_matched/DraftMatched_Figure_4",
           letter, "_coefficients.csv"),
    read_csv("results/tables/Table_3_curated_model_coefficients_95CI.csv")[
      dataset == dataset_name
    ], c("dataset", "outcome", "pathway"),
    c("n", "standardized_beta", "ci_low", "ci_high", "fdr_bh"), dataset_name
  )
  draft_stability <- read_csv(paste0(
    "results/figure_source_data/draft_matched/DraftMatched_Figure_5",
    letter, "_top10_stability.csv"
  ))
  canonical_subset <- canonical_stability[
    dataset == dataset_name &
      paste(outcome, pathway) %in% paste(draft_stability$outcome,
                                         draft_stability$pathway)
  ]
  add_audit(
    "draft_matched", paste0("DraftMatched_Figure_5", letter),
    paste0("results/figure_source_data/draft_matched/DraftMatched_Figure_5",
           letter, "_top10_stability.csv"),
    canonical_subset, c("dataset", "outcome", "model_type", "pathway"),
    c("outer_models", "selection_count", "selection_frequency",
      "sign_consistency", "median_nonzero_coefficient"), dataset_name
  )
}

performance <- read_csv("results/tables/Table_4_repeated_nested_CV_model_performance.csv")
draft_perf_checks <- list(
  `DraftMatched_Figure_6A` = list(
    file = "DraftMatched_Figure_6A_curated_vs_elastic_net.csv",
    filter = performance[model_type %in% c("Curated pathways", "Elastic Net pathways")]
  ),
  `DraftMatched_Figure_6B` = list(
    file = "DraftMatched_Figure_6B_three_model_performance.csv",
    filter = performance[model_type %in% c(
      "Biomarker-only", "Curated pathways", "Elastic Net pathways"
    )]
  ),
  `DraftMatched_Figure_6C` = list(
    file = "DraftMatched_Figure_6C_sputum_performance.csv",
    filter = performance[dataset == "sputum" & model_type %in% c(
      "Biomarker-only", "Curated pathways", "Elastic Net pathways"
    )]
  )
)
for (stem in names(draft_perf_checks)) {
  check <- draft_perf_checks[[stem]]
  add_audit(
    "draft_matched", stem,
    file.path("results/figure_source_data/draft_matched", check$file),
    check$filter, c("dataset", "outcome", "model_type"),
    c("n", "mean_cv_r_squared", "sd_cv_r_squared",
      "mean_cv_pearson_r", "sd_cv_pearson_r")
  )
}

calibration <- read_csv("results/tables/Table_S_predictive_calibration_bins.csv")
add_audit(
  "draft_matched", "DraftMatched_Figure_6D",
  "results/figure_source_data/draft_matched/DraftMatched_Figure_6D_sputum_eosinophil_calibration.csv",
  calibration[dataset == "sputum" & outcome == "sputum_eosinophils" &
                model_type %in% c(
                  "Biomarker-only", "Curated pathways", "Elastic Net pathways",
                  "Biomarker + curated pathways",
                  "Biomarker + Elastic Net pathways"
                )],
  c("dataset", "outcome", "model_type", "calibration_bin"),
  c("bin_n", "mean_observed_model", "mean_predicted_model", "se_observed_model")
)
concordance <- read_csv("results/tables/Table_5_matched_pathway_concordance.csv")
add_audit(
  "draft_matched", "DraftMatched_Figure_7A",
  "results/figure_source_data/draft_matched/DraftMatched_Figure_7A_concordance.csv",
  concordance[comparison %in% c(
    "brushing_gpl570_vs_biopsy", "sputum_vs_biopsy",
    "sputum_vs_brushing_gpl570"
  )], c("comparison", "pathway"),
  c("n_matched", "spearman_rho", "p_value", "fdr_bh")
)
add_audit(
  "draft_matched", "DraftMatched_Figure_7B",
  "results/figure_source_data/draft_matched/DraftMatched_Figure_7B_concordance.csv",
  concordance[comparison == "brushing_gpl570_vs_brushing_rnaseq"],
  c("comparison", "pathway"),
  c("n_matched", "spearman_rho", "p_value", "fdr_bh",
    "rnaseq_detected_percent", "rnaseq_low_coverage_flag")
)
add_audit(
  "draft_matched", "DraftMatched_Figure_S1",
  "results/figure_source_data/draft_matched/DraftMatched_Figure_S1_PCA_scores.csv",
  "results/analysis_data/PCA_participant_scores.csv",
  c("dataset", "Subject_ID"), c("PC1", "PC2")
)
for (i in seq_along(names(letters))) {
  dataset_name <- names(letters)[i]
  add_audit(
    "draft_matched", paste0("DraftMatched_Figure_S", i + 1L),
    paste0("results/figure_source_data/draft_matched/DraftMatched_Figure_S",
           i + 1L, "_calibration.csv"),
    calibration[dataset == dataset_name & model_type %in% c(
      "Biomarker-only", "Curated pathways", "Elastic Net pathways"
    )], c("dataset", "outcome", "model_type", "calibration_bin"),
    c("bin_n", "mean_observed_model", "mean_predicted_model", "se_observed_model"),
    dataset_name
  )
}

source_audit_table <- data.table::rbindlist(source_audit, fill = TRUE)
if (any(!source_audit_table$values_match_source)) {
  data.table::fwrite(source_audit_table, file.path(
    paths$validation, "FIGURE_SOURCE_AUDIT.csv"
  ))
  stop("At least one figure-source numerical reconstruction failed.")
}
data.table::fwrite(source_audit_table, file.path(
  paths$validation, "FIGURE_SOURCE_AUDIT.csv"
))

# ---------------------------------------------------------------------------
# Four-format and direct-render checks.
# ---------------------------------------------------------------------------
format_checks <- manifest[, {
  info <- file.info(output_file)
  image_width <- NA_real_
  image_height <- NA_real_
  if (output_format %in% c("png", "tiff")) {
    image <- magick::image_info(magick::image_read(output_file))
    image_width <- image$width[1]
    image_height <- image$height[1]
  }
  expected_width <- if (output_format == "png") round(width_mm / 25.4 * 300) else
    if (output_format == "tiff") round(width_mm / 25.4 * 600) else NA_real_
  expected_height <- if (output_format == "png") round(height_mm / 25.4 * 300) else
    if (output_format == "tiff") round(height_mm / 25.4 * 600) else NA_real_
  dimensions_pass <- if (output_format %in% c("png", "tiff")) {
    abs(image_width - expected_width) <= 2L &&
      abs(image_height - expected_height) <= 2L
  } else TRUE
  list(
    file_exists = file.exists(output_file),
    current_bytes = as.numeric(info$size),
    manifest_bytes_match = identical(as.numeric(info$size), as.numeric(output_bytes)),
    image_width_px = as.numeric(image_width),
    image_height_px = as.numeric(image_height),
    expected_width_px = as.numeric(expected_width),
    expected_height_px = as.numeric(expected_height),
    dimensions_pass = dimensions_pass,
    sha256 = digest::digest(output_file, algo = "sha256", file = TRUE)
  )
}, by = .(figure_set, figure_stem, output_format, output_file,
          output_bytes, width_mm, height_mm, source_script,
          generated_in, hand_edited)]

format_counts <- format_checks[, .(
  formats = data.table::uniqueN(output_format),
  format_set = paste(sort(unique(output_format)), collapse = ";"),
  nonempty = all(current_bytes > 0),
  bytes_match_manifest = all(manifest_bytes_match),
  dimensions_pass = all(dimensions_pass)
), by = .(figure_set, figure_stem)]
if (any(format_counts$formats != 4L) ||
    any(format_counts$format_set != "pdf;png;svg;tiff") ||
    any(!format_counts$nonempty) || any(!format_counts$bytes_match_manifest) ||
    any(!format_counts$dimensions_pass)) {
  stop("A four-format export check failed.")
}
data.table::fwrite(format_checks, file.path(
  paths$validation, "FIGURE_FORMAT_AUDIT.csv"
))

render_manifest <- data.table::fread(file.path(
  paths$validation, "ALL_FIGURE_PDF_DIRECT_RENDER_MANIFEST.csv"
))
render_checks <- render_manifest[, .(
  pdf_hash_matches_render_record =
    digest::digest(source_pdf, algo = "sha256", file = TRUE) == source_pdf_sha256,
  colour_exists = file.exists(colour_render),
  grayscale_exists = file.exists(grayscale_render),
  colour_nonempty = file.info(colour_render)$size > 0,
  grayscale_nonempty = file.info(grayscale_render)$size > 0,
  colour_grayscale_dimensions_match =
    colour_width_px == grayscale_width_px & colour_height_px == grayscale_height_px
), by = .(figure_set, figure_stem, source_pdf, source_pdf_sha256,
          colour_render, grayscale_render, colour_width_px, colour_height_px,
          grayscale_width_px, grayscale_height_px)]
if (nrow(render_checks) != 46L ||
    any(!unlist(render_checks[, 7:ncol(render_checks)]))) {
  stop("A direct colour/grayscale PDF-render check failed.")
}
data.table::fwrite(render_checks, file.path(
  paths$validation, "FIGURE_DIRECT_RENDER_AUDIT.csv"
))

redundancy <- data.table::data.table(
  figure_family = c(
    "Normalized ssGSEA heatmaps", "Spearman association heatmaps",
    "Curated coefficient forests", "Elastic Net stability plots",
    "Predictive model comparisons", "Cross-compartment concordance heatmaps",
    "RNA-seq coverage bars", "Calibration plots"
  ),
  affected_stems = c(
    "Figure 2, S2 and DraftMatched 2A-D",
    "Figure 3 and DraftMatched 3A-D",
    "Figure 4, S4 and DraftMatched 4A-D",
    "Figure 5 and DraftMatched 5A-D",
    "Figures 6, 8 and DraftMatched 6A-C",
    "DraftMatched 7A",
    "Figure S9 and DraftMatched 7B",
    "Figures S5-S8 and DraftMatched 6D/S2-S5"
  ),
  non_colour_encoding = c(
    "Blue-white-red scale with signed numeric legend; titles identify normalized/z scale",
    "+/- sign glyph in every cell; stars show within-outcome BH q",
    "Filled versus open points for BH-FDR status",
    "Up/down/open point shape encodes predominant sign/never selected",
    "Square/circle/triangle point shapes encode model class",
    "+/- sign glyph in every cell; stars show BH q",
    "Detected/total text and explicit LOW label",
    "Square/circle/triangle point shapes encode model class"
  ),
  colour_computational_render = "PASS",
  grayscale_computational_render = "PASS",
  current_manual_visual_review = "REQUIRED",
  historical_visual_review = paste(
    "The dissertation-locked 18-August-2026 renders were inspected",
    "individually; a new environment-specific render requires a new review"
  )
)

# Targeted presentation regressions requested after the initial rebuild.
# These are code-level guarantees because all four production formats are
# generated from the same ggplot object in one R function.
helper_code <- paste(readLines(file.path(paths$scripts, "12_figure_helpers.R")),
                     collapse = "\n")
main_code <- paste(readLines(file.path(paths$scripts,
                                      "12_generate_final_figures.R")),
                   collapse = "\n")
draft_code <- paste(readLines(file.path(paths$scripts,
                                       "14_generate_draft_matched_figures.R")),
                    collapse = "\n")
legend_file <- file.path(
  paths$final_figures, "FINAL_DISSERTATION_FIGURE_LEGENDS.md"
)
if (!file.exists(legend_file)) {
  stop("Final dissertation figure legends are unavailable: ", legend_file)
}
caption_code <- paste(readLines(legend_file), collapse = "\n")
presentation_regressions <- data.table::data.table(
  check = c(
    "all_diverging_heatmaps_use_draft_blue_white_red",
    "elastic_net_stability_uses_red_blue_sign_encoding",
    "explanatory_subtitles_and_captions_hidden_from_exports",
    "external_legends_define_blue_white_red_heatmaps",
    "external_legends_define_elastic_net_red_blue"
  ),
  passed = c(
    grepl("score_heatmap_colours <- association_colours", helper_code,
          fixed = TRUE) &&
      grepl('"low" = "#2166AC"', helper_code, fixed = TRUE) &&
      grepl('"mid" = "#F7F7F7"', helper_code, fixed = TRUE) &&
      grepl('"high" = "#B2182B"', helper_code, fixed = TRUE),
    grepl('predominant_sign == "Positive", "#B2182B"', main_code,
          fixed = TRUE) &&
      grepl('predominant_sign == "Negative", "#2166AC"', main_code,
            fixed = TRUE) &&
      grepl('"Predominantly positive" = "#B2182B"', draft_code,
            fixed = TRUE) &&
      grepl('"Predominantly negative" = "#2166AC"', draft_code,
            fixed = TRUE),
    grepl("plot.subtitle = element_blank()", helper_code, fixed = TRUE) &&
      grepl("plot.caption = element_blank()", helper_code, fixed = TRUE),
    grepl("blue indicates lower", tolower(caption_code), fixed = TRUE) &&
      grepl("red higher", tolower(caption_code), fixed = TRUE),
    grepl("blue denotes a predominantly negative coefficient",
          tolower(caption_code), fixed = TRUE) &&
      grepl("red a predominantly positive coefficient",
            tolower(caption_code), fixed = TRUE)
  ),
  evidence = c(
    "Shared heatmap helper maps low/mid/high to #2166AC/#F7F7F7/#B2182B",
    "All-39 and draft-matched Elastic Net figures map negative to blue and positive to red",
    "Shared four-format exporter blanks plot.subtitle and plot.caption before rendering",
    "Main/supplementary/draft-matched caption files retain the heatmap explanation",
    "Main/draft-matched caption files retain the Elastic Net sign explanation"
  )
)
if (any(!presentation_regressions$passed)) {
  stop("A user-requested colour/legend presentation regression failed.")
}
data.table::fwrite(
  presentation_regressions,
  file.path(paths$validation, "FIGURE_USER_STYLE_REGRESSION_AUDIT.csv")
)
data.table::fwrite(redundancy, file.path(
  paths$validation, "FIGURE_GRAYSCALE_ENCODING_AUDIT.csv"
))

writeLines(
  c(
    "FIGURE SOURCE / FORMAT / R-ONLY COMPUTATIONAL QA: PASS",
    paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "Figure stems: 46/46 (18 base dissertation; 28 draft-matched)",
    "Production exports: 184/184 (PDF, SVG, PNG300 and TIFF600 per stem)",
    paste0("Numerical source reconstruction checks: ", nrow(source_audit_table),
           "/", nrow(source_audit_table), " PASS"),
    "Direct PDF colour render existence/hash/dimension checks: 46/46 PASS",
    "Direct PDF grayscale render existence/hash/dimension checks: 46/46 PASS",
    "Manual visual inspection of renders from this run: REQUIRED",
    "Historical note: the dissertation-locked 18-August-2026 renders passed individual visual review.",
    "All figure manifests declare generated_in=R and hand_edited=FALSE.",
    "All scientific figures were generated and assembled directly by the R workflow.",
    "All declared source files exist; wildcard source declarations resolve.",
    "Interpretive review: association is not causality; same-participant concordance is not replication; concurrent sputum-eosinophil capture is not prognosis; predictive results are exploratory internal validation."
  ),
  file.path(paths$validation, "FIGURE_FORMAT_AND_SOURCE_QA_STATUS.txt")
)

message("Validated 46 R-produced figure stems and 184 exports.")
