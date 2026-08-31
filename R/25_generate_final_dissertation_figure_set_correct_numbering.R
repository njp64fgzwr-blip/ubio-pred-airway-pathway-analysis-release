# Build the final, ordered dissertation figure set entirely in R.
#
# This script does not edit the dissertation and does not hand-edit images.
# It reuses validated R plot objects, assigns the final dissertation figure
# numbers/titles in R, and renders fresh PDF, SVG, PNG300 and TIFF600 exports.

source(file.path("R", "00_common.R"), local = TRUE)
assert_packages(c(
  "data.table", "ggplot2", "patchwork", "svglite", "ragg", "scales"
))
message_rule("25: Generating the correctly renumbered dissertation figure set in R")

# Recreate the validated clean and draft-shaped R plot objects. These scripts
# read only the saved analysis matrices/tables and export directly from R.
source(file.path("R", "12_generate_final_figures.R"), local = TRUE)
source(file.path("R", "14_generate_draft_matched_figures.R"), local = TRUE)

final_root <- paths$final_figures
final_main <- file.path(final_root, "01_MAIN_FIGURES")
final_supp <- file.path(final_root, "02_SUPPLEMENTARY_FIGURES")
final_qa <- file.path(final_root, "03_QA_RENDERS")

# This directory contains generated deliverables only. Clearing it prevents a
# superseded title or number from surviving a reproducible rerun.
if (dir.exists(final_root)) unlink(final_root, recursive = TRUE, force = TRUE)
dir.create(final_main, recursive = TRUE, showWarnings = FALSE)
dir.create(final_supp, recursive = TRUE, showWarnings = FALSE)
dir.create(final_qa, recursive = TRUE, showWarnings = FALSE)

# The helper's manifest is reset so this manifest contains only the final set.
figure_manifest_rows <- list()
final_script <-
  "R/25_generate_final_dissertation_figure_set_correct_numbering.R"

set_final_title <- function(plot, title) {
  # Align titles to the full figure rather than the plotting panel. This keeps
  # long dataset-specific titles clear of the right boundary on 180-mm pages.
  title_size <- if (nchar(title) >= 78L) {
    9.4
  } else if (nchar(title) >= 68L) {
    10.0
  } else {
    11.0
  }
  final_title_theme <- theme(
    plot.title = element_text(
      size = title_size, face = "bold", hjust = 0,
      margin = margin(b = 5)
    ),
    plot.title.position = "plot"
  )
  if (inherits(plot, "patchwork")) {
    plot + patchwork::plot_annotation(
      title = title,
      subtitle = NULL,
      caption = NULL,
      theme = figure_title_theme + final_title_theme
    )
  } else {
    plot +
      ggplot2::labs(title = title, subtitle = NULL, caption = NULL) +
      final_title_theme
  }
}

export_final <- function(plot, directory, stem, title, width_mm, height_mm,
                         sources) {
  final_plot <- set_final_title(plot, title)
  save_figure_four_formats(
    final_plot, directory, stem, width_mm, height_mm, sources,
    source_script = final_script
  )
}

dataset_slug <- c(
  sputum = "sputum_gpl570",
  brushing_gpl570 = "bronchial_brushing_gpl570",
  biopsy = "bronchial_biopsy_gpl570",
  brushing_rnaseq = "bronchial_brushing_rnaseq"
)
dataset_title <- c(
  sputum = "Sputum GPL570",
  brushing_gpl570 = "Bronchial brushing GPL570",
  biopsy = "Bronchial biopsy GPL570",
  brushing_rnaseq = "Bronchial brushing RNA-seq"
)

specifications <- list()
add_spec <- function(section, order, figure_id, title, stem, plot, width_mm,
                     height_mm, source_plot_basis, sources) {
  specifications[[length(specifications) + 1L]] <<- list(
    section = section,
    section_order = order,
    figure_id = figure_id,
    title = title,
    stem = stem,
    plot = plot,
    width_mm = width_mm,
    height_mm = height_mm,
    source_plot_basis = source_plot_basis,
    sources = sources
  )
}

# ---------------------------------------------------------------------------
# Main figures: nine numbered figure groups. Cohort size, outcome availability
# and overlap are separate Figures 1-3. Four-dataset results retain A-D panels
# in the fixed anatomical/platform order.
# ---------------------------------------------------------------------------
add_spec(
  "Main", 1L, "Figure 1",
  "Figure 1. Airway transcriptomic dataset sizes",
  "01_Figure_1_airway_transcriptomic_dataset_sizes",
  p1a, 180, 92, "Clean separate R Figure 1",
  "results/tables/Table_S_sample_availability.csv"
)
add_spec(
  "Main", 2L, "Figure 2",
  "Figure 2. Clinical outcome availability by airway dataset",
  "02_Figure_2_clinical_outcome_availability_by_airway_dataset",
  p1b, 180, 118, "Clean separate R Figure 2",
  "results/tables/Table_S_outcome_missingness.csv"
)
add_spec(
  "Main", 3L, "Figure 3",
  "Figure 3. Exact participant-overlap patterns across airway datasets",
  "03_Figure_3_exact_participant_overlap_patterns",
  p1c, 180, 142, "Clean separate R Figure 3",
  "results/tables/Table_S_exact_dataset_overlap.csv"
)

for (i in seq_along(dataset_order)) {
  dataset_name <- dataset_order[i]
  letter <- dataset_letters[dataset_name]
  add_spec(
    "Main", 3L + i, paste0("Figure 4", letter),
    paste0(
      "Figure 4", letter, ". Normalized ssGSEA pathway activity: ",
      dataset_title[dataset_name]
    ),
    paste0(
      sprintf("%02d", 3L + i), "_Figure_4", letter, "_",
      dataset_slug[dataset_name], "_normalized_ssgsea_pathway_activity"
    ),
    draft_figure_plot_objects[[paste0("Figure_2", letter)]],
    180, 135, "Validated separate R heatmap",
    paste0(
      "results/ssgsea/", dataset_name,
      "_ssgsea_normalized_scores.rds"
    )
  )
}

for (i in seq_along(dataset_order)) {
  dataset_name <- dataset_order[i]
  letter <- dataset_letters[dataset_name]
  add_spec(
    "Main", 7L + i, paste0("Figure 5", letter),
    paste0(
      "Figure 5", letter, ". Pathway-clinical Spearman associations: ",
      dataset_title[dataset_name]
    ),
    paste0(
      sprintf("%02d", 7L + i), "_Figure_5", letter, "_",
      dataset_slug[dataset_name], "_pathway_clinical_spearman_associations"
    ),
    draft_figure_plot_objects[[paste0("Figure_3", letter)]],
    180, 178, "Validated separate R Spearman heatmap",
    "results/tables/Table_2_pathway_clinical_spearman_associations.csv"
  )
}

for (i in seq_along(dataset_order)) {
  dataset_name <- dataset_order[i]
  letter <- dataset_letters[dataset_name]
  add_spec(
    "Main", 11L + i, paste0("Figure 6", letter),
    paste0(
      "Figure 6", letter, ". Prespecified curated pathway models: ",
      dataset_title[dataset_name]
    ),
    paste0(
      sprintf("%02d", 11L + i), "_Figure_6", letter, "_",
      dataset_slug[dataset_name], "_prespecified_curated_pathway_models"
    ),
    draft_figure_plot_objects[[paste0("Figure_4", letter)]],
    180, 205, "Validated separate R curated-model forest plot",
    c(
      "results/tables/Table_3_curated_model_coefficients_95CI.csv",
      "results/tables/Table_S_curated_model_formulas.csv"
    )
  )
}

for (i in seq_along(dataset_order)) {
  dataset_name <- dataset_order[i]
  letter <- dataset_letters[dataset_name]
  add_spec(
    "Main", 15L + i, paste0("Figure 7", letter),
    paste0(
      "Figure 7", letter, ". Top-10 Elastic Net pathway-selection stability: ",
      dataset_title[dataset_name]
    ),
    paste0(
      sprintf("%02d", 15L + i), "_Figure_7", letter, "_",
      dataset_slug[dataset_name], "_top10_elastic_net_stability"
    ),
    draft_figure_plot_objects[[paste0("Figure_5", letter)]],
    180, 240, "Validated separate R top-10 Elastic Net stability plot",
    c(
      "results/tables/Table_5_elastic_net_selection_stability.csv",
      "results/tables/Table_4_repeated_nested_CV_model_performance.csv"
    )
  )
}

add_spec(
  "Main", 20L, "Figure 8",
  "Figure 8. Internally held-out predictive performance",
  "20_Figure_8_internally_held_out_predictive_performance",
  figure6, 180, 225, "Clean R predictive-performance figure",
  c(
    "results/tables/Table_7_pathway_vs_biomarker_model_comparison.csv",
    "results/tables/Table_S_combined_vs_biomarker_delta_inference.csv",
    "results/analysis_data/nested_cv_averaged_predictions.csv.gz"
  )
)

add_spec(
  "Main", 21L, "Figure 9A",
  "Figure 9A. Same-platform cross-compartment pathway concordance",
  "21_Figure_9A_same_platform_cross_compartment_pathway_concordance",
  draft_figure_plot_objects[["Figure_7A"]], 180, 175,
  "Validated separate R cross-compartment concordance heatmap",
  "results/tables/Table_5_matched_pathway_concordance.csv"
)
add_spec(
  "Main", 22L, "Figure 9B",
  "Figure 9B. Brushing cross-platform pathway concordance",
  "22_Figure_9B_brushing_cross_platform_pathway_concordance",
  draft_figure_plot_objects[["Figure_7B"]], 180, 185,
  "Validated separate R cross-platform concordance plot",
  c(
    "results/tables/Table_5_matched_pathway_concordance.csv",
    "results/tables/Table_S_signature_gene_coverage_by_dataset.csv"
  )
)

# ---------------------------------------------------------------------------
# Supplementary figures: final S1-S11 numbering agreed for the dissertation.
# ---------------------------------------------------------------------------
add_spec(
  "Supplementary", 1L, "Figure S1",
  "Supplementary Figure S1. PCA of z-standardised ssGSEA scores",
  "01_Figure_S1_pathway_PCA",
  draft_figure_plot_objects[["Figure_S1"]], 180, 150,
  "Validated R PCA",
  c(
    "results/analysis_data/PCA_participant_scores.csv",
    "results/tables/Table_S_PCA_variance_explained.csv",
    "results/analysis_data/*_pathway_zscores_for_PCA.csv"
  )
)

calibration_titles <- c(
  sputum = "Sputum GPL570",
  brushing_gpl570 = "Bronchial brushing GPL570",
  biopsy = "Bronchial biopsy GPL570",
  brushing_rnaseq = "Bronchial brushing RNA-seq"
)
for (i in seq_along(dataset_order)) {
  dataset_name <- dataset_order[i]
  final_number <- i + 1L
  add_spec(
    "Supplementary", final_number, paste0("Figure S", final_number),
    paste0(
      "Supplementary Figure S", final_number, ". Binned calibration: ",
      calibration_titles[dataset_name]
    ),
    paste0(
      sprintf("%02d", final_number), "_Figure_S", final_number, "_",
      dataset_slug[dataset_name], "_binned_calibration"
    ),
    draft_figure_plot_objects[[paste0("Figure_S", final_number)]],
    180, 150, "Validated dataset-specific R calibration plot",
    c(
      "results/tables/Table_S_predictive_calibration_bins.csv",
      "results/analysis_data/nested_cv_averaged_predictions.csv.gz"
    )
  )
}

add_spec(
  "Supplementary", 6L, "Figure S6",
  "Supplementary Figure S6. Normalized versus z-standardised ssGSEA heatmaps",
  "06_Figure_S6_normalized_vs_zstandardised_ssgsea_heatmaps",
  figure_s2, 180, 300, "Clean R score-scaling comparison",
  unlist(lapply(dataset_order, function(dataset_name) c(
    paste0(
      "results/ssgsea/", dataset_name,
      "_ssgsea_normalized_scores.rds"
    ),
    paste0(
      "results/ssgsea/", dataset_name,
      "_ssgsea_pathway_zscores.rds"
    )
  )))
)
add_spec(
  "Supplementary", 7L, "Figure S7",
  "Supplementary Figure S7. Between-participant variation in normalized ssGSEA scores",
  "07_Figure_S7_between_participant_ssgsea_score_variation",
  p_s3, 180, 180, "Clean R pathway-variation sensitivity",
  "results/tables/Table_S_pathway_score_variation_SD_IQR.csv"
)
add_spec(
  "Supplementary", 8L, "Figure S8",
  "Supplementary Figure S8. Covariate-adjusted curated-model sensitivity",
  "08_Figure_S8_covariate_adjusted_curated_model_sensitivity",
  figure_s4, 180, 245, "Clean R covariate-adjusted sensitivity",
  "results/tables/Table_S_covariate_adjusted_curated_coefficients.csv"
)
add_spec(
  "Supplementary", 9L, "Figure S9",
  "Supplementary Figure S9. Elastic Net stability across all 39 pathways",
  "09_Figure_S9_complete_all39_elastic_net_stability",
  p5, 180, 245, "Clean R all-39 Elastic Net stability plot",
  "results/tables/Table_5_elastic_net_selection_stability.csv"
)
add_spec(
  "Supplementary", 10L, "Figure S10",
  "Supplementary Figure S10. RNA-seq signature-gene coverage",
  "10_Figure_S10_rnaseq_signature_gene_coverage",
  p_s9, 180, 185, "Clean R RNA-seq coverage audit",
  "results/tables/Table_S_signature_gene_coverage_by_dataset.csv"
)
add_spec(
  "Supplementary", 11L, "Figure S11",
  "Supplementary Figure S11. Below-zero prediction sensitivity",
  "11_Figure_S11_below_zero_prediction_sensitivity",
  p_s10, 180, 160, "Clean R below-zero prediction sensitivity",
  "results/tables/Table_S_below_zero_prediction_clipping_sensitivity.csv"
)

# Export all 33 final figure files (22 main panel files + 11 supplementary).
for (specification in specifications) {
  destination <- if (specification$section == "Main") final_main else final_supp
  export_final(
    specification$plot, destination, specification$stem,
    specification$title, specification$width_mm, specification$height_mm,
    specification$sources
  )
}

order_table <- data.table::rbindlist(lapply(specifications, function(x) {
  data.table::data.table(
    section = x$section,
    section_order = x$section_order,
    figure_id = x$figure_id,
    final_title = x$title,
    figure_stem = x$stem,
    source_plot_basis = x$source_plot_basis,
    width_mm = x$width_mm,
    height_mm = x$height_mm,
    source_data_files = paste(x$sources, collapse = "; ")
  )
}))
order_table[, section_rank := match(section, c("Main", "Supplementary"))]
data.table::setorder(order_table, section_rank, section_order)
order_table[, section_rank := NULL]

final_export_manifest <- data.table::rbindlist(figure_manifest_rows, fill = TRUE)
final_export_manifest <- merge(
  final_export_manifest,
  order_table[, .(
    figure_stem, section, section_order, figure_id, final_title,
    source_plot_basis
  )],
  by = "figure_stem", all.x = TRUE, sort = FALSE
)
final_export_manifest[, `:=`(
  section_rank = match(section, c("Main", "Supplementary")),
  format_rank = match(output_format, c("pdf", "svg", "png", "tiff"))
)]
data.table::setorder(
  final_export_manifest, section_rank, section_order, format_rank
)
final_export_manifest[, c("section_rank", "format_rank") := NULL]
data.table::fwrite(
  order_table, file.path(final_root, "FINAL_DISSERTATION_FIGURE_ORDER.csv")
)
data.table::fwrite(
  final_export_manifest,
  file.path(final_root, "FINAL_DISSERTATION_FIGURE_EXPORT_MANIFEST.csv")
)

main_order_lines <- c(
  "1. Figure 1 — Airway transcriptomic dataset sizes",
  "2. Figure 2 — Clinical outcome availability by airway dataset",
  "3. Figure 3 — Exact participant-overlap patterns across airway datasets",
  "4. Figures 4A-D — Normalized ssGSEA pathway activity in the four airway datasets",
  "5. Figures 5A-D — Pathway-clinical Spearman associations",
  "6. Figures 6A-D — Prespecified curated pathway models",
  "7. Figures 7A-D — Top-10 Elastic Net pathway-selection stability",
  "8. Figure 8 — Internally held-out predictive performance",
  "9. Figures 9A-B — Cross-compartment and cross-platform pathway concordance"
)
supp_order_lines <- paste0(
  seq_len(11), ". ",
  c(
    "Figure S1 — PCA of z-standardised ssGSEA scores",
    "Figure S2 — Binned calibration: sputum GPL570",
    "Figure S3 — Binned calibration: bronchial brushing GPL570",
    "Figure S4 — Binned calibration: bronchial biopsy GPL570",
    "Figure S5 — Binned calibration: bronchial brushing RNA-seq",
    "Figure S6 — Normalized versus z-standardised ssGSEA heatmaps",
    "Figure S7 — Between-participant variation in normalized ssGSEA scores",
    "Figure S8 — Covariate-adjusted curated-model sensitivity",
    "Figure S9 — Elastic Net stability across all 39 pathways",
    "Figure S10 — RNA-seq signature-gene coverage",
    "Figure S11 — Below-zero prediction sensitivity"
  )
)
writeLines(
  c(
    "# Final dissertation figure order",
    "",
    "All scientific graphics, titles, panel assembly and exports were produced directly in R. No image was manually edited. The dissertation itself was not modified.",
    "",
    "## Main figures",
    "",
    main_order_lines,
    "",
    "## Supplementary figures",
    "",
    supp_order_lines
  ),
  file.path(final_root, "FINAL_DISSERTATION_FIGURE_ORDER.md")
)

legend_lines <- c(
  "# Final dissertation figure legends",
  "",
  "## Figure 1. Airway transcriptomic dataset sizes",
  "Bars show the numbers of adult U-BIOPRED transcriptomic profiles in sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570 and bronchial brushing RNA-seq. Bar height and the displayed label give n.",
  "",
  "## Figure 2. Clinical outcome availability by airway dataset",
  "Bars report the participants available for FEV1 % predicted, FeNO, blood eosinophils, sputum eosinophils and previous-12-month exacerbations relative to each transcriptomic dataset total. Differences in completeness explain outcome-specific analysis denominators.",
  "",
  "## Figure 3. Exact participant-overlap patterns across airway datasets",
  "Bars show mutually exclusive exact transcriptomic-dataset overlap patterns. The 495 profiles represented 243 unique participants; the four datasets partially overlap and must not be treated as independent cohorts. Baseline characteristics are reported separately in Table 3 and sample/omics preprocessing in Table 1.",
  "",
  "## Figure 4. Normalized ssGSEA pathway activity",
  "Panels A-D show the 39 fixed asthma-relevant pathway scores across participants in sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570 and bronchial brushing RNA-seq, respectively. Blue indicates lower, white approximately central and red higher saved normalized ssGSEA scores on common -1 to +1 display limits. Hierarchical ordering is descriptive. The panels show graded variation across continuous-valued scores but do not test for, prove or exclude discrete latent classes or endotypes.",
  "",
  "## Figure 5. Pathway-clinical Spearman associations",
  "Panels A-D show marginal pathway-outcome Spearman correlations in the four airway datasets. Blue indicates negative and red positive rho; plus/minus glyphs preserve direction in grayscale. Asterisks mark Benjamini-Hochberg false-discovery-rate support within each dataset-outcome family (* q<0.05, ** q<0.01, *** q<0.001). These are univariable associations, not conditional effects, causal estimates or predictive performance.",
  "",
  "## Figure 6. Prespecified curated pathway models",
  "Panels A-D show standardised coefficients and residual-degrees-of-freedom t-based 95% confidence intervals for biologically prespecified outcome-specific pathway subsets in each dataset. Filled red symbols indicate within-family BH q<0.05 and open grey symbols q>=0.05. The plots do not show all 39 pathways. Coefficients are conditional on other pathways in the same model; correlated signatures can produce collinearity or suppression.",
  "",
  "## Figure 7. Elastic Net pathway-selection stability",
  "Panels A-D show the ten pathways with the highest pathway-only Elastic Net selection frequency for each outcome in each dataset. Selection frequency is the proportion of 50 outer-training fits (10 repeats x 5 folds) with a non-zero coefficient and is not a p value. Red upward triangles indicate predominantly positive coefficients, blue downward triangles predominantly negative coefficients, and point size is a coefficient-based relative share rather than variance decomposition or causal importance.",
  "",
  "## Figure 8. Internally held-out predictive performance",
  "Panel A compares mean model-scale cross-validated R-squared for biomarker-only, curated-pathway-only and Elastic-Net-pathway-only models on identical common cases and outer folds. Panel B shows paired change in model-scale R-squared when curated or Elastic Net pathways are added to the biomarker baseline. Intervals are participant-bootstrap intervals conditional on fixed averaged repeated out-of-fold predictions and do not include refitting, tuning or fold-generation uncertainty. The analysis is exploratory internal validation; concurrent sputum-eosinophil estimation is same-sample phenotype capture rather than prognosis.",
  "",
  "## Figure 9. Compartment and platform concordance",
  "Panel A shows same-pathway matched-participant Spearman concordance among GPL570 airway compartments. Panel B shows same-compartment concordance between bronchial brushing GPL570 and RNA-seq in 118 matched participants, with low RNA-seq signature coverage flagged. These are same-participant concordance analyses, not independent replication or external validation; matched sputum subsets are small.",
  "",
  "## Supplementary Figure S1. PCA of z-standardised ssGSEA scores",
  "PCA was run separately in each airway dataset after pathway-wise z-standardisation. PC1 and PC2 percentages are reported on the axes. The projection is descriptive; axes are dataset-specific and do not test for, prove or exclude latent classes or endotypes.",
  "",
  "## Supplementary Figures S2-S5. Binned internal calibration",
  "Figures S2-S5 show calibration for sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570 and bronchial brushing RNA-seq, respectively. Five equal-frequency bins summarise observed versus predicted values for biomarker-only, curated-pathway and Elastic-Net-pathway models. Each participant's ten held-out predictions were averaged before binning. Bars are standard errors of observed bin means on the modelling scale, not uncertainty from model refitting.",
  "",
  "## Supplementary Figure S6. Normalized versus z-standardised ssGSEA heatmaps",
  "Paired panels use identical participant and pathway ordering. Normalized panels retain between-pathway score location and dispersion; z-standardised panels centre each pathway to mean zero and standard deviation one within its dataset. Blue indicates lower, white central and red higher values on the labelled scales.",
  "",
  "## Supplementary Figure S7. Between-participant variation in normalized ssGSEA scores",
  "Points report the between-participant standard deviation of each normalized ssGSEA pathway score before z-standardisation. Larger values indicate greater participant-to-participant spread on the saved score scale. This is descriptive and is not a test of biological importance or cross-dataset calibration.",
  "",
  "## Supplementary Figure S8. Covariate-adjusted curated-model sensitivity",
  "Primary prespecified pathway subsets were refitted with age, sex, smoking status, BMI and current oral-corticosteroid exposure where complete. Symbols and intervals follow Figure 6. Changes can reflect both covariate adjustment and differences in the covariate-complete sample.",
  "",
  "## Supplementary Figure S9. Elastic Net stability across all 39 pathways",
  "Selection frequency is the proportion of 50 pathway-only Elastic Net outer-training fits in which a pathway coefficient was non-zero. Blue denotes a predominantly negative coefficient and red a predominantly positive coefficient; colour intensity gives selection frequency and plus/minus glyphs redundantly show sign. Selection frequency is a stability descriptor, not a p value.",
  "",
  "## Supplementary Figure S10. RNA-seq signature-gene coverage",
  "Bars show the percentage of each fixed signature detected in the processed bronchial brushing RNA-seq matrix. Labels give detected/total harmonised genes and LOW marks the prespecified low-coverage flag. All retained signatures met the ssGSEA minimum of three genes, but low proportional coverage requires pathway-specific caution.",
  "",
  "## Supplementary Figure S11. Below-zero prediction sensitivity",
  "Bars show counts of below-zero participant-level averaged predictions for clinical outcomes that cannot be negative; FEV1 is excluded. Unconstrained predictions remain primary. Clipping was evaluated only as a sensitivity for original-scale error metrics, while model-scale R-squared was never clipped."
)
writeLines(
  legend_lines,
  file.path(final_root, "FINAL_DISSERTATION_FIGURE_LEGENDS.md")
)

writeLines(
  c(
    "# Final dissertation figures — R-only set",
    "",
    "This folder contains the final dissertation figures in their final numbering and order.",
    "",
    "- Scientific calculations, plotting, titles, panel assembly and exports were produced directly in R.",
    "- Every scientific image was generated directly by the R workflow without manual editing.",
    "- Every figure is supplied as PDF, SVG, PNG at 300 dpi and TIFF at 600 dpi.",
    "- Cohort size, outcome availability and participant overlap are separate Figures 1, 2 and 3.",
    "- Every multi-compartment display uses the fixed order: sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570, bronchial brushing RNA-seq.",
    "- Explanatory detail is in FINAL_DISSERTATION_FIGURE_LEGENDS.md rather than small grey text inside the plots.",
    "- The dissertation Word document was not edited.",
    "",
    "Use PDF for the clearest vector insertion into Word where supported; use PNG300 if Word compatibility is preferred."
  ),
  file.path(final_root, "README.md")
)

# Structural title/name/order checks. Visual rendering is completed by the
# companion QA script after all files have been exported.
expected_stems <- vapply(specifications, `[[`, character(1), "stem")
if (length(expected_stems) != 33L || anyDuplicated(expected_stems)) {
  stop("The final set must contain 33 unique figure files.")
}
if (!setequal(unique(final_export_manifest$figure_stem), expected_stems)) {
  stop("The final export manifest does not match the agreed figure set.")
}
if (nrow(final_export_manifest) != 132L) {
  stop("Expected 33 figures x four formats = 132 exports.")
}
if (any(final_export_manifest$output_bytes <= 0L)) {
  stop("At least one final dissertation figure export is empty.")
}
if (any(final_export_manifest$generated_in != "R") ||
    any(final_export_manifest$hand_edited)) {
  stop("Every final dissertation figure must be an unedited R export.")
}

writeLines(
  c(
    "FINAL DISSERTATION FIGURE GENERATION: COMPUTATION PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "Main figure groups: 9",
    "Main panel files: 22",
    "Supplementary figures: 11",
    "Unique final figure stems: 33",
    "Production exports: 132 (PDF, SVG, PNG300 and TIFF600)",
    "Figures 1-3 were exported as three separately numbered full-width figures.",
    "Canonical multi-compartment order: sputum GPL570; bronchial brushing GPL570; bronchial biopsy GPL570; bronchial brushing RNA-seq.",
    "Titles and final numbering were applied in R.",
    "Every scientific image was generated directly by the R workflow without manual editing.",
    "The dissertation Word document was not modified."
  ),
  file.path(final_root, "FINAL_DISSERTATION_FIGURE_GENERATION_STATUS.txt")
)

message("Generated 33 final dissertation figures in 132 R-produced formats.")
