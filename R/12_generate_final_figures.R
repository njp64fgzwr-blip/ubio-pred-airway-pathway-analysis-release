# Generate the dissertation figure set directly in R from validated sources.

source(file.path("R", "00_common.R"), local = TRUE)
assert_packages(c(
  "data.table", "ggplot2", "patchwork", "svglite", "ragg", "scales"
))
source(file.path("R", "12_figure_helpers.R"), local = TRUE)
message_rule("12: Generating base dissertation figures in R")

# Figure 1: sample availability, outcome completeness and exact overlap.
availability <- data.table::fread(file.path(
  paths$tables, "Table_S_sample_availability.csv"
))
missingness <- data.table::fread(file.path(
  paths$tables, "Table_S_outcome_missingness.csv"
))
overlap <- data.table::fread(file.path(
  paths$tables, "Table_S_exact_dataset_overlap.csv"
))

availability[, dataset_factor := factor(dataset, levels = rev(dataset_order))]
availability[, display_label := wrap_dataset_label(dataset_label)]
availability_labels <- setNames(
  availability$display_label, as.character(availability$dataset_factor)
)
p1a <- ggplot(availability, aes(transcriptomic_profiles, dataset_factor)) +
  geom_col(width = 0.68, fill = "#3C78A8") +
  geom_text(
    aes(label = transcriptomic_profiles), hjust = 1.25,
    colour = "white", fontface = "bold", size = 3.1
  ) +
  scale_y_discrete(labels = availability_labels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(title = "A  Airway transcriptomic datasets", x = "Participants", y = NULL) +
  theme(panel.grid.major.y = element_blank())

missingness[, dataset_factor := factor(dataset, levels = dataset_order)]
missingness[, outcome_factor := factor(outcome, levels = rev(outcome_order))]
missingness[, text_colour := ifelse(missing_percent > 35, "white", "black")]
p1b <- ggplot(missingness, aes(dataset_factor, outcome_factor,
                               fill = missing_percent)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = paste0(available_n, "\n/", total_n), colour = text_colour),
            size = 2.35, lineheight = 0.85) +
  scale_colour_identity() +
  scale_fill_gradient(
    low = "#F7FBFF", high = "#08306B", limits = c(0, 60),
    name = "Missing (%)"
  ) +
  scale_x_discrete(labels = setNames(
    c(
      "Sputum\nGPL570", "Brush\nGPL570", "Biopsy\nGPL570",
      "Brush\nRNA-seq"
    ), dataset_order
  )) +
  scale_y_discrete(labels = setNames(outcome_labels, outcome_order)) +
  labs(
    title = "B  Outcome availability",
    subtitle = "Cell labels are available n / transcriptomic n",
    x = NULL, y = NULL
  ) +
  theme(axis.text.x = element_text(size = 7), panel.grid = element_blank())

pattern_labels <- c(
  sputum = "Sputum only",
  brushing_gpl570 = "Brushing GPL570 only",
  biopsy = "Biopsy only",
  `brushing_gpl570+biopsy` = "Brushing GPL570 + biopsy",
  `brushing_gpl570+brushing_rnaseq` = "Both brushing platforms",
  `brushing_gpl570+biopsy+brushing_rnaseq` = "Both brushings + biopsy",
  `sputum+brushing_gpl570` = "Sputum + brushing GPL570",
  `sputum+biopsy` = "Sputum + biopsy",
  `sputum+brushing_gpl570+biopsy` = "Sputum + brushing GPL570 + biopsy",
  `sputum+brushing_gpl570+brushing_rnaseq` = "Sputum + both brushing platforms",
  `sputum+brushing_gpl570+biopsy+brushing_rnaseq` = "All four datasets"
)
overlap[, pattern_label := unname(pattern_labels[exact_pattern])]
if (any(is.na(overlap$pattern_label))) stop("An overlap pattern lacks a label.")
data.table::setorder(overlap, participants)
overlap[, pattern_factor := factor(pattern_label, levels = pattern_label)]
p1c <- ggplot(overlap, aes(participants, pattern_factor)) +
  geom_col(width = 0.65, fill = "#6A3D9A") +
  geom_text(aes(label = participants), hjust = -0.25, size = 2.8) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "C  Exact participant overlap patterns",
    subtitle = "Each participant appears in one mutually exclusive row",
    x = "Participants", y = NULL
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7)
  )

figure1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    title = "Figure 1. Cohort structure and outcome completeness",
    subtitle = paste0(
      "All transcriptomic sample IDs matched clinical IDs. Missingness explains why\n",
      "complete-case n differs by outcome and dataset; baseline characteristics are in Table 1."
    ),
    caption = paste0(
      "The four datasets include 243 unique participants. Exact-overlap rows are descriptive\n",
      "and should not be summed as independent cohorts."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure1, paths$figures_main, "Figure_1_cohort_structure_and_missingness",
  180, 170,
  c(
    "results/tables/Table_S_sample_availability.csv",
    "results/tables/Table_S_outcome_missingness.csv",
    "results/tables/Table_S_exact_dataset_overlap.csv",
    "results/tables/Table_1_clinical_characteristics.csv"
  )
)
write_figure_source(availability, "Figure_1A_sample_availability.csv")
write_figure_source(missingness, "Figure_1B_outcome_availability.csv")
write_figure_source(overlap, "Figure_1C_exact_overlap_patterns.csv")

# Figure 2: normalized ssGSEA values, deliberately not pathway-wise z scores.
raw_heatmap_panels <- list()
raw_heatmap_source <- character()
for (dataset_index in seq_along(dataset_order)) {
  dataset_name <- dataset_order[dataset_index]
  matrix <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")
  ))
  long <- make_long_matrix(matrix, dataset_name, "Normalized ssGSEA")
  raw_heatmap_panels[[dataset_name]] <- make_heatmap_panel(
    long, c(-1, 1), datasets[[dataset_name]]$label,
    paste0(
      datasets[[dataset_name]]$compartment, " | ",
      datasets[[dataset_name]]$platform, " | n=", ncol(matrix)
    ),
    show_y_labels = dataset_index %in% c(1L, 3L),
    legend_title = "Normalized ssGSEA",
    show_legend = FALSE
  )
  raw_heatmap_source <- c(
    raw_heatmap_source,
    paste0("results/ssgsea/", dataset_name,
           "_ssgsea_normalized_scores.rds")
  )
}
figure2 <- (raw_heatmap_panels[[1]] | raw_heatmap_panels[[2]]) /
  (raw_heatmap_panels[[3]] | raw_heatmap_panels[[4]]) +
  plot_annotation(
    title = "Figure 2. Participant heterogeneity in normalized ssGSEA pathway scores",
    subtitle = paste0(
      "Rows are the fixed 39 signatures; columns are participants. All panels use the same\n",
      "-1 to +1 display limits and saved normalized ssGSEA values (blue lower; white central; red higher)."
    ),
    caption = paste0(
      "Participants showed graded variation across multiple continuous-valued pathway scores.\n",
      "Hierarchical ordering is descriptive and does not test for, prove or exclude discrete latent classes/endotypes."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure2, paths$figures_main, "Figure_2_normalized_ssgsea_pathway_heatmaps",
  180, 210, raw_heatmap_source
)

# Figure 3: all 780 pathway-clinical Spearman associations.
associations <- data.table::fread(file.path(
  paths$tables, "Table_2_pathway_clinical_spearman_associations.csv"
))
associations[, dataset_factor := factor(dataset, levels = dataset_order)]
associations[, outcome_factor := factor(outcome, levels = outcome_order)]
associations[, pathway_factor := factor(pathway, levels = rev(names(pathway_labels)))]
associations[, significance_marker := data.table::fcase(
  fdr_bh < 0.001, "***",
  fdr_bh < 0.01, "**",
  fdr_bh < 0.05, "*",
  default = ""
)]
associations[, direction_marker := paste0(
  ifelse(spearman_rho >= 0, "+", "-"), significance_marker
)]
associations[, marker_colour := ifelse(
  abs(spearman_rho) >= 0.55, "white", "black"
)]
n_labels <- associations[, .(n = unique(n)), by = .(dataset, outcome)]
n_labels[, display := paste0(outcome_labels_compact[outcome], "\nn=", n)]
axis_labels_by_dataset <- lapply(dataset_order, function(dataset_name) {
  x <- n_labels[dataset == dataset_name]
  setNames(x$display, x$outcome)
})
names(axis_labels_by_dataset) <- dataset_order

association_panels <- lapply(seq_along(dataset_order), function(i) {
  dataset_name <- dataset_order[i]
  x <- associations[dataset == dataset_name]
  ggplot(x, aes(outcome_factor, pathway_factor, fill = spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = direction_marker, colour = marker_colour),
      size = 1.75, fontface = "bold"
    ) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = association_colours["low"], mid = association_colours["mid"],
      high = association_colours["high"], midpoint = 0, limits = c(-1, 1),
      name = "Spearman rho"
    ) +
    scale_x_discrete(labels = axis_labels_by_dataset[[dataset_name]]) +
    scale_y_discrete(
      labels = if (i %in% c(1L, 3L)) pathway_labels else NULL,
      drop = FALSE
    ) +
    labs(
      title = datasets[[dataset_name]]$label,
      subtitle = paste(
        datasets[[dataset_name]]$compartment,
        datasets[[dataset_name]]$platform, sep = " | "
      ),
      x = NULL, y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 6.4),
      axis.text.y = element_text(size = 5.25),
      plot.title = element_text(size = 9, face = "bold"),
      plot.subtitle = element_text(size = 7),
      legend.position = "bottom"
    )
})
figure3 <- ((association_panels[[1]] | association_panels[[2]]) /
  (association_panels[[3]] | association_panels[[4]])) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Figure 3. Pathway-clinical Spearman associations by airway dataset",
    subtitle = paste0(
      "Colour is Spearman rho; +/- preserves direction in grayscale. Asterisks mark BH-FDR\n",
      "within each dataset-outcome family (* q<0.05, ** q<0.01, *** q<0.001)."
    ),
    caption = paste0(
      "These 780 tests describe marginal univariable associations. They are not conditional effects,\n",
      "causal estimates or predictive performance."
    ),
    theme = figure_title_theme
  ) & theme(legend.position = "bottom")
save_figure_four_formats(
  figure3, paths$figures_main, "Figure_3_pathway_clinical_spearman_associations",
  180, 220,
  "results/tables/Table_2_pathway_clinical_spearman_associations.csv"
)
write_figure_source(associations, "Figure_3_spearman_associations.csv")

# Figure 4: curated prespecified multivariable pathway models.
curated <- data.table::fread(file.path(
  paths$tables, "Table_3_curated_model_coefficients_95CI.csv"
))
curated[, outcome_factor := factor(outcome, levels = outcome_order)]
curated[, pathway_factor := factor(
  pathway_label, levels = rev(unique(pathway_label))
)]
curated[, fdr_supported := fdr_bh < 0.05]
outcome_strip_labels <- c(
  FEV1 = "FEV1",
  FeNO = "FeNO",
  blood_eosinophils = "Blood eosinophils",
  sputum_eosinophils = "Sputum eosinophils",
  exacerbation_frequency = "Prior-year exacerbations"
)
make_curated_forest_stack <- function(table, x_axis_label) {
  panels <- lapply(seq_along(outcome_order), function(outcome_index) {
    outcome_name <- outcome_order[outcome_index]
    x <- data.table::copy(table[outcome == outcome_name])
    n_by_dataset <- unique(x[, .(dataset, n)])
    data.table::setorder(n_by_dataset, dataset)
    n_by_dataset <- n_by_dataset[match(dataset_order, dataset)]
    n_labels_local <- setNames(
      paste0(
        dataset_labels_compact[n_by_dataset$dataset], "\nn=", n_by_dataset$n
      ),
      n_by_dataset$dataset
    )
    x[, dataset_n_factor := factor(
      dataset, levels = dataset_order, labels = n_labels_local[dataset_order]
    )]
    ggplot(x, aes(standardized_beta, pathway_factor)) +
      geom_vline(
        xintercept = 0, linetype = "dashed", linewidth = 0.35,
        colour = "#666666"
      ) +
      geom_errorbar(
        aes(xmin = ci_low, xmax = ci_high, colour = fdr_supported),
        orientation = "y", width = 0, linewidth = 0.62
      ) +
      geom_point(
        aes(colour = fdr_supported, shape = fdr_supported),
        size = 1.95, stroke = 0.7, fill = "white"
      ) +
      facet_grid(. ~ dataset_n_factor) +
      scale_colour_manual(
        values = c(`FALSE` = "#555555", `TRUE` = "#B2182B"),
        breaks = c(FALSE, TRUE),
        labels = c("q>=0.05", "BH-FDR q<0.05"),
        drop = FALSE,
        name = NULL
      ) +
      scale_shape_manual(
        values = c(`FALSE` = 21, `TRUE` = 16),
        breaks = c(FALSE, TRUE),
        labels = c("q>=0.05", "BH-FDR q<0.05"),
        drop = FALSE,
        name = NULL,
        guide = "none"
      ) +
      scale_x_continuous(
        limits = c(-1.08, 1.23), breaks = c(-1, -0.5, 0, 0.5, 1),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      labs(
        title = outcome_strip_labels[outcome_name],
        x = if (outcome_index == length(outcome_order)) x_axis_label else NULL,
        y = NULL
      ) +
      theme(
        plot.title = element_text(size = 8.5, face = "bold", margin = margin(b = 2)),
        strip.text = element_text(size = 6.5, face = "bold", lineheight = 0.9),
        axis.text.y = element_text(size = 5.7),
        axis.text.x = element_text(size = 5.8),
        axis.title.x = element_text(size = 8),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        panel.spacing.x = grid::unit(3.0, "mm"),
        plot.margin = margin(1, 2, 1, 2)
      )
  })
  wrap_plots(panels, ncol = 1, guides = "collect") &
    theme(legend.position = "none") &
    guides(colour = "none", shape = "none")
}

p4 <- make_curated_forest_stack(
  curated, "Standardised beta (95% residual-df t CI)"
)
figure4 <- p4 +
  plot_annotation(
    title = "Figure 4. Prespecified curated multivariable pathway models",
    subtitle = paste0(
      "Only each outcome's biologically prespecified pathway subset is shown; these were not selected\n",
      "as the strongest results from all 39 signatures. Filled red = family BH-FDR q<0.05; open grey = q>=0.05."
    ),
    caption = paste0(
      "Coefficients are conditional on the other displayed pathways. Correlated signatures can cause\n",
      "collinearity or suppression; the negative adjusted IL-5 coefficient is not evidence of biological protection."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure4, paths$figures_main, "Figure_4_curated_pathway_model_forest_plots",
  180, 245,
  c(
    "results/tables/Table_3_curated_model_coefficients_95CI.csv",
    "results/tables/Table_S_curated_model_formulas.csv"
  )
)
write_figure_source(curated, "Figure_4_curated_coefficients.csv")

# Figure 5: complete pathway-only Elastic Net selection frequencies.
stability <- data.table::fread(file.path(
  paths$tables, "Table_5_elastic_net_selection_stability.csv"
))
stability <- stability[model_type == "Elastic Net pathways"]
stability[, outcome_factor := factor(outcome, levels = outcome_order)]
stability[, pathway_factor := factor(pathway, levels = rev(names(pathway_labels)))]
stability[, sign_glyph := data.table::fcase(
  selection_count == 0L, "",
  predominant_sign == "Positive", "+",
  predominant_sign == "Negative", "-",
  default = "+/-"
)]
stability[, sign_colour := data.table::fcase(
  predominant_sign == "Positive", "#B2182B",
  predominant_sign == "Negative", "#2166AC",
  predominant_sign == "Mixed", "#762A83",
  default = "white"
)]
stability_panels <- lapply(seq_along(dataset_order), function(i) {
  dataset_name <- dataset_order[i]
  x <- stability[dataset == dataset_name]
  ggplot(x, aes(outcome_factor, pathway_factor, alpha = selection_frequency)) +
    geom_tile(aes(fill = sign_colour), colour = "#D9D9D9", linewidth = 0.22) +
    geom_text(
      aes(label = sign_glyph), colour = "black", size = 1.75,
      fontface = "bold"
    ) +
    scale_fill_identity() +
    scale_alpha_continuous(
      range = c(0.05, 1), limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1),
      name = "Selection frequency"
    ) +
    scale_x_discrete(labels = outcome_labels_compact) +
    scale_y_discrete(
      labels = if (i %in% c(1L, 3L)) pathway_labels else NULL,
      drop = FALSE
    ) +
    labs(title = datasets[[dataset_name]]$label, x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 6.1),
      axis.text.y = element_text(size = 5.15),
      plot.title = element_text(size = 9, face = "bold"),
      legend.position = "bottom"
    )
})
p5 <- (stability_panels[[1]] | stability_panels[[2]]) /
  (stability_panels[[3]] | stability_panels[[4]]) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Figure 5. Elastic Net selection stability across all 39 pathways",
    subtitle = paste0(
      "Selection frequency is the proportion of 50 outer models (10 repeats x 5 folds) with a non-zero\n",
      "coefficient. +/- shows the predominant selected-coefficient sign; +/- in one cell denotes a mixed sign."
    ),
    caption = paste0(
      "Selection frequency is a stability descriptor, not a p value. A frequently selected pathway can\n",
      "occur in a model with poor or negative held-out R-squared; predictive performance is shown in Figure 6."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  p5, paths$figures_main, "Figure_5_elastic_net_pathway_selection_stability",
  180, 245,
  "results/tables/Table_5_elastic_net_selection_stability.csv"
)
write_figure_source(stability, "Figure_5_all_pathway_selection_stability.csv")

# Figure 6 and required Figure 8: predictive performance versus biomarkers.
comparison <- data.table::fread(file.path(
  paths$tables, "Table_7_pathway_vs_biomarker_model_comparison.csv"
))
delta <- data.table::fread(file.path(
  paths$tables, "Table_S_combined_vs_biomarker_delta_inference.csv"
))
comparison[, dataset_factor := factor(dataset, levels = dataset_order)]
comparison[, outcome_factor := factor(outcome, levels = outcome_order)]
comparison[, dataset_factor_compact := factor(
  dataset, levels = dataset_order, labels = dataset_labels_compact[dataset_order]
)]
comparison[, outcome_factor_compact := factor(
  outcome, levels = outcome_order, labels = gsub("\n", " ", outcome_labels_compact[outcome_order])
)]
comparison[, model_factor := factor(
  model_type,
  levels = c("Biomarker-only", "Curated pathways", "Elastic Net pathways")
)]
p6a <- ggplot(
  comparison,
  aes(
    model_factor, mean_cross_validated_r_squared,
    colour = model_factor, shape = model_factor
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 2.1) +
  facet_grid(
    dataset_factor_compact ~ outcome_factor_compact
  ) +
  scale_colour_manual(values = model_colours, guide = "none") +
  scale_shape_manual(
    values = c(
      "Biomarker-only" = 15,
      "Curated pathways" = 16,
      "Elastic Net pathways" = 17
    ),
    guide = "none"
  ) +
  scale_x_discrete(labels = c(
    "Biomarker-only" = "Biomarker",
    "Curated pathways" = "Curated",
    "Elastic Net pathways" = "Elastic Net"
  )) +
  labs(
    title = "A  Three standalone model classes",
    subtitle = "Mean R-squared across 10 repeated outer validations",
    x = NULL, y = "Mean cross-validated R-squared"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
    axis.text.y = element_text(size = 6),
    strip.text.x = element_text(size = 6.0),
    strip.text.y = element_text(size = 5.8, angle = 0),
    panel.grid.minor = element_blank()
  )

delta[, augmented_short := factor(
  augmented_model,
  levels = c(
    "Biomarker + curated pathways",
    "Biomarker + Elastic Net pathways"
  ),
  labels = c("Biomarker + curated", "Biomarker + Elastic Net")
)]
delta[, outcome_factor := factor(outcome, levels = outcome_order)]
delta[, dataset_factor := factor(dataset, levels = dataset_order)]
delta[, dataset_factor_compact := factor(
  dataset, levels = dataset_order, labels = dataset_labels_compact[dataset_order]
)]
p6b <- ggplot(
  delta,
  aes(delta_model_r_squared_point, outcome_factor,
      colour = augmented_short, shape = augmented_short)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_errorbar(
    aes(
      xmin = delta_model_r_squared_lower_95,
      xmax = delta_model_r_squared_upper_95
    ),
    orientation = "y", width = 0.16, linewidth = 0.5,
    position = position_dodge(width = 0.35)
  ) +
  geom_point(size = 2.0, position = position_dodge(width = 0.35)) +
  facet_wrap(
    ~ dataset_factor_compact, ncol = 2
  ) +
  scale_y_discrete(labels = outcome_labels) +
  scale_colour_manual(
    values = c(
      "Biomarker + curated" =
        unname(model_colours["Biomarker + curated pathways"]),
      "Biomarker + Elastic Net" =
        unname(model_colours["Biomarker + Elastic Net pathways"])
    ),
    name = NULL
  ) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(
    title = "B  Added pathways above the biomarker baseline",
    subtitle = paste0(
      "Delta R-squared from fixed averaged repeated held-out predictions; bars are\n",
      "conditional 95% participant-bootstrap intervals"
    ),
    x = "Delta R-squared versus biomarker-only", y = NULL
  ) +
  theme(
    axis.text = element_text(size = 6.3),
    strip.text = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  )

figure6 <- p6a / p6b +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = "Figure 6. Internally held-out performance is outcome- and model-specific",
    subtitle = paste0(
      "All model classes use identical common cases and outer folds within each comparison;\n",
      "finite negative R-squared values are retained."
    ),
    caption = paste0(
      "Bootstrap intervals are conditional on fixed averages of ten held-out predictions and exclude refitting, tuning and fold-generation\n",
      "uncertainty. This is exploratory internal validation. Sputum-eosinophil prediction is concurrent phenotype capture, not prognosis."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure6, paths$figures_main, "Figure_6_internal_predictive_performance",
  180, 225,
  c(
    "results/tables/Table_7_pathway_vs_biomarker_model_comparison.csv",
    "results/tables/Table_S_combined_vs_biomarker_delta_inference.csv",
    "results/analysis_data/nested_cv_averaged_predictions.csv.gz"
  )
)

p8 <- ggplot(
  comparison,
  aes(
    model_factor, mean_cross_validated_r_squared,
    colour = model_factor, shape = model_factor
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 2.5) +
  facet_grid(
    dataset_factor_compact ~ outcome_factor_compact
  ) +
  scale_colour_manual(values = model_colours, name = "Model") +
  scale_shape_manual(
    values = c(
      "Biomarker-only" = 15,
      "Curated pathways" = 16,
      "Elastic Net pathways" = 17
    ),
    name = "Model"
  ) +
  scale_x_discrete(labels = c(
    "Biomarker-only" = "Biomarker",
    "Curated pathways" = "Curated",
    "Elastic Net pathways" = "Elastic Net"
  )) +
  labs(
    title = "Figure 8. Biomarker-only versus standalone pathway models",
    subtitle = paste0(
      "Mean repeated nested-CV R-squared on identical common cases and outer folds;\n",
      "no model class wins for every outcome or dataset."
    ),
    caption = paste0(
      "Negative R-squared means performance worse than the held-out mean benchmark. This is exploratory\n",
      "internal validation and does not establish external transportability or clinical utility."
    ),
    x = NULL, y = "Mean cross-validated R-squared"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.7),
    strip.text.x = element_text(size = 6.0),
    strip.text.y = element_text(size = 5.8, angle = 0),
    legend.position = "bottom"
  ) + figure_title_theme
save_figure_four_formats(
  p8, paths$figures_main, "Figure_8_pathway_vs_biomarker_model_performance",
  180, 155,
  c(
    "results/tables/Table_6_biomarker_only_model_performance.csv",
    "results/tables/Table_7_pathway_vs_biomarker_model_comparison.csv"
  )
)
write_figure_source(
  comparison, "Figure_6A_and_Figure_8_model_performance.csv"
)
write_figure_source(delta, "Figure_6B_delta_performance.csv")

# Figure 7: matched-participant compartment/platform concordance.
concordance <- data.table::fread(file.path(
  paths$tables, "Table_5_matched_pathway_concordance.csv"
))
concordance_summary <- data.table::fread(file.path(
  paths$tables, "Table_S_matched_concordance_summary.csv"
))
comparison_labels <- c(
  brushing_gpl570_vs_biopsy = "Brushing GPL570 vs biopsy\n(n=99)",
  sputum_vs_biopsy = "Sputum vs biopsy\n(n=23)",
  sputum_vs_brushing_gpl570 = "Sputum vs brushing GPL570\n(n=31)",
  brushing_gpl570_vs_brushing_rnaseq =
    "Brushing GPL570 vs RNA-seq\n(n=118)"
)
main_comparisons <- names(comparison_labels)
concordance_main <- concordance[comparison %in% main_comparisons]
concordance_main[, comparison_factor := factor(
  comparison, levels = main_comparisons
)]
concordance_main[, low_coverage_display := data.table::fcase(
  comparison == "brushing_gpl570_vs_brushing_rnaseq" &
    rnaseq_low_coverage_flag,
  "RNA-seq low coverage",
  default = "Standard coverage"
)]
p7a <- ggplot(concordance_main, aes(spearman_rho, comparison_factor)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_jitter(
    aes(shape = low_coverage_display, fill = low_coverage_display),
    width = 0, height = 0.12, size = 1.65, alpha = 0.8,
    colour = "black", stroke = 0.3
  ) +
  scale_y_discrete(labels = comparison_labels) +
  scale_shape_manual(
    values = c("Standard coverage" = 21, "RNA-seq low coverage" = 24),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Standard coverage" = "#3C78A8",
               "RNA-seq low coverage" = "#E69F00"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(-0.65, 1)) +
  labs(
    title = "A  Pathway-level matched-participant concordance",
    x = "Spearman rho across matched participants", y = NULL
  ) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

summary_main <- concordance_summary[comparison %in% main_comparisons]
summary_main[, comparison_factor := factor(comparison, levels = main_comparisons)]
p7b <- ggplot(summary_main, aes(median_spearman_rho, comparison_factor)) +
  geom_col(width = 0.58, fill = "#4D4D4D") +
  geom_text(
    aes(label = sprintf(
      "median %.2f\n%d/39 q<0.05",
      median_spearman_rho, fdr_supported_pathways
    )),
    hjust = -0.08, size = 2.45
  ) +
  scale_y_discrete(labels = comparison_labels) +
  scale_x_continuous(limits = c(0, 1.03), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "B  Summary across the fixed 39 pathways",
    x = "Median pathway concordance", y = NULL
  ) +
  theme(panel.grid.major.y = element_blank())

figure7 <- p7a / p7b +
  plot_annotation(
    title = "Figure 7. Matched pathway concordance across airway datasets",
    subtitle = paste0(
      "Same-compartment brushing platform concordance exceeds most cross-compartment concordance.\n",
      "All correlations compare the same pathway between matched participants."
    ),
    caption = paste0(
      "These are same-participant concordance analyses, not independent replication or external validation.\n",
      "The smaller sputum pairings are less precise than the n=118 brushing-platform comparison."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure7, paths$figures_main, "Figure_7_matched_pathway_concordance",
  180, 150,
  c(
    "results/tables/Table_5_matched_pathway_concordance.csv",
    "results/tables/Table_S_matched_concordance_summary.csv"
  )
)
write_figure_source(concordance_main, "Figure_7_pathway_concordance.csv")

# Supplementary Figure S1: PCA on pathway-wise z-standardised scores.
pca_scores <- data.table::fread(file.path(
  paths$analysis, "PCA_participant_scores.csv"
))
pca_variance <- data.table::fread(file.path(
  paths$tables, "Table_S_PCA_variance_explained.csv"
))
pca_axis_labels <- lapply(dataset_order, function(dataset_name) {
  x <- pca_variance[
    dataset == dataset_name & component %in% c("PC1", "PC2")
  ]
  setNames(sprintf("%s (%.1f%%)", x$component, x$variance_percent), x$component)
})
names(pca_axis_labels) <- dataset_order
pca_panels <- lapply(dataset_order, function(dataset_name) {
  x <- pca_scores[dataset == dataset_name]
  labels_local <- pca_axis_labels[[dataset_name]]
  ggplot(x, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#DDDDDD") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "#DDDDDD") +
    geom_point(size = 1.65, alpha = 0.72, colour = "#3C78A8") +
    labs(
      title = paste0(datasets[[dataset_name]]$label, " (n=", nrow(x), ")"),
      x = labels_local["PC1"], y = labels_local["PC2"]
    ) +
    theme(plot.title = element_text(size = 9, face = "bold"))
})
figure_s1 <- wrap_plots(pca_panels, ncol = 2) +
  plot_annotation(
    title = "Supplementary Figure S1. PCA of z-standardised ssGSEA scores",
    subtitle = paste0(
      "PCA was performed separately for each dataset; axes and loadings are dataset-specific\n",
      "and cannot be treated as a shared coordinate system."
    ),
    caption = paste0(
      "This is a descriptive two-dimensional projection. It does not test for, prove or exclude\n",
      "discrete latent classes/endotypes."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure_s1, paths$figures_supp, "Figure_S1_pathway_PCA",
  180, 150,
  c(
    "results/analysis_data/PCA_participant_scores.csv",
    "results/tables/Table_S_PCA_variance_explained.csv",
    "results/analysis_data/*_pathway_zscores_for_PCA.csv"
  )
)
write_figure_source(pca_scores, "Figure_S1_PCA_scores.csv")

# Supplementary Figure S2: normalized versus z-standardised scores side by side.
raw_vs_z_panels <- list()
raw_vs_z_sources <- character()
for (dataset_name in dataset_order) {
  raw_matrix <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")
  ))
  z_matrix <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_pathway_zscores.rds")
  ))
  raw_long <- make_long_matrix(raw_matrix, dataset_name, "Normalized ssGSEA")
  raw_row_order <- as.character(unique(raw_long[order(pathway_order), pathway]))
  raw_subject_order <- as.character(unique(
    raw_long[order(participant_order), Subject_ID]
  ))
  z_ordered <- z_matrix[raw_row_order, raw_subject_order, drop = FALSE]
  z_long <- data.table::as.data.table(as.table(z_ordered))
  data.table::setnames(z_long, c("pathway", "Subject_ID", "score"))
  z_long[, `:=`(
    pathway_label = unname(pathway_labels[as.character(pathway)]),
    participant_order = match(as.character(Subject_ID), colnames(z_ordered)),
    pathway_order = match(as.character(pathway), rownames(z_ordered))
  )]
  raw_panel <- make_heatmap_panel(
    raw_long, c(-1, 1),
    paste0(datasets[[dataset_name]]$label, " - normalized"),
    "Saved normalized ssGSEA; between-pathway magnitudes retained",
    TRUE, "Normalized ssGSEA", FALSE
  )
  z_panel <- make_heatmap_panel(
    z_long, c(-3, 3),
    paste0(datasets[[dataset_name]]$label, " - z-standardised"),
    "Each pathway centred to mean 0 and SD 1 within this dataset",
    FALSE, "Within-pathway z", FALSE
  )
  raw_vs_z_panels[[dataset_name]] <- raw_panel | z_panel
  raw_vs_z_sources <- c(
    raw_vs_z_sources,
    paste0("results/ssgsea/", dataset_name,
           "_ssgsea_normalized_scores.rds"),
    paste0("results/ssgsea/", dataset_name,
           "_ssgsea_pathway_zscores.rds")
  )
}
figure_s2 <- wrap_plots(raw_vs_z_panels, ncol = 1) +
  plot_annotation(
    title = "Supplementary Figure S2. Normalized and z-standardised ssGSEA heatmaps",
    subtitle = paste0(
      "Within each dataset, participant and pathway order is identical in the paired panels. Blue = lower;\n",
      "white = central; red = higher. Normalized limits are -1 to +1; within-pathway z-score limits are -3 to +3."
    ),
    caption = paste0(
      "Z-standardisation supports relative within-pathway patterning, PCA and standardised modelling, but it forces\n",
      "every pathway to mean 0 and SD 1, removing between-pathway differences in mean and variance."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure_s2, paths$figures_supp,
  "Figure_S2_normalized_vs_zstandardised_heatmaps",
  180, 300, raw_vs_z_sources
)

# Supplementary Figure S3: normalized-score between-participant variability.
variation <- data.table::fread(file.path(
  paths$tables, "Table_S_pathway_score_variation_SD_IQR.csv"
))
variation <- variation[score_scale == "normalized_ssgsea"]
variation[, dataset_factor := factor(dataset, levels = dataset_order)]
variation[, pathway_factor := factor(pathway, levels = rev(names(pathway_labels)))]
p_s3 <- ggplot(
  variation,
  aes(sd, pathway_factor, colour = dataset_factor, shape = dataset_factor)
) +
  geom_point(size = 1.65, alpha = 0.85,
             position = position_dodge(width = 0.5)) +
  scale_y_discrete(labels = pathway_labels) +
  scale_colour_manual(
    values = c("#0072B2", "#009E73", "#D55E00", "#CC79A7"),
    labels = dataset_labels_compact, name = "Dataset"
  ) +
  scale_shape_manual(
    values = c(16, 17, 15, 18), labels = dataset_labels_compact,
    name = "Dataset"
  ) +
  labs(
    title = "Supplementary Figure S3. Variation in normalized ssGSEA scores",
    subtitle = paste0(
      "SD is calculated before pathway-wise z-standardisation; larger values indicate greater\n",
      "participant-to-participant spread on the saved normalized ssGSEA scale."
    ),
    caption = paste0(
      "The plot preserves score-dispersion differences hidden after every pathway is rescaled to SD 1.\n",
      "It is descriptive, not a test of biological importance or cross-dataset calibration."
    ),
    x = "SD of normalized ssGSEA score across participants", y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 6.2),
    panel.grid.major.y = element_blank()
  ) + figure_title_theme
save_figure_four_formats(
  p_s3, paths$figures_supp,
  "Figure_S3_normalized_score_between_subject_variation",
  180, 180,
  "results/tables/Table_S_pathway_score_variation_SD_IQR.csv"
)
write_figure_source(variation, "Figure_S3_pathway_variation.csv")

# Supplementary Figure S4: covariate-adjusted curated sensitivity.
adjusted <- data.table::fread(file.path(
  paths$tables, "Table_S_covariate_adjusted_curated_coefficients.csv"
))
adjusted[, outcome_factor := factor(outcome, levels = outcome_order)]
adjusted[, pathway_factor := factor(
  pathway_label, levels = rev(unique(pathway_label))
)]
adjusted[, fdr_supported := !is.na(fdr_bh) & fdr_bh < 0.05]
p_s4 <- make_curated_forest_stack(
  adjusted, "Adjusted standardised beta (95% residual-df t CI)"
)
figure_s4 <- p_s4 +
  plot_annotation(
    title = "Supplementary Figure S4. Covariate-adjusted curated-model sensitivity",
    subtitle = paste0(
      "Models add age, sex, smoking status, BMI and current oral corticosteroid exposure\n",
      "to the same outcome-specific, prespecified pathway subsets. Filled red = BH-FDR q<0.05; open grey = q>=0.05."
    ),
    caption = paste0(
      "Adjusted analyses use covariate-complete cases and may differ from primary models because\n",
      "both covariate adjustment and sample composition changed."
    ),
    theme = figure_title_theme
  )
save_figure_four_formats(
  figure_s4, paths$figures_supp,
  "Figure_S4_covariate_adjusted_curated_models",
  180, 245,
  "results/tables/Table_S_covariate_adjusted_curated_coefficients.csv"
)
write_figure_source(adjusted, "Figure_S4_adjusted_coefficients.csv")

# Supplementary Figures S5-S8: complete calibration by dataset.
calibration_bins <- data.table::fread(file.path(
  paths$tables, "Table_S_predictive_calibration_bins.csv"
))
calibration_models <- c(
  "Biomarker-only", "Curated pathways", "Elastic Net pathways"
)
for (dataset_index in seq_along(dataset_order)) {
  dataset_name <- dataset_order[dataset_index]
  x <- calibration_bins[
    dataset == dataset_name & model_type %in% calibration_models
  ]
  x[, outcome_factor := factor(outcome, levels = outcome_order)]
  x[, model_factor := factor(model_type, levels = calibration_models)]
  p <- ggplot(
    x,
    aes(
      mean_predicted_model, mean_observed_model,
      colour = model_factor, shape = model_factor, group = model_factor
    )
  ) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed",
      linewidth = 0.4, colour = "#555555"
    ) +
    geom_errorbar(
      aes(
        ymin = mean_observed_model - se_observed_model,
        ymax = mean_observed_model + se_observed_model
      ),
      width = 0, linewidth = 0.45
    ) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 2) +
    facet_wrap(
      ~ outcome_factor, scales = "free", ncol = 3,
      labeller = as_labeller(outcome_labels)
    ) +
    scale_colour_manual(values = model_colours, name = "Model") +
    scale_shape_manual(values = c(16, 17, 15), name = "Model") +
    labs(
      title = paste0(
        "Supplementary Figure S", dataset_index + 4L,
        ". Binned calibration - ", dataset_labels_compact[dataset_name]
      ),
      subtitle = paste0(
        "Five equal-frequency bins of each participant's average of ten held-out predictions;\n",
        "bars are SE of observed bin means on the modelling scale."
      ),
      caption = paste0(
        "The dashed line is ideal agreement. Calibration is internal and conditional on fixed repeated-CV\n",
        "predictions; it does not establish external transportability."
      ),
      x = "Mean predicted outcome (model scale)",
      y = "Mean observed outcome (model scale)"
    ) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    figure_title_theme
  save_figure_four_formats(
    p, paths$figures_supp,
    paste0("Figure_S", dataset_index + 4L, "_calibration_", dataset_name),
    180, 150,
    c(
      "results/tables/Table_S_predictive_calibration_bins.csv",
      "results/analysis_data/nested_cv_averaged_predictions.csv.gz"
    )
  )
}
write_figure_source(
  calibration_bins, "Figures_S5_to_S8_calibration_bins.csv"
)

# Supplementary Figure S9: RNA-seq signature coverage.
coverage <- data.table::fread(file.path(
  paths$tables, "Table_S_signature_gene_coverage_by_dataset.csv"
))
coverage_rnaseq <- coverage[dataset == "brushing_rnaseq"]
coverage_rnaseq[, pathway_factor := factor(
  pathway, levels = coverage_rnaseq[order(detected_percent), pathway]
)]
coverage_rnaseq[, coverage_code := ifelse(low_coverage_flag, "LOW", "")]
p_s9 <- ggplot(
  coverage_rnaseq,
  aes(detected_percent, pathway_factor, fill = low_coverage_flag)
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = paste0(detected_gene_count, "/", harmonised_unique_genes)),
    hjust = 1.1, colour = "white", size = 2.35, fontface = "bold"
  ) +
  geom_text(
    aes(label = coverage_code), hjust = -0.15, colour = "black",
    size = 2.1, fontface = "bold"
  ) +
  scale_y_discrete(labels = pathway_labels) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%"), limits = c(0, 105),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    values = c(`FALSE` = "#3C78A8", `TRUE` = "#E69F00"),
    labels = c(`FALSE` = "Standard coverage", `TRUE` = "Low coverage flag"),
    name = NULL
  ) +
  labs(
    title = "Supplementary Figure S9. RNA-seq signature-gene coverage",
    subtitle = "Inside-bar labels are detected / signature genes; LOW is a redundant low-coverage flag",
    caption = paste0(
      "Low coverage was flagged rather than removing signatures after examining results.\n",
      "All retained scores met the ssGSEA minimum of three detected genes."
    ),
    x = "Detected signature genes (%)", y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid.major.y = element_blank()
  ) + figure_title_theme
save_figure_four_formats(
  p_s9, paths$figures_supp, "Figure_S9_rnaseq_signature_gene_coverage",
  180, 185,
  "results/tables/Table_S_signature_gene_coverage_by_dataset.csv"
)
write_figure_source(coverage_rnaseq, "Figure_S9_RNAseq_coverage.csv")

# Supplementary Figure S10: below-zero original-scale clipping sensitivity.
clipping <- data.table::fread(file.path(
  paths$tables, "Table_S_below_zero_prediction_clipping_sensitivity.csv"
))
clip_average <- clipping[prediction_level == "Averaged repeated OOF"]
clip_average[, dataset_factor := factor(dataset, levels = dataset_order)]
clip_average[, outcome_factor := factor(outcome, levels = outcome_order)]
clip_average[, model_factor := factor(model_type, levels = model_types)]
clip_average[, dataset_factor_compact := factor(
  dataset, levels = dataset_order, labels = dataset_labels_compact[dataset_order]
)]
clip_average[, outcome_factor_compact := factor(
  outcome, levels = outcome_order, labels = gsub("\n", " ", outcome_labels_compact[outcome_order])
)]
p_s10 <- ggplot(
  clip_average,
  aes(model_factor, below_zero_prediction_rows, fill = model_factor)
) +
  geom_col(width = 0.7) +
  facet_grid(
    dataset_factor_compact ~ outcome_factor_compact,
    scales = "free_y"
  ) +
  scale_fill_manual(values = model_colours, guide = "none") +
  scale_x_discrete(labels = c(
    "Biomarker-only" = "Bio",
    "Curated pathways" = "Cur",
    "Elastic Net pathways" = "EN",
    "Biomarker + curated pathways" = "B+Cur",
    "Biomarker + Elastic Net pathways" = "B+EN"
  )) +
  labs(
    title = "Supplementary Figure S10. Below-zero prediction sensitivity",
    subtitle = paste0(
      "Counts use participant-level averages of ten held-out predictions for non-negative\n",
      "outcomes; FEV1 is excluded from clipping."
    ),
    caption = paste0(
      "Unconstrained predictions remain primary. Zero-clipping is evaluated only for original-scale\n",
      "RMSE/MAE; model-scale R-squared is never clipped."
    ),
    x = NULL, y = "Below-zero averaged predictions"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    strip.text.x = element_text(size = 5.9),
    strip.text.y = element_text(size = 5.6, angle = 0)
  ) + figure_title_theme
save_figure_four_formats(
  p_s10, paths$figures_supp,
  "Figure_S10_below_zero_prediction_sensitivity",
  180, 160,
  "results/tables/Table_S_below_zero_prediction_clipping_sensitivity.csv"
)
write_figure_source(clip_average, "Figure_S10_clipping_counts.csv")

# Final figure export manifest and structural checks.
figure_manifest <- data.table::rbindlist(figure_manifest_rows, fill = TRUE)
data.table::fwrite(
  figure_manifest,
  file.path(paths$validation, "FIGURE_EXPORT_MANIFEST.csv")
)
expected_stems <- c(
  "Figure_1_cohort_structure_and_missingness",
  "Figure_2_normalized_ssgsea_pathway_heatmaps",
  "Figure_3_pathway_clinical_spearman_associations",
  "Figure_4_curated_pathway_model_forest_plots",
  "Figure_5_elastic_net_pathway_selection_stability",
  "Figure_6_internal_predictive_performance",
  "Figure_7_matched_pathway_concordance",
  "Figure_8_pathway_vs_biomarker_model_performance",
  "Figure_S1_pathway_PCA",
  "Figure_S2_normalized_vs_zstandardised_heatmaps",
  "Figure_S3_normalized_score_between_subject_variation",
  "Figure_S4_covariate_adjusted_curated_models",
  paste0("Figure_S", 5:8, "_calibration_", dataset_order),
  "Figure_S9_rnaseq_signature_gene_coverage",
  "Figure_S10_below_zero_prediction_sensitivity"
)
if (!setequal(unique(figure_manifest$figure_stem), expected_stems)) {
  stop("Figure stem set differs from the expected main/supplementary set.")
}
if (any(figure_manifest$output_bytes <= 0)) stop("A figure export is empty.")

writeLines(
  c(
    "R-ONLY FIGURE GENERATION: COMPUTATION PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Unique figure stems:", data.table::uniqueN(figure_manifest$figure_stem)),
    paste("Production exports:", nrow(figure_manifest)),
    "Every scientific figure was rendered directly in R from named source data.",
    "Every stem was exported as PDF, SVG, 300-dpi PNG and 600-dpi TIFF.",
    "All exported figures were generated directly by the R workflow without manual image editing.",
    "Figure 2 uses normalized ssGSEA values on fixed -1 to +1 display limits.",
    "Supplementary Figure S2 presents normalized and pathway-wise z-standardised scores side by side with identical ordering.",
    "PCA uses pathway-wise z-standardised scores and is descriptive only.",
    "Selection frequency is the fraction of 50 outer models with a non-zero coefficient and is not a p value.",
    "Predictive comparisons are exploratory internal validation and retain negative R-squared values."
  ),
  file.path(paths$validation, "FIGURE_GENERATION_STATUS.txt")
)
message("All dissertation figures generated directly in R.")
