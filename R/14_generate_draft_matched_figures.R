# Recreate the scientific figure types used in the corrected dissertation draft.
#
# IMPORTANT:
# - All scientific values, plots, panel assembly and exports are produced in R.
# - These figures use the validated analysis tables and score matrices.
# - No pre-existing image or externally edited graphic is reused as scientific
#   content.
# - The draft-matched set is separate from the base dissertation figure set.

source(file.path("R", "00_common.R"), local = TRUE)
assert_packages(c(
  "data.table", "ggplot2", "patchwork", "svglite", "ragg", "scales"
))
source(file.path("R", "12_figure_helpers.R"), local = TRUE)
message_rule("14: Generating draft-matched figures directly in R")

draft_main <- file.path(paths$root, "figures", "draft_matched", "main")
draft_supp <- file.path(
  paths$root, "figures", "draft_matched", "supplementary"
)
draft_sources <- file.path(
  paths$results, "figure_source_data", "draft_matched"
)
dir.create(draft_main, recursive = TRUE, showWarnings = FALSE)
dir.create(draft_supp, recursive = TRUE, showWarnings = FALSE)
dir.create(draft_sources, recursive = TRUE, showWarnings = FALSE)

# Retain the R plot objects so the separately numbered dissertation set can
# reuse the exact validated plotting constructions without reading or editing
# any exported image file. The dissertation exporter changes titles and output
# names in R, then renders fresh PDF/SVG/PNG/TIFF files.
draft_figure_plot_objects <- list()

draft_script <- "R/14_generate_draft_matched_figures.R"
draft_export <- function(plot, directory, stem, width_mm, height_mm, sources) {
  save_figure_four_formats(
    plot, directory, stem, width_mm, height_mm, sources,
    source_script = draft_script
  )
}
write_draft_source <- function(x, filename) {
  data.table::fwrite(x, file.path(draft_sources, filename))
}

dataset_letters <- setNames(c("A", "B", "C", "D"), dataset_order)
dataset_colours <- c(
  sputum = "#3C78A8",
  brushing_gpl570 = "#4D9221",
  biopsy = "#8C6BB1",
  brushing_rnaseq = "#E66101"
)
outcome_display <- c(
  FEV1 = "FEV1 % predicted",
  FeNO = "FeNO",
  blood_eosinophils = "Blood eosinophils",
  sputum_eosinophils = "Sputum eosinophils",
  exacerbation_frequency = "Prior-year exacerbations"
)

# ---------------------------------------------------------------------------
# Draft-matched Figure 1: sample availability and outcome missingness.
# Table 1 remains a table, rather than being relabelled as Figure 1C.
# ---------------------------------------------------------------------------
availability <- data.table::fread(file.path(
  paths$tables, "Table_S_sample_availability.csv"
))
missingness <- data.table::fread(file.path(
  paths$tables, "Table_S_outcome_missingness.csv"
))
availability[, dataset_factor := factor(dataset, levels = dataset_order)]
p1a <- ggplot(
  availability,
  aes(dataset_factor, transcriptomic_profiles, fill = dataset_factor)
) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = paste0("n=", transcriptomic_profiles)),
    vjust = -0.35, size = 3.0, fontface = "bold"
  ) +
  scale_fill_manual(values = dataset_colours, guide = "none") +
  scale_x_discrete(labels = c(
    sputum = "Sputum\nGPL570",
    brushing_gpl570 = "Brushing\nGPL570",
    biopsy = "Biopsy\nGPL570",
    brushing_rnaseq = "Brushing\nRNA-seq"
  )) +
  scale_y_continuous(
    limits = c(0, 165), breaks = seq(0, 160, 40),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    title = "A  Airway transcriptomic sample availability",
    subtitle = "Transcriptomic profiles define the maximum analysis n",
    x = NULL, y = "Transcriptomic profiles"
  ) +
  theme(panel.grid.major.x = element_blank())

missingness[, dataset_factor := factor(dataset, levels = dataset_order)]
missingness[, outcome_factor := factor(outcome, levels = rev(outcome_order))]
p1b <- ggplot(missingness, aes(missing_percent, outcome_factor)) +
  geom_col(width = 0.64, fill = "#9DB7D5") +
  geom_text(
    aes(label = paste0(available_n, "/", total_n, " available")),
    hjust = -0.05, size = 2.55
  ) +
  facet_wrap(
    ~ dataset_factor, ncol = 2,
    labeller = as_labeller(dataset_labels_compact)
  ) +
  scale_y_discrete(labels = outcome_display) +
  scale_x_continuous(
    limits = c(0, 75), breaks = seq(0, 60, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B  Clinical outcome missingness",
    subtitle = "Labels show complete outcome records / transcriptomic profiles",
    x = "Missing clinical outcome data", y = NULL
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 7.3, face = "bold"),
    axis.text.y = element_text(size = 6.5)
  )

draft_figure1 <- p1a / p1b +
  plot_layout(heights = c(0.72, 1.28)) +
  plot_annotation(
    title = "Draft-matched Figure 1. Sample and clinical-data availability",
    subtitle = paste0(
      "These panels define the analysed samples and explain outcome-specific differences in complete-case n.\n",
      "Validated baseline summaries remain in Table 1."
    ),
    theme = figure_title_theme
  )
draft_export(
  draft_figure1, draft_main,
  "DraftMatched_Figure_1_sample_availability_and_missingness",
  180, 205,
  c(
    "results/tables/Table_S_sample_availability.csv",
    "results/tables/Table_S_outcome_missingness.csv",
    "results/tables/Table_1_clinical_characteristics.csv"
  )
)
write_draft_source(availability, "DraftMatched_Figure_1A_sample_availability.csv")
write_draft_source(missingness, "DraftMatched_Figure_1B_missingness.csv")

# ---------------------------------------------------------------------------
# Draft-matched Figure 2A-D: one normalized ssGSEA heatmap per dataset.
# ---------------------------------------------------------------------------
for (dataset_name in dataset_order) {
  matrix <- readRDS(file.path(
    paths$scores, paste0(dataset_name, "_ssgsea_normalized_scores.rds")
  ))
  long <- make_long_matrix(
    matrix, dataset_name, "Normalized ssGSEA",
    cluster_rows = TRUE, cluster_columns = TRUE
  )
  coverage_note <- if (dataset_name == "brushing_rnaseq") {
    "Mast-cell and mucus signatures carry a low detected-gene-coverage flag"
  } else {
    "All 39 fixed signatures; hierarchical ordering is descriptive"
  }
  p <- make_heatmap_panel(
    long, c(-1, 1),
    paste0(
      "Draft-matched Figure 2", dataset_letters[dataset_name], ". ",
      datasets[[dataset_name]]$label, " pathway activity"
    ),
    paste0(
      datasets[[dataset_name]]$compartment, " | ",
      datasets[[dataset_name]]$platform, " | n=", ncol(matrix),
      " | ", coverage_note
    ),
    show_y_labels = TRUE,
    legend_title = "Normalized ssGSEA",
    show_legend = TRUE
  ) +
    labs(
      caption = paste0(
        "Saved normalized ssGSEA scores use fixed -1 to +1 limits; blue is lower, white central and red higher.\n",
        "The display does not test for, ",
        "prove or exclude discrete latent classes/endotypes."
      )
    ) +
    theme(plot.caption = element_text(size = 7, hjust = 0))
  stem <- paste0(
    "DraftMatched_Figure_2", dataset_letters[dataset_name], "_",
    dataset_name, "_normalized_ssgsea_heatmap"
  )
  draft_figure_plot_objects[[paste0(
    "Figure_2", dataset_letters[dataset_name]
  )]] <- p
  source_file <- paste0(
    "results/ssgsea/", dataset_name, "_ssgsea_normalized_scores.rds"
  )
  draft_export(p, draft_main, stem, 180, 135, source_file)
  write_draft_source(
    long,
    paste0("DraftMatched_Figure_2", dataset_letters[dataset_name], "_heatmap.csv")
  )
}

# ---------------------------------------------------------------------------
# Draft-matched Supplementary Figure S1: PCA.
# ---------------------------------------------------------------------------
pca_scores <- data.table::fread(file.path(
  paths$analysis, "PCA_participant_scores.csv"
))
pca_variance <- data.table::fread(file.path(
  paths$tables, "Table_S_PCA_variance_explained.csv"
))
pca_panels <- lapply(dataset_order, function(dataset_name) {
  x <- pca_scores[dataset == dataset_name]
  v <- pca_variance[
    dataset == dataset_name & component %in% c("PC1", "PC2")
  ]
  v <- setNames(v$variance_percent, v$component)
  ggplot(x, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "#DDDDDD", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "#DDDDDD", linewidth = 0.25) +
    geom_point(size = 1.7, alpha = 0.72, colour = dataset_colours[dataset_name]) +
    labs(
      title = paste0(dataset_labels_compact[dataset_name], " (n=", nrow(x), ")"),
      x = sprintf("PC1 (%.1f%%)", v["PC1"]),
      y = sprintf("PC2 (%.1f%%)", v["PC2"])
    ) +
    theme(plot.title = element_text(size = 9, face = "bold"))
})
draft_s1 <- wrap_plots(pca_panels, ncol = 2) +
  plot_annotation(
    title = "Draft-matched Supplementary Figure S1. PCA of pathway z-scores",
    subtitle = paste0(
      "PCA was performed separately in each dataset using pathway-wise z-standardised ssGSEA scores.\n",
      "Axes are dataset-specific; this descriptive display does not test latent classes."
    ),
    theme = figure_title_theme
  )
draft_figure_plot_objects[["Figure_S1"]] <- draft_s1
draft_export(
  draft_s1, draft_supp, "DraftMatched_Figure_S1_pathway_PCA",
  180, 150,
  c(
    "results/analysis_data/PCA_participant_scores.csv",
    "results/tables/Table_S_PCA_variance_explained.csv",
    "results/analysis_data/*_pathway_zscores_for_PCA.csv"
  )
)
write_draft_source(pca_scores, "DraftMatched_Figure_S1_PCA_scores.csv")

# ---------------------------------------------------------------------------
# Draft-matched Figure 3A-D: one 39 x 5 association heatmap per dataset.
# ---------------------------------------------------------------------------
associations <- data.table::fread(file.path(
  paths$tables, "Table_2_pathway_clinical_spearman_associations.csv"
))
associations[, significance_marker := data.table::fcase(
  fdr_bh < 0.001, "***",
  fdr_bh < 0.01, "**",
  fdr_bh < 0.05, "*",
  default = ""
)]
associations[, direction_marker := paste0(
  ifelse(spearman_rho >= 0, "+", "-"), significance_marker
)]
associations[, marker_colour := ifelse(abs(spearman_rho) >= 0.55, "white", "black")]

for (dataset_name in dataset_order) {
  x <- data.table::copy(associations[dataset == dataset_name])
  matrix <- data.table::dcast(
    x, pathway ~ outcome, value.var = "spearman_rho"
  )
  matrix_values <- as.matrix(matrix[, ..outcome_order])
  rownames(matrix_values) <- matrix$pathway
  pathway_order_local <- rownames(matrix_values)[
    stats::hclust(stats::dist(matrix_values))$order
  ]
  x[, pathway_factor := factor(pathway, levels = rev(pathway_order_local))]
  x[, outcome_factor := factor(outcome, levels = outcome_order)]
  n_axis <- unique(x[, .(outcome, n)])
  n_axis[, display := paste0(outcome_display[outcome], "\n(n=", n, ")")]
  n_labels <- setNames(n_axis$display, n_axis$outcome)
  p <- ggplot(x, aes(outcome_factor, pathway_factor, fill = spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = direction_marker, colour = marker_colour),
      size = 2.1, fontface = "bold"
    ) +
    scale_colour_identity() +
    scale_fill_gradient2(
      low = association_colours["low"], mid = association_colours["mid"],
      high = association_colours["high"], midpoint = 0, limits = c(-1, 1),
      name = "Spearman rho"
    ) +
    scale_x_discrete(labels = n_labels) +
    scale_y_discrete(labels = pathway_labels, drop = FALSE) +
    labs(
      title = paste0(
        "Figure 3", dataset_letters[dataset_name], ". ",
        dataset_labels_compact[dataset_name], " Spearman associations"
      ),
      subtitle = paste0(
        "+/- preserves direction in grayscale; * q<0.05, ** q<0.01, *** q<0.001 ",
        "after BH correction within each outcome"
      ),
      caption = paste0(
        "All 39 pathways are shown. These are marginal univariable associations,\n",
        "not conditional effects, causal estimates or predictive performance."
      ),
      x = "Clinical outcome", y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 6.8, angle = 15, hjust = 1),
      axis.text.y = element_text(size = 6.1),
      legend.position = "right"
    ) + figure_title_theme
  stem <- paste0(
    "DraftMatched_Figure_3", dataset_letters[dataset_name], "_",
    dataset_name, "_spearman_heatmap"
  )
  draft_figure_plot_objects[[paste0(
    "Figure_3", dataset_letters[dataset_name]
  )]] <- p
  draft_export(
    p, draft_main, stem, 180, 178,
    "results/tables/Table_2_pathway_clinical_spearman_associations.csv"
  )
  write_draft_source(
    x,
    paste0("DraftMatched_Figure_3", dataset_letters[dataset_name], "_spearman.csv")
  )
}

# ---------------------------------------------------------------------------
# Draft-matched Figure 4A-D: separate curated forests by dataset.
# ---------------------------------------------------------------------------
curated <- data.table::fread(file.path(
  paths$tables, "Table_3_curated_model_coefficients_95CI.csv"
))
curated[, fdr_supported := !is.na(fdr_bh) & fdr_bh < 0.05]

make_draft_curated_dataset <- function(dataset_name) {
  dataset_table <- data.table::copy(curated[dataset == dataset_name])
  panels <- lapply(seq_along(outcome_order), function(i) {
    outcome_name <- outcome_order[i]
    x <- data.table::copy(dataset_table[outcome == outcome_name])
    x[, pathway_factor := factor(
      pathway_label, levels = rev(pathway_label)
    )]
    ggplot(x, aes(standardized_beta, pathway_factor)) +
      geom_vline(
        xintercept = 0, linetype = "dashed", linewidth = 0.35,
        colour = "#666666"
      ) +
      geom_errorbar(
        aes(
          xmin = ci_low, xmax = ci_high,
          colour = fdr_supported
        ),
        orientation = "y", width = 0, linewidth = 0.62
      ) +
      geom_point(
        aes(colour = fdr_supported, shape = fdr_supported),
        size = 2.1, stroke = 0.75, fill = "white"
      ) +
      scale_colour_manual(
        values = c(`TRUE` = "#B2182B", `FALSE` = "#555555"),
        labels = c(`TRUE` = "BH-FDR q<0.05", `FALSE` = "q>=0.05"),
        name = NULL,
        guide = "none"
      ) +
      scale_shape_manual(
        values = c(`TRUE` = 16, `FALSE` = 21),
        name = NULL, guide = "none"
      ) +
      coord_cartesian(xlim = c(-1.1, 1.25), clip = "off") +
      labs(
        title = paste0(outcome_display[outcome_name], " (n=", unique(x$n), ")"),
        x = if (i >= 4L) "Standardised beta (95% residual-df t CI)" else NULL,
        y = NULL
      ) +
      theme(
        plot.title = element_text(size = 8.2, face = "bold"),
        axis.text.y = element_text(size = 6.7),
        axis.text.x = element_text(size = 6.5),
        panel.grid.major.y = element_blank(),
        legend.position = "none"
      )
  })
  panels[[6]] <- plot_spacer()
  wrap_plots(panels, ncol = 2, guides = "collect") +
    plot_annotation(
      title = paste0(
        "Draft-matched Figure 4", dataset_letters[dataset_name], ". ",
        datasets[[dataset_name]]$label, " curated pathway models"
      ),
      subtitle = paste0(
        "Only biologically prespecified pathway subsets are shown; these are not top ",
        "pathways selected from all 39 signatures.\n",
        "Filled red = BH-FDR q<0.05; open grey = q>=0.05."
      ),
      caption = paste0(
        "Coefficients are conditional on the other displayed pathways. Correlated signatures can cause\n",
        "collinearity or suppression; a negative IL-5 coefficient ",
        "does not indicate biological protection."
      ),
      theme = figure_title_theme
    ) & theme(legend.position = "none")
}

for (dataset_name in dataset_order) {
  p <- make_draft_curated_dataset(dataset_name)
  draft_figure_plot_objects[[paste0(
    "Figure_4", dataset_letters[dataset_name]
  )]] <- p
  stem <- paste0(
    "DraftMatched_Figure_4", dataset_letters[dataset_name], "_",
    dataset_name, "_curated_model_forests"
  )
  draft_export(
    p, draft_main, stem, 180, 205,
    c(
      "results/tables/Table_3_curated_model_coefficients_95CI.csv",
      "results/tables/Table_S_curated_model_formulas.csv"
    )
  )
  write_draft_source(
    curated[dataset == dataset_name],
    paste0("DraftMatched_Figure_4", dataset_letters[dataset_name], "_coefficients.csv")
  )
}

# ---------------------------------------------------------------------------
# Draft-matched Figure 5A-D: top-ten Elastic Net stability dot plots.
# The old 'final coefficient' claim is replaced with predominant outer-fit sign.
# ---------------------------------------------------------------------------
stability <- data.table::fread(file.path(
  paths$tables, "Table_5_elastic_net_selection_stability.csv"
))[model_type == "Elastic Net pathways"]
performance <- data.table::fread(file.path(
  paths$tables, "Table_4_repeated_nested_CV_model_performance.csv"
))
predictive_n <- unique(performance[
  model_type == "Elastic Net pathways", .(dataset, outcome, n)
])
stability <- merge(
  stability, predictive_n, by = c("dataset", "outcome"), all.x = TRUE
)
stability[, coefficient_weight := data.table::fifelse(
  is.finite(median_nonzero_coefficient), abs(median_nonzero_coefficient), 0
)]
stability[, coefficient_share_percent := {
  denominator <- sum(coefficient_weight)
  if (denominator > 0) 100 * coefficient_weight / denominator else 0
}, by = .(dataset, outcome)]
data.table::setorder(
  stability, dataset, outcome, -selection_frequency,
  -coefficient_weight, pathway
)
stability[, stability_rank := seq_len(.N), by = .(dataset, outcome)]
stability_top <- stability[stability_rank <= 10L]
stability_top[, direction_display := data.table::fcase(
  predominant_sign == "Positive", "Predominantly positive",
  predominant_sign == "Negative", "Predominantly negative",
  predominant_sign == "Mixed", "Mixed sign",
  default = "Never selected"
)]
direction_shapes <- c(
  "Predominantly positive" = 24,
  "Predominantly negative" = 25,
  "Mixed sign" = 23,
  "Never selected" = 21
)
direction_fills <- c(
  "Predominantly positive" = "#B2182B",
  "Predominantly negative" = "#2166AC",
  "Mixed sign" = "#762A83",
  "Never selected" = "white"
)

make_stability_panel <- function(dataset_name, outcome_name) {
  x <- data.table::copy(stability_top[
    dataset == dataset_name & outcome == outcome_name
  ])
  data.table::setorder(x, selection_frequency, coefficient_weight)
  x[, pathway_factor := factor(pathway_label, levels = pathway_label)]
  ggplot(
    x,
    aes(
      selection_frequency, pathway_factor,
      size = coefficient_share_percent,
      shape = direction_display,
      fill = direction_display
    )
  ) +
    geom_point(colour = "black", stroke = 0.45, alpha = 0.9) +
    scale_x_continuous(
      limits = c(0, 1.04), breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_shape_manual(values = direction_shapes, name = "Outer-fit sign") +
    scale_fill_manual(values = direction_fills, name = "Outer-fit sign") +
    scale_size_continuous(
      range = c(1.8, 6.2), limits = c(0, 40), breaks = c(5, 15, 30),
      name = "Median-coefficient\nshare (%)"
    ) +
    labs(
      title = paste0(outcome_display[outcome_name], " (n=", unique(x$n), ")"),
      x = "Selection frequency across 50 outer models", y = NULL
    ) +
    theme(
      plot.title = element_text(size = 8.3, face = "bold"),
      axis.text.y = element_text(size = 6.25),
      axis.text.x = element_text(size = 6.2),
      axis.title.x = element_text(size = 7.2),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    ) +
    guides(shape = "none", fill = "none", size = "none")
}

for (dataset_name in dataset_order) {
  panels <- lapply(outcome_order, function(outcome_name) {
    make_stability_panel(dataset_name, outcome_name)
  })
  panels[[6]] <- plot_spacer()
  p <- wrap_plots(panels, ncol = 2, guides = "collect") +
    plot_annotation(
      title = paste0(
        "Draft-matched Figure 5", dataset_letters[dataset_name], ". ",
        datasets[[dataset_name]]$label, "\nElastic Net pathway stability"
      ),
      subtitle = paste0(
        "Top 10 pathways per outcome. Selection frequency is the fraction of 50 ",
        "outer models with a non-zero coefficient; it is not a p value.\n",
        "Up triangle = positive, down triangle = negative, open circle = never selected; point size = coefficient share."
      ),
      caption = paste0(
        "Shape and fill show the predominant selected-coefficient sign across outer fits.\n",
        "Point size is a coefficient-based share derived from median non-zero coefficients; ",
        "it is not variance decomposition or causal importance."
      ),
      theme = figure_title_theme
    ) & theme(legend.position = "none")
  draft_figure_plot_objects[[paste0(
    "Figure_5", dataset_letters[dataset_name]
  )]] <- p
  stem <- paste0(
    "DraftMatched_Figure_5", dataset_letters[dataset_name], "_",
    dataset_name, "_top10_elastic_net_stability"
  )
  draft_export(
    p, draft_main, stem, 180, 240,
    c(
      "results/tables/Table_5_elastic_net_selection_stability.csv",
      "results/tables/Table_4_repeated_nested_CV_model_performance.csv"
    )
  )
  write_draft_source(
    stability_top[dataset == dataset_name],
    paste0("DraftMatched_Figure_5", dataset_letters[dataset_name], "_top10_stability.csv")
  )
}

# ---------------------------------------------------------------------------
# Draft-matched Figure 6A-D: the predictive-performance forms in the draft.
# ---------------------------------------------------------------------------
standalone_models <- c(
  "Biomarker-only", "Curated pathways", "Elastic Net pathways"
)
pathway_models <- c("Curated pathways", "Elastic Net pathways")
performance[, dataset_factor := factor(
  dataset, levels = dataset_order, labels = dataset_labels_compact[dataset_order]
)]
performance[, outcome_factor := factor(outcome, levels = rev(outcome_order))]
performance[, model_factor := factor(model_type, levels = standalone_models)]

make_metric_plot <- function(x, mean_column, sd_column, title, x_label,
                             x_limits, models) {
  plot_data <- data.table::copy(x[model_type %in% models])
  plot_data[, metric_mean := get(mean_column)]
  plot_data[, metric_sd := get(sd_column)]
  plot_data[, outcome_n := paste0(outcome_display[outcome], "\n(n=", n, ")")]
  outcome_n_labels <- unique(plot_data[, .(outcome, outcome_n)])
  outcome_n_labels <- setNames(outcome_n_labels$outcome_n, outcome_n_labels$outcome)
  plot_data[, outcome_factor := factor(outcome, levels = rev(outcome_order))]
  ggplot(
    plot_data,
    aes(metric_mean, outcome_factor, colour = model_type, shape = model_type)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = metric_mean - metric_sd, xmax = metric_mean + metric_sd),
      orientation = "y", width = 0.12, linewidth = 0.48,
      position = position_dodge(width = 0.38)
    ) +
    geom_point(size = 1.9, position = position_dodge(width = 0.38)) +
    facet_wrap(~ dataset_factor, nrow = 1) +
    scale_y_discrete(labels = outcome_n_labels) +
    scale_colour_manual(values = model_colours, name = NULL) +
    scale_shape_manual(
      values = c(
        "Biomarker-only" = 15,
        "Curated pathways" = 16,
        "Elastic Net pathways" = 17
      ),
      name = NULL
    ) +
    coord_cartesian(xlim = x_limits) +
    labs(title = title, x = x_label, y = NULL) +
    theme(
      strip.text = element_text(size = 7.2, face = "bold"),
      axis.text.y = element_text(size = 5.8),
      axis.text.x = element_text(size = 6.2),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom"
    )
}

p6a_r2 <- make_metric_plot(
  performance,
  "mean_cv_r_squared", "sd_cv_r_squared",
  "Cross-validated R-squared", "Mean R-squared +/- repeat SD",
  c(-0.48, 0.92), pathway_models
)
p6a_r <- make_metric_plot(
  performance,
  "mean_cv_pearson_r", "sd_cv_pearson_r",
  "Pearson correlation", "Mean Pearson r +/- repeat SD",
  c(-0.35, 1.02), pathway_models
)
p6a <- p6a_r2 / p6a_r +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Draft-matched Figure 6A. Curated versus Elastic Net performance",
    subtitle = paste0(
      "Means +/- one SD across ten repeated outer validations on identical common cases and folds.\n",
      "Finite negative R-squared values are retained; Pearson r is association, not calibration."
    ),
    theme = figure_title_theme
  ) & theme(legend.position = "bottom")
draft_export(
  p6a, draft_main,
  "DraftMatched_Figure_6A_curated_vs_elastic_net_performance",
  180, 210,
  "results/tables/Table_4_repeated_nested_CV_model_performance.csv"
)
write_draft_source(
  performance[model_type %in% pathway_models],
  "DraftMatched_Figure_6A_curated_vs_elastic_net.csv"
)

p6b <- make_metric_plot(
  performance,
  "mean_cv_r_squared", "sd_cv_r_squared",
  "Biomarker-only versus pathway-model R-squared",
  "Mean cross-validated R-squared +/- repeat SD",
  c(-0.48, 0.92), standalone_models
) +
  plot_annotation(
    title = "Draft-matched Figure 6B. Biomarker-only versus pathway models",
    subtitle = paste0(
      "All three standalone models use identical common cases and outer folds within each\n",
      "dataset-outcome comparison; performance is outcome- and dataset-specific."
    ),
    caption = paste0(
      "This is exploratory internal validation. Negative R-squared indicates performance\n",
      "worse than the held-out mean benchmark and is not truncated to zero."
    ),
    theme = figure_title_theme
  )
draft_export(
  p6b, draft_main,
  "DraftMatched_Figure_6B_biomarker_vs_pathway_performance",
  180, 150,
  "results/tables/Table_4_repeated_nested_CV_model_performance.csv"
)
write_draft_source(
  performance[model_type %in% standalone_models],
  "DraftMatched_Figure_6B_three_model_performance.csv"
)

sputum_performance <- performance[
  dataset == "sputum" & model_type %in% standalone_models
]
p6c_r2 <- make_metric_plot(
  sputum_performance,
  "mean_cv_r_squared", "sd_cv_r_squared",
  "Sputum R-squared", "Mean R-squared +/- repeat SD",
  c(-0.15, 0.92), standalone_models
) + theme(strip.text = element_blank())
p6c_r <- make_metric_plot(
  sputum_performance,
  "mean_cv_pearson_r", "sd_cv_pearson_r",
  "Sputum Pearson correlation", "Mean Pearson r +/- repeat SD",
  c(-0.15, 1.02), standalone_models
) + theme(strip.text = element_blank())
p6c <- (p6c_r2 | p6c_r) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Draft-matched Figure 6C. Detailed sputum model performance",
    subtitle = paste0(
      "Sputum is shown in detail because concurrent sputum-eosinophil estimation was the strongest\n",
      "evaluated local phenotype result; it is same-sample phenotype capture, not prognosis."
    ),
    theme = figure_title_theme
  ) & theme(legend.position = "bottom")
draft_export(
  p6c, draft_main,
  "DraftMatched_Figure_6C_sputum_detailed_performance",
  180, 145,
  "results/tables/Table_4_repeated_nested_CV_model_performance.csv"
)
write_draft_source(
  sputum_performance,
  "DraftMatched_Figure_6C_sputum_performance.csv"
)

calibration_bins <- data.table::fread(file.path(
  paths$tables, "Table_S_predictive_calibration_bins.csv"
))
calibration_models <- standalone_models
make_calibration_plot <- function(x, title, subtitle = NULL) {
  x <- data.table::copy(x[model_type %in% calibration_models])
  x[, model_factor := factor(model_type, levels = calibration_models)]
  ggplot(
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
    geom_line(linewidth = 0.5) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = model_colours, name = "Model") +
    scale_shape_manual(values = c(16, 17, 15), name = "Model") +
    labs(
      title = title, subtitle = subtitle,
      x = "Mean held-out prediction within bin (model scale)",
      y = "Mean observed outcome within bin (model scale)"
    ) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

sputum_eos_calibration <- calibration_bins[
  dataset == "sputum" & outcome == "sputum_eosinophils"
]
p6d <- make_calibration_plot(
  sputum_eos_calibration,
  "Draft-matched Figure 6D. Representative sputum-eosinophil calibration",
  "Five equal-frequency bins of each participant's average of ten held-out predictions"
) +
  labs(
    caption = paste0(
      "Bars are SE of observed bin means; the diagonal is ideal agreement.\n",
      "internal calibration conditional on fixed repeated-CV predictions, not external validation."
    )
  ) + figure_title_theme
draft_export(
  p6d, draft_main,
  "DraftMatched_Figure_6D_sputum_eosinophil_calibration",
  180, 125,
  c(
    "results/tables/Table_S_predictive_calibration_bins.csv",
    "results/analysis_data/nested_cv_averaged_predictions.csv.gz"
  )
)
write_draft_source(
  sputum_eos_calibration,
  "DraftMatched_Figure_6D_sputum_eosinophil_calibration.csv"
)

# ---------------------------------------------------------------------------
# Draft-matched Figure 7A-B: concordance heatmap and brushing-platform bars.
# ---------------------------------------------------------------------------
concordance <- data.table::fread(file.path(
  paths$tables, "Table_5_matched_pathway_concordance.csv"
))
cross_compartment <- c(
  "brushing_gpl570_vs_biopsy",
  "sputum_vs_biopsy",
  "sputum_vs_brushing_gpl570"
)
comparison_labels <- c(
  brushing_gpl570_vs_biopsy = "Brushing GPL570 vs biopsy\n(n=99)",
  sputum_vs_biopsy = "Sputum vs biopsy\n(n=23)",
  sputum_vs_brushing_gpl570 = "Sputum vs brushing GPL570\n(n=31)"
)
x7a <- data.table::copy(concordance[comparison %in% cross_compartment])
x7a[, comparison_factor := factor(comparison, levels = cross_compartment)]
matrix7 <- data.table::dcast(
  x7a, pathway ~ comparison, value.var = "spearman_rho"
)
matrix7_values <- as.matrix(matrix7[, ..cross_compartment])
rownames(matrix7_values) <- matrix7$pathway
pathway_order7 <- rownames(matrix7_values)[
  stats::hclust(stats::dist(matrix7_values))$order
]
x7a[, pathway_factor := factor(pathway, levels = rev(pathway_order7))]
x7a[, direction_marker := paste0(
  ifelse(spearman_rho >= 0, "+", "-"),
  ifelse(!is.na(fdr_bh) & fdr_bh < 0.05, "*", "")
)]
x7a[, marker_colour := ifelse(abs(spearman_rho) >= 0.55, "white", "black")]
p7a <- ggplot(x7a, aes(comparison_factor, pathway_factor, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.28) +
  geom_text(
    aes(label = direction_marker, colour = marker_colour),
    size = 2.1, fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = association_colours["low"], mid = association_colours["mid"],
    high = association_colours["high"], midpoint = 0, limits = c(-1, 1),
    name = "Spearman rho"
  ) +
  scale_x_discrete(labels = comparison_labels) +
  scale_y_discrete(labels = pathway_labels, drop = FALSE) +
  labs(
    title = paste0(
      "Draft-matched Figure 7A. Same-platform\n",
      "cross-compartment concordance"
    ),
    subtitle = paste0(
      "All comparisons use GPL570; +/- preserves direction in grayscale and * marks BH-FDR q<0.05.\n",
      "RNA-seq is excluded to hold platform constant."
    ),
    caption = paste0(
      "Each cell is a same-pathway Spearman correlation across matched participants.\n",
      "These are concordance analyses, not independent replication."
    ),
    x = NULL, y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7.2, angle = 10, hjust = 1),
    axis.text.y = element_text(size = 6.0),
    legend.position = "right"
  ) + figure_title_theme
draft_figure_plot_objects[["Figure_7A"]] <- p7a
draft_export(
  p7a, draft_main,
  "DraftMatched_Figure_7A_cross_compartment_concordance_heatmap",
  180, 175,
  "results/tables/Table_5_matched_pathway_concordance.csv"
)
write_draft_source(x7a, "DraftMatched_Figure_7A_concordance.csv")

x7b <- data.table::copy(concordance[
  comparison == "brushing_gpl570_vs_brushing_rnaseq"
])
data.table::setorder(x7b, spearman_rho)
x7b[, pathway_factor := factor(pathway_label, levels = pathway_label)]
x7b[, coverage_display := ifelse(
  rnaseq_low_coverage_flag, "Low coverage", "Standard coverage"
)]
x7b[, coverage_code := ifelse(rnaseq_low_coverage_flag, "LOW", "")]
p7b <- ggplot(x7b, aes(spearman_rho, pathway_factor, fill = coverage_display)) +
  geom_col(width = 0.68) +
  geom_text(
    aes(label = coverage_code), hjust = -0.12,
    size = 2.2, fontface = "bold", colour = "black"
  ) +
  scale_fill_manual(
    values = c("Standard coverage" = "#3C78A8", "Low coverage" = "#E69F00"),
    name = "RNA-seq signature coverage"
  ) +
  scale_x_continuous(
    limits = c(0, 1.03), breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = paste0(
      "Draft-matched Figure 7B. Brushing GPL570 versus RNA-seq\n",
      "same-participant concordance"
    ),
    subtitle = paste0(
      "Same compartment and 118 matched participants; LOW is a redundant\n",
      "detected-gene-coverage flag that remains interpretable in grayscale."
    ),
    caption = paste0(
      "This is same-participant cross-platform concordance, not an external cohort or independent validation.\n",
      "Coverage informs pathway-specific caution but does not alone determine agreement."
    ),
    x = "Spearman rho", y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 6.2),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom"
  ) + figure_title_theme
draft_figure_plot_objects[["Figure_7B"]] <- p7b
draft_export(
  p7b, draft_main,
  "DraftMatched_Figure_7B_brushing_cross_platform_concordance",
  180, 185,
  c(
    "results/tables/Table_5_matched_pathway_concordance.csv",
    "results/tables/Table_S_signature_gene_coverage_by_dataset.csv"
  )
)
write_draft_source(x7b, "DraftMatched_Figure_7B_concordance.csv")

# ---------------------------------------------------------------------------
# Draft-matched Supplementary Figures S2-S5: complete binned calibration.
# ---------------------------------------------------------------------------
for (dataset_index in seq_along(dataset_order)) {
  dataset_name <- dataset_order[dataset_index]
  x <- calibration_bins[
    dataset == dataset_name & model_type %in% calibration_models
  ]
  x[, outcome_factor := factor(outcome, levels = outcome_order)]
  p <- make_calibration_plot(
    x,
    paste0(
      "Draft-matched Supplementary Figure S", dataset_index + 1L,
      ". Calibration - ", dataset_labels_compact[dataset_name]
    ),
    paste0(
      "Each participant's ten held-out predictions were averaged before five equal-frequency bins were formed;\n",
      "bars are SE of observed bin means."
    )
  ) +
    facet_wrap(
      ~ outcome_factor, scales = "free", ncol = 3,
      labeller = as_labeller(outcome_display)
    ) +
    labs(
      caption = paste0(
        "The dashed line is ideal agreement. Calibration is internal and conditional\n",
        "on fixed repeated-CV predictions; it does not establish external transportability."
      )
    ) + figure_title_theme
  draft_figure_plot_objects[[paste0(
    "Figure_S", dataset_index + 1L
  )]] <- p
  stem <- paste0(
    "DraftMatched_Figure_S", dataset_index + 1L, "_calibration_",
    dataset_name
  )
  draft_export(
    p, draft_supp, stem, 180, 150,
    c(
      "results/tables/Table_S_predictive_calibration_bins.csv",
      "results/analysis_data/nested_cv_averaged_predictions.csv.gz"
    )
  )
  write_draft_source(
    x,
    paste0("DraftMatched_Figure_S", dataset_index + 1L, "_calibration.csv")
  )
}

# ---------------------------------------------------------------------------
# Manifest and structural checks for the separate draft-matched set.
# ---------------------------------------------------------------------------
draft_manifest <- data.table::rbindlist(figure_manifest_rows, fill = TRUE)
expected_stems <- c(
  "DraftMatched_Figure_1_sample_availability_and_missingness",
  paste0(
    "DraftMatched_Figure_2", dataset_letters, "_", dataset_order,
    "_normalized_ssgsea_heatmap"
  ),
  paste0(
    "DraftMatched_Figure_3", dataset_letters, "_", dataset_order,
    "_spearman_heatmap"
  ),
  paste0(
    "DraftMatched_Figure_4", dataset_letters, "_", dataset_order,
    "_curated_model_forests"
  ),
  paste0(
    "DraftMatched_Figure_5", dataset_letters, "_", dataset_order,
    "_top10_elastic_net_stability"
  ),
  "DraftMatched_Figure_6A_curated_vs_elastic_net_performance",
  "DraftMatched_Figure_6B_biomarker_vs_pathway_performance",
  "DraftMatched_Figure_6C_sputum_detailed_performance",
  "DraftMatched_Figure_6D_sputum_eosinophil_calibration",
  "DraftMatched_Figure_7A_cross_compartment_concordance_heatmap",
  "DraftMatched_Figure_7B_brushing_cross_platform_concordance",
  "DraftMatched_Figure_S1_pathway_PCA",
  paste0(
    "DraftMatched_Figure_S", 2:5, "_calibration_", dataset_order
  )
)
if (!setequal(unique(draft_manifest$figure_stem), expected_stems)) {
  missing <- setdiff(expected_stems, unique(draft_manifest$figure_stem))
  unexpected <- setdiff(unique(draft_manifest$figure_stem), expected_stems)
  stop(
    "Draft-matched stem mismatch. Missing: ", paste(missing, collapse = ", "),
    "; unexpected: ", paste(unexpected, collapse = ", ")
  )
}
if (nrow(draft_manifest) != length(expected_stems) * 4L) {
  stop("Draft-matched figure manifest does not contain four formats per stem.")
}
if (any(draft_manifest$output_bytes <= 0)) {
  stop("At least one draft-matched figure export is empty.")
}
data.table::fwrite(
  draft_manifest,
  file.path(paths$validation, "DRAFT_MATCHED_FIGURE_EXPORT_MANIFEST.csv")
)
writeLines(
  c(
    "DRAFT-MATCHED R-ONLY FIGURE GENERATION: COMPUTATION PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Unique draft-matched figure stems:", length(expected_stems)),
    paste("Production exports:", nrow(draft_manifest)),
    "Every scientific value and plot was generated directly in R from the current rerun outputs.",
    "No pre-existing or manually edited scientific graphic was reused.",
    "The separate set recreates the draft's figure types using the current validated results and final labels.",
    "Table 1 remains a table and is not relabelled as Figure 1C."
  ),
  file.path(paths$validation, "DRAFT_MATCHED_FIGURE_GENERATION_STATUS.txt")
)

message(
  "Generated ", length(expected_stems),
  " draft-matched figure stems in four R-produced formats."
)
