# Shared plotting and export helpers for the dissertation figure set.
# Every production export is created directly by R from validated source data.

library(ggplot2)
library(patchwork)

theme_set(theme_minimal(base_family = "Helvetica", base_size = 9))

dataset_order <- canonical_dataset_order
dataset_labels <- vapply(datasets, `[[`, character(1), "label")
names(dataset_labels) <- dataset_order
outcome_order <- names(outcome_specs)
outcome_labels <- vapply(outcome_specs, `[[`, character(1), "short_label")
names(outcome_labels) <- outcome_order
outcome_labels_compact <- c(
  FEV1 = "FEV1",
  FeNO = "FeNO",
  blood_eosinophils = "Blood eos",
  sputum_eosinophils = "Sputum eos",
  exacerbation_frequency = "Prior-year\nexacerbations"
)
dataset_labels_compact <- c(
  sputum = "Sputum GPL570",
  brushing_gpl570 = "Brushing GPL570",
  biopsy = "Biopsy GPL570",
  brushing_rnaseq = "Brushing RNA-seq"
)

model_colours <- c(
  "Biomarker-only" = "#555555",
  "Curated pathways" = "#0072B2",
  "Elastic Net pathways" = "#D55E00",
  "Biomarker + curated pathways" = "#009E73",
  "Biomarker + Elastic Net pathways" = "#CC79A7"
)
sign_colours <- c("Negative" = "#2166AC", "Positive" = "#B2182B")
association_colours <- c(
  "low" = "#2166AC", "mid" = "#F7F7F7", "high" = "#B2182B"
)
# Preserve the blue-white-red diverging palette requested by the user and used
# in the corrected dissertation draft. Blue is lower/negative, white is the
# midpoint and red is higher/positive.
score_heatmap_colours <- association_colours

figure_manifest_rows <- list()

figure_title_theme <- theme(
  plot.title = element_text(size = 13, face = "bold", margin = margin(b = 5)),
  plot.subtitle = element_text(
    size = 8, colour = "#333333", lineheight = 1.08
  ),
  plot.caption = element_text(
    size = 7.5, colour = "#444444", hjust = 0,
    lineheight = 1.1, margin = margin(t = 6)
  ),
  strip.text = element_text(size = 8.5, face = "bold"),
  panel.grid.minor = element_blank(),
  legend.position = "bottom",
  legend.title = element_text(face = "bold")
)

wrap_dataset_label <- function(x) {
  data.table::fcase(
    x == "Sputum GPL570", "Sputum\nGPL570",
    x == "Bronchial brushing GPL570", "Brushing\nGPL570",
    x == "Bronchial biopsy GPL570", "Biopsy\nGPL570",
    x == "Bronchial brushing RNA-seq", "Brushing\nRNA-seq",
    default = x
  )
}

save_figure_four_formats <- function(plot, directory, stem, width_mm,
                                     height_mm, source_files,
                                     source_script =
                                       "R/12_generate_final_figures.R") {
  # Explanatory prose belongs in the external figure legend, not inside the
  # plotting area. Retain figure/panel titles, axes, n labels and required keys.
  plot <- plot & theme(
    plot.subtitle = element_blank(),
    plot.caption = element_blank()
  )
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  output_files <- c(
    pdf = file.path(directory, paste0(stem, ".pdf")),
    svg = file.path(directory, paste0(stem, ".svg")),
    png = file.path(directory, paste0(stem, "_PNG300.png")),
    tiff = file.path(directory, paste0(stem, "_TIFF600.tiff"))
  )

  # Use R's native vector PDF device. The local R installation does not have
  # XQuartz/Cairo libraries, and scientific content must not depend on them.
  native_pdf <- function(filename, ...) {
    grDevices::pdf(file = filename, useDingbats = FALSE, ...)
  }
  ggplot2::ggsave(
    output_files["pdf"], plot = plot, width = width_in, height = height_in,
    units = "in", device = native_pdf, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    output_files["svg"], plot = plot, width = width_in, height = height_in,
    units = "in", device = svglite::svglite, bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    output_files["png"], plot = plot, width = width_in, height = height_in,
    units = "in", device = ragg::agg_png, dpi = 300, bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    output_files["tiff"], plot = plot, width = width_in, height = height_in,
    units = "in", device = ragg::agg_tiff, dpi = 600,
    compression = "lzw", bg = "white", limitsize = FALSE
  )

  for (format_name in names(output_files)) {
    output <- unname(output_files[format_name])
    if (!file.exists(output) || file.info(output)$size <= 0) {
      stop("Figure export failed: ", output)
    }
    figure_manifest_rows[[length(figure_manifest_rows) + 1L]] <<-
      data.table::data.table(
        figure_stem = stem,
        output_format = format_name,
        output_file = normalizePath(output, mustWork = TRUE),
        output_bytes = file.info(output)$size,
        width_mm = width_mm,
        height_mm = height_mm,
        source_data_files = paste(source_files, collapse = "; "),
        source_script = source_script,
        generated_in = "R",
        hand_edited = FALSE
      )
  }
  invisible(output_files)
}

make_long_matrix <- function(matrix, dataset_name, score_scale,
                             cluster_rows = TRUE, cluster_columns = TRUE) {
  row_order <- seq_len(nrow(matrix))
  column_order <- seq_len(ncol(matrix))
  if (cluster_rows && nrow(matrix) > 1L) {
    row_order <- stats::hclust(stats::dist(matrix))$order
  }
  if (cluster_columns && ncol(matrix) > 1L) {
    column_order <- stats::hclust(stats::dist(t(matrix)))$order
  }
  ordered <- matrix[row_order, column_order, drop = FALSE]
  long <- data.table::as.data.table(as.table(ordered))
  data.table::setnames(long, c("pathway", "Subject_ID", "score"))
  long[, `:=`(
    pathway_label = unname(pathway_labels[as.character(pathway)]),
    participant_order = match(as.character(Subject_ID), colnames(ordered)),
    pathway_order = match(as.character(pathway), rownames(ordered)),
    dataset = dataset_name,
    dataset_label = datasets[[dataset_name]]$label,
    score_scale = score_scale
  )]
  long
}

make_heatmap_panel <- function(long, fill_limits, title, subtitle,
                               show_y_labels = TRUE, legend_title,
                               show_legend = TRUE) {
  ggplot(long, aes(participant_order, pathway_order, fill = score)) +
    geom_raster() +
    scale_fill_gradient2(
      low = score_heatmap_colours["low"], mid = score_heatmap_colours["mid"],
      high = score_heatmap_colours["high"], midpoint = 0,
      limits = fill_limits, oob = scales::squish, name = legend_title
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_reverse(
      breaks = if (show_y_labels) unique(long$pathway_order) else NULL,
      labels = if (show_y_labels) unique(long$pathway_label) else NULL,
      expand = c(0, 0)
    ) +
    labs(title = title, subtitle = subtitle, x = "Participants", y = NULL) +
    theme_void(base_family = "Helvetica", base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.2, colour = "#333333", hjust = 0),
      axis.title.x = element_text(size = 7.5, margin = margin(t = 3)),
      axis.text.y = element_text(size = 5.5, hjust = 1),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text = element_text(size = 6.5),
      legend.key.width = grid::unit(12, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

write_figure_source <- function(table, filename) {
  figure_source_dir <- file.path(paths$results, "figure_source_data")
  dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(table, file.path(figure_source_dir, filename))
}
