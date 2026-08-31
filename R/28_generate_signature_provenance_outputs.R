# Generate the public-safe provenance outputs and Supplementary Figure S12.
#
# The frozen signature-membership file and audited provenance CSV are supplied
# locally by an authorised user. Exact gene membership is checked in memory but
# is never written by this script.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "ggplot2", "stringr", "ragg", "svglite"))
message_rule("28: Generating public-safe signature provenance outputs")

if (is.na(SIGNATURE_FILE) || !nzchar(SIGNATURE_FILE) ||
    !file.exists(SIGNATURE_FILE)) {
  stop(
    "The frozen signature file is unavailable. Set AIRWAY_SIGNATURE_FILE or ",
    "inputs.signature_file in config/analysis_parameters.yml."
  )
}
require_signature_provenance_input()
if (sha256_file(SIGNATURE_PROVENANCE_FILE) !=
    signature_settings$provenance_sha256) {
  stop("The provenance input does not match the locked non-membership SHA-256 value.")
}

signature_lines <- readLines(SIGNATURE_FILE, warn = FALSE, encoding = "UTF-8")
signature_parts <- lapply(
  signature_lines,
  function(line) trimws(strsplit(line, "\t", fixed = TRUE)[[1L]])
)
signature_parts <- lapply(signature_parts, function(x) x[nzchar(x)])
signature_ids <- vapply(signature_parts, `[`, character(1), 1L)
signature_counts <- vapply(
  signature_parts,
  function(x) length(unique(x[-1L])),
  integer(1)
)
names(signature_counts) <- signature_ids

if (length(signature_ids) != signature_settings$expected_count ||
    anyDuplicated(signature_ids)) {
  stop("The local signature file is not the expected unique 39-signature library.")
}
if (sha256_file(SIGNATURE_FILE) != signature_settings$frozen_sha256) {
  stop("The local signature file does not match the frozen SHA-256 value.")
}

provenance_input <- data.table::fread(
  SIGNATURE_PROVENANCE_FILE, check.names = FALSE, showProgress = FALSE
)
required_columns <- c(
  "signature_id", "reader_label", "genes_n", "biological_family",
  "derivation_group", "source_material", "source_contrast_or_method",
  "direction", "asthma_specificity", "primary_reference",
  "traceability_grade", "interpretation_limit"
)
missing_columns <- setdiff(required_columns, names(provenance_input))
if (length(missing_columns)) {
  stop(
    "The provenance input is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
if (nrow(provenance_input) != signature_settings$expected_count ||
    anyDuplicated(provenance_input$signature_id) ||
    !setequal(provenance_input$signature_id, signature_ids)) {
  stop("The provenance input does not match the frozen signature identifiers.")
}

provenance <- provenance_input[
  match(signature_ids, provenance_input$signature_id),
  ..required_columns
]
provenance[, genes_n := as.integer(genes_n)]
if (any(provenance$genes_n != signature_counts[provenance$signature_id])) {
  stop("Provenance gene-entry counts do not match the frozen signature file.")
}
if (!all(provenance$traceability_grade %in% c("A", "B", "C"))) {
  stop("Traceability grades must be A, B or C.")
}

# Whitelisted public tables: no gene-symbol or participant-level column is
# selected from the locally supplied source.
public_provenance <- provenance[, .(
  signature_id, reader_label, genes_n, biological_family, derivation_group,
  source_material, source_contrast_or_method, direction, asthma_specificity,
  primary_reference, traceability_grade, interpretation_limit
)]
public_counts <- public_provenance[, .(
  signature_id, reader_label, genes_n, biological_family, derivation_group,
  traceability_grade
)]
data.table::fwrite(
  public_provenance,
  file.path(paths$tables, "Table_S_signature_provenance_public.csv")
)
data.table::fwrite(
  public_counts,
  file.path(paths$tables, "Table_S_signature_counts_public.csv")
)

figure_source_dir <- file.path(paths$results, "figure_source_data")
dir.create(figure_source_dir, recursive = TRUE, showWarnings = FALSE)
figure_source_file <- file.path(
  figure_source_dir,
  "Figure_S12_signature_scope_and_provenance_public.csv"
)
data.table::fwrite(public_counts, figure_source_file)

family_levels <- c(
  "Type 2 / eosinophil / epithelium",
  "Innate / neutrophil / inflammasome",
  "Mast cell / IgE / IL-33",
  "IL-17 / ILC3",
  "Macrophage / regulatory T cell",
  "Structural / metabolic / cell stress"
)
plot_data <- data.table::copy(public_counts)
plot_data[, biological_family := factor(
  biological_family, levels = family_levels
)]
if (any(is.na(plot_data$biological_family))) {
  stop("A provenance biological-family label is outside the fixed display set.")
}
plot_data[, label := sprintf("%s  (%d)", reader_label, genes_n)]
plot_data[, label := stringr::str_wrap(label, width = 39)]
plot_data[, row_in_family := seq_len(.N), by = biological_family]
plot_data[, panel_y := .N - row_in_family + 1L, by = biological_family]

derivation_colours <- c(
  "Asthma-cohort-derived" = "#F2B6BC",
  "Asthma/airway-derived" = "#F5C5B8",
  "Airway/cell perturbation" = "#F9D9C4",
  "Immune-cell reference" = "#BFD9EA",
  "Canonical/general pathway" = "#B4C9E6",
  "Non-asthma disease-derived" = "#D8C5E5",
  "Curated marker panel" = "#D1D1D1"
)
if (any(!as.character(plot_data$derivation_group) %in%
        names(derivation_colours))) {
  stop("A provenance derivation-group label is outside the fixed display set.")
}
trace_shapes <- c(A = 21, B = 22, C = 24)

figure_s12 <- ggplot2::ggplot(
  plot_data, ggplot2::aes(x = 0.5, y = panel_y)
) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = derivation_group), width = 0.96, height = 0.82,
    colour = "white", linewidth = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = label), x = 0.035, hjust = 0, size = 2.35,
    lineheight = 0.88, colour = "#111111"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = traceability_grade), x = 0.955, size = 2.5,
    fill = "white", colour = "#111111", stroke = 0.7
  ) +
  ggplot2::facet_wrap(~ biological_family, ncol = 2, scales = "free_y") +
  ggplot2::scale_fill_manual(
    values = derivation_colours, name = "Derivation context"
  ) +
  ggplot2::scale_shape_manual(
    values = trace_shapes, name = "Traceability",
    labels = c(
      A = "A: exact list/source reconstructed",
      B = "B: source supported; list partly reconstructed",
      C = "C: material ambiguity"
    )
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") +
  ggplot2::labs(
    title = paste(
      "Supplementary Figure S12. Scope and provenance of 39 airway",
      "pathway signatures"
    ),
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 9, base_family = "sans") +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold", size = 9, colour = "#163A5F"
    ),
    strip.background = ggplot2::element_rect(fill = "#E8EEF5", colour = NA),
    plot.title = ggplot2::element_text(
      face = "bold", size = 12, margin = ggplot2::margin(b = 9)
    ),
    plot.title.position = "plot",
    panel.spacing.y = grid::unit(5, "mm"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = ggplot2::element_text(size = 7.4),
    legend.title = ggplot2::element_text(size = 8, face = "bold"),
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      nrow = 2, byrow = TRUE, order = 1,
      title.position = "top", title.hjust = 0.5
    ),
    shape = ggplot2::guide_legend(
      nrow = 1, order = 2, title.position = "top", title.hjust = 0.5
    )
  )

supplementary_dir <- file.path(paths$final_figures, "02_SUPPLEMENTARY_FIGURES")
dir.create(supplementary_dir, recursive = TRUE, showWarnings = FALSE)
stem <- "12_Figure_S12_pathway_signature_scope_and_provenance"
width_mm <- 180
height_mm <- 255
output_files <- c(
  pdf = file.path(supplementary_dir, paste0(stem, ".pdf")),
  svg = file.path(supplementary_dir, paste0(stem, ".svg")),
  png = file.path(supplementary_dir, paste0(stem, "_PNG300.png")),
  tiff = file.path(supplementary_dir, paste0(stem, "_TIFF600.tiff"))
)
ggplot2::ggsave(
  output_files[["pdf"]], figure_s12, width = width_mm, height = height_mm,
  units = "mm", device = grDevices::pdf
)
ggplot2::ggsave(
  output_files[["svg"]], figure_s12, width = width_mm, height = height_mm,
  units = "mm", device = svglite::svglite
)
ggplot2::ggsave(
  output_files[["png"]], figure_s12, width = width_mm, height = height_mm,
  units = "mm", dpi = 300, device = ragg::agg_png, background = "white"
)
ggplot2::ggsave(
  output_files[["tiff"]], figure_s12, width = width_mm, height = height_mm,
  units = "mm", dpi = 600, device = ragg::agg_tiff, compression = "lzw",
  background = "white"
)
if (any(!file.exists(output_files)) || any(file.info(output_files)$size <= 0)) {
  stop("One or more Figure S12 exports are missing or empty.")
}

relative_source <- file.path(
  "results", "figure_source_data",
  "Figure_S12_signature_scope_and_provenance_public.csv"
)
order_file <- file.path(
  paths$final_figures, "FINAL_DISSERTATION_FIGURE_ORDER.csv"
)
manifest_file <- file.path(
  paths$final_figures, "FINAL_DISSERTATION_FIGURE_EXPORT_MANIFEST.csv"
)
if (!file.exists(order_file) || !file.exists(manifest_file)) {
  stop("Run exporter 25 before adding Figure S12.")
}

order_table <- data.table::fread(order_file)
order_table <- order_table[figure_id != "Figure S12"]
s12_order <- data.table::data.table(
  section = "Supplementary",
  section_order = 12L,
  figure_id = "Figure S12",
  final_title = paste(
    "Supplementary Figure S12. Scope and provenance of 39 airway",
    "pathway signatures"
  ),
  figure_stem = stem,
  source_plot_basis = "R-generated public-safe signature provenance map",
  width_mm = width_mm,
  height_mm = height_mm,
  source_data_files = relative_source
)
order_table <- data.table::rbindlist(list(order_table, s12_order), fill = TRUE)
data.table::fwrite(order_table, order_file)

export_manifest <- data.table::fread(manifest_file)
export_manifest <- export_manifest[figure_id != "Figure S12"]
s12_manifest <- data.table::rbindlist(lapply(names(output_files), function(fmt) {
  data.table::data.table(
    figure_stem = stem,
    output_format = fmt,
    output_file = normalizePath(output_files[[fmt]], mustWork = TRUE),
    output_bytes = as.numeric(file.info(output_files[[fmt]])$size),
    width_mm = width_mm,
    height_mm = height_mm,
    source_data_files = relative_source,
    source_script = "R/28_generate_signature_provenance_outputs.R",
    generated_in = "R",
    hand_edited = FALSE,
    section = "Supplementary",
    section_order = 12L,
    figure_id = "Figure S12",
    final_title = s12_order$final_title,
    source_plot_basis = s12_order$source_plot_basis
  )
}))
export_manifest <- data.table::rbindlist(
  list(export_manifest, s12_manifest), fill = TRUE
)
data.table::fwrite(export_manifest, manifest_file)

order_markdown_file <- file.path(
  paths$final_figures, "FINAL_DISSERTATION_FIGURE_ORDER.md"
)
if (file.exists(order_markdown_file)) {
  order_lines <- readLines(order_markdown_file, warn = FALSE)
  if (!any(grepl("Figure S12", order_lines, fixed = TRUE))) {
    writeLines(
      c(order_lines, "12. Figure S12 — Signature scope and provenance"),
      order_markdown_file
    )
  }
}

legend_file <- file.path(
  paths$final_figures, "FINAL_DISSERTATION_FIGURE_LEGENDS.md"
)
if (file.exists(legend_file)) {
  legend_lines <- readLines(legend_file, warn = FALSE)
  if (!any(grepl("Supplementary Figure S12", legend_lines, fixed = TRUE))) {
    writeLines(
      c(
        legend_lines, "",
        "## Supplementary Figure S12. Biological scope and derivation context of the 39 prespecified airway signatures",
        paste(
          "Each tile represents one ssGSEA signature and reports its reader-facing label and frozen source-list gene-entry count.",
          "Tile colour identifies derivation context; shape gives provenance traceability: A, exact analysed list and source reconstructed;",
          "B, source biology or contrast supported but exact-list construction partly reconstructed; C, material, label, source, symbol or construction ambiguity.",
          "Asthma-relevant means selected to represent biology implicated in asthma or airway disease; not every list was derived in asthma.",
          "This is a conceptual provenance map, not a causal pathway network."
        )
      ),
      legend_file
    )
  }
}

writeLines(
  c(
    "SIGNATURE PROVENANCE OUTPUTS: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Signatures verified:", nrow(public_provenance)),
    paste("Frozen signature SHA-256:", sha256_file(SIGNATURE_FILE)),
    "Exact gene membership was checked in memory and was not exported.",
    "Figure S12 exports: PDF, SVG, PNG300 and TIFF600."
  ),
  file.path(paths$validation, "SIGNATURE_PROVENANCE_OUTPUT_STATUS.txt")
)
message("Generated public-safe provenance tables and Figure S12.")
