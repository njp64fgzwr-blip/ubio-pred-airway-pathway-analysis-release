# Regenerate dissertation Table 1 with its correct embedded number in R.
#
# Numerical content is reconstructed from the validated pipeline tables. The
# layout and values are unchanged; only the dissertation table number is fixed.

source(file.path("R", "00_common.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(ragg)
  library(svglite)
})

sample_source_csv <- file.path(
  paths$tables, "Table_S_sample_availability.csv"
)
preprocessing_source_csv <- file.path(
  paths$tables, "Table_S_expression_preparation_and_mapping_QC.csv"
)
required_sources <- c(sample_source_csv, preprocessing_source_csv)
if (any(!file.exists(required_sources))) {
  stop(
    "Table 1 source files are missing: ",
    paste(required_sources[!file.exists(required_sources)], collapse = "; ")
  )
}
output_dir <- paths$final_tables
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_source <- fread(sample_source_csv)
preprocessing_source <- fread(preprocessing_source_csv)
table_source <- merge(
  sample_source,
  preprocessing_source[, .(dataset, final_genes, normalization, technology)],
  by = "dataset", all = FALSE, sort = FALSE
)
canonical_order <- canonical_dataset_order
table_source[, dataset_order := match(dataset, canonical_order)]
setorder(table_source, dataset_order)
if (!identical(table_source$dataset, canonical_order)) {
  stop("The canonical four-dataset order is not present in the source table.")
}

wrap_for_table <- function(x, width) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

table_source[, display_preprocessing := fifelse(
  dataset == "brushing_rnaseq",
  paste0(
    "Sparse counts; restore zeros; aggregate symbols; retain genes with ",
    "CPM >=1 in >=12 samples; log2(CPM+1)"
  ),
  paste0(
    "Supplied LOG2E; complete gene symbols; harmonise aliases; collapse ",
    "duplicate probes by the median"
  )
)]

overview_display <- table_source[, .(
  Dataset = wrap_for_table(dataset_label, 18),
  Specimen = wrap_for_table(compartment, 14),
  Platform = platform,
  n = transcriptomic_profiles,
  `Final genes` = final_genes,
  `Clinical IDs` = paste0(matched_clinical_ids, "/", transcriptomic_profiles)
)]
preprocessing_display <- table_source[, .(
  Dataset = wrap_for_table(dataset_label, 19),
  Preprocessing = wrap_for_table(display_preprocessing, 68)
)]

word_table_theme <- function(base_size) {
  ttheme_minimal(
    base_size = base_size,
    padding = unit(c(3.8, 4.6), "mm"),
    core = list(
      fg_params = list(hjust = 0, x = 0.025, col = "#111111"),
      bg_params = list(
        fill = rep(c("#FFFFFF", "#EEF3F7"), 4), col = "#D9E1E8"
      )
    ),
    colhead = list(
      fg_params = list(
        col = "white", fontface = "bold", hjust = 0.5, x = 0.5
      ),
      bg_params = list(fill = "#1F4E78", col = "#1F4E78")
    )
  )
}

overview_grob <- tableGrob(
  overview_display, rows = NULL, theme = word_table_theme(9.8)
)
overview_grob$widths <- unit(
  c(1.65, 1.35, 0.85, 0.45, 0.90, 1.00), "null"
)
preprocessing_grob <- tableGrob(
  preprocessing_display, rows = NULL, theme = word_table_theme(9.8)
)
preprocessing_grob$widths <- unit(c(1.65, 5.35), "null")

layout_grob <- arrangeGrob(
  textGrob(
    "Table 1. Airway transcriptomic datasets and preprocessing",
    x = 0, hjust = 0,
    gp = gpar(fontsize = 13, fontface = "bold", col = "#111111")
  ),
  textGrob(
    "A  Dataset and sample overview (all four inputs are transcriptomic)",
    x = 0, hjust = 0,
    gp = gpar(fontsize = 12, fontface = "bold", col = "#111111")
  ),
  overview_grob,
  textGrob(
    "B  Preprocessing applied before ssGSEA",
    x = 0, hjust = 0,
    gp = gpar(fontsize = 12, fontface = "bold", col = "#111111")
  ),
  preprocessing_grob,
  ncol = 1,
  heights = c(0.07, 0.07, 0.34, 0.07, 0.45),
  padding = unit(2, "mm")
)
table_grob <- grobTree(
  rectGrob(gp = gpar(fill = "white", col = NA)), layout_grob
)

stem <- "01_Table_1_sample_omics_and_preprocessing"
width_mm <- 180
height_mm <- 138
ggsave(
  file.path(output_dir, paste0(stem, ".pdf")), table_grob,
  width = width_mm, height = height_mm, units = "mm", device = grDevices::pdf
)
ggsave(
  file.path(output_dir, paste0(stem, ".svg")), table_grob,
  width = width_mm, height = height_mm, units = "mm", device = svglite
)
ggsave(
  file.path(output_dir, paste0(stem, "_PNG300.png")), table_grob,
  width = width_mm, height = height_mm, units = "mm", dpi = 300,
  device = agg_png, background = "white"
)
ggsave(
  file.path(output_dir, paste0(stem, "_TIFF600.tiff")), table_grob,
  width = width_mm, height = height_mm, units = "mm", dpi = 600,
  device = agg_tiff, compression = "lzw", background = "white"
)

writeLines(
  c(
    "TABLE 1 CORRECT NUMBERING: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "Generated directly in R from validated sample-availability and preprocessing tables.",
    "Embedded title: Table 1. Airway transcriptomic datasets and preprocessing",
    paste("Rows:", nrow(table_source)),
    paste("Canonical dataset order:", paste(table_source$dataset, collapse = "; ")),
    "Exports: PDF, SVG, PNG300 and TIFF600"
  ),
  file.path(output_dir, "TABLE_1_R_GENERATION_STATUS.txt")
)

message("Regenerated correctly numbered Table 1 in four R-produced formats.")
