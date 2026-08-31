# Render every clean and draft-matched production PDF for direct visual QA.
# This script is run from R/RStudio. It does not alter scientific plots.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "magick", "digest"))
message_rule("15: Rendering all clean and draft-matched PDFs for QA")

clean_manifest_file <- file.path(
  paths$validation, "FIGURE_EXPORT_MANIFEST.csv"
)
draft_manifest_file <- file.path(
  paths$validation, "DRAFT_MATCHED_FIGURE_EXPORT_MANIFEST.csv"
)
clean_manifest <- data.table::fread(clean_manifest_file)
draft_manifest <- data.table::fread(draft_manifest_file)
clean_manifest[, figure_set := "base_dissertation"]
draft_manifest[, figure_set := "draft_matched"]
manifest <- data.table::rbindlist(
  list(clean_manifest, draft_manifest), fill = TRUE
)
pdf_manifest <- manifest[output_format == "pdf"]

expected_counts <- c(base_dissertation = 18L, draft_matched = 28L)
actual_counts <- pdf_manifest[, .N, by = figure_set]
for (set_name in names(expected_counts)) {
  if (actual_counts[figure_set == set_name, N] != expected_counts[set_name]) {
    stop("Unexpected PDF count for figure set: ", set_name)
  }
}
if (data.table::uniqueN(pdf_manifest$figure_stem) != sum(expected_counts)) {
  stop("Figure stems are not unique across the combined production manifests.")
}

run_id <- format(Sys.time(), "run_%Y%m%d_%H%M%S")
run_dir <- file.path(paths$figures_qa, paste0("all_figures_", run_id))
pdftoppm <- Sys.getenv("PDFTOPPM_BIN", unset = Sys.which("pdftoppm"))
if (!nzchar(pdftoppm) || !file.exists(pdftoppm)) {
  stop("Poppler pdftoppm is unavailable. Install Poppler or set PDFTOPPM_BIN.")
}

render_rows <- vector("list", nrow(pdf_manifest))
for (i in seq_len(nrow(pdf_manifest))) {
  row <- pdf_manifest[i]
  set_dir <- file.path(run_dir, row$figure_set)
  colour_dir <- file.path(set_dir, "colour")
  grayscale_dir <- file.path(set_dir, "grayscale")
  dir.create(colour_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(grayscale_dir, recursive = TRUE, showWarnings = FALSE)

  colour_file <- file.path(colour_dir, paste0(row$figure_stem, ".png"))
  gray_file <- file.path(grayscale_dir, paste0(row$figure_stem, ".png"))
  colour_prefix <- sub("[.]png$", "", colour_file)
  poppler_output <- system2(
    pdftoppm,
    args = c(
      "-png", "-r", "180", "-singlefile",
      row$output_file, colour_prefix
    ),
    stdout = TRUE, stderr = TRUE
  )
  exit_status <- attr(poppler_output, "status")
  if (!is.null(exit_status) && exit_status != 0L) {
    stop(
      "Poppler render failed for ", row$output_file, ": ",
      paste(poppler_output, collapse = "\n")
    )
  }
  if (!file.exists(colour_file)) {
    stop("Poppler did not create: ", colour_file)
  }
  colour_image <- magick::image_read(colour_file)
  if (length(colour_image) != 1L) {
    stop("Production figure PDF is not single-page: ", row$output_file)
  }
  gray_image <- magick::image_convert(colour_image, colorspace = "gray")
  magick::image_write(gray_image, path = gray_file, format = "png")

  colour_info <- magick::image_info(colour_image)
  gray_info <- magick::image_info(magick::image_read(gray_file))
  render_rows[[i]] <- data.table::data.table(
    figure_set = row$figure_set,
    figure_stem = row$figure_stem,
    source_pdf = row$output_file,
    source_pdf_sha256 = digest::digest(
      row$output_file, algo = "sha256", file = TRUE
    ),
    colour_render = normalizePath(colour_file, mustWork = TRUE),
    grayscale_render = normalizePath(gray_file, mustWork = TRUE),
    colour_width_px = colour_info$width,
    colour_height_px = colour_info$height,
    grayscale_width_px = gray_info$width,
    grayscale_height_px = gray_info$height,
    colour_bytes = file.info(colour_file)$size,
    grayscale_bytes = file.info(gray_file)$size,
    render_density_dpi = 180L,
    rendered_by = "R system2 Poppler; grayscale via R magick"
  )
}

render_manifest <- data.table::rbindlist(render_rows)
if (any(render_manifest$colour_bytes <= 0L) ||
    any(render_manifest$grayscale_bytes <= 0L) ||
    any(render_manifest$colour_width_px <= 0L) ||
    any(render_manifest$colour_height_px <= 0L)) {
  stop("At least one combined QA render is empty or invalid.")
}
data.table::fwrite(
  render_manifest,
  file.path(paths$validation, "ALL_FIGURE_PDF_DIRECT_RENDER_MANIFEST.csv")
)
writeLines(
  normalizePath(run_dir),
  file.path(paths$validation, "LATEST_ALL_FIGURE_QA_RENDER_DIR.txt")
)
writeLines(
  c(
    "ALL FIGURE PDF DIRECT RENDER: COMPUTATIONAL PASS",
    paste0("Run directory: ", normalizePath(run_dir)),
    paste0("Clean PDFs rendered in colour: 18/18"),
    paste0("Clean PDFs rendered in grayscale: 18/18"),
    paste0("Draft-matched PDFs rendered in colour: 28/28"),
    paste0("Draft-matched PDFs rendered in grayscale: 28/28"),
    "Total direct PDF renders: 92 (46 colour + 46 grayscale).",
    "Density: 180 dpi at each figure's native 180-mm publication width.",
    "Visual inspection status is recorded separately after every render is reviewed."
  ),
  file.path(paths$validation, "ALL_FIGURE_QA_RENDER_STATUS.txt")
)

message(
  "Rendered ", nrow(render_manifest),
  " clean/draft-matched PDFs in colour and grayscale: ", run_dir
)
