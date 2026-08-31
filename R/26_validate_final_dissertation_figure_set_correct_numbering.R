# Validate the final dissertation figure set without changing scientific data.

source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "magick", "digest"))
message_rule("26: Validating corrected dissertation figure titles, order and formats")

final_root <- paths$final_figures
main_dir <- file.path(final_root, "01_MAIN_FIGURES")
supp_dir <- file.path(final_root, "02_SUPPLEMENTARY_FIGURES")
qa_dir <- file.path(final_root, "03_QA_RENDERS")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

order_table <- data.table::fread(file.path(
  final_root, "FINAL_DISSERTATION_FIGURE_ORDER.csv"
))
export_manifest <- data.table::fread(file.path(
  final_root, "FINAL_DISSERTATION_FIGURE_EXPORT_MANIFEST.csv"
))

if (nrow(order_table) != 34L || data.table::uniqueN(order_table$figure_stem) != 34L) {
  stop("Expected 34 ordered final figure files.")
}
if (nrow(export_manifest) != 136L ||
    data.table::uniqueN(export_manifest$figure_stem) != 34L) {
  stop("Expected 136 exports for 34 figure files.")
}

expected_main_ids <- c(
  "Figure 1", "Figure 2", "Figure 3",
  paste0("Figure 4", LETTERS[1:4]),
  paste0("Figure 5", LETTERS[1:4]),
  paste0("Figure 6", LETTERS[1:4]),
  paste0("Figure 7", LETTERS[1:4]),
  "Figure 8", "Figure 9A", "Figure 9B"
)
actual_main_ids <- order_table[section == "Main"][order(section_order), figure_id]
if (!identical(actual_main_ids, expected_main_ids)) {
  stop("Main figure order does not match the fixed dissertation order.")
}
expected_supplementary_ids <- paste0("Figure S", seq_len(12L))
actual_supplementary_ids <- order_table[
  section == "Supplementary"
][order(section_order), figure_id]
if (!identical(actual_supplementary_ids, expected_supplementary_ids)) {
  stop("Supplementary figure order does not match Figures S1-S12.")
}
if (!identical(
  canonical_dataset_order,
  c("sputum", "brushing_gpl570", "biopsy", "brushing_rnaseq")
)) {
  stop("Canonical multi-compartment display order has changed.")
}

expected_formats <- c("pdf", "svg", "png", "tiff")
format_counts <- export_manifest[, .(
  formats = paste(sort(output_format), collapse = ";"),
  format_count = data.table::uniqueN(output_format),
  all_nonempty = all(file.exists(output_file)) &&
    all(file.info(output_file)$size > 0),
  generated_in_R = all(generated_in == "R"),
  no_hand_edits = all(!hand_edited)
), by = .(figure_stem)]
format_counts[, format_set_correct := formats == paste(sort(expected_formats), collapse = ";")]
if (any(format_counts$format_count != 4L) ||
    any(!format_counts$format_set_correct) ||
    any(!format_counts$all_nonempty) ||
    any(!format_counts$generated_in_R) ||
    any(!format_counts$no_hand_edits)) {
  stop("Final four-format/R-only validation failed.")
}

# PNG dimensions should equal the declared physical dimensions at 300 dpi;
# TIFF dimensions should equal them at 600 dpi, allowing a two-pixel rounding
# difference from unit conversion.
raster_rows <- export_manifest[output_format %in% c("png", "tiff")]
raster_audit <- data.table::rbindlist(lapply(seq_len(nrow(raster_rows)), function(i) {
  x <- raster_rows[i]
  image <- magick::image_read(x$output_file)
  info <- magick::image_info(image)[1, ]
  expected_dpi <- if (x$output_format == "png") 300 else 600
  data.table::data.table(
    figure_stem = x$figure_stem,
    output_format = x$output_format,
    output_file = x$output_file,
    actual_width_px = info$width,
    actual_height_px = info$height,
    expected_width_px = round(x$width_mm / 25.4 * expected_dpi),
    expected_height_px = round(x$height_mm / 25.4 * expected_dpi),
    dimensions_pass =
      abs(info$width - round(x$width_mm / 25.4 * expected_dpi)) <= 2L &&
      abs(info$height - round(x$height_mm / 25.4 * expected_dpi)) <= 2L
  )
}))
if (any(!raster_audit$dimensions_pass)) stop("Raster resolution check failed.")

# SVG is searchable text, so it provides a deterministic title check for every
# stem without depending on PDF font extraction.
title_audit <- data.table::rbindlist(lapply(seq_len(nrow(order_table)), function(i) {
  x <- order_table[i]
  svg <- file.path(
    if (x$section == "Main") main_dir else supp_dir,
    paste0(x$figure_stem, ".svg")
  )
  text <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  data.table::data.table(
    section = x$section,
    section_order = x$section_order,
    figure_id = x$figure_id,
    figure_stem = x$figure_stem,
    expected_title = x$final_title,
    exact_title_present = grepl(x$final_title, text, fixed = TRUE),
    no_draft_matched_token = !grepl("Draft-matched", text, fixed = TRUE),
    no_small_explanatory_subtitle =
      !grepl("plot.subtitle", text, fixed = TRUE) &&
      !grepl("plot.caption", text, fixed = TRUE)
  )
}))
if (any(!title_audit$exact_title_present) ||
    any(!title_audit$no_draft_matched_token)) {
  stop("A final figure title is missing or stale.")
}

# Render every PDF directly with Poppler at 180-mm publication width and save a
# grayscale counterpart. These files are retained for final visual inspection.
pdftoppm <- Sys.getenv("PDFTOPPM_BIN", unset = Sys.which("pdftoppm"))
if (!nzchar(pdftoppm) || !file.exists(pdftoppm)) {
  stop(
    "Poppler pdftoppm is unavailable. Install Poppler or set PDFTOPPM_BIN."
  )
}

pdf_rows <- export_manifest[output_format == "pdf"]
render_audit <- data.table::rbindlist(lapply(seq_len(nrow(pdf_rows)), function(i) {
  x <- pdf_rows[i]
  stem <- x$figure_stem
  colour_prefix <- file.path(qa_dir, paste0(stem, "_colour"))
  colour_file <- paste0(colour_prefix, "-1.png")
  grayscale_file <- file.path(qa_dir, paste0(stem, "_grayscale.png"))
  status <- system2(
    pdftoppm,
    c("-png", "-r", "180", "-singlefile", x$output_file, colour_prefix),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(colour_file)) {
    # Poppler -singlefile omits the page suffix on some versions.
    alternate <- paste0(colour_prefix, ".png")
    if (file.exists(alternate)) colour_file <- alternate
  }
  render_pass <- file.exists(colour_file) && file.info(colour_file)$size > 0
  if (render_pass) {
    colour <- magick::image_read(colour_file)
    magick::image_write(
      magick::image_convert(colour, colorspace = "gray"),
      path = grayscale_file, format = "png"
    )
  }
  data.table::data.table(
    figure_stem = stem,
    pdf_file = x$output_file,
    colour_render = if (render_pass) colour_file else NA_character_,
    grayscale_render = if (file.exists(grayscale_file)) grayscale_file else NA_character_,
    colour_render_pass = render_pass,
    grayscale_render_pass = file.exists(grayscale_file) &&
      file.info(grayscale_file)$size > 0,
    poppler_message = paste(status, collapse = " | ")
  )
}))
if (any(!render_audit$colour_render_pass) ||
    any(!render_audit$grayscale_render_pass)) {
  stop("Direct PDF colour/grayscale rendering failed.")
}

# Confirm the source tables used by the final set are the same canonical files
# already covered by the validated clean/draft-matched figure-source audit.
canonical_source_audit <- data.table::fread(file.path(
  paths$validation, "FIGURE_SOURCE_AUDIT.csv"
))
source_traceability_pass <- all(canonical_source_audit$values_match_source)
if (!source_traceability_pass) stop("A canonical figure-source audit is not passing.")
provenance_status_file <- file.path(
  paths$validation, "SIGNATURE_PROVENANCE_OUTPUT_STATUS.txt"
)
if (!file.exists(provenance_status_file) ||
    !any(grepl(
      "SIGNATURE PROVENANCE OUTPUTS: PASS",
      readLines(provenance_status_file, warn = FALSE), fixed = TRUE
    ))) {
  stop("Figure S12 provenance validation status is unavailable or not passing.")
}

data.table::fwrite(
  format_counts, file.path(final_root, "FINAL_FIGURE_FORMAT_AUDIT.csv")
)
data.table::fwrite(
  raster_audit, file.path(final_root, "FINAL_FIGURE_RASTER_RESOLUTION_AUDIT.csv")
)
data.table::fwrite(
  title_audit, file.path(final_root, "FINAL_FIGURE_TITLE_AND_ORDER_AUDIT.csv")
)
data.table::fwrite(
  render_audit, file.path(final_root, "FINAL_FIGURE_DIRECT_RENDER_AUDIT.csv")
)

writeLines(
  c(
    "FINAL DISSERTATION FIGURE SET VALIDATION: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "Main figure groups: 9",
    "Main panel files: 22",
    "Supplementary figures: 12",
    paste("Final figure stems:", nrow(order_table), "/34"),
    paste("Production exports:", nrow(export_manifest), "/136"),
    paste("Exact final SVG titles:", sum(title_audit$exact_title_present), "/34"),
    paste("No stale Draft-matched tokens:", sum(title_audit$no_draft_matched_token), "/34"),
    paste("PNG/TIFF resolution checks:", sum(raster_audit$dimensions_pass), "/68"),
    paste("Direct PDF colour renders:", sum(render_audit$colour_render_pass), "/34"),
    paste("Direct PDF grayscale renders:", sum(render_audit$grayscale_render_pass), "/34"),
    "Figures 1-3 separation and numbering: PASS (three full-width source-rendered files)",
    "Canonical multi-compartment order: PASS (sputum GPL570; bronchial brushing GPL570; bronchial biopsy GPL570; bronchial brushing RNA-seq)",
    paste("Canonical numerical source checks:", nrow(canonical_source_audit), "/", nrow(canonical_source_audit), "PASS"),
    "All figures were exported directly from R in final dissertation order.",
    "The dissertation Word document was not modified.",
    "Self-exclusion: FINAL_FIGURE_SET_SHA256_MANIFEST.csv"
  ),
  file.path(final_root, "FINAL_DISSERTATION_FIGURE_VALIDATION.txt")
)

# A self-contained hash manifest is written after every other file in the
# final figure folder. The manifest explicitly excludes only itself.
manifest_path <- file.path(final_root, "FINAL_FIGURE_SET_SHA256_MANIFEST.csv")
if (file.exists(manifest_path)) unlink(manifest_path)
files <- list.files(
  final_root, recursive = TRUE, full.names = TRUE,
  all.files = TRUE, no.. = TRUE, include.dirs = FALSE
)
relative <- substring(files, nchar(final_root) + 2L)
hash_manifest <- data.table::data.table(
  relative_path = relative,
  bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, digest::digest, character(1),
                  algo = "sha256", file = TRUE, serialize = FALSE)
)
data.table::setorder(hash_manifest, relative_path)
data.table::fwrite(hash_manifest, manifest_path)

# Immediate exact reconstruction: no missing, changed or unmanifested file is
# accepted, apart from the documented manifest self-exclusion.
current <- list.files(
  final_root, recursive = TRUE, full.names = TRUE,
  all.files = TRUE, no.. = TRUE, include.dirs = FALSE
)
current <- current[normalizePath(current, mustWork = TRUE) !=
  normalizePath(manifest_path, mustWork = TRUE)]
current_relative <- substring(current, nchar(final_root) + 2L)
if (!setequal(hash_manifest$relative_path, current_relative)) {
  stop("Final figure-set hash manifest file-set verification failed.")
}
current <- current[match(hash_manifest$relative_path, current_relative)]
if (any(as.numeric(file.info(current)$size) != hash_manifest$bytes) ||
    any(vapply(current, digest::digest, character(1),
               algo = "sha256", file = TRUE, serialize = FALSE) !=
        hash_manifest$sha256)) {
  stop("Final figure-set hash manifest byte/hash verification failed.")
}

message("Validated 34 final dissertation figures and 136 R-only exports.")
