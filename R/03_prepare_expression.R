# Prepare one gene-by-participant expression matrix for each airway dataset.
#
# Microarray: use supplied LOG2E values (verified as log2(VALUE)), fill blank
# symbols from hgu133plus2.db, harmonise unambiguous aliases, and collapse probes
# mapping to the same gene by their median within each participant.
#
# RNA-seq: use sparse raw VALUE counts, map Ensembl IDs with org.Hs.eg.db,
# supplement unresolved IDs with supplied symbols and unambiguous aliases,
# reconstruct omitted zeroes, sum features mapping to one gene, retain genes
# with CPM >= 1 in at least 10% of participants, and calculate log2(CPM + 1).

source(file.path("R", "00_common.R"))
assert_packages(c(
  "data.table", "AnnotationDbi", "hgu133plus2.db", "org.Hs.eg.db"
))
message_rule("03: Preparing gene-level airway expression matrices")

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

first_nonmissing <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[1L] else NA_character_
}

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("", "na", "nan", "null")] <- NA_character_
  # The export occasionally uses R-style duplicate suffixes such as SYNE2..1.
  sub("[.][.][0-9]+$", "", x)
}

build_symbol_harmonisation <- function(symbols) {
  source_symbol <- unique(clean_symbols(symbols))
  source_symbol <- source_symbol[!is.na(source_symbol)]
  official_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")
  alias_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ALIAS")

  output <- data.table::data.table(
    source_symbol = source_symbol,
    canonical_symbol = source_symbol,
    mapping_status = ifelse(
      source_symbol %in% official_keys, "official_symbol", "unresolved"
    ),
    candidate_count = ifelse(source_symbol %in% official_keys, 1L, 0L),
    candidate_symbols = ifelse(source_symbol %in% official_keys, source_symbol, "")
  )

  aliases <- intersect(
    output[mapping_status == "unresolved", source_symbol], alias_keys
  )
  if (length(aliases)) {
    mapping <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db,
      keys = aliases,
      keytype = "ALIAS",
      columns = "SYMBOL"
    ))
    mapping <- data.table::as.data.table(mapping)[
      !is.na(SYMBOL) & nzchar(SYMBOL), .(
        candidates = list(sort(unique(SYMBOL))),
        candidate_count = data.table::uniqueN(SYMBOL),
        candidate_symbols = paste(sort(unique(SYMBOL)), collapse = ";")
      ), by = ALIAS
    ]
    unique_map <- mapping[candidate_count == 1L]
    if (nrow(unique_map)) {
      replacement <- vapply(unique_map$candidates, `[`, character(1), 1L)
      output[match(unique_map$ALIAS, source_symbol), `:=`(
        canonical_symbol = replacement,
        mapping_status = "unambiguous_alias",
        candidate_count = 1L,
        candidate_symbols = replacement
      )]
    }
    ambiguous_map <- mapping[candidate_count > 1L]
    if (nrow(ambiguous_map)) {
      output[match(ambiguous_map$ALIAS, source_symbol), `:=`(
        mapping_status = "ambiguous_alias_retained",
        candidate_count = ambiguous_map$candidate_count,
        candidate_symbols = ambiguous_map$candidate_symbols
      )]
    }
  }
  output
}

all_exported_symbols <- list()
for (dataset_name in names(datasets)) {
  all_exported_symbols[[dataset_name]] <- data.table::fread(
    datasets[[dataset_name]]$expression_file,
    select = "GENE SYMBOL",
    na.strings = c("", "NA", "NaN", "null", "NULL"),
    showProgress = FALSE
  )[[1L]]
}

signature_sets_source <- readRDS(
  file.path(paths$analysis, "airway_39_signatures_analysis.rds")
)
all_symbols <- c(
  unlist(signature_sets_source, use.names = FALSE),
  unlist(all_exported_symbols, use.names = FALSE)
)
symbol_harmonisation <- build_symbol_harmonisation(all_symbols)
symbol_lookup <- stats::setNames(
  symbol_harmonisation$canonical_symbol,
  symbol_harmonisation$source_symbol
)

harmonise_symbols <- function(x) {
  cleaned <- clean_symbols(x)
  mapped <- unname(symbol_lookup[cleaned])
  mapped[is.na(mapped) & !is.na(cleaned)] <- cleaned[is.na(mapped) & !is.na(cleaned)]
  mapped
}

# Harmonise identifiers for analysis without changing the frozen source file.
signature_mapping_rows <- list()
gene_sets <- lapply(names(signature_sets_source), function(pathway) {
  source_genes <- signature_sets_source[[pathway]]
  analysis_genes <- harmonise_symbols(source_genes)
  signature_mapping_rows[[pathway]] <<- data.table::data.table(
    pathway = pathway,
    source_gene = source_genes,
    analysis_gene = analysis_genes,
    changed_by_annotation = source_genes != analysis_genes
  )
  unique(analysis_genes[!is.na(analysis_genes) & nzchar(analysis_genes)])
})
names(gene_sets) <- names(signature_sets_source)

saveRDS(
  gene_sets,
  file.path(paths$analysis, "airway_39_signatures_harmonised.rds"),
  compress = "gzip"
)
data.table::fwrite(
  data.table::rbindlist(signature_mapping_rows),
  file.path(paths$validation, "SIGNATURE_SYMBOL_HARMONISATION_AUDIT.csv")
)
data.table::fwrite(
  symbol_harmonisation,
  file.path(paths$validation, "GLOBAL_SYMBOL_HARMONISATION_AUDIT.csv")
)
rm(all_exported_symbols, all_symbols)
invisible(gc())

prepare_microarray <- function(dataset_name, specification) {
  samples <- data.table::fread(
    specification$samples_file, check.names = FALSE, showProgress = FALSE
  )
  samples[, `Assay ID` := as.character(`Assay ID`)]
  samples[, `Subject ID` := as.character(`Subject ID`)]

  expression <- data.table::fread(
    specification$expression_file,
    select = c("Assay ID", "LOG2E", "PROBE", "GENE SYMBOL"),
    na.strings = c("", "NA", "NaN", "null", "NULL"),
    showProgress = TRUE
  )
  data.table::setnames(
    expression,
    c("Assay ID", "LOG2E", "PROBE", "GENE SYMBOL"),
    c("assay_id", "expression", "probe", "exported_symbol")
  )
  expression[, assay_id := as.character(assay_id)]

  annotation <- unique(expression[, .(probe, exported_symbol)])
  annotation <- annotation[, .(
    exported_symbol = first_nonmissing(exported_symbol)
  ), by = probe]
  annotation[, annotation_probe := sub("_PM_", "_", probe, fixed = TRUE)]

  database_symbol <- AnnotationDbi::mapIds(
    hgu133plus2.db::hgu133plus2.db,
    keys = annotation$annotation_probe,
    keytype = "PROBEID",
    column = "SYMBOL",
    multiVals = "first"
  )
  annotation[, database_symbol := unname(database_symbol)]
  annotation[, selected_source_symbol := data.table::fcoalesce(
    clean_symbols(exported_symbol), clean_symbols(database_symbol)
  )]
  annotation[, gene_symbol := harmonise_symbols(selected_source_symbol)]
  annotation[, mapping_source := data.table::fcase(
    !is.na(clean_symbols(exported_symbol)), "exported_symbol",
    !is.na(clean_symbols(database_symbol)), "hgu133plus2.db_fill",
    default = "unmapped"
  )]
  annotation[, exported_database_conflict :=
    !is.na(clean_symbols(exported_symbol)) &
    !is.na(clean_symbols(database_symbol)) &
    harmonise_symbols(exported_symbol) != harmonise_symbols(database_symbol)]

  gene_by_probe <- stats::setNames(annotation$gene_symbol, annotation$probe)
  expression[, gene_symbol := unname(gene_by_probe[probe])]
  input_rows <- nrow(expression)
  expression <- expression[
    is.finite(expression) & !is.na(gene_symbol) & nzchar(gene_symbol)
  ]

  collapsed <- expression[, .(
    expression = stats::median(expression)
  ), by = .(gene_symbol, assay_id)]
  wide <- data.table::dcast(
    collapsed, gene_symbol ~ assay_id, value.var = "expression"
  )
  genes <- wide$gene_symbol
  wide[, gene_symbol := NULL]
  matrix_object <- as.matrix(wide)
  storage.mode(matrix_object) <- "double"
  rownames(matrix_object) <- genes

  expected_assays <- samples[["Assay ID"]]
  if (!setequal(colnames(matrix_object), expected_assays)) {
    stop(dataset_name, ": matrix/sample assay mismatch.")
  }
  matrix_object <- matrix_object[, expected_assays, drop = FALSE]
  colnames(matrix_object) <- samples[["Subject ID"]]
  if (anyNA(matrix_object) || any(!is.finite(matrix_object))) {
    stop(dataset_name, ": non-finite values after preparation.")
  }

  qc <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    technology = "microarray",
    input_rows = input_rows,
    input_features = nrow(annotation),
    exported_symbol_features = sum(!is.na(clean_symbols(annotation$exported_symbol))),
    database_filled_features = sum(annotation$mapping_source == "hgu133plus2.db_fill"),
    unmapped_features = sum(annotation$mapping_source == "unmapped"),
    exported_database_conflicts = sum(annotation$exported_database_conflict, na.rm = TRUE),
    final_genes = nrow(matrix_object),
    final_samples = ncol(matrix_object),
    minimum_expression = min(matrix_object),
    maximum_expression = max(matrix_object),
    normalization = paste(
      "Supplied LOG2E; blank symbols filled from hgu133plus2.db;",
      "unambiguous aliases harmonised; duplicate probes collapsed by median"
    )
  )
  list(matrix = matrix_object, qc = qc, annotation = annotation)
}

prepare_rnaseq <- function(dataset_name, specification) {
  samples <- data.table::fread(
    specification$samples_file, check.names = FALSE, showProgress = FALSE
  )
  samples[, `Assay ID` := as.character(`Assay ID`)]
  samples[, `Subject ID` := as.character(`Subject ID`)]

  expression <- data.table::fread(
    specification$expression_file,
    select = c("Assay ID", "VALUE", "PROBE", "GENE SYMBOL"),
    na.strings = c("", "NA", "NaN", "null", "NULL"),
    showProgress = TRUE
  )
  data.table::setnames(
    expression,
    c("Assay ID", "VALUE", "PROBE", "GENE SYMBOL"),
    c("assay_id", "count", "probe", "exported_symbol")
  )
  expression[, assay_id := as.character(assay_id)]
  if (any(expression$count < 0, na.rm = TRUE)) {
    stop(dataset_name, ": negative RNA-seq counts.")
  }

  # Library sizes include all exported features, including those not mapped.
  library_sizes <- expression[, .(library_size = sum(count)), by = assay_id]
  annotation <- unique(expression[, .(probe, exported_symbol)])
  annotation <- annotation[, .(
    exported_symbol = first_nonmissing(exported_symbol)
  ), by = probe]
  annotation[, ensembl_id := sub("[.].*$", "", probe)]

  ensembl_map <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(annotation$ensembl_id),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
  ))
  ensembl_map <- data.table::as.data.table(ensembl_map)[
    !is.na(SYMBOL) & nzchar(SYMBOL), .(
      candidate_count = data.table::uniqueN(SYMBOL),
      ensembl_candidates = paste(sort(unique(SYMBOL)), collapse = ";"),
      unique_ensembl_symbol = if (data.table::uniqueN(SYMBOL) == 1L) {
        sort(unique(SYMBOL))[1L]
      } else NA_character_
    ), by = ENSEMBL
  ]
  annotation <- merge(
    annotation, ensembl_map,
    by.x = "ensembl_id", by.y = "ENSEMBL", all.x = TRUE, sort = FALSE
  )
  annotation[, exported_canonical := harmonise_symbols(exported_symbol)]

  # Prefer a unique Ensembl mapping. If Ensembl is ambiguous, use a supplied
  # symbol only when it is one of the Ensembl candidates. Otherwise retain an
  # unambiguous supplied symbol as a documented fallback.
  annotation[, gene_symbol := unique_ensembl_symbol]
  annotation[
    is.na(gene_symbol) & candidate_count > 1L &
      !is.na(exported_canonical) &
      mapply(
        function(symbol, candidates) symbol %in% strsplit(candidates, ";", fixed = TRUE)[[1]],
        exported_canonical, ensembl_candidates
      ),
    gene_symbol := exported_canonical
  ]
  annotation[
    is.na(gene_symbol) & (is.na(candidate_count) | candidate_count == 0L) &
      !is.na(exported_canonical),
    gene_symbol := exported_canonical
  ]
  annotation[, gene_symbol := harmonise_symbols(gene_symbol)]
  annotation[, mapping_status := data.table::fcase(
    !is.na(unique_ensembl_symbol), "unique_ensembl",
    candidate_count > 1L & !is.na(gene_symbol), "ambiguous_ensembl_resolved_by_export",
    (is.na(candidate_count) | candidate_count == 0L) & !is.na(gene_symbol), "exported_symbol_fallback",
    candidate_count > 1L & is.na(gene_symbol), "ambiguous_ensembl_unresolved",
    default = "unmapped"
  )]
  annotation[, export_ensembl_conflict :=
    !is.na(unique_ensembl_symbol) & !is.na(exported_canonical) &
    harmonise_symbols(unique_ensembl_symbol) != exported_canonical]

  gene_by_probe <- stats::setNames(annotation$gene_symbol, annotation$probe)
  expression[, gene_symbol := unname(gene_by_probe[probe])]
  input_rows <- nrow(expression)
  mapped_expression <- expression[
    is.finite(count) & !is.na(gene_symbol) & nzchar(gene_symbol)
  ]
  mapped_counts_by_assay <- mapped_expression[, .(
    mapped_count = sum(count)
  ), by = assay_id]

  collapsed <- mapped_expression[, .(count = sum(count)), by = .(
    gene_symbol, assay_id
  )]
  wide <- data.table::dcast(
    collapsed, gene_symbol ~ assay_id, value.var = "count", fill = 0
  )
  genes <- wide$gene_symbol
  wide[, gene_symbol := NULL]
  count_matrix <- as.matrix(wide)
  storage.mode(count_matrix) <- "double"
  rownames(count_matrix) <- genes

  expected_assays <- samples[["Assay ID"]]
  missing_assays <- setdiff(expected_assays, colnames(count_matrix))
  if (length(missing_assays)) {
    stop(dataset_name, ": mapped matrix missing assays: ",
         paste(missing_assays, collapse = ", "))
  }
  count_matrix <- count_matrix[, expected_assays, drop = FALSE]
  library_sizes <- library_sizes[match(expected_assays, assay_id), library_size]
  if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop(dataset_name, ": invalid total library size.")
  }

  cpm <- sweep(count_matrix, 2L, library_sizes / 1e6, "/")
  minimum_samples <- ceiling(
    preprocessing_settings$rnaseq_minimum_sample_fraction * ncol(cpm)
  )
  keep <- rowSums(cpm >= preprocessing_settings$rnaseq_cpm_threshold) >=
    minimum_samples
  matrix_object <- log2(cpm[keep, , drop = FALSE] + 1)
  colnames(matrix_object) <- samples[["Subject ID"]]
  if (anyNA(matrix_object) || any(!is.finite(matrix_object))) {
    stop(dataset_name, ": non-finite normalized RNA-seq values.")
  }

  mapped_fraction <- merge(
    data.table::data.table(assay_id = expected_assays, library_size = library_sizes),
    mapped_counts_by_assay, by = "assay_id", all.x = TRUE
  )
  mapped_fraction[, mapped_count_fraction := mapped_count / library_size]

  qc <- data.table::data.table(
    dataset = dataset_name,
    dataset_label = specification$label,
    technology = "rnaseq",
    input_rows = input_rows,
    input_features = nrow(annotation),
    exported_symbol_features = sum(!is.na(clean_symbols(annotation$exported_symbol))),
    unique_ensembl_features = sum(annotation$mapping_status == "unique_ensembl"),
    exported_fallback_features = sum(annotation$mapping_status == "exported_symbol_fallback"),
    ambiguous_ensembl_resolved = sum(
      annotation$mapping_status == "ambiguous_ensembl_resolved_by_export"
    ),
    ambiguous_ensembl_unresolved = sum(
      annotation$mapping_status == "ambiguous_ensembl_unresolved"
    ),
    unmapped_features = sum(annotation$mapping_status == "unmapped"),
    export_ensembl_conflicts = sum(annotation$export_ensembl_conflict, na.rm = TRUE),
    genes_before_cpm_filter = nrow(count_matrix),
    final_genes = nrow(matrix_object),
    final_samples = ncol(matrix_object),
    cpm_threshold = preprocessing_settings$rnaseq_cpm_threshold,
    minimum_samples_for_cpm = minimum_samples,
    minimum_mapped_count_fraction = min(mapped_fraction$mapped_count_fraction, na.rm = TRUE),
    median_mapped_count_fraction = stats::median(mapped_fraction$mapped_count_fraction, na.rm = TRUE),
    normalization = paste0(
      "Sparse raw counts; total-library CPM; CPM>=",
      preprocessing_settings$rnaseq_cpm_threshold, " in >=", minimum_samples,
      " participants; log2(CPM+1)"
    )
  )

  list(
    matrix = matrix_object,
    qc = qc,
    annotation = annotation,
    library_audit = mapped_fraction
  )
}

qc_rows <- list()
for (dataset_name in names(datasets)) {
  specification <- datasets[[dataset_name]]
  message("Preparing ", specification$label, "...")
  prepared <- if (specification$technology == "microarray") {
    prepare_microarray(dataset_name, specification)
  } else {
    prepare_rnaseq(dataset_name, specification)
  }

  saveRDS(
    prepared$matrix,
    file.path(paths$processed, paste0(dataset_name, "_expression_gene_level.rds")),
    compress = "gzip"
  )
  data.table::fwrite(
    prepared$annotation,
    file.path(paths$validation, paste0(dataset_name, "_feature_mapping_audit.csv.gz"))
  )
  if (!is.null(prepared$library_audit)) {
    data.table::fwrite(
      prepared$library_audit,
      file.path(paths$validation, paste0(dataset_name, "_library_size_mapping_audit.csv"))
    )
  }
  qc_rows[[dataset_name]] <- prepared$qc
  rm(prepared)
  invisible(gc())
}

expression_qc <- data.table::rbindlist(qc_rows, fill = TRUE)
data.table::fwrite(
  expression_qc,
  file.path(paths$tables, "Table_S_expression_preparation_and_mapping_QC.csv")
)

expected_n <- dataset_metadata_table()[, .(dataset, expected_n)]
observed_n <- expression_qc[, .(dataset, observed_n = final_samples)]
n_check <- merge(expected_n, observed_n, by = "dataset", all = TRUE)
if (any(n_check$expected_n != n_check$observed_n)) {
  stop("Final expression sample count differs from source metadata.")
}

writeLines(
  c(
    "EXPRESSION PREPARATION: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste("Datasets prepared:", nrow(expression_qc)),
    paste("Final genes:", paste(
      paste0(expression_qc$dataset, "=", expression_qc$final_genes),
      collapse = "; "
    )),
    "Microarray input: supplied LOG2E; duplicate probes collapsed by median.",
    "RNA-seq input: sparse raw counts; zeros reconstructed; total-library CPM; log2(CPM+1).",
    "Source signature file remained unchanged; annotation harmonisation is fully audited."
  ),
  file.path(paths$validation, "EXPRESSION_PREPARATION_STATUS.txt")
)
message("Expression preparation completed successfully.")
