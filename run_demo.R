# Run a data-free software demonstration on wholly synthetic inputs.
#
# This is a smoke test of file handling, ssGSEA scoring and fold-contained
# prediction. It has no biological meaning and cannot reproduce dissertation
# estimates.

find_entrypoint_root <- function(start = getwd()) {
  candidates <- c(start)
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument)) {
    candidates <- c(
      dirname(normalizePath(sub("^--file=", "", file_argument[1L]))),
      candidates
    )
  }
  for (candidate in unique(candidates)) {
    current <- normalizePath(candidate, mustWork = TRUE)
    repeat {
      if (file.exists(file.path(
        current, "ubio-pred-airway-pathway-analysis.Rproj"
      ))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not locate the repository root.")
}

setwd(find_entrypoint_root())
source(file.path("R", "00_common.R"))
assert_packages(c("data.table", "GSVA", "BiocParallel", "digest"))
source(file.path("R", "09_nested_cv_helpers.R"))

message(paste(
  "SYNTHETIC SOFTWARE DEMONSTRATION ONLY:",
  "no U-BIOPRED value, identifier or signature membership is used."
))

fixture_dir <- file.path(PROJECT_ROOT, "tests", "synthetic")
fixture_files <- file.path(
  fixture_dir,
  c("clinical.tsv", "samples.tsv", "expression.tsv", "signatures.tsv")
)
if (any(!file.exists(fixture_files))) {
  sys.source(
    file.path(fixture_dir, "generate_fixture.R"),
    envir = new.env(parent = globalenv())
  )
}

clinical <- data.table::fread(file.path(fixture_dir, "clinical.tsv"))
samples <- data.table::fread(file.path(fixture_dir, "samples.tsv"))
expression_table <- data.table::fread(file.path(fixture_dir, "expression.tsv"))
signature_table <- data.table::fread(file.path(fixture_dir, "signatures.tsv"))

if (data.table::uniqueN(signature_table$signature_id) != 39L) {
  stop("Synthetic fixture must contain 39 demonstration signatures.")
}
if (anyDuplicated(clinical$subject_id) || anyDuplicated(samples$sample_id)) {
  stop("Synthetic clinical or sample identifiers are not unique.")
}
if (!setequal(samples$subject_id, clinical$subject_id)) {
  stop("Synthetic sample-to-clinical linkage failed.")
}

expression <- as.matrix(expression_table[, -"gene_symbol"])
storage.mode(expression) <- "double"
rownames(expression) <- expression_table$gene_symbol
gene_sets <- split(signature_table$gene_symbol, signature_table$signature_id)
gene_sets <- lapply(gene_sets, unique)

demo_rows <- list()
for (dataset_name in unique(samples$dataset)) {
  dataset_samples <- samples[dataset == dataset_name]
  dataset_expression <- expression[, dataset_samples$sample_id, drop = FALSE]
  parameter <- GSVA::ssgseaParam(
    exprData = dataset_expression,
    geneSets = gene_sets,
    minSize = 3L,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE,
    checkNA = "yes",
    use = "everything",
    verbose = FALSE
  )
  scores <- GSVA::gsva(
    parameter, verbose = FALSE,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  )
  score_table <- data.table::as.data.table(
    data.frame(sample_id = colnames(scores), t(scores), check.names = FALSE)
  )
  score_table <- merge(score_table, dataset_samples, by = "sample_id")
  score_table <- merge(score_table, clinical, by = "subject_id")

  predictor_names <- names(gene_sets)[seq_len(3L)]
  x <- as.matrix(score_table[, ..predictor_names])
  y <- score_table$fev1_percent_predicted
  fold_id <- make_stratified_folds(y, k = 5L, seed = 20260820L)
  prediction <- rep(NA_real_, length(y))
  for (fold in seq_len(5L)) {
    fit <- fit_outer_ols(
      x_train = x[fold_id != fold, , drop = FALSE],
      y_train = y[fold_id != fold],
      x_test = x[fold_id == fold, , drop = FALSE]
    )
    prediction[fold_id == fold] <- fit$prediction
  }
  metrics <- calculate_metrics(y, prediction)
  demo_rows[[dataset_name]] <- data.table::data.table(
    dataset = dataset_name,
    participants = length(y),
    signatures_scored = nrow(scores),
    demonstration_model = "Three synthetic signatures; fivefold OLS",
    cross_validated_r_squared = unname(metrics[["r_squared"]]),
    rmse = unname(metrics[["rmse"]]),
    mae = unname(metrics[["mae"]])
  )
}

demo_results <- data.table::rbindlist(demo_rows)
demo_output_dir <- file.path(PROJECT_ROOT, "work", "demo")
dir.create(demo_output_dir, recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(
  demo_results, file.path(demo_output_dir, "synthetic_demo_metrics.csv")
)
writeLines(
  c(
    "SYNTHETIC SOFTWARE DEMONSTRATION: PASS",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "Inputs are deterministic and wholly synthetic.",
    "The demonstration used 39 invented signatures across four invented datasets.",
    "Reported metrics are smoke-test values without scientific interpretation.",
    "This run does not reproduce or validate dissertation estimates."
  ),
  file.path(demo_output_dir, "SYNTHETIC_DEMO_STATUS.txt")
)
print(demo_results)
message("\nSYNTHETIC SOFTWARE DEMONSTRATION: PASS")
