# Generate a deterministic, wholly synthetic fixture for software testing.
# No value or identifier is copied or derived from U-BIOPRED.

set.seed(20260820L)

output_dir <- file.path("tests", "synthetic")
if (!dir.exists(output_dir)) {
  stop("Run this script from the repository root.")
}

n_participants <- 20L
dataset_keys <- c(
  "sputum", "brushing_gpl570", "biopsy", "brushing_rnaseq"
)
participant_ids <- sprintf("DEMO_P%03d", seq_len(n_participants))

latent_type2 <- stats::rnorm(n_participants)
latent_innate <- stats::rnorm(n_participants)

clinical <- data.frame(
  subject_id = participant_ids,
  age = round(stats::runif(n_participants, 20, 70), 1),
  sex = rep(c("Female", "Male"), length.out = n_participants),
  smoking_status = rep(c("Never", "Former", "Current", "Never"),
                       length.out = n_participants),
  bmi = round(stats::rnorm(n_participants, 27, 4), 1),
  fev1_percent_predicted = round(78 - 7 * latent_innate +
                                  stats::rnorm(n_participants, 0, 6), 1),
  feno_ppb = round(exp(3.0 + 0.35 * latent_type2 +
                       stats::rnorm(n_participants, 0, 0.25)) - 1, 1),
  blood_eosinophils_x10_3_uL = round(exp(-1.4 + 0.25 * latent_type2 +
                                          stats::rnorm(n_participants, 0, 0.2)),
                                        3),
  sputum_eosinophils_percent = round(exp(1.2 + 0.50 * latent_type2 +
                                          stats::rnorm(n_participants, 0, 0.3)) - 1,
                                        1),
  exacerbations_previous_12m = pmax(
    0L, stats::rpois(n_participants, exp(0.5 + 0.20 * latent_innate))
  ),
  stringsAsFactors = FALSE
)

samples <- do.call(rbind, lapply(seq_along(dataset_keys), function(dataset_i) {
  data.frame(
    sample_id = sprintf("DEMO_S%03d", (dataset_i - 1L) * n_participants +
                          seq_len(n_participants)),
    subject_id = participant_ids,
    dataset = dataset_keys[dataset_i],
    stringsAsFactors = FALSE
  )
}))

n_genes <- 120L
genes <- sprintf("DEMO_G%03d", seq_len(n_genes))
expression <- matrix(
  stats::rnorm(n_genes * nrow(samples), mean = 0, sd = 1),
  nrow = n_genes,
  dimnames = list(genes, samples$sample_id)
)

participant_index <- match(samples$subject_id, participant_ids)
dataset_effect <- match(samples$dataset, dataset_keys) * 0.08
expression[1:20, ] <- expression[1:20, ] +
  rep(latent_type2[participant_index], each = 20L) * 0.45
expression[21:40, ] <- expression[21:40, ] +
  rep(latent_innate[participant_index], each = 20L) * 0.45
expression <- sweep(expression, 2L, dataset_effect, "+")

expression_table <- data.frame(
  gene_symbol = rownames(expression), expression,
  check.names = FALSE, stringsAsFactors = FALSE
)

signature_rows <- do.call(rbind, lapply(seq_len(39L), function(signature_i) {
  start <- ((signature_i - 1L) * 3L) %% (n_genes - 7L) + 1L
  data.frame(
    signature_id = sprintf("DEMO_SIGNATURE_%02d", signature_i),
    gene_order = seq_len(8L),
    gene_symbol = genes[start:(start + 7L)],
    stringsAsFactors = FALSE
  )
}))

write.table(
  clinical, file.path(output_dir, "clinical.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE, na = ""
)
write.table(
  samples, file.path(output_dir, "samples.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE, na = ""
)
write.table(
  expression_table, file.path(output_dir, "expression.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE, na = ""
)
write.table(
  signature_rows, file.path(output_dir, "signatures.tsv"), sep = "\t",
  row.names = FALSE, quote = FALSE, na = ""
)

message("Synthetic fixture written to ", output_dir)
