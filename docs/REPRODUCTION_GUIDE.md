# Reproduction guide

## Purpose

This guide supports two distinct tasks:

- a **data-free repository audit**, which runs on any machine with the required R environment; and
- an **authorised full reproduction**, which additionally requires controlled U-BIOPRED inputs and the frozen 39-signature file.

Do not interpret successful data-free checks as independent scientific replication. Do not describe a successful full rerun on the supplied cohort as external validation.

## 1. Software environment

The validated rebuild used:

- R 4.5.2;
- GSVA 2.4.9;
- `hgu133plus2.db` 3.13.0;
- `org.Hs.eg.db` 3.22.0; and
- `glmnet` 5.0.

Full figure generation and visual QA also require Poppler (`pdftoppm`) and the
system libraries used by the R `magick` package. On macOS, compiling some
locked packages may additionally require the current GNU Fortran toolchain.
The data-free synthetic validation does not render the controlled figures.

All remaining package versions are recorded in `renv.lock`. From the repository root:

```sh
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")'
Rscript -e 'renv::restore(prompt = FALSE)'
```

Package restoration may require access to CRAN and Bioconductor repositories. A version mismatch should be resolved before comparing regenerated hashes or numerical results.

## 2. Data-free repository audit

Run:

```sh
Rscript validate_results.R --public
```

This command includes the wholly synthetic demonstration. To run only that component, use `Rscript run_demo.R`.

The expected sequence is:

1. Generate or load wholly synthetic inputs.
2. Exercise the portable analysis components on those invented records.
3. Run data-free tests.
4. Verify hashes for the tracked schemas, synthetic fixtures, provenance map
   and legends, and validate the withheld-output checksum inventories.
5. Scan the repository for forbidden controlled-data paths, participant-level public tables, personal absolute paths, unsafe file types and unexpectedly large files.
6. End with an explicit pass/fail result.

If the demonstration passes, it establishes that the exercised synthetic code
paths execute in the restored environment. It does not exercise every
preprocessing, nested-validation, inference or figure-generation stage, and it
does not show that synthetic estimates agree with the dissertation.

## 3. Configure an authorised full run

Keep all controlled files outside the cloned repository.

### Environment-variable route

```sh
export UBIOPRED_SOURCE_ROOT="/authorised/location/processed-adult-export"
export AIRWAY_SIGNATURE_FILE="/authorised/location/frozen-39-signatures.tsv"
# Optional: override the tracked, gene-membership-free Figure S12 provenance.
# export AIRWAY_SIGNATURE_PROVENANCE_FILE="/authorised/location/provenance.csv"

# Only needed if data_clinical.tsv is not directly beneath the source root:
export UBIOPRED_CLINICAL_FILE="/authorised/location/data_clinical.tsv"
```

### Local configuration route

```sh
cp config/config.example.yml config/config.yml
```

Edit the local input path fields in `config/config.yml`, then activate it:

```sh
export UBIOPRED_ANALYSIS_CONFIG=config/config.yml
```

`config/config.yml` is ignored by Git. Do not commit it because local paths may reveal controlled-data locations.

The tracked `config/analysis_parameters.yml` records the locked scientific settings. Do not alter the dataset order, expected sample counts, signature checksum, preprocessing, ssGSEA or validation settings when testing exact reproduction.

## 4. Validate controlled inputs

The first part of the full analysis checks:

- all 13 required TSV files and the signature file are present;
- the four expected dataset directories are correctly resolved;
- expression, sample and platform files have compatible schemas and orientation;
- expected profile counts are 120, 149, 108 and 118 in the locked dataset order;
- participant and sample links are internally consistent;
- the signature file contains 39 entries; and
- the signature file SHA-256 matches the locked value.

An input failure is intentional. Resolve the source file or configuration mismatch; do not weaken the check to force execution.

## 5. Run the scientific analysis

```sh
Rscript run_analysis.R
```

This entry point performs, in order:

1. controlled-input audit;
2. GPL570 and RNA-seq expression preparation;
3. ssGSEA scoring of the 39 frozen signatures;
4. clinical-data merging and cohort summaries;
5. descriptive, PCA and marginal Spearman analyses;
6. prespecified curated regressions and sensitivity analyses;
7. matched-participant compartment/platform concordance;
8. five-strategy repeated nested internal validation;
9. independent nested-validation leakage and reconstruction checks; and
10. predictive inference, calibration and sensitivity analyses.

The nested-validation stage is computationally expensive. Runtime depends on processor, memory, package build and filesystem performance.

Participant-level derived files are written under the configured output root, normally `work/full_analysis/`. That directory is ignored by Git.

## 6. Generate figures and tables

After a successful full analysis:

```sh
Rscript run_outputs.R
```

This entry point regenerates the correctly numbered dissertation figures and
R-rendered Table 1 from the validated local analysis outputs. Figures 1–9 and
S1–S11 are produced by the RStudio rebuild exporters. Figure S12 and safe
signature-provenance summaries are regenerated from the authorised local
signature membership and the tracked provenance map. Other aggregate analysis
tables are produced during `run_analysis.R`. These controlled artifacts remain
in the ignored local output root.

Figures are generated from explicit plotting-source tables before export. The output stage also records manifests and checks that the final visual sources correspond to the dissertation numbering and source lock.

## 7. Validate the complete run

```sh
Rscript validate_results.R --full
```

With authorised local outputs available, validation should check:

- package and input lineage;
- exact signature identity and coverage;
- expected dataset, outcome and model combinations;
- 10 repeats and five outer folds where required;
- absence of preprocessing or tuning leakage;
- retention of negative cross-validated R² values;
- reconstructability of reported metrics from saved out-of-fold results;
- numerical correspondence between plot-source tables and figures;
- expected figure/table inventory and format.

The run is complete only when required checks end in `PASS`. A failed check
should remain visible and be investigated; do not replace a failed output with
an older archived value. `Rscript validate_results.R --public` is a separate
data-free software and repository check; it does not validate the controlled
scientific results.

## 8. Comparing with the dissertation lock

Use `docs/FIGURE_TABLE_PROVENANCE.md` to map each dissertation object to its
generating stage. Compare regenerated files with authorised local copies and
the checksum inventories rather than relying on filenames alone.

A byte-level difference does not automatically indicate a scientific difference: PDF metadata and graphics-device versions can change bytes. When hashes differ, compare the underlying plotting data, reported values, dimensions and rendered appearance. Conversely, a visually similar figure is not sufficient when its numerical source differs.

## 9. Expected and non-error conditions

- Negative cross-validated R² values are valid and must not be deleted.
- Missing output combinations should be assessed against complete-case availability rather than silently imputed.
- Low proportional RNA-seq coverage for a signature is an interpretation warning, not permission to redefine the frozen set.
- Cross-platform score concordance concerns participant ranking, not equality of absolute score scales.
- The primary heatmap limit of −1 to +1 is a display convention for saved normalised scores; z-standardised figures answer a different descriptive question.

## 10. Minimum reporting record

For a rerun, retain locally:

- Git commit identifier;
- active configuration with sensitive paths redacted;
- R session information and `renv.lock` state;
- source-input and signature checksums;
- script checksums;
- validation logs;
- controlled-output manifest; and
- a note explaining any environment-dependent non-scientific rendering differences.

Do not publish the active configuration, participant-level logs or controlled-data manifests without a separate disclosure review.
