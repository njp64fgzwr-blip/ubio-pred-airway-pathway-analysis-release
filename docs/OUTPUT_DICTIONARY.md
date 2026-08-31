# Controlled output dictionary

## Scope

Full authorised runs write results beneath the configured output root,
normally `work/full_analysis/`. This directory is ignored by Git. No real
U-BIOPRED-derived result file is committed to this repository.

`outputs/README.md` records this exclusion. The CSV inventories in
`validation/reference_manifests/` contain only logical filenames, dimensions,
byte sizes, SHA-256 values and release classifications for locked outputs.
They are reconciliation evidence, not the outputs themselves.

## Local output groups

| Local location | Contents | Git status |
|---|---|---|
| `work/full_analysis/results/processed_expression/` | prepared gene-level expression matrices | ignored; controlled |
| `work/full_analysis/results/ssgsea/` | normalised ssGSEA scores and z-standardised derivatives | ignored; controlled |
| `work/full_analysis/results/analysis_data/` | merged clinical/pathway data, PCA, fold records and predictions | ignored; controlled |
| `work/full_analysis/results/tables/` | aggregate scientific and diagnostic tables | ignored pending approval |
| `work/full_analysis/results/models/` | fitted model objects | ignored; controlled |
| `work/full_analysis/figures/` | working R-generated figures | ignored pending approval |
| `work/full_analysis/final_dissertation_figures/` | correctly numbered Figures 1–9 and S1–S12 | ignored pending approval |
| `work/full_analysis/final_dissertation_tables/` | R-generated main Table 1 | ignored pending approval |
| `work/full_analysis/validation/` | run-specific PASS/FAIL records and source checks | ignored; may contain local paths |

## Main result families

The analysis creates the following controlled local table families:

- marginal pathway–outcome Spearman estimates and BH-adjusted q values;
- prespecified curated-model coefficients, confidence intervals and
  diagnostics;
- five-strategy repeated nested-cross-validation performance;
- Elastic Net selection-frequency and tuning summaries;
- participant-averaged out-of-fold calibration;
- exploratory paired incremental-performance comparisons;
- matched-participant compartment and platform concordance;
- RNA-seq signature coverage and pathway-score variation summaries; and
- sample, clinical-outcome and preprocessing summaries for the main tables.

The authoritative mapping from these families to dissertation Figures 1–9,
Supplementary Figures S1–S12 and Tables 1–3 is in
`docs/FIGURE_TABLE_PROVENANCE.md`.

## Material that must remain outside Git

- expression matrices and clinical exports;
- participant-level ssGSEA scores or PCA coordinates;
- complete-case identifiers and sample links;
- cross-validation assignments and participant-level predictions;
- fitted models and serialized R objects;
- exact 39-signature gene membership; and
- any derived output not expressly approved for its intended audience.

An institutional release decision can later authorise selected aggregate
artifacts, but such files should be added in a separate reviewed change rather
than silently mixed into the code archive.
