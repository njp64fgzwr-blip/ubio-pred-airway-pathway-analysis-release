# Quantifying airway pathway activity and disease heterogeneity in severe asthma

This repository contains the R/RStudio analysis code and reproducibility materials for an MSc dissertation examining continuous airway pathway activity in adult U-BIOPRED transcriptomic datasets.

The analysis scores **39 prespecified airway- and asthma-relevant signatures with ssGSEA** and relates those scores to five clinical outcomes across four airway datasets. It is a continuous pathway-activity analysis, not primarily a clustering or endotype-discovery analysis.

## Reproducibility boundary

This repository **does not contain original or processed U-BIOPRED data, machine-readable participant records, or U-BIOPRED-derived figures and result tables**. Raw CEL, FASTQ, BAM and CRAM files were not available for the dissertation. The reproducible scientific workflow therefore begins with separately authorised processed U-BIOPRED expression, sample, platform and clinical TSV exports.

There are two supported review routes:

1. **Review without controlled data:** inspect the complete code, fixed parameters, figure legends, source/output checksum inventories and documentation, then run the wholly synthetic demonstration.
2. **Full reproduction with authorised data:** connect a separately authorised copy of the processed U-BIOPRED exports and the frozen 39-signature file, then rerun the complete analysis and validation workflow.

The internal repeated nested cross-validation reported here is **not external validation**. The repository does not claim clinical utility, prognosis or universal superiority of pathway-based models.

The real dissertation figures and result tables remain only in the authorised local project. They can be added later if institutional and U-BIOPRED release approval is documented.

## Analysis at a glance

### Airway datasets

| Key | Dataset | Platform | Profiles |
|---|---|---:|---:|
| `sputum` | Sputum | GPL570 microarray | 120 |
| `brushing_gpl570` | Bronchial brushing | GPL570 microarray | 149 |
| `biopsy` | Bronchial biopsy | GPL570 microarray | 108 |
| `brushing_rnaseq` | Bronchial brushing | RNA sequencing | 118 |

The 495 profiles represent 243 partially overlapping participants. Datasets are analysed separately except for prespecified matched-participant concordance analyses.

### Clinical outcomes

- FEV1 % predicted, analysed without transformation.
- FeNO, analysed as `log1p`.
- Blood eosinophils, analysed as `log1p`.
- Sputum eosinophils, analysed as `log1p`.
- Exacerbations during the 12 months before baseline, analysed as `log1p`.

The exacerbation measure is retrospective burden, not a prospective-risk endpoint.

### Five model strategies

- Biomarker-only.
- Curated-pathway-only.
- Pathway-only Elastic Net.
- Biomarker-plus-curated.
- Biomarker-plus-Elastic-Net.

Models use common complete cases and matched outer test splits within each dataset–outcome comparison. Validation comprises 10 repetitions of five-fold outer cross-validation, with five inner folds for Elastic Net tuning. All transformations, scaling, tuning and model selection are confined to training data. Negative cross-validated R² values are retained.

## Repository map

| Location | Purpose |
|---|---|
| `R/` | Numbered scientific analysis and output-generation scripts |
| `config/analysis_parameters.yml` | Tracked, fixed dissertation parameters |
| `config/config.example.yml` | Local path template for authorised users |
| `data/` | Input schemas and non-membership signature provenance; no U-BIOPRED records |
| `outputs/` | Explanation of locally generated, Git-ignored controlled outputs |
| `validation/` | Source locks, checksum inventories and repository safety checks |
| `docs/` | Reproduction, methods, provenance, output and data-access guidance |
| `tests/` | Data-free tests and wholly synthetic fixtures |

## Quick start without U-BIOPRED access

Use R 4.5.2. From the repository root:

```sh
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv"); renv::restore(prompt = FALSE)'
Rscript validate_results.R --public
```

The data-free validation command checks R syntax, runs the synthetic demonstration, validates the safe tracked manifests and applies the repository privacy gate. To run only the demonstration, use `Rscript run_demo.R`. It uses invented records and verifies workflow mechanics only; it cannot reproduce the dissertation's scientific estimates.

`Rscript validate_results.R` without a flag also defaults to this data-free public mode. Full scientific validation must be requested explicitly with `--full`.

The mapping from dissertation objects to generating scripts and local outputs is given in [`docs/FIGURE_TABLE_PROVENANCE.md`](docs/FIGURE_TABLE_PROVENANCE.md); complete scientific legends are in [`docs/FIGURE_LEGENDS.md`](docs/FIGURE_LEGENDS.md).

## Full run for an authorised user

Do not copy controlled data into this repository. Either set the input locations with environment variables:

```sh
export UBIOPRED_SOURCE_ROOT="/authorised/location/processed-adult-export"
export AIRWAY_SIGNATURE_FILE="/authorised/location/frozen-39-signatures.tsv"
# Optional: override the tracked, gene-membership-free provenance metadata.
# export AIRWAY_SIGNATURE_PROVENANCE_FILE="/authorised/location/byte-identical-provenance.csv"
# Optional if the clinical file is not beneath the source root:
export UBIOPRED_CLINICAL_FILE="/authorised/location/data_clinical.tsv"

Rscript run_analysis.R
Rscript run_outputs.R
Rscript validate_results.R --full
```

or create a local ignored configuration:

```sh
cp config/config.example.yml config/config.yml
# Edit only the local paths in config/config.yml.
export UBIOPRED_ANALYSIS_CONFIG=config/config.yml

Rscript run_analysis.R
Rscript run_outputs.R
Rscript validate_results.R --full
```

The full run is computationally intensive because it includes repeated nested cross-validation. Derived participant-level files are written only beneath the ignored local output directory.

See [`docs/REPRODUCTION_GUIDE.md`](docs/REPRODUCTION_GUIDE.md) for the complete procedure and expected checks.

## Scientific interpretation

- Primary heatmaps show saved normalised ssGSEA scores with common display limits of −1 to +1.
- Complementary z-standardised heatmaps show within-pathway relative variation; they are not absolute pathway activity.
- Spearman analyses are marginal associations, not causal or predictive effects.
- Curated-model coefficients are conditional on the other prespecified pathways in the same model.
- Elastic Net selection frequency is a stability descriptor, not a p value or causal-importance measure.
- Matched-sample concordance is not independent replication.
- Pearson correlation, cross-validated R² and calibration answer different questions and are reported separately.

Detailed fixed parameters are recorded in [`docs/METHODS_AND_PARAMETERS.md`](docs/METHODS_AND_PARAMETERS.md) and `config/analysis_parameters.yml`.

## Data access and licensing

The U-BIOPRED data, participant-level derivatives and exact frozen signature membership are not distributed here. See [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md), [`data/README.md`](data/README.md) and [`LICENSES_AND_DATA_RESTRICTIONS.md`](LICENSES_AND_DATA_RESTRICTIONS.md).

Original repository code and documentation are released under the MIT License. That licence does not grant rights to U-BIOPRED data, third-party gene-set content, publications or other controlled material.

## Citation

Please use the metadata in [`CITATION.cff`](CITATION.cff) when citing this software archive. When referring to a scientific method or source signature, also cite the corresponding primary publication listed in the dissertation.
