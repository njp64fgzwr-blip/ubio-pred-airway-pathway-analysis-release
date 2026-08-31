# Data directory

This directory contains only non-sensitive input schemas and
gene-membership-free signature provenance. It does not contain U-BIOPRED
participant data. Wholly synthetic demonstration fixtures are stored under
`tests/synthetic/`, and validation manifests are stored under
`validation/reference_manifests/`.

## Controlled inputs used by the full workflow

The full analysis expects an authorised local source root with the following relative structure:

```text
processed-adult-export/
├── data_clinical.tsv
└── OMICS+DATA/
    └── Transcriptomics/
        └── Baseline+Visit/
            ├── Sputum/
            │   └── GPL570_SPUTUM/
            │       ├── data_mrna.tsv
            │       ├── samples.tsv
            │       └── platform.tsv
            ├── Bronchial+Brushings/
            │   ├── GPL570_BRBRUSHING/
            │   │   ├── data_mrna.tsv
            │   │   ├── samples.tsv
            │   │   └── platform.tsv
            │   └── BRONCHIAL_BRUSHING_RNASEQ/
            │       ├── data_mrna.tsv
            │       ├── samples.tsv
            │       └── platform.tsv
            └── Bronchial+Biopsy/
                └── GPL570_BIOPSY/
                    ├── data_mrna.tsv
                    ├── samples.tsv
                    └── platform.tsv
```

The frozen 39-signature tab-delimited file is configured separately. It is not
stored here because redistribution rights for every inherited gene set have
not been established. Figure S12 uses the tracked audited provenance CSV under
`data/signatures/`, which contains no gene membership; its checksum is fixed in
the analysis configuration.

Input paths are supplied through environment variables or a local ignored `config/config.yml`; see `config/config.example.yml`. Controlled files should never be copied into the repository.

## Supporting locations

- `schemas/` describes expected input structure and required fields without containing participant records.
- `../tests/synthetic/` contains invented values generated solely to exercise the data-free demonstration.
- `../validation/reference_manifests/` contains non-identifying file inventories, checksums and repository safety records.

Synthetic records are not transformed, sampled or perturbed U-BIOPRED participants. Results produced from them are software checks and have no biological interpretation.

## Expected datasets and checks

The input audit expects four datasets in the fixed order below:

| Dataset key | Profiles | Expression technology |
|---|---:|---|
| `sputum` | 120 | GPL570 microarray |
| `brushing_gpl570` | 149 | GPL570 microarray |
| `biopsy` | 108 | GPL570 microarray |
| `brushing_rnaseq` | 118 | RNA sequencing |

The audit checks file presence, schema, sample counts, identifier relationships, expression orientation and the checksum and count of the frozen signature input before analysis proceeds.

## Files that must remain outside Git

- all expression and clinical exports;
- all sample, assay and participant identifiers;
- participant-level pathway scores and PCA coordinates;
- complete-case membership and cross-validation fold assignments;
- individual observed or predicted outcome values;
- fitted models and serialized R objects; and
- any derived table from which participant-level information could be reconstructed.

The repository's `.gitignore` and `validation/check_public_release.R` provide technical safeguards, but the person publishing the repository remains responsible for compliance with the applicable data-use agreement and institutional guidance.
