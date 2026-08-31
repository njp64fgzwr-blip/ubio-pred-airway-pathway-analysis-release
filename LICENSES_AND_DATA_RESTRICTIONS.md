# Licences and data restrictions

## Repository code and documentation

The original source code and documentation written for this repository are available under the MIT License in `LICENSE`.

The MIT License applies only to material for which the repository author is entitled to grant that licence. It does not override third-party rights, controlled-data agreements, research-governance requirements or the licences of R packages used by the workflow.

## U-BIOPRED material

No original or processed U-BIOPRED data, machine-readable participant records, or U-BIOPRED-derived figures and result tables are distributed through this repository. This exclusion covers, among other material:

- processed expression matrices;
- clinical and sample metadata;
- participant, sample and assay identifiers;
- participant-level ssGSEA scores and PCA coordinates;
- cross-validation folds and participant-level predictions;
- fitted model objects that may retain controlled information; and
- intermediate files from which individual measurements could be recovered.

Possession of this code does not confer permission to obtain, use or redistribute U-BIOPRED data. Users must obtain access independently through the applicable U-BIOPRED data-governance route and comply with their own approval and data-use terms.

## Dissertation outputs

Real dissertation figures and result tables remain in the authorised local analysis environment and are ignored by Git. Some locked heatmaps and PCA graphics visibly encode de-identified participant-level patterns, and some aggregate tables contain small cells. Distribution requires prior approval through the applicable institutional and U-BIOPRED governance routes.

## Pathway-signature membership

The analysis used one frozen local file containing the exact membership of 39 curated signatures. Some lists originate from, reproduce or derive from third-party publications and pathway resources. Redistribution permission for all exact memberships has not been established; consequently, the controlled membership file is not included.

The repository instead provides:

- the expected signature count;
- the SHA-256 checksum of the exact frozen input;
- reader-facing labels, sizes and provenance summaries where redistribution is appropriate;
- per-dataset detected-gene coverage without individual gene membership; and
- a configuration field through which an authorised user can supply the exact local file.

This separation preserves computational identity without asserting redistribution rights.

## Third-party software and publications

R, Bioconductor, CRAN packages and cited publications remain subject to their own licences and terms. No licence to reproduce article text, supplementary files or proprietary database content is granted here.

## Questions before public release

Before adding any U-BIOPRED-derived output or exact signature membership, the owner should obtain approval through the applicable institutional and U-BIOPRED governance routes. The automated repository check is an additional technical safeguard; it is not a substitute for institutional approval.
