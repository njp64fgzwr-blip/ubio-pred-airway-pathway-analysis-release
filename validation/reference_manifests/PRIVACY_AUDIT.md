# Repository privacy and disclosure audit

Audit updated: 31 August 2026

## Result

**PASS for a code-only repository.** Release remains subject to normal
institutional review, but no original, processed or derived U-BIOPRED result
file is tracked by Git.

## Tracked content

- portable R analysis, output-generation and validation code;
- fixed non-secret parameters and an R dependency lock;
- input schemas and a wholly synthetic fixture using only `DEMO_*` records;
- a 39-row provenance map without exact gene membership;
- scientific figure legends and dissertation-output mappings; and
- non-content inventories containing logical paths, dimensions, byte sizes,
  SHA-256 values and release classifications.

## Explicit exclusions

- clinical, expression, sample and platform exports;
- participant, sample or assay identifiers;
- processed expression and participant-level ssGSEA matrices;
- PCA coordinates, complete-case lists, fold assignments and predictions;
- fitted models and R serialisation files;
- real dissertation figures and result tables;
- dissertation drafts, review files, EndNote libraries and literature;
- exact frozen or harmonised gene membership for the 39 signatures; and
- active configurations containing controlled local paths.

## Safeguards

- `.gitignore` blocks controlled filenames, sensitive formats, exact-membership
  patterns and all generated outputs except `outputs/README.md`.
- The output resolver permits in-repository full-run outputs only below the
  ignored `work/` directory.
- `validation/check_public_release.R` scans the tracked Git set for forbidden
  paths, formats, symbolic links, personal absolute paths and oversized files.
- `validate_results.R --public` operates without controlled data and verifies
  tracked safe-file hashes plus the structure of withheld-output inventories.

These controls are safeguards, not permission to add a U-BIOPRED-derived
artifact. Any later addition of a real figure, table or signature membership
requires a separate, documented institutional and data-owner decision.
