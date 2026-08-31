# Data access and availability

## Repository boundary

This GitHub repository is deliberately code-only. It does not contain original
or processed participant-level U-BIOPRED data, identifiers, expression
matrices, clinical records, pathway scores, fold assignments, predictions,
fitted models, or U-BIOPRED-derived figures and result tables.

These materials originate from a controlled human-research resource.
De-identification alone does not establish permission to redistribute them,
and some figures or small aggregate cells may still carry disclosure risk.
Private GitHub storage is not treated as a substitute for U-BIOPRED and
institutional approval.

## Reproducibility starting point

Raw CEL, FASTQ, BAM and CRAM files were unavailable for the dissertation. The
earliest inputs were validated processed adult U-BIOPRED TSV exports. Exact
scientific reproduction therefore begins at that processed-data boundary; the
work should not be described as an instrument-level raw-data reanalysis.

## What can be verified without controlled data

An independent reviewer can:

- inspect every scientific analysis, output-generation and validation script;
- inspect the fixed parameters and model definitions;
- restore the recorded R environment;
- run a wholly synthetic demonstration of the packaged workflow;
- run syntax, integrity and privacy checks;
- inspect input schemas, scientific figure legends and the mapping from every
  dissertation output to its generating stage;
- inspect SHA-256 inventories describing the locked controlled inputs and
  withheld dissertation outputs; and
- verify that no controlled study data are tracked by Git.

The synthetic demonstration tests software mechanics and expected object
structure. It does not reproduce or validate the dissertation's biological
findings.

## What is required for exact rerunning

An authorised user must independently supply:

1. the processed adult U-BIOPRED source tree containing the 13 locked clinical,
   expression, sample and platform TSV files;
2. the exact frozen 39-signature input whose SHA-256 is recorded in the
   configuration; and
3. an R environment restored from `renv.lock`.

The repository includes gene-membership-free signature provenance metadata
needed for Figure S12. An authorised alternative can be supplied only if its
checksum and fields agree with the source lock.

The user must also have permission to use the controlled files under the
applicable U-BIOPRED and institutional governance arrangements. This
repository neither grants access nor changes existing data-use conditions.

## Withheld outputs

Real dissertation figures and result tables remain in the authorised local
project and are ignored by Git. Their filenames, sizes and SHA-256 values are
recorded as non-content inventories so an authorised rerun can be reconciled
with the dissertation lock. The inventories do not contain the underlying
results.

The exact 39-signature membership is also excluded because the library
incorporates lists from multiple publications and pathway resources, and
redistribution rights for every inherited list have not been established.

## Suggested data-availability statement

> Analysis code, dependency metadata, fixed parameters, input schemas,
> checksums and a wholly synthetic software demonstration are available in the
> accompanying repository. Individual-level U-BIOPRED clinical and
> transcriptomic data, linked derivatives and dissertation result artifacts
> are not redistributed because they are subject to controlled-access
> governance. Exact rerunning begins from the validated processed U-BIOPRED
> exports and requires independently authorised access to those files and the
> frozen 39-signature input. Raw CEL and sequencing-read files were unavailable
> for this analysis.

The final wording should follow the applicable institutional and U-BIOPRED
governance requirements.
