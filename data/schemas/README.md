# Controlled-input schemas

This directory documents the structure expected from the authorised U-BIOPRED
exports without distributing any study data. It contains field names, data
types and validation rules only.

The full analysis starts from the processed expression and clinical exports
available to the project. It does not start from CEL, FASTQ, BAM or CRAM files,
which were not available.

The schema files are deliberately separated from `tests/synthetic/`. The
synthetic fixture uses invented identifiers and simulated values and is only a
software demonstration. It must not be interpreted as U-BIOPRED data.

The fixed 39-signature membership is not included here pending confirmation of
redistribution rights. Only the signature-file structure, public pathway names,
signature sizes and provenance metadata may be released at this stage.
