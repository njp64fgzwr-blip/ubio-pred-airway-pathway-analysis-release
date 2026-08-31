# Fully synthetic demonstration fixture

Every identifier and every numeric value in this directory is generated and
invented. Nothing was sampled, perturbed, transformed or copied from U-BIOPRED.

Run `Rscript tests/synthetic/generate_fixture.R` from the repository root to
recreate the four fixture files deterministically:

- `clinical.tsv`: 20 invented participants and simulated clinical variables;
- `samples.tsv`: 80 invented samples across four demonstration datasets;
- `expression.tsv`: a simulated 120-gene by 80-sample matrix;
- `signatures.tsv`: 39 demonstration signatures using synthetic gene names.

The fixture is intentionally small and has no scientific meaning. It is
provided to exercise input validation and software flow without controlled
study data. The demonstration signature memberships are not the memberships
used in the dissertation.
