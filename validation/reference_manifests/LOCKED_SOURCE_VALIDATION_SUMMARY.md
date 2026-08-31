# Locked-source validation summary

Source decision date: 18 August 2026

The repository figure package was assembled only from
`AUTHORITATIVE_DISSERTATION_FIGURES_FROM_DRAFT`. That source records:

- 22 main figure stems and 12 supplementary figure stems;
- four validated production formats per source stem;
- 136/136 source exports matching their audited R sources;
- 43/43 canonical numerical-source checks passing;
- 780/780 Spearman values, 104/104 curated-model coefficient rows,
  780/780 selection-stability rows, 60/60 performance rows and 156/156
  concordance rows passing keyed comparisons;
- no explanatory subtitles, teaching notes or presentation footers in the
  locked scientific figures.

During local assembly, one PDF and one 300-dpi PNG per stem were compared
byte-for-byte with their locked source. Their sizes and SHA-256 values are
retained in `AUTHORITATIVE_FIGURE_PUBLIC_MANIFEST.csv`; the image files are not
committed to GitHub.

The audited aggregate result-table candidates came from the validated RStudio
rebuild. Three reduced-disclosure projections were also checked mechanically:

1. the signature provenance table omits `analysed_gene_symbols`;
2. the signature coverage table omits `detected_genes`;
3. the pathway-variation table omits per-pathway minima and maxima.

Their metadata remain in `PUBLIC_TABLE_MANIFEST.csv`, but none of these real
tables is committed. No numerical value in the locked figures was changed
during repository assembly.
