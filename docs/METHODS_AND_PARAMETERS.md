# Methods and locked analysis parameters

This document records the computational specification represented by the dissertation and the validated RStudio rebuild. `config/analysis_parameters.yml` is the machine-readable companion. Changes to these parameters constitute a different analysis and should not be described as exact reproduction.

## Study boundary

The study is a cross-sectional secondary analysis of de-identified adult U-BIOPRED baseline data. Raw microarray CEL files and raw sequencing reads were unavailable. The workflow begins from validated processed expression, sample, platform and clinical TSV exports.

Four partially overlapping datasets are analysed separately:

| Key | Compartment and platform | n | Final genes |
|---|---|---:|---:|
| `sputum` | Sputum GPL570 microarray | 120 | 22,090 |
| `brushing_gpl570` | Bronchial brushing GPL570 microarray | 149 | 22,085 |
| `biopsy` | Bronchial biopsy GPL570 microarray | 108 | 22,085 |
| `brushing_rnaseq` | Bronchial brushing RNA-seq | 118 | 16,691 |

The 495 profiles represent 243 unique participants. Separate dataset analyses avoid treating repeated profiles from the same participant as independent cross-compartment replication. Matched participants enter only the prespecified concordance analyses.

## Expression preparation

### GPL570 microarrays

- Each supplied matrix contains 54,715 probe features.
- Supplied non-blank symbols are retained as the starting annotation.
- Of 13,355 probes with blank symbols, `hgu133plus2.db` 3.13.0 supplies symbols for 4,362; the remaining 8,993 unmapped probes are excluded.
- Unambiguous aliases are harmonised.
- Probes assigned to the same final symbol are collapsed by the within-participant median.
- The supplied expression scale is retained; no cross-dataset batch correction or merged expression matrix is created.

### Bronchial-brushing RNA sequencing

- Zeros omitted from the sparse export are restored.
- Ensembl identifiers are mapped with `org.Hs.eg.db` 3.22.0, preferring a unique Ensembl-to-symbol mapping and using a compatible supplied symbol where the Ensembl mapping is ambiguous or absent.
- Counts assigned to the same final symbol are summed.
- A gene is retained when counts per million are at least 1 in at least 12 of 118 participants.
- Retained counts are transformed as `log2(CPM + 1)`.

## Frozen signature library and ssGSEA

- The feature library contains exactly 39 prespecified airway- and asthma-relevant signatures.
- Membership was frozen before clinical-outcome analysis and was not changed in response to the results.
- The expected SHA-256 of the exact local membership file is `00aaf4a958729605fd2005294b36f944699067f10cc3e160809f21ac12c8a3d8`.
- Exact membership is not redistributed in this repository; see `LICENSES_AND_DATA_RESTRICTIONS.md`.
- Detected-gene coverage is checked for each dataset, with at least three detected genes required for scoring.

ssGSEA is implemented in GSVA 2.4.9 under R 4.5.2 using `ssgseaParam()` and `gsva()` with:

| Parameter | Locked value |
|---|---|
| minimum gene-set size | 3 |
| maximum gene-set size | infinite |
| alpha | 0.25 |
| normalisation | enabled |
| parallel backend | serial |

The method should be described scientifically as **ssGSEA pathway scoring**, not as a generic GSVA analysis.

Primary heatmaps display the saved normalised ssGSEA values with common colour limits of −1 to +1. The source values are not rescaled merely to fill those limits. Complementary z-standardised displays centre and scale each pathway across participants within a dataset; they show within-pathway relative variation and are not absolute pathway activity.

PCA uses the 39 within-dataset z-standardised scores and is descriptive. It does not test for, prove or exclude discrete endotypes.

## Clinical outcomes

| Key | Interpretation | Unit | Modelling transform |
|---|---|---|---|
| `FEV1` | supplied FEV1 % predicted | % predicted | identity |
| `FeNO` | fractional exhaled nitric oxide | ppb | `log1p` |
| `blood_eosinophils` | circulating eosinophils | ×10³/µL | `log1p` |
| `sputum_eosinophils` | non-squamous differential | % | `log1p` |
| `exacerbation_frequency` | events during the 12 months before baseline | events/year | `log1p` |

No outcome or covariate is imputed. The supplied exacerbation count represents retrospective baseline burden rather than future risk. The supplied FEV1 percentage is analysed as provided because its reference equation and pre/post-bronchodilator status could not be reconstructed confidently from the export.

## Marginal associations

- Spearman correlations use pairwise-complete observations within each dataset and pathway–outcome pair.
- Benjamini–Hochberg adjustment is performed across 39 pathways within each dataset–outcome family.
- `q < 0.05` denotes false-discovery-rate support in that defined family.
- These are marginal associations, not conditional, causal or predictive effects.

## Prespecified curated models

Ordinary least-squares curated models use complete cases and fixed outcome-specific subsets:

| Outcome | Prespecified pathways |
|---|---|
| FEV1 | IL-13 signalling; Eosinophils; IL-17 signalling; Neutrophils; ECM organisation |
| FeNO | IL-13 signalling; IL-5 signalling; Eosinophils; Mast cells; IL-33 signalling |
| Blood eosinophils | Eosinophils; IL-5 signalling; TAC1; ILC2; severe-asthma eosinophil signature |
| Sputum eosinophils | Eosinophils; TAC1; IL-5 signalling; IL-13 signalling; Mast cells |
| Prior-year exacerbations | IL-13 signalling; IL-17 signalling; Neutrophils; GM-CSF/TNF-alpha macrophage; ECM organisation; Inflammasome |

Reported coefficients and residual-degrees-of-freedom t-based 95% confidence intervals are standardised as `beta × SD(X) / SD(Y)`. Benjamini–Hochberg adjustment is applied within each dataset–outcome curated model.

Sensitivity models additionally include age, sex, smoking category, body mass index and current/ongoing oral-corticosteroid exposure. Curated coefficients are conditional on the other variables in the same model; signature overlap and correlation can produce multicollinearity or suppression.

## Five predictive modelling strategies

The same complete-case participants and outer test splits are used for:

1. Biomarker-only.
2. Curated-pathway-only.
3. Pathway-only Elastic Net using all 39 pathway scores.
4. Biomarker-plus-curated.
5. Biomarker-plus-Elastic-Net, with biomarkers unpenalised.

Biomarker-only and combined-model predictors exclude the target outcome:

| Target | Biomarker predictors |
|---|---|
| FEV1 | FeNO; blood eosinophils; sputum eosinophils |
| FeNO | blood eosinophils; sputum eosinophils |
| Blood eosinophils | FeNO; sputum eosinophils |
| Sputum eosinophils | FeNO; blood eosinophils |
| Prior-year exacerbations | FeNO; blood eosinophils; sputum eosinophils |

Biomarkers are `log1p`-transformed where appropriate. Elastic Net is fitted with `glmnet` 5.0.

## Repeated nested internal validation

The design contains:

- 10 outer repetitions;
- 5 outcome-stratified outer folds in each repetition;
- 5 inner folds for Elastic Net tuning;
- 50 outer-training fits per dataset–outcome–model strategy;
- alpha values 0.10, 0.25, 0.50, 0.75 and 1.00;
- a fixed 100-point relative-lambda grid from 1 to 10^-4;
- seed 20260813; and
- the one-standard-error rule for lambda selection, with deterministic tie-breaking.

Within every inner split, preprocessing and scaling parameters are learned from inner-training participants and applied to inner-validation participants. Following tuning, preprocessing and the selected model are refitted on the complete outer-training set and applied once to the untouched outer-test participants. Only Elastic Net requires inner-loop tuning; biomarker and curated ordinary least-squares models are fitted directly inside each outer-training set.

Selection frequency is the proportion of 50 outer-training fits with a non-zero pathway coefficient. It is a stability descriptor, not a p value, formal variable-importance estimate or causal measure.

## Performance and calibration

Metrics are calculated from out-of-fold predictions after pooling the five outer-test folds within each repetition, then summarised across the 10 repetitions.

For observed values `y` and predictions `p`:

```text
cross-validated R² = 1 - sum((y - p)^2) / sum((y - mean(y))^2)
RMSE               = sqrt(mean((y - p)^2))
MAE                = mean(abs(y - p))
Pearson r           = cor(y, p)
```

Negative R² values are retained because they indicate performance worse than the corresponding mean-only benchmark. Pearson correlation measures ranking, not calibration, and must not be reported as R².

Each participant's 10 held-out predictions are averaged for calibration and paired incremental comparisons. Calibration intercept and slope are obtained from observed-on-predicted linear regression on the model scale. Five equal-frequency bins are used for calibration displays.

Conditional 95% intervals for fixed averaged predictions use 2,000 participant bootstraps with seed 20260814. These intervals do not propagate model refitting, tuning or cross-validation resampling uncertainty. Paired incremental ΔR² p values are Benjamini–Hochberg-adjusted across 40 exploratory combined-versus-biomarker comparisons.

AUC-ROC is not calculated because all five outcomes are continuous.

## Matched-participant concordance

Two-sided Spearman correlations assess same-pathway rank concordance in matched participants:

| Comparison | Matched n | Purpose |
|---|---:|---|
| Brushing GPL570–biopsy GPL570 | 99 | cross-compartment, same platform |
| Sputum GPL570–biopsy GPL570 | 23 | cross-compartment, same platform |
| Sputum GPL570–brushing GPL570 | 31 | cross-compartment, same platform |
| Brushing GPL570–brushing RNA-seq | 118 | same compartment, cross-platform |

Benjamini–Hochberg adjustment is performed across 39 signatures separately within each comparison. Rank concordance is not absolute agreement, independent replication or external validation.

## Interpretation limits

- The design is observational and cannot establish causality.
- Bulk-transcriptome scores can reflect cellular composition, within-cell transcriptional state or both.
- Dataset-specific case mix and missingness affect comparisons.
- Signature overlap and correlated predictors affect coefficient interpretation.
- Same-participant comparisons are not independent validation cohorts.
- The predictive analysis is internally held out but externally unvalidated.
- Concurrent sputum transcriptomic and eosinophil measurements represent same-visit phenotype capture, not prognosis.
- No model class is assumed to perform best for every outcome or compartment.

## Key method citations

The dissertation should be consulted for the complete reference list. Core methods include Subramanian et al. (2005) for GSEA, Barbie et al. (2009) and Hänzelmann et al. (2013) for ssGSEA/GSVA, Zou and Hastie (2005) for Elastic Net, Friedman et al. (2010) for `glmnet`, Varma and Simon (2006) for nested validation, and Van Calster et al. (2019) for predictive-performance interpretation and calibration.
