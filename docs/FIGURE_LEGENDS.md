# Final dissertation figure legends

## Figure 1. Airway transcriptomic dataset sizes
Bars show the numbers of adult U-BIOPRED transcriptomic profiles in sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570 and bronchial brushing RNA-seq. Bar height and the displayed label give n.

## Figure 2. Clinical outcome availability by airway dataset
Bars report the participants available for FEV1 % predicted, FeNO, blood eosinophils, sputum eosinophils and previous-12-month exacerbations relative to each transcriptomic dataset total. Differences in completeness explain outcome-specific analysis denominators.

## Figure 3. Exact participant-overlap patterns across airway datasets
Bars show mutually exclusive exact transcriptomic-dataset overlap patterns. The 495 profiles represented 243 unique participants; the four datasets partially overlap and must not be treated as independent cohorts. Baseline characteristics are reported separately in Table 3 and sample/omics preprocessing in Table 1.

## Figure 4. Normalized ssGSEA pathway activity
Panels A-D show the 39 fixed asthma-relevant pathway scores across participants in sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570 and bronchial brushing RNA-seq, respectively. Blue indicates lower, white approximately central and red higher saved normalized ssGSEA scores on common -1 to +1 display limits. Hierarchical ordering is descriptive. The panels show graded variation across continuous-valued scores but do not test for, prove or exclude discrete latent classes or endotypes.

## Figure 5. Pathway-clinical Spearman associations
Panels A-D show marginal pathway-outcome Spearman correlations in the four airway datasets. Blue indicates negative and red positive rho; plus/minus glyphs preserve direction in grayscale. Asterisks mark Benjamini-Hochberg false-discovery-rate support within each dataset-outcome family (* q<0.05, ** q<0.01, *** q<0.001). These are univariable associations, not conditional effects, causal estimates or predictive performance.

## Figure 6. Prespecified curated pathway models
Panels A-D show standardised coefficients and residual-degrees-of-freedom t-based 95% confidence intervals for biologically prespecified outcome-specific pathway subsets in each dataset. Filled red symbols indicate within-family BH q<0.05 and open grey symbols q>=0.05. The plots do not show all 39 pathways. Coefficients are conditional on other pathways in the same model; correlated signatures can produce collinearity or suppression.

## Figure 7. Elastic Net pathway-selection stability
Panels A-D show the ten pathways with the highest pathway-only Elastic Net selection frequency for each outcome in each dataset. Selection frequency is the proportion of 50 outer-training fits (10 repeats x 5 folds) with a non-zero coefficient and is not a p value. Red upward triangles indicate predominantly positive coefficients, blue downward triangles predominantly negative coefficients, and point size is a coefficient-based relative share rather than variance decomposition or causal importance.

## Figure 8. Internally held-out predictive performance
Panel A compares mean model-scale cross-validated R-squared for biomarker-only, curated-pathway-only and Elastic-Net-pathway-only models on identical common cases and outer folds. Panel B shows paired change in model-scale R-squared when curated or Elastic Net pathways are added to the biomarker baseline. Intervals are participant-bootstrap intervals conditional on fixed averaged repeated out-of-fold predictions and do not include refitting, tuning or fold-generation uncertainty. The analysis is exploratory internal validation; concurrent sputum-eosinophil estimation is same-sample phenotype capture rather than prognosis.

## Figure 9. Compartment and platform concordance
Panel A shows same-pathway matched-participant Spearman concordance among GPL570 airway compartments. Panel B shows same-compartment concordance between bronchial brushing GPL570 and RNA-seq in 118 matched participants, with low RNA-seq signature coverage flagged. These are same-participant concordance analyses, not independent replication or external validation; matched sputum subsets are small.

## Supplementary Figure S1. PCA of z-standardised ssGSEA scores
PCA was run separately in each airway dataset after pathway-wise z-standardisation. PC1 and PC2 percentages are reported on the axes. The projection is descriptive; axes are dataset-specific and do not test for, prove or exclude latent classes or endotypes.

## Supplementary Figures S2-S5. Binned internal calibration
Figures S2-S5 show calibration for sputum GPL570, bronchial brushing GPL570, bronchial biopsy GPL570 and bronchial brushing RNA-seq, respectively. Five equal-frequency bins summarise observed versus predicted values for biomarker-only, curated-pathway and Elastic-Net-pathway models. Each participant's ten held-out predictions were averaged before binning. Bars are standard errors of observed bin means on the modelling scale, not uncertainty from model refitting.

## Supplementary Figure S6. Normalized versus z-standardised ssGSEA heatmaps
Paired panels use identical participant and pathway ordering. Normalized panels retain between-pathway score location and dispersion; z-standardised panels centre each pathway to mean zero and standard deviation one within its dataset. Blue indicates lower, white central and red higher values on the labelled scales.

## Supplementary Figure S7. Between-participant variation in normalized ssGSEA scores
Points report the between-participant standard deviation of each normalized ssGSEA pathway score before z-standardisation. Larger values indicate greater participant-to-participant spread on the saved score scale. This is descriptive and is not a test of biological importance or cross-dataset calibration.

## Supplementary Figure S8. Covariate-adjusted curated-model sensitivity
Primary prespecified pathway subsets were refitted with age, sex, smoking status, BMI and current oral-corticosteroid exposure where complete. Symbols and intervals follow Figure 6. Changes can reflect both covariate adjustment and differences in the covariate-complete sample.

## Supplementary Figure S9. Elastic Net stability across all 39 pathways
Selection frequency is the proportion of 50 pathway-only Elastic Net outer-training fits in which a pathway coefficient was non-zero. Blue denotes a predominantly negative coefficient and red a predominantly positive coefficient; colour intensity gives selection frequency and plus/minus glyphs redundantly show sign. Selection frequency is a stability descriptor, not a p value.

## Supplementary Figure S10. RNA-seq signature-gene coverage
Bars show the percentage of each fixed signature detected in the processed bronchial brushing RNA-seq matrix. Labels give detected/total harmonised genes and LOW marks the prespecified low-coverage flag. All retained signatures met the ssGSEA minimum of three genes, but low proportional coverage requires pathway-specific caution.

## Supplementary Figure S11. Below-zero prediction sensitivity
Bars show counts of below-zero participant-level averaged predictions for clinical outcomes that cannot be negative; FEV1 is excluded. Unconstrained predictions remain primary. Clipping was evaluated only as a sensitivity for original-scale error metrics, while model-scale R-squared was never clipped.

## Supplementary Figure S12. Biological scope and derivation context of the 39 prespecified airway signatures
Each tile represents one ssGSEA signature and reports its reader-facing label and frozen source-list gene-entry count. Tile colour identifies derivation context; shape gives provenance traceability: A, exact analysed list and source reconstructed; B, source biology or contrast supported but exact-list construction partly reconstructed; C, material, label, source, symbol or construction ambiguity. Asthma-relevant means selected to represent biology implicated in asthma or airway disease; not every list was derived in asthma. This is a conceptual provenance map, not a causal pathway network.
