# beef-exudate-metabolomics
Unravelling novel metabolites in beef exudate to ensure meat freshness and safety. An untargeted LC-MS/MS metabolomics pipeline for detecting microbial spoilage in beef, comparing a conventional destructive matrix (muscle biopsy tissue) against a non-destructive one (meat exudate / drip loss) over 21 days of cold storage.

## 👥 Research Team
Bioinformatics team - Horizon Science Communication LLC
https://www.horizon-sci-comm.us/bioinformatics-services

- Hesham Abdullah, PhD 
- Abeer Farag, PhD    
- Ahmed Abdelmaksoud
- Mai Mohamed Salah


**Overview**

Early detection of microbial contamination in meat currently depends on destructive
sampling, which limits how often and how widely a carcass or cut can be tested. This
project asks whether meat exudate — the fluid that collects passively during cold
storage — carries a metabolic signature of contamination that is at least as
informative as the tissue itself.

Beef samples were inoculated with E. coli or Salmonella and profiled alongside
uninoculated controls across a three-week storage window (days 1, 7, 14, 21), with
36 samples per comparison. The analysis tracks temporal trajectories of
contamination-associated metabolites in both matrices and evaluates which matrix
yields the earlier and cleaner discrimination.

Principal finding: exudate shows more pronounced treatment-specific divergence
and clearer biomarker signals than tissue, supporting non-destructive sampling as a
viable — and possibly more sensitive — route to spoilage surveillance.

**Project Overview*
<p align="center">
  <img src="Metabolite_analysis_pipeline.png" alt="RNA-seq pipeline" width="900">
</p>

**Pipeline**

1. Experimental design. Paired tissue and exudate matrices; contaminated
(E. coli, Salmonella) vs. control groups; 36 samples per comparison; four
storage timepoints spanning 21 days.

2. Preprocessing and QC. Features with >40% zeros are filtered; remaining zeros
receive half-minimum imputation. Log2 transformation, per-sample median
normalization, and Pareto scaling bring raw MS intensities from heavily
right-skewed (skewness 17.29) to approximately symmetric (-0.16).

3. Univariate analysis. Volcano plots at p < 0.05 and |log2FC| > 0.5, with
Cohen's d effect sizes to separate statistically significant features from
biologically meaningful ones. The significant set is dominated by large (>0.5) and
very large (>0.8) effects.

4. Temporal pattern discovery. K-means clustering resolves four distinct
metabolic trajectories across the storage window; heatmaps of the top 50
discriminating features reveal treatment-specific clusters, most sharply in exudate.

5. Multivariate modeling. PCA (unsupervised) alongside PLS-DA (supervised, VIP
≥ 1.0, 100-iteration validation). Random forest classifiers are evaluated under
leave-one-day-out validation, so models are tested on storage days they never
saw during training and cannot succeed by memorizing the temporal pattern.

6. Annotation confidence. All reported features carry an MSI confidence level
assigned by R/msi_classifier.R (see attached R code).

**Metabolite Identification workflow*
<p align="center">
  <img src="Metabolite_Identification_steps.png" alt="Metabolite Identification workflow" width="900">
</p>











