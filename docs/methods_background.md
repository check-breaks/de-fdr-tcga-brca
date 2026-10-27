## Methods and background

This project asks one question: **can differential expression (DE) tools for RNA-seq be trusted, and when do they fail?** It answers it in two ways: a simulation where the truth is known, and a real-data case study (TCGA-BRCA) where it is not.

### 1. Multiple testing and the false discovery rate

A DE analysis tests about 20,000 genes at once. At a p-value cutoff of 0.05, roughly 1,000 genes would pass by chance even if nothing were truly different. Instead of controlling the error rate per test, the analysis controls the **false discovery rate (FDR)**: the expected proportion of called genes that are false positives.

`src/scratch/multiple_testing.py` implements the **Benjamini-Hochberg (BH)** procedure. Sort the m p-values, scale the i-th smallest by m/i, then take a reverse cumulative minimum so that adjusted values never decrease with rank. BH controls the FDR under independence or positive dependence of the tests. Genes are correlated, so this is an assumption rather than a guarantee, which is one reason the simulation is worth running.

The implementation is validated against `scipy.stats.false_discovery_control`.

### 2. Permutation tests

`src/scratch/permutation.py` implements a two-sided permutation test. If condition labels carry no information, shuffling them should not change the data's distribution. The observed statistic is compared with its distribution under many random label shuffles. The p-value uses the +1 correction, `(1 + #extreme) / (1 + n_perm)`, so it can never be exactly zero.

A permutation test has a hard floor on its p-values. With 3 vs 3 samples there are only 20 possible label assignments, so the smallest two-sided p-value is 2/20 = 0.1. The rank-based Wilcoxon test has the same limit, which is why it shows zero power at 3 samples per group in the simulation.

### 3. The negative binomial GLM

RNA-seq counts are non-negative integers whose variance exceeds their mean, because biological replicates genuinely differ. The negative binomial (NB) distribution captures this with `Var = mu + alpha * mu^2`, where `alpha` is the dispersion.

`src/scratch/nbglm.py` fits, for each gene, `log(mu_ij) = X_i * beta + log(s_j)`:

- **Log link:** the condition coefficient is a natural-log fold change.
- **Offset `log(s_j)`:** corrects for sequencing depth, so a deeper sample is not read as higher expression.
- **IRLS:** each iteration solves a weighted least-squares problem whose weights depend on the current fitted means, until the coefficients converge.
- **Dispersion:** estimated per gene by maximum likelihood, alternating with the IRLS fit.
- **Tests:** a Wald test (coefficient / standard error against a normal) and a likelihood-ratio test.

One design detail matters. In the likelihood-ratio test, dispersion is estimated once from the full model and held fixed in the reduced model. Re-estimating it under the reduced model lets `alpha` absorb the effect being tested and makes the test conservative. DESeq2 and edgeR fix dispersion across both models for the same reason.

**Validation:** with a two-group design the maximum-likelihood coefficients equal the log group means, which is checked in closed form. Fits with an offset and a covariate are checked against a brute-force BFGS optimization of the likelihood. A cross-check against `statsmodels` runs when statsmodels is installed.

### 4. Why dispersion shrinkage matters

With 3 samples per group, a gene's dispersion is estimated from very little information, and maximum likelihood tends to underestimate it, much as dividing by n instead of n-1 underestimates variance. A dispersion that is too small gives a standard error that is too small, a z-score that is too large, and false positives.

`simulation/calibration_check.py` isolates this effect on pure-null data (mean 100, true `alpha` = 0.1). The false-positive rate at a 0.05 cutoff:

| samples per group | true alpha known | alpha estimated per gene |
|---|---|---|
| 3 | 0.056 | 0.185 |
| 5 | 0.050 | 0.123 |
| 10 | 0.045 | 0.076 |
| 20 | 0.043 | 0.060 |

The model itself is calibrated. The inflation comes from the dispersion estimate. DESeq2, edgeR and limma fix this by **borrowing information across genes**: genes of similar expression have similar dispersion, so each noisy gene-wise estimate is pulled toward a fitted trend (empirical Bayes shrinkage).

### 5. The tools compared

| Method | Data scale | Distribution | Dispersion / variance handling |
|---|---|---|---|
| t-test, Wilcoxon | log-CPM | normal / rank-based | none |
| Scratch NB GLM | raw counts | negative binomial | gene-wise ML, no shrinkage |
| DESeq2 | raw counts | negative binomial | shrunk toward a fitted trend, Wald test |
| edgeR (QL) | raw counts | negative binomial | trended dispersion, quasi-likelihood F-test |
| limma-voom | log-CPM with precision weights | normal | mean-variance weights plus moderated variances |

All three established tools share one idea: with few replicates, pool information across genes.

### 6. Simulation design

Real data cannot reveal the true FDR because the truly differential genes are unknown. `simulation/simulate.py` generates NB counts with 3,000 genes, of which 10% are truly DE. Base means are lognormal, dispersion follows a DESeq2-style trend `phi(mu) = 0.05 + 1/mu` with lognormal gene-level scatter, and sequencing depth varies by sample. Each tool's p-values are adjusted with the scratch BH implementation and scored against the truth:

- **Realized FDR** = false calls / total calls.
- **Power** = true DE genes called / all true DE genes.

The two must be read together: calling every gene gives high power and a terrible FDR.

**Results (5 simulated datasets per sample size, nominal FDR 5%):**

| Method | FDR at 3 / 5 / 10 per group | Power at 3 / 5 / 10 per group |
|---|---|---|
| DESeq2 | 0.061 / 0.084 / 0.072 | 0.43 / 0.62 / 0.80 |
| edgeR (QL) | 0.036 / 0.047 / 0.047 | 0.38 / 0.59 / 0.79 |
| limma-voom | 0.033 / 0.042 / 0.045 | 0.37 / 0.57 / 0.77 |
| Scratch NB GLM | 0.503 / 0.325 / 0.285 | 0.63 / 0.72 / 0.84 |
| t-test | 0.000 / 0.037 / 0.127 | 0.01 / 0.36 / 0.74 |
| Wilcoxon | 0.000 / 0.000 / 0.087 | 0.00 / 0.00 / 0.70 |

Interpretation:

- DESeq2 is slightly liberal and slightly more powerful. edgeR (QL) and limma-voom are slightly conservative. Differences of a few percentage points are within the replicate-to-replicate spread (SD of 1 to 2 points), so they should be treated as tentative.
- The scratch NB GLM has high power but 28 to 50% FDR. Since DESeq2 uses the same model and reaches 6 to 8%, the difference is attributable to dispersion shrinkage.
- The t-test and Wilcoxon lose almost all power at small n. An FDR of 0 means no genes were called, not good calibration.
- Power rises from roughly 40% at 3 per group to roughly 80% at 10 per group, which matters for the small mutant group in the TCGA analysis.

**Limitations:** the simulation draws from an NB model, which favors NB-based methods. Real data contain outliers, correlated genes and other structure it does not capture. DESeq2 was run with raw p-values and without its independent filtering or outlier handling, so its power may differ from a standard run. Five replicate datasets per sample size is thin.

### 7. TCGA-BRCA: BRCA1/BRCA2 and confounding

**Data.** Open-access TCGA-BRCA from the GDC: raw STAR gene counts for primary tumors (one per patient), somatic mutation calls, PAM50 subtype and clinical data. A tumor is labelled BRCA-mutant if it carries a damaging somatic variant (nonsense, frame-shift, splice-site or translation-start) in BRCA1 or BRCA2.

**Confounding.** A confounder is associated with both the exposure and the outcome. Here the exposure is BRCA1/2 mutation, the outcome is gene expression, and the confounder is **tumor subtype**: BRCA1-mutant tumors are predominantly basal-like, and basal-like tumors differ from other subtypes in expression regardless of mutation. A naive comparison (`~ mutation`) therefore largely rediscovers subtype. Adding subtype to the model (`~ subtype + mutation`) estimates the mutation effect within subtypes. The analysis compares the number of DE genes under the two designs across DESeq2, edgeR and limma-voom.

Adjustment has limits. If nearly all mutants fall in one subtype, there is little information to separate the two effects, and the adjusted estimate becomes noisy.

**Caveats.**
- Open-access MAFs contain somatic mutations only. Germline BRCA carriers require controlled access (dbGaP).
- The mutant group is small, so power is limited. The simulation's power curves indicate what to expect.
- BRCA1 can also be silenced by promoter methylation, so some wild-type tumors may be functionally BRCA1-deficient.
- GDC and TCGAbiolinks field names change between releases and were verified at the time of analysis.

**Results.** *To be filled in after the analysis runs: mutant counts per subtype, number of DE genes under each design and tool, tool agreement, and BRCA1/2 expression and survival by subtype.*

### 8. Reproducibility

```
conda env create -f environment.yml && conda activate de-fdr-tcga-brca
pytest -q
python simulation/simulate.py --reps 5 --n_genes 3000
python simulation/run_python_methods.py results/sim
Rscript simulation/run_r_methods.R results/sim
python simulation/evaluate.py results/sim
python simulation/calibration_check.py
Rscript tcga/01_download_prepare.R
Rscript tcga/02_confounded_de.R
Rscript tcga/03_survival_expression.R
```
