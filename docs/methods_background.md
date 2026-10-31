# Methods, results and limitations

This project asks one question: **can RNA-seq differential expression (DE) tools be trusted, and when do they fail?**
It answers it in three ways: statistics implemented from scratch, simulations where the truth is known, and a real-data
case study on TCGA-BRCA (BRCA1/BRCA2 mutation) where it is not.

## 1. Multiple testing and the false discovery rate

A DE analysis tests about 20,000 genes at once. At a p-value cutoff of 0.05, roughly 1,000 genes would pass by chance even if
nothing differed. Instead of controlling the error rate per test, the analysis controls the **false discovery rate (FDR)**: the
expected proportion of called genes that are false positives.

`src/scratch/multiple_testing.py` implements **Benjamini-Hochberg (BH)**: sort the m p-values, scale the i-th smallest by m/i,
then take a reverse cumulative minimum so adjusted values never decrease with rank. BH controls the FDR under independence or
positive dependence of the tests. Genes are correlated, so this is an assumption, not a guarantee, which is one reason the
simulations are worth running. The implementation is validated against `scipy.stats.false_discovery_control`.

## 2. Permutation tests

`src/scratch/permutation.py` implements a two-sided permutation test. If labels carry no information, shuffling them should not
change the data's distribution, so the observed statistic is compared with its distribution under many random shuffles. The
p-value uses the +1 correction, `(1 + #extreme) / (1 + n_perm)`, so it can never be exactly zero.

A permutation test has a hard floor on its p-values. With 3 vs 3 samples there are 20 possible assignments, so the smallest
two-sided p-value is 2/20 = 0.1. The rank-based Wilcoxon test has the same floor, which is why it has zero power at 3 samples per
group (and at 5 per group, where the smallest p is about 0.008, too large to survive a correction over thousands of genes).

## 3. The negative binomial GLM

RNA-seq counts are non-negative integers whose variance exceeds their mean, because biological replicates genuinely differ.
The negative binomial (NB) distribution captures this with `Var = mu + alpha * mu^2`, where `alpha` is the dispersion.
`src/scratch/nbglm.py` fits, per gene, `log(mu_ij) = X_i * beta + log(s_j)`:

- **Log link:** the condition coefficient is a natural-log fold change.
- **Offset `log(s_j)`:** corrects for sequencing depth.
- **IRLS:** each iteration solves a weighted least-squares problem whose weights depend on the current fitted means.
- **Dispersion:** estimated per gene by maximum likelihood, alternating with the IRLS fit.
- **Tests:** a Wald test and a likelihood-ratio test.

In the likelihood-ratio test, dispersion is estimated once from the full model and held fixed in the reduced model. Re-estimating
it under the reduced model lets `alpha` absorb the effect being tested and makes the test conservative (a bug found and fixed
during development; DESeq2 and edgeR also hold dispersion fixed across the two models).

**Validation** (`tests/test_scratch.py`, 10 tests): with a two-group design the maximum-likelihood coefficients equal the log group
means (closed form); fits with an offset and a covariate match a brute-force BFGS optimization of the likelihood; BH matches
scipy; permutation p-values are uniform under the null. A cross-check against `statsmodels` runs when statsmodels is installed.

## 4. Why dispersion shrinkage matters

With few samples per group, a gene's dispersion is estimated from very little information, and maximum likelihood tends to
underestimate it. A dispersion that is too small gives a standard error that is too small and false positives.
`simulation/calibration_check.py` isolates this on pure-null data (mean 100, true `alpha` = 0.1, 1,500 genes per setting).
False-positive rate at a 0.05 cutoff:

| samples per group | true alpha known | alpha estimated per gene |
|---|---|---|
| 3 | 0.056 | 0.185 |
| 5 | 0.050 | 0.123 |
| 10 | 0.045 | 0.076 |
| 20 | 0.043 | 0.060 |

The model itself is calibrated. The inflation comes from the dispersion estimate. DESeq2, edgeR and limma address it by
borrowing information across genes: genes of similar expression have similar dispersion, so each noisy gene-wise estimate is
pulled toward a fitted trend (empirical Bayes shrinkage).

## 5. The methods compared

| Method | Data scale | Distribution | Dispersion / variance handling |
|---|---|---|---|
| t-test, Wilcoxon | log-CPM | normal / rank-based | none |
| Scratch NB GLM | raw counts | negative binomial | gene-wise ML, no shrinkage |
| DESeq2 | raw counts | negative binomial | shrunk toward a fitted trend, Wald test |
| edgeR (QL) | raw counts | negative binomial | trended dispersion, quasi-likelihood F-test |
| limma-voom | log-CPM with precision weights | normal | mean-variance weights plus moderated variances |

All three established tools share one idea: with few replicates, pool information across genes.

In the simulations, the t-test and Wilcoxon use log2 CPM with **median-of-ratios size factors**. An earlier version used
library-size CPM only, which inflated their FDR whenever many genes shifted between groups (see section 10).
On real data, the Wilcoxon uses TMM-normalized log-CPM; in the "adjusted" design the log-CPM values are first regressed on
subtype and the residuals are tested (an approximation to a subtype-stratified rank test).

## 6. Balanced simulation

Real data cannot reveal the true FDR because the truly differential genes are unknown. `simulation/simulate.py` generates NB
counts for 3,000 genes, 10% truly DE (log2 fold change at least 0.6). Base means are lognormal, dispersion follows a DESeq2-style
trend `phi(mu) = 0.05 + 1/mu` with lognormal gene-level scatter, and depth varies by sample. There are 5 datasets for each of
3, 5 and 10 samples per group. p-values are adjusted with the scratch BH and scored against the truth.

- **Realized FDR** = false calls / total calls. **Power** = true DE genes called / all true DE genes. Read them together:
  calling every gene gives high power and a terrible FDR.

| Method | FDR at 3 / 5 / 10 per group | Power at 3 / 5 / 10 per group |
|---|---|---|
| DESeq2 | 0.061 / 0.084 / 0.072 | 0.43 / 0.62 / 0.80 |
| edgeR (QL) | 0.036 / 0.047 / 0.047 | 0.38 / 0.59 / 0.79 |
| limma-voom | 0.033 / 0.042 / 0.045 | 0.37 / 0.57 / 0.77 |
| Scratch NB GLM | 0.477 / 0.297 / 0.140 | 0.63 / 0.72 / 0.83 |
| t-test | 0.000 / 0.028 / 0.036 | 0.01 / 0.36 / 0.74 |
| Wilcoxon | 0.000 / 0.000 / 0.023 | 0.00 / 0.00 / 0.70 |

The standard deviation of FDR across the 5 datasets is roughly 1 to 4 percentage points, so differences of a few points between
DESeq2, edgeR and limma-voom are tentative.

- DESeq2 is slightly liberal and slightly more powerful; edgeR (QL) and limma-voom are slightly conservative.
- The scratch NB GLM has high power but 14% to 48% FDR at small n. DESeq2 uses the same model and reaches 6% to 8%, consistent
  with the calibration check above (the remaining difference is dispersion shrinkage).
- An FDR of 0 means no genes were called (t-test and Wilcoxon at the smallest sample sizes), not good calibration.
- Power rises from about 40% at 3 per group to about 80% at 10 per group for the established tools.

## 7. Unbalanced and confounded simulation

TCGA mutants are few and unevenly spread across subtypes. `simulation/simulate_unbalanced.py` mimics this: 7 or 11 cases
against 300 controls, with three latent "subtypes". Subtypes 1 and 2 differ from subtype 0 in 30% of genes, independent of case
status. 10% of genes have a true case effect (log2 fold change at least 0.6). Two scenarios:

- **mixed:** cases come from the same subtype mix as controls (subtype is not a confounder).
- **confounded:** 64% of cases come from one subtype, like 7 of 11 BRCA1-mutant tumors being basal-like.

Genes driven by subtype count as false positives, because they have no true case effect. Each tool is run with a naive design
(`~ condition`) and an adjusted design (`~ subtype + condition`). 20 simulated datasets per cell.

**Realized FDR (nominal 5%):**

| Scenario | Cases | Design | DESeq2 | edgeR | limma-voom | Wilcoxon |
|---|---|---|---|---|---|---|
| mixed | 7 | naive | 0.040 | 0.034 | 0.037 | 0.021 |
| mixed | 7 | adjusted | 0.045 | 0.043 | 0.047 | 0.027 |
| mixed | 11 | naive | 0.055 | 0.046 | 0.048 | 0.032 |
| mixed | 11 | adjusted | 0.050 | 0.046 | 0.047 | 0.034 |
| confounded | 7 | naive | 0.219 | 0.212 | 0.224 | 0.138 |
| confounded | 7 | adjusted | 0.051 | 0.042 | 0.048 | 0.028 |
| confounded | 11 | naive | 0.364 | 0.357 | 0.362 | 0.279 |
| confounded | 11 | adjusted | 0.040 | 0.037 | 0.046 | 0.029 |

**Power (adjusted design):**

| Scenario | Cases | DESeq2 | edgeR | limma-voom | Wilcoxon |
|---|---|---|---|---|---|
| mixed | 7 | 0.847 | 0.846 | 0.823 | 0.811 |
| mixed | 11 | 0.897 | 0.894 | 0.879 | 0.889 |
| confounded | 7 | 0.837 | 0.835 | 0.816 | 0.792 |
| confounded | 11 | 0.886 | 0.886 | 0.871 | 0.879 |

- In this simulation, **all four methods are calibrated once subtype is in the model**, including with 7 to 11 cases against 300
  controls. Adjusting costs nothing when there is no confounding.
- **Confounding is large when it is present.** The naive FDR reaches 0.22 to 0.36 for the NB-based tools. The standard deviation
  of the naive FDR (0.17 to 0.21) is as large as its mean, because the damage depends on how many cases happen to fall in the
  enriched subtype in a given draw (for example 3 of 7 in one dataset and 8 of 11 in another).
- A clean NB simulation with a heavily unbalanced design does **not** reproduce the inflation DESeq2 and edgeR show on real data
  (section 9.2).

## 8. TCGA-BRCA data

Open-access data from the GDC, accessed with TCGAbiolinks: STAR gene counts for primary tumors (1,111 files, 1,095 patients
after keeping one tumor per patient), somatic mutation calls, PAM50 subtype, and clinical data (vital status, days to death
or follow-up, age). A tumor is labelled BRCA-mutant if it carries a damaging somatic variant (nonsense, frame-shift, splice-site
or translation-start) in BRCA1 or BRCA2. Twelve tumors have no PAM50 label and are excluded from the DE analyses, leaving 1,083.

| PAM50 subtype | tumors | BRCA1-mutant | BRCA2-mutant |
|---|---|---|---|
| Basal | 190 | 7 | 1 |
| Her2 | 82 | 0 | 1 |
| LumA | 562 | 2 | 4 |
| LumB | 209 | 2 | 1 |
| Normal | 40 | 0 | 0 |
| total | 1,083 | 11 | 7 |

No tumor carries both a BRCA1 and a BRCA2 mutation. Nine tumors carry a BRCA2 mutation in total, but two lack a PAM50 label.
Seven of the 11 BRCA1-mutant tumors (64%) are basal-like, against 18% of the cohort, so subtype is a strong confounder for BRCA1.
BRCA2 mutants are spread more evenly (mostly LumA).

## 9. TCGA-BRCA results

### 9.1 Mutant vs wild-type DE, naive and subtype-adjusted

Controls are a random sample of about 300 wild-type-for-the-gene tumors, stratified by subtype, so that memory use stays low.
Genes called at q < 0.05. The number of genes tested differs between rows (the expression filter depends on group sizes), so
compare each naive result with its own adjusted result, not across rows.

| Comparison | Design | DESeq2 | edgeR | limma-voom | all three |
|---|---|---|---|---|---|
| Positive control: 100 Basal vs 100 LumA | none | 15,151 | 14,960 | 14,545 | 13,304 |
| BRCA1 (11 vs 302) | naive | 2,210 | 1,195 | 293 | 174 |
| BRCA1 | subtype-adjusted | 264 | 202 | 2 | 0 |
| BRCA2 (7 vs 302) | naive | 134 | 100 | 10 | 3 |
| BRCA2 | subtype-adjusted | 87 | 119 | 21 | 6 |
| BRCA1 or BRCA2 (18 vs 303) | naive | 1,638 | 1,043 | 31 | 18 |
| BRCA1 or BRCA2 | subtype-adjusted | 297 | 244 | 2 | 2 |
| BRCA1 within Basal only (7 vs 183) | none | 55 | 105 | 0 | 0 |

- **The positive control works.** The three tools largely agree on a very strong effect (88% to 92% of each tool's calls are shared).
  About 70% of tested genes are "significant", which reflects 100 samples per group on the largest axis of breast-cancer
  variation, not 15,000 distinct pieces of biology.
- **Most of the naive BRCA1 signal was subtype.** Adjustment removes about 88% of DESeq2's calls and 83% of edgeR's. BRCA2
  barely changes, consistent with its more even subtype distribution. The pooled result resembles BRCA1.
- **The tools disagree.** After adjustment, the BRCA1 overlap between DESeq2 and edgeR is 152 genes, and limma-voom shares
  none with DESeq2. In the basal-only comparison DESeq2 and edgeR share 27 genes and limma-voom calls none.
- No gene-level result here should be read as BRCA biology until it is compared with the null in the next section.

### 9.2 Permutation null on real data

To measure false-positive behavior in exactly this design, the real BRCA1 mutants were replaced by **pseudo-mutants**: 11
tumors drawn at random from the 1,065 tumors wild-type for both genes, with the same subtype mix as the real mutants (7 Basal,
2 LumA, 2 LumB), against 150 stratified controls, with the subtype-adjusted design. Ten draws.

| Method | mean calls | median | max | draws with any call | mean % of raw p < 0.05 (expected 5) |
|---|---|---|---|---|---|
| DESeq2 | 513 | 473 | 1,179 | 100% | 12.8 |
| edgeR | 389 | 319 | 824 | 100% | 11.3 |
| limma-voom | 8.4 | 1.5 | 62 | 70% | 6.1 |
| Wilcoxon | 0 | 0 | 0 | 0% | 4.7 |

- DESeq2 and edgeR call hundreds of genes on random labels, in every draw, with raw p-values about 2.3 to 2.5 times the expected 5%.
- limma-voom is usually well controlled but not always (one draw reached 62 calls, with about 11% of raw p-values below 0.05).
  The Wilcoxon never called a gene.
- **The real BRCA1 result (264 DESeq2, 202 edgeR, 2 voom, 0 Wilcoxon calls, adjusted design) lies within these null ranges**,
  so it cannot be distinguished from what random labels produce in this design. The real run used 300 controls and the null 150,
  so the comparison is approximate.
- This does not match the simulations, where a correctly specified NB model is calibrated. The cause was not tested. A
  plausible candidate is real-data structure the NB model does not capture (outliers, sample-level factors such as tumor purity
  or proliferation, gene-gene correlation). A partial test with surrogate variables is in section 9.6.
- The number of calls varies widely between draws (DESeq2 from about 180 to 1,179), so a single run's count is not a stable summary.

### 9.3 Wilcoxon check: power and null

| Check | Result |
|---|---|
| Positive control (100 Basal vs 100 LumA) | 14,698 of 21,139 genes at q < 0.05 (71.6% raw p < 0.05) |
| Null: 50 pseudo-mutant draws, adjusted | 0 calls in every draw; 5.0% of raw p < 0.05 |
| Real BRCA1 (11), naive and adjusted | 0 calls in all 10 control draws |
| Real BRCA2 (7), naive and adjusted | 0 calls in all 10 control draws |
| Real BRCA1 or BRCA2 (18) | adjusted: 0 in all draws; naive: 2 calls in 2 of 10 draws, otherwise 0 |

The Wilcoxon is calibrated on the null and has power on a strong effect, and it finds nothing in the real mutants. In the
simulation, with effects of at least 0.6 log2 in 10% of genes, it had 79% to 89% power at 7 to 11 cases, so effects that large
and widespread would probably have been detected, *if* the simulation is a fair guide. Smaller or sparser effects would not
have been. The accurate statement is "no effect detectable at this sample size", not "no effect".

### 9.4 The mutated genes themselves

| Analysis | log2FC of the gene | q (DESeq2 / edgeR / voom) |
|---|---|---|
| BRCA1 mutants, naive | +0.14 | 0.82 / 0.85 / 0.83 |
| BRCA1 mutants, adjusted | +0.19 | 0.87 / 0.89 / 1.00 |
| BRCA1 mutants, basal only | +0.22 | 0.93 / 0.94 / 1.00 |
| BRCA2 mutants, naive | -1.38 | 0.044 / 0.39 / 0.999 |
| BRCA2 mutants, adjusted | -1.20 | 0.051 / 0.22 / 0.77 |

BRCA1 expression is not lower in tumors with somatic BRCA1 mutations. This is weak evidence about the labels (a truncating
mutation lowers its transcript only if the mutant allele is degraded and the tumor does not compensate, and BRCA1 expression
also tracks proliferation) and is reported as an open caveat, not as validation. BRCA2 is about 2.3 to 2.6 times lower in
mutants, the expected direction, but only DESeq2 reaches significance, and only just.

### 9.5 Survival

Median BRCA1 expression by subtype (approximate, read from the box plot, log2 CPM): LumB highest at about 5, Basal about 4.3,
Her2 about 4.2, LumA about 4.1, Normal-like lowest at about 3.5.

A Cox model of overall survival on BRCA1 expression (log2 CPM), subtype and age, with 1,036 complete cases and 144 deaths:

| Term | Hazard ratio | 95% CI | p |
|---|---|---|---|
| BRCA1 expression (per log2 unit) | 1.15 | 0.96 to 1.38 | 0.13 |
| Age (per year) | 1.038 | 1.024 to 1.052 | 1.3e-7 |
| Her2 vs Basal | 1.86 | 1.00 to 3.48 | 0.051 |
| LumA vs Basal | 0.73 | 0.46 to 1.15 | 0.18 |
| LumB vs Basal | 1.11 | 0.65 to 1.88 | 0.70 |
| Normal vs Basal | 1.52 | 0.65 to 3.60 | 0.34 |

Concordance is 0.68. Age is the only clear predictor. BRCA1 expression shows no detectable association; with 144 events the
interval is wide, so the result is inconclusive, not evidence of no effect. The proportional-hazards test found no significant
violation for any term (BRCA1 p = 0.25, subtype 0.13, age 0.31), though the global test was borderline (p = 0.058).
Fifty-nine patients were dropped for missing data (48 missing age, 12 missing subtype, 1 missing time).
Twenty-one patients have a recorded time of zero or below (19 alive with a 0-day follow-up, 2 deaths recorded at day 0, and one
with an impossible follow-up of -7 days); they contribute no usable follow-up. Median follow-up is about 2.3 years, which
limits the number of events. A sensitivity analysis (`tcga/03c_cox_sensitivity.R`) leaves the BRCA1 estimate unchanged:
excluding the 21 patients (20 of whom were in the model, leaving 1,016 patients and 142 deaths) gives HR 1.15 (0.96 to 1.38,
p = 0.14), and additionally stratifying by subtype, so that each subtype has its own baseline hazard, gives HR 1.14
(0.94 to 1.36, p = 0.18).
The Kaplan-Meier plot with a median split gave log-rank p = 0.37 (1,094 patients, 152 deaths). Only 19 and 21 patients remain at
risk at 10 years, and 11 and 8 at 12 years, so late steps rest on a handful of people.
Only 11 patients carry a BRCA1 mutation, too few for a survival comparison, so this analysis concerns expression, not mutation.

### 9.6 Surrogate variables on the permutation null

To test whether hidden sample-level structure explains the false calls in section 9.2, five of the ten null draws were rerun with
surrogate variables (`tcga/06_sva_null.R`): the same pseudo-mutants, the same 150 controls and the same gene filter, with
`svaseq` estimating the variables from the 10,000 most variable genes while keeping subtype in both models and protecting the
case/control label. `num.sv` returned 10 in every draw, which is the cap that was set, so the data may call for more.
Mean over the five draws:

| Method | calls without SVs | calls with 10 SVs | raw p < 0.05 without (%) | raw p < 0.05 with (%) |
|---|---|---|---|---|
| DESeq2 | 575 | 402 | 13.1 | 11.5 |
| edgeR | 497 | 161 | 11.7 | 9.7 |
| limma-voom | 3.2 | 12.0 | 6.2 | 7.2 |
| Wilcoxon | 0 | 0 | 4.4 | 5.3 |

- Surrogate variables **reduce but do not remove** the inflation. edgeR's false calls fall by about two thirds and DESeq2's by
  about 30%, yet both still call hundreds of genes and their raw p-values stay at about twice the expected 5%. The DESeq2 effect
  is not consistent: it rose in one draw (262 to 298 calls) and fell in another (561 to 346).
- They made the two well-behaved methods **slightly worse** (limma-voom 3.2 to 12.0 calls; raw p-values 6.2% to 7.2% for
  limma-voom and 4.4% to 5.3% for Wilcoxon). One possible reason, not tested, is that variables estimated from the same data can
  correlate with a random label by chance and inflate test statistics. SVA is therefore not a free fix here.
- The variables relate loosely to known biology. The strongest Spearman correlation of any variable with a marker gene,
  averaged over draws, was 0.56 for PTPRC (immune cells), 0.49 for MKI67 (proliferation) and 0.42 for KRT19 (epithelial
  content). This suggests a mixture of immune infiltration, proliferation and tumor content rather than one factor. The marker
  genes were probably among those used to estimate the variables, so these correlations are somewhat inflated and are not
  causal evidence.
- Hidden structure therefore explains part of the inflation, more for edgeR than for DESeq2. What remains is not explained;
  outliers or heavy-tailed counts are candidates and were not examined.
- A comparison on the real BRCA1 mutants with surrogate variables was attempted but not completed: the first draw had run for
  over two hours, against 30 to 40 minutes for a null draw, and was stopped. No real-mutant result with surrogate variables is
  reported.

### 9.7 PCA and quality control

`tcga/07_pca_qc.R` ran a principal component analysis (PCA) of all 1,095 tumors: genes with CPM above 1 in at least 20 tumors
(21,915 of 60,660), log2-CPM, then the 500 most variable genes. The table gives the share of variance in each component and the
R-squared of a linear model of that component on PAM50 subtype and on tissue source site (the two-letter code in the patient ID).

| Component | % of variance | R-squared, subtype | R-squared, site |
|---|---|---|---|
| PC1 | 19.6 | 0.67 | 0.05 |
| PC2 | 11.3 | 0.30 | 0.06 |
| PC3 | 5.3 | 0.30 | 0.12 |
| PC4 | 3.9 | 0.31 | 0.04 |
| PC5 | 2.3 | 0.05 | 0.03 |

- Subtype is the largest single source of variation: basal tumors separate from luminal tumors along PC1. The BRCA1 mutants fall
  mostly inside the basal cluster, which is the confounding of section 9.1 shown as a picture (`results/figures/pca_subtype_site.png`).
- Tissue source site shows no visible pattern. Its R-squared values are small, but a label with many small groups gets a few
  percent by chance, so they should not be read as a real site effect.
- The sample names in the stored object are patient IDs, so **sequencing plate could not be checked**. No site structure does not
  rule out a plate effect.
- The PCA shows structure in the data. It does not show that this structure causes the inflation in section 9.2, because
  shuffled labels are unrelated to subtype.

### 9.8 Does group size or subtype mixing explain the inflation?

Under the null, about 5% of genes should have a raw p-value below 0.05. This fraction is steadier than the number of calls, which
varied from 407 to 1,179 for DESeq2 across draws. `tcga/04b_permutation_groupsize.R` changes the number of pseudo-mutants and
controls, and `tcga/04d_permutation_single_subtype.R` draws both groups from one subtype (LumA, wild-type tumors only, so no
subtype term is needed). The controls are 43 rather than 40 because the script rounds up within each subtype.

| Setting (subtype-adjusted design) | Draws | DESeq2 | edgeR | limma-voom | Wilcoxon |
|---|---|---|---|---|---|
| 11 cases, 150 controls (section 9.2) | 10 | 0.128 | 0.113 | 0.061 | 0.047 |
| 11 cases, 43 controls | 3 | 0.125 | 0.113 | 0.063 | 0.049 |
| 22 cases, 153 controls | 3 | 0.134 | 0.126 | 0.059 | 0.040 |
| LumA only, 11 cases, 150 controls (no subtype term) | 3 | 0.109 | 0.093 | 0.048 | 0.049 |

Mean number of calls at 5% FDR: DESeq2 513, 539, 522 and 217; edgeR 389, 366, 604 and 234; limma-voom 8, 0.3, 0 and 1.7;
Wilcoxon 0 in all four.

- Doubling the cases or cutting the controls to about a quarter barely changes the raw p-value fractions. The inflation is
  therefore **not a small-sample effect of the kind seen in section 6**, which should shrink as groups grow.
- Holding subtype constant lowered the DESeq2 and edgeR fractions somewhat (0.128 to 0.109 and 0.113 to 0.093), but they stay
  about twice the expected rate. Subtype mixing may contribute; it does not explain the whole effect.
- Each extra setting has only 3 draws and individual draws differ a lot (LumA DESeq2: 0.145, 0.090, 0.092), so these are first
  looks. A 400-control run was planned and not done.

### 9.9 Per-gene p-value diagnostics on one shuffled draw

`tcga/04c_pvalue_diagnostics.R` repeats draw 1 of section 9.2 (same seed; its raw p < 0.05 fractions of 0.106, 0.095, 0.042 and
0.040 match the original run) and saves every gene's p-value, expression level and DESeq2 dispersion. Raw p < 0.05 by expression
quintile (Q1 lowest):

| Quintile | DESeq2 | edgeR | limma-voom | Wilcoxon |
|---|---|---|---|---|
| Q1 | 0.103 | 0.116 | 0.046 | 0.054 |
| Q2 | 0.128 | 0.113 | 0.044 | 0.043 |
| Q3 | 0.138 | 0.116 | 0.042 | 0.038 |
| Q4 | 0.087 | 0.072 | 0.040 | 0.033 |
| Q5 | 0.072 | 0.059 | 0.036 | 0.031 |

- The excess is not confined to the lowest-count genes. It is largest in the low and middle quintiles and smaller, though
  still present for DESeq2, in the most highly expressed genes.
- Median DESeq2 dispersion (final estimate) was higher for genes with raw p < 0.05 than for the others within every quintile:
  1.24 against 1.04 (Q1), 1.09 against 0.67 (Q2), 0.94 against 0.58 (Q3), 0.28 against 0.245 (Q4) and 0.21 against 0.17 (Q5).
  The false calls therefore come disproportionately from genes whose counts vary more between tumors.
- With about 160 samples, DESeq2's final dispersion is almost identical to the per-gene estimate, so shrinkage toward the trend
  does little here. This argues against the section 4 explanation (unreliable per-gene dispersion at small n) for this inflation.
- DESeq2 and edgeR rank genes almost identically (Spearman correlation of p-values 0.96); limma-voom and Wilcoxon agree more with
  each other (0.81) than with either (0.52 to 0.60).
- This is **one draw, and an association, not a demonstrated mechanism**. One explanation is that the single dispersion per gene
  shared by cases and controls misjudges the variance of a small, mostly basal case group; another is that noisy genes catch
  extreme samples by chance. Neither was tested. Draw 1 was also milder than average (262 DESeq2 calls against a mean of 513).

### 9.10 Positive control: basal against LumA

To check that the pipeline finds real signal, `tcga/02_confounded_de.R` compared 100 basal with 100 LumA tumors (basal = case,
no confounder), and `tcga/08_volcano_upset.R` plots the saved table (`results/figures/volcano_basal_vs_lumA.png` and
`upset_basal_vs_lumA.png`). Wilcoxon results are not in this table, so the overlap covers three tools.

- At q < 0.05, DESeq2, edgeR and limma-voom called 15,151, 14,960 and 14,545 genes; 13,304 genes were called by all three. The
  other groups were: DESeq2 and edgeR only 1,438; limma-voom only 951; DESeq2 and limma-voom only 235; DESeq2 only 174;
  edgeR only 163; edgeR and limma-voom only 55.
- Known markers have the expected direction (positive log fold change = higher in basal): FOXC1 +4.0, SOX10 +4.5, KRT17 +2.6,
  KRT5 +2.7 and KRT14 +2.2 are higher in basal; ESR1 -6.0, PGR -6.5, FOXA1 -5.4, GATA3 -3.6 and XBP1 -3.2 are higher in LumA.
  All have DESeq2 q-values below 1e-10.
- The over-calling on shuffled labels is therefore a calibration problem, not a failure to detect large real differences.
- About half of all tested genes are called, so "q < 0.05" says little about which tool is more liberal; the permutation
  results are the better evidence for that. The log fold change column is used as saved in the file, and I did not check which
  tool it comes from.

## 10. Limitations

- **Somatic mutations only.** Germline BRCA carriers need controlled-access data (dbGaP) and are missing. "Wild-type" tumors may
  be functionally BRCA-deficient (germline carriage, BRCA1 promoter methylation), which weakens any comparison.
- **Few mutants** (11 BRCA1, 7 BRCA2 with subtype labels). Power is low, so absence of calls is not absence of effect.
- **Simulation realism.** The simulations draw from an NB model, favor NB-based tools, and use true effects of at least 0.6 log2
  in 10% of genes, so power figures do not transfer to real BRCA effects.
- **The cause of the real-data inflation is only partly explained:** surrogate variables reduced the false calls but did not remove them (section 9.6), holding subtype constant lowered them slightly (section 9.8), and group size does not explain them. Outlier samples and tumor purity were not examined.
- **Few draws in the follow-up tests.** The group-size and single-subtype settings use 3 draws each and the per-gene diagnostics use one, so they are first looks.
- **PCA colours by site, not plate.** The stored sample names are patient IDs (section 9.7).
- **Simulation parameters.** Dispersion and effect size were not varied in the simulations, only sample size and case counts.
- **Additive adjustment.** `subtype + mutation` assumes the mutation effect is the same in every subtype.
- **Control sampling.** Real-data DE runs used one random draw of about 300 controls; the null used 150. DESeq2 was run with its
  default independent filtering for adjusted p-values.
- **Follow-up data quality.** Median follow-up is short (about 2.3 years), 21 patients have a recorded time of zero or below, and
  follow-up times come from the latest GDC follow-up record, which may not equal the true last contact.
- **Replicate counts are modest** (5 balanced and 20 unbalanced simulated datasets per cell, 10 permutation draws).
- Not done: a per-dataset Snakemake workflow, adding the newer diagnostic scripts (sections 9.7 to 9.10) to the workflow, a
  test of outliers and tumor purity as causes of the inflation, a sweep over dispersion and effect size, and pathway enrichment.

## 11. Lessons from the analysis

- **Normalization matters as much as the test.** The t-test and Wilcoxon were first run on library-size-normalized log-CPM, and
  their FDR looked inflated (for example, Wilcoxon adjusted in the mixed 11-case scenario: 0.19). Switching to median-of-ratios
  size factors brought it to 0.034 with unchanged power. Part of the scratch NB GLM's apparent inflation also came from this.
- **A silent data bug.** The first survival model reported 144 patients and 143 events out of 1,095 tumors. That ratio showed that
  patients who were alive had been dropped: follow-up times were missing because the GDC stores them in a different field
  (`follow_ups`) than the one first queried. Fixing the query gave 1,036 patients with 144 events, and the sign of the BRCA1
  hazard ratio flipped (0.88 to 1.15), illustrating the bias of analyzing only the patients who died.
- **Naive designs can look fine by chance.** With few cases, the naive FDR depends on how the cases fall across subtypes in a given draw.
- **Test an explanation before accepting it.** The small-sample dispersion story from the simulations (section 4) was a natural explanation for the real-data inflation. Changing the group sizes did not change the inflation, so that explanation does not fit, and the cause remains open.
- **Use the raw p-value fraction when call counts are noisy.** The share of genes with p < 0.05 on shuffled labels (expected 0.05) was far steadier across draws than the number of FDR calls.

## 12. Reproducibility

The same steps are available as a Snakemake workflow (`workflow/Snakefile`; see the README). The commands below run the scripts directly.

```
conda env create -f environment.yml && conda activate de-fdr-tcga-brca
pytest -q tests/

# balanced simulation
python simulation/simulate.py --reps 5 --n_genes 3000
python simulation/run_python_methods.py results/sim
Rscript simulation/run_r_methods.R results/sim
python simulation/evaluate.py results/sim
python simulation/calibration_check.py

# unbalanced / confounded simulation
python simulation/simulate_unbalanced.py --reps 20
Rscript simulation/run_r_unbalanced.R results/sim_unbal
python simulation/run_fast_unbalanced.py results/sim_unbal
python simulation/evaluate_unbalanced.py results/sim_unbal

# TCGA-BRCA (needs about 5 GB and a few hours of downloading and compute)
Rscript tcga/01_download_prepare.R     # includes follow-up times (01b_fix_clinical.R only repairs older objects)
Rscript tcga/02_confounded_de.R
Rscript tcga/04_permutation_null.R 10 adjusted BRCA1_mut 150
Rscript tcga/05_wilcoxon_check.R 50 10 150
Rscript tcga/06_sva_null.R 5 0 150            # about 35 minutes per draw
Rscript tcga/03_survival_expression.R
Rscript tcga/03b_km_plot.R
Rscript tcga/07_pca_qc.R                                              # PCA by subtype and tissue source site
Rscript tcga/04b_permutation_groupsize.R 3 adjusted BRCA1_mut 150 2   # args: draws, design, target, controls, case-scale
Rscript tcga/04c_pvalue_diagnostics.R 1                               # per-gene p-values of one shuffled draw
Rscript tcga/04d_permutation_single_subtype.R 3 LumA 11 150           # args: draws, subtype, cases, controls
Rscript tcga/08_volcano_upset.R                                       # volcano and overlap plots, basal vs LumA

# figures
python simulation/make_figures.py
```
