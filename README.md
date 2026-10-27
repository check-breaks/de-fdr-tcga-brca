# Do DE tools control FDR? Simulation + TCGA-BRCA case study (BRCA1/BRCA2)

One repo, one story: (1) implement the core statistics from scratch, (2) test whether RNA-seq DE methods
control FDR in simulation, (3) apply the lessons to TCGA-BRCA, where **tumor subtype confounds BRCA1/2 mutation status**.

## Status
| Part | Status |
|---|---|
| `src/scratch/` BH, permutation test, NB GLM (IRLS) | **Done, tested** (`python tests/test_scratch.py`) |
| `simulation/` simulator, t-test/Wilcoxon/scratch-NB, evaluator | **Done, run** (figure in `results/figures/fdr_power.png`) |
| `simulation/run_r_methods.R` DESeq2 / edgeR / limma-voom | Written, **not yet run** (needs R) |
| `tcga/` download, confounded DE, survival | Written, **not yet run** (needs R + internet) |
| `workflow/` Snakemake | Day 6 |

## Run order
```
conda env create -f environment.yml && conda activate de-fdr-tcga-brca
python tests/test_scratch.py                     # or: pytest -q  (also runs the statsmodels cross-check)
python simulation/simulate.py --reps 5 --n_genes 3000
python simulation/run_python_methods.py results/sim
Rscript simulation/run_r_methods.R results/sim
python simulation/evaluate.py results/sim        # FDR/power table + figure, now with the R tools
python simulation/calibration_check.py
Rscript tcga/01_download_prepare.R               # GDC download (large; start early)
Rscript tcga/02_confounded_de.R
Rscript tcga/03_survival_expression.R
```

## Findings so far (Python-side only)
- Scratch NB GLM is correctly calibrated with the *true* dispersion (FPR ~5% at every n) but anti-conservative with
  gene-wise ML dispersion at small n (FPR 18.5% at n=3/group, 7.6% at n=10). This motivates DESeq2/edgeR dispersion shrinkage.
- Wilcoxon has zero power at 3 vs 3 (minimum two-sided p = 0.1). Naive tests lose power at small n.
- Reported FDR is 0 when no genes are called; read it together with power.
- Simulation uses 4 replicate datasets per n and 2,000 genes: noisy. Raise `--reps` for the final figures.

## TCGA caveats (state these in the write-up)
- Open-access MAFs contain **somatic** mutations only; germline BRCA carriers need dbGaP-controlled access.
- Few mutant tumors (expect tens, not hundreds), so power is low. Link this to the power curves above.
- BRCA1 is also silenced by promoter methylation; "wild-type" tumors may be functionally BRCA1-deficient.
- Verify GDC/TCGAbiolinks column names (PAM50, clinical fields); they change between releases.

## Methods and background
See [docs/methods_background.md](docs/methods_background.md) for the statistical background, simulation design and results.
