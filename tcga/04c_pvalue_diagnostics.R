#!/usr/bin/env Rscript
# Per-gene p-value diagnostics on ONE shuffled-label draw (BRCA1 pseudo-mutants, adjusted design).
# Question: are the false calls from DESeq2/edgeR concentrated in low-count genes, or spread across all genes?
# Usage: Rscript tcga/04c_pvalue_diagnostics.R [iter=1]    (iter 1 uses the same draw as iter01 of script 04)
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(edgeR); library(limma)})
args <- commandArgs(TRUE)
iter <- if (length(args) >= 1) as.integer(args[1]) else 1
target <- "BRCA1_mut"; n_ctrl <- 150
dir.create("results/tcga/diag", recursive = TRUE, showWarnings = FALSE)

se <- readRDS("data/tcga/brca_se.rds")
se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
cd <- as.data.frame(colData(se)); cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
counts_all <- assay(se, "unstranded")
case_n <- table(cd$subtype[cd[[target]]])
pool <- which(!cd$BRCA_mut)

# same shuffling code as script 04
set.seed(1000 + iter)
cases <- unlist(lapply(names(case_n), function(s) {
  if (case_n[[s]] == 0) return(NULL)
  p <- pool[cd$subtype[pool] == s]
  p[sample.int(length(p), case_n[[s]])]
}))
rest <- setdiff(pool, cases)
per  <- split(rest, cd$subtype[rest]); per <- per[lengths(per) > 0]
ctrl <- unlist(lapply(per, function(ix)
  ix[sample.int(length(ix), min(length(ix), ceiling(n_ctrl * length(ix) / length(rest))))]))
idx <- unname(c(cases, ctrl))
m <- data.frame(subtype = droplevels(cd$subtype[idx]),
                grp = factor(rep(c("case", "ctrl"), c(length(cases), length(ctrl))), levels = c("ctrl", "case")))
cat("cases:", length(cases), "controls:", length(ctrl), "\n")

fml <- ~ subtype + grp
des <- model.matrix(fml, m)
cm  <- counts_all[, idx]
keep <- filterByExpr(DGEList(cm), design = des); cm <- cm[keep, ]
cat("genes tested:", nrow(cm), "\n")

# DESeq2 (also keep its dispersion estimates)
dds <- DESeq(DESeqDataSetFromMatrix(cm, m, fml), quiet = TRUE)
d <- results(dds, contrast = c("grp", "case", "ctrl"))
disp_gene <- mcols(dds)$dispGeneEst; disp_fit <- mcols(dds)$dispFit; disp_final <- mcols(dds)$dispersion
rm(dds); gc()

# edgeR quasi-likelihood
y <- calcNormFactors(DGEList(cm)); yd <- estimateDisp(y, des)
e <- glmQLFTest(glmQLFit(yd, des), coef = "grpcase")$table

# limma-voom
l <- topTable(eBayes(lmFit(voom(yd, des), des)), coef = "grpcase", number = Inf, sort.by = "none")

# Wilcoxon on log-CPM residuals after removing subtype
lc <- cpm(y, log = TRUE, prior.count = 1)
lc <- t(lm.fit(model.matrix(~ subtype, m), t(lc))$residuals)
is_case <- m$grp == "case"; n1 <- sum(is_case); n2 <- sum(!is_case)
rk <- t(apply(lc, 1, rank))
U  <- rowSums(rk[, is_case, drop = FALSE]) - n1 * (n1 + 1) / 2
z  <- (U - n1 * n2 / 2) / sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
pw <- 2 * pnorm(-abs(z))

df <- data.frame(gene = rownames(cm), logCPM = e$logCPM, baseMean = d$baseMean,
                 p_deseq2 = d$pvalue, p_edger = e$PValue, p_limmavoom = l$P.Value, p_wilcoxon = pw,
                 disp_gene = disp_gene, disp_fit = disp_fit, disp_final = disp_final)
write.csv(df, sprintf("results/tcga/diag/pvalues_iter%02d.csv", iter), row.names = FALSE)

# Fraction of genes with raw p < 0.05, by expression quintile (Q1 = lowest expressed)
tools <- c(deseq2 = "p_deseq2", edgeR = "p_edger", voom = "p_limmavoom", wilcoxon = "p_wilcoxon")
bins <- cut(df$logCPM, quantile(df$logCPM, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))
mat <- sapply(tools, function(cn) as.vector(tapply(df[[cn]] < 0.05, bins, mean, na.rm = TRUE)))
rownames(mat) <- levels(bins)
cat("\nFraction of genes with raw p < 0.05, by expression quintile (expect ~0.05):\n"); print(round(mat, 3))
cat("\nOverall fraction raw p < 0.05:\n"); print(round(sapply(tools, function(cn) mean(df[[cn]] < 0.05, na.rm = TRUE)), 3))
cat("\nSpearman correlation between tools' p-values:\n")
print(round(cor(df[, tools], method = "spearman", use = "pairwise.complete.obs"), 2))
write.csv(mat, sprintf("results/tcga/diag/frac_p05_by_expression_iter%02d.csv", iter))

png(sprintf("results/figures/pvalue_diagnostics_iter%02d.png", iter), width = 1600, height = 800, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
matplot(1:5, mat, type = "b", pch = 19, lty = 1, col = c("red", "orange", "blue", "darkgreen"),
        xaxt = "n", ylim = c(0, max(mat) * 1.1), xlab = "Expression quintile (1 = lowest)",
        ylab = "Fraction of genes with raw p < 0.05", main = "Shuffled labels: by expression level")
axis(1, 1:5, levels(bins)); abline(h = 0.05, lty = 2)
legend("topleft", legend = names(tools), col = c("red", "orange", "blue", "darkgreen"), pch = 19, lty = 1, cex = 0.8)
hist(df$p_deseq2, breaks = 20, col = "grey80", main = "DESeq2 raw p-values", xlab = "p-value")
dev.off()
cat("\nsaved results/tcga/diag/ and results/figures/pvalue_diagnostics_iter", sprintf("%02d", iter), ".png\n", sep = "")
