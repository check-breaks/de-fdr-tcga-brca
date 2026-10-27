#!/usr/bin/env Rscript
# DESeq2, edgeR (QL F-test), limma-voom on every simulated dataset.
# NOTE: written without an R runtime available when authored -- run it and fix any
# version-specific issues. Usage: Rscript simulation/run_r_methods.R results/sim
suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(limma)})
args <- commandArgs(TRUE); root <- if (length(args)) args[1] else "results/sim"
for (d in list.dirs(root, recursive = FALSE)) {
  counts <- as.matrix(read.csv(file.path(d, "counts.csv"), row.names = 1, check.names = FALSE))
  meta <- read.csv(file.path(d, "meta.csv")); meta$condition <- factor(meta$condition, levels = c("A", "B"))
  keep <- rowSums(counts >= 5) >= 2                       # same filter as the Python methods
  cm <- counts[keep, ]
  res <- data.frame(gene = rownames(counts), p_deseq2 = NA_real_, p_edger = NA_real_, p_limmavoom = NA_real_)
  rownames(res) <- res$gene

  dds <- DESeqDataSetFromMatrix(cm, meta, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res[rownames(cm), "p_deseq2"] <- results(dds, name = "condition_B_vs_A")$pvalue   # keep DESeq2 independent filtering OFF-by-NA? see below
  # DESeq2 returns NA padj for independent-filtered genes but pvalue is untouched, so p-values are comparable.

  y <- DGEList(cm); y <- calcNormFactors(y)
  des <- model.matrix(~ condition, meta)
  y <- estimateDisp(y, des); fit <- glmQLFit(y, des)
  res[rownames(cm), "p_edger"] <- glmQLFTest(fit, coef = 2)$table$PValue

  v <- voom(y, des); vf <- eBayes(lmFit(v, des))
  res[rownames(cm), "p_limmavoom"] <- topTable(vf, coef = 2, number = Inf, sort.by = "none")$P.Value

  write.csv(res, file.path(d, "pvals_r.csv"), row.names = FALSE)
  cat(basename(d), "done\n")
}
