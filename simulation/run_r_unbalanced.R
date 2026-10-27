#!/usr/bin/env Rscript
# DESeq2 / edgeR (QL) / limma-voom on unbalanced sims, naive (~condition) and adjusted (~subtype + condition).
# Written without an R runtime available; run on one folder first. Usage: Rscript simulation/run_r_unbalanced.R results/sim_unbal
suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(limma)})
args <- commandArgs(TRUE); root <- if (length(args)) args[1] else "results/sim_unbal"
dirs <- list.dirs(root, recursive = FALSE)
if (!length(dirs)) stop("No dataset folders in ", root, " - run simulate_unbalanced.py first")
for (d in dirs) {
  out <- file.path(d, "pvals_r.csv"); if (file.exists(out)) next        # resumable
  counts <- as.matrix(read.csv(file.path(d, "counts.csv"), row.names = 1, check.names = FALSE))
  meta <- read.csv(file.path(d, "meta.csv"))
  meta$condition <- factor(meta$condition, levels = c("ctrl", "case")); meta$subtype <- factor(meta$subtype)
  keep <- filterByExpr(DGEList(counts), design = model.matrix(~ subtype + condition, meta))   # same genes for both designs
  cm <- counts[keep, ]
  res <- data.frame(gene = rownames(counts)); rownames(res) <- res$gene
  for (design in c("naive", "adjusted")) {
    fml <- if (design == "naive") ~ condition else ~ subtype + condition
    des <- model.matrix(fml, meta)
    dds <- DESeq(DESeqDataSetFromMatrix(cm, meta, fml), quiet = TRUE)
    res[rownames(cm), paste0("p_deseq2_", design)] <- results(dds, contrast = c("condition", "case", "ctrl"))$pvalue
    y <- calcNormFactors(DGEList(cm)); y <- estimateDisp(y, des)
    res[rownames(cm), paste0("p_edger_", design)] <- glmQLFTest(glmQLFit(y, des), coef = "conditioncase")$table$PValue
    v <- eBayes(lmFit(voom(y, des), des))
    res[rownames(cm), paste0("p_limmavoom_", design)] <- topTable(v, coef = "conditioncase", number = Inf, sort.by = "none")$P.Value
  }
  write.csv(res, out, row.names = FALSE); cat(basename(d), "done\n")
}
