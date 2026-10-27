#!/usr/bin/env Rscript
# BRCA-mutant vs wild-type DE: naive design vs subtype-adjusted, across DESeq2 / edgeR / limma-voom.
# Real-data confounding case study: subtype explains much of the "mutation signature".
suppressPackageStartupMessages({library(SummarizedExperiment); library(DESeq2); library(edgeR); library(limma)
  library(dplyr); library(ggplot2)})
dir.create("results/tcga", recursive = TRUE, showWarnings = FALSE)
se <- readRDS("data/tcga/brca_se.rds")
se <- se[, !is.na(se$BRCA_Subtype_PAM50)]
counts <- assay(se, "unstranded"); cd <- as.data.frame(colData(se))
cd$mut <- factor(ifelse(cd$BRCA_mut, "mut", "wt"), levels = c("wt", "mut"))
cd$subtype <- factor(cd$BRCA_Subtype_PAM50)
keep <- filterByExpr(DGEList(counts), group = cd$mut); counts <- counts[keep, ]

run_all <- function(design_formula, tag) {
  des <- model.matrix(design_formula, cd); coef <- "mutmut"
  dds <- DESeq(DESeqDataSetFromMatrix(counts, cd, design_formula), quiet = TRUE)
  d <- results(dds, name = coef)
  y <- calcNormFactors(DGEList(counts)); y <- estimateDisp(y, des)
  e <- glmQLFTest(glmQLFit(y, des), coef = coef)$table
  v <- eBayes(lmFit(voom(y, des), des)); l <- topTable(v, coef = coef, number = Inf, sort.by = "none")
  tab <- data.frame(gene = rownames(counts), deseq2 = d$padj, edger = p.adjust(e$PValue, "BH"),
                    limmavoom = p.adjust(l$P.Value, "BH"), lfc = d$log2FoldChange, design = tag)
  write.csv(tab, sprintf("results/tcga/de_%s.csv", tag), row.names = FALSE); tab
}
naive <- run_all(~ mut, "naive")
adj   <- run_all(~ subtype + mut, "subtype_adjusted")

n_sig <- function(t) sapply(c("deseq2", "edger", "limmavoom"), function(m) sum(t[[m]] < 0.05, na.rm = TRUE))
summary_tab <- rbind(naive = n_sig(naive), subtype_adjusted = n_sig(adj))
print(summary_tab); write.csv(summary_tab, "results/tcga/n_sig_naive_vs_adjusted.csv")
cat("Mutant n =", sum(cd$mut == "mut"), "; WT n =", sum(cd$mut == "wt"), "\n")
# Agreement between tools (feed to an UpSet plot): sets of genes with q<0.05 in each tool
